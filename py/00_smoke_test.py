"""
00 — Prueba de humo. Ejecútalo ANTES que nada.

Comprueba en ~2 minutos que los seis servicios responden y que tu cuenta
de Copernicus autentica correctamente. Es preferible descubrir aquí que
un endpoint cambió de nombre, y no a los 40 minutos de una descarga.

No descarga datos: solo pide metadatos y una tesela mínima.
"""

from __future__ import annotations

import sys

from common import (
    URL_AQ_API,
    URL_AQ_STATIONS,
    URL_NUTS,
    URL_IED,
    URL_EUROSTAT,
    URL_OPENEO,
    SESSION,
    TIMEOUT,
    log,
)

OK, FALLO = "  [OK]  ", "  [FALLO]"
resultados: dict[str, bool] = {}


def prueba(nombre: str, fn) -> None:
    try:
        detalle = fn()
        print(f"{OK} {nombre}: {detalle}")
        resultados[nombre] = True
    except Exception as exc:  # noqa: BLE001
        print(f"{FALLO} {nombre}: {type(exc).__name__}: {exc}")
        resultados[nombre] = False


def _arcgis(url: str) -> str:
    r = SESSION.get(url, params={"f": "pjson"}, timeout=TIMEOUT)
    r.raise_for_status()
    j = r.json()
    capas = j.get("layers", [])
    return (
        f"{len(capas)} capa(s), maxRecordCount={j.get('maxRecordCount')}, "
        f"v{j.get('currentVersion')}"
    )


def _aq_api() -> str:
    # /Country devuelve la lista de países disponibles: la petición más
    # barata que confirma que el servicio de descarga está vivo.
    r = SESSION.get(f"{URL_AQ_API}/Country", timeout=TIMEOUT)
    r.raise_for_status()
    datos = r.json()
    codigos = [d.get("countryCode") for d in datos] if isinstance(datos, list) else []
    if "ES" not in codigos:
        raise RuntimeError(f"'ES' no aparece entre los países: {codigos[:10]}")
    return f"{len(codigos)} países, 'ES' disponible"


def _eurostat() -> str:
    r = SESSION.get(
        f"{URL_EUROSTAT}/demo_r_d3dens",
        params={"format": "JSON", "lang": "EN", "time": "2023", "geo": "ES300"},
        timeout=TIMEOUT,
    )
    r.raise_for_status()
    j = r.json()
    v = list(j.get("value", {}).values())
    return f"Madrid (ES300) densidad {v[0] if v else 'sin dato'} hab/km2"


def _openeo() -> str:
    import openeo  # import local: solo se necesita aquí

    con = openeo.connect(URL_OPENEO)
    con.authenticate_oidc()
    cols = {c["id"] for c in con.list_collections()}
    if "SENTINEL_5P_L2" not in cols:
        candidatos = [c for c in cols if "5P" in c or "5p" in c]
        raise RuntimeError(
            f"No veo SENTINEL_5P_L2. Colecciones parecidas: {candidatos}"
        )
    return f"autenticado, {len(cols)} colecciones, SENTINEL_5P_L2 disponible"


def main() -> None:
    print("\n=== Prueba de humo: servicios y credenciales ===\n")
    prueba("EEA — API de descarga de calidad del aire", _aq_api)
    prueba("EEA — estaciones (ArcGIS)", lambda: _arcgis(URL_AQ_STATIONS))
    prueba("EEA — NUTS 2021 (ArcGIS)", lambda: _arcgis(URL_NUTS))
    prueba("EEA — IED SiteMap (ArcGIS)", lambda: _arcgis(URL_IED))
    prueba("Eurostat — demo_r_d3dens", _eurostat)
    prueba("Copernicus — openEO (abre el navegador)", _openeo)

    fallidos = [k for k, v in resultados.items() if not v]
    print("\n" + "=" * 60)
    if fallidos:
        print(f"{len(fallidos)} servicio(s) con problemas:")
        for f in fallidos:
            print(f"   - {f}")
        print(
            "\nSi falla un ArcGIS, suele ser mantenimiento del EEA: reintenta "
            "en una hora.\nSi falla openEO, borra el token cacheado en "
            "~/.local/share/openeo-python-client/ y repite."
        )
        sys.exit(1)
    print("Todo correcto. Puedes lanzar 01_nuts.py .. 07_corine.py")


if __name__ == "__main__":
    main()
