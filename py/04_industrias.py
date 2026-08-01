"""
04 — Instalaciones industriales que reportan bajo IED / E-PRTR en España.

Fuente : EEA, Air/IED_SiteMap (MapServer) — ubicación y sector
         Portal industry.eea.europa.eu — masas emitidas de NOx
Salida : data/raw/ied_es.geojson       (EPSG:3035)
         data/raw/ied_campos.json

Este es el patrón de puntos del Objetivo 3. Dos advertencias que determinan
la validez del modelo de intensidad:

1) SESGO DE UMBRAL DE NOTIFICACIÓN. E-PRTR solo obliga a declarar si se
   superan umbrales por contaminante (100 t/año de NOx). El patrón NO es el
   de "todas las fuentes industriales", sino el de "fuentes por encima del
   umbral". Es un proceso puntual *truncado por la marca*. Hay que decirlo
   explícitamente en el informe: sesga hacia grandes instalaciones de
   combustión y subrepresenta polígonos de pymes.

2) VENTANA DE OBSERVACIÓN. La ventana del proceso puntual debe ser la
   frontera peninsular, no el bounding box. Con bounding box, la estimación
   de intensidad asigna masa al mar y al territorio francés/portugués, y
   la K de Ripley queda sesgada por efecto borde mal corregido.
   La corrección de borde se aplica en R (spatstat), con la ventana real.

Sobre las marcas
----------------
La masa de NOx emitida por instalación se une aquí como MARCA del patrón.
Un patrón de puntos sin marcas trata igual a una central térmica de 20 000 t
y a una fábrica de cerámica de 110 t. Ponderar por masa es lo que hace
posible el Objetivo 4 (superficie de presión de emisión).
"""

from __future__ import annotations

import json

from common import (
    URL_IED,
    RAW,
    PAIS,
    arcgis_service_info,
    arcgis_layer_fields,
    arcgis_query_paged,
    guardar_geojson,
    log,
)


def main() -> None:
    info = arcgis_service_info(URL_IED)
    capas = {c["id"]: c["name"] for c in info.get("layers", [])}
    log(f"Capas del servicio IED_SiteMap: {capas}")
    log(f"MaxRecordCount del servicio: {info.get('maxRecordCount')}")

    inventario = {}
    for lid, nombre in capas.items():
        campos = arcgis_layer_fields(URL_IED, lid)
        inventario[f"{lid}:{nombre}"] = campos
        log(f"\nCapa {lid} — {nombre}")
        for f in campos:
            log(f"    {f['name']:<38} {f['type']:<22} {f['alias']}")

    (RAW / "ied_campos.json").write_text(
        json.dumps(inventario, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    # Elegimos la capa de sitios/instalaciones (la de mayor granularidad
    # con geometría de punto).
    layer = min(capas) if capas else 0
    nombres = {f["name"] for f in inventario[f"{layer}:{capas[layer]}"]}

    campo_pais = next(
        (
            c
            for c in (
                "countryCode",
                "CountryCode",
                "countryName",
                "CNTR_CODE",
                "country",
            )
            if c in nombres
        ),
        None,
    )
    if campo_pais is None:
        raise SystemExit(
            "No identifiqué el campo de país en la capa IED. "
            "Consulta data/raw/ied_campos.json y fíjalo a mano."
        )

    where = f"({campo_pais}='{PAIS}' OR {campo_pais}='Spain')"
    log(f"\nDescargando instalaciones con where = {where}")
    feats = list(arcgis_query_paged(URL_IED, layer, where=where, page_size=1000))
    guardar_geojson(feats, RAW / "ied_es.geojson")

    if feats:
        props = feats[0]["properties"]
        log(f"Atributos disponibles: {sorted(props)}")
        # Resumen por sector: comprueba que hay un campo de actividad NACE/EPRTR
        for clave in props:
            if any(k in clave.lower() for k in ("sector", "activity", "nace")):
                valores = {}
                for f in feats:
                    v = f["properties"].get(clave)
                    valores[v] = valores.get(v, 0) + 1
                top = sorted(valores.items(), key=lambda kv: -kv[1])[:12]
                log(f"\n  Top sectores según '{clave}':")
                for v, n in top:
                    log(f"    {n:>5}  {v}")

    log(
        "\nSIGUIENTE PASO MANUAL: descarga la base tabular de emisiones desde\n"
        "  https://industry.eea.europa.eu/industrial-emissions/dataset\n"
        "y guárdala como data/raw/eprtr_emisiones.csv. La unión con este\n"
        "GeoJSON se hace por `facilityInspireId` en R/10_prepare.R.\n"
        "Filtra: pollutant = 'Nitrogen oxides (NOX/NO2)', medium = 'AIR', "
        f"reportingYear = 2023."
    )


if __name__ == "__main__":
    main()
