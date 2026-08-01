"""
01 — Geometría de regiones administrativas NUTS 2021 para España.

Fuente : EEA, servicio NUTS_RG_01M_2021_3035_EU27 (MapServer)
Salida : data/raw/nuts3_es.geojson  (EPSG:3035)

ESTRUCTURA REAL DEL SERVICIO (verificada)
-----------------------------------------
Capa 0 : NUTS_RG_01M_2021_3035_EU27_NUTS3_cntr  -> polígonos NUTS-3
Capa 1 : NUTS_RG_01M_2021_3035_EU27_NUTS0_cntr  -> polígonos NUTS-0 (países)

Es decir, el servicio del EEA solo publica los niveles 3 y 0. NO hay
NUTS-2. Los nombres de campo van en minúsculas: cntr_code, levl_code,
nuts_id, nuts_name.

CONSTRUCCIÓN DE NUTS-2
----------------------
NUTS-2 se obtiene disolviendo NUTS-3 por los cuatro primeros caracteres
del código: ES300 (Madrid provincia) -> ES30 (Comunidad de Madrid). La
jerarquía NUTS está codificada en el propio identificador, así que la
agregación es exacta, no aproximada. Lo necesitamos para el chequeo de
MAUP del Objetivo 2.
"""

import json

from common import (
    URL_NUTS,
    RAW,
    PAIS,
    arcgis_service_info,
    arcgis_layer_fields,
    arcgis_query_paged,
    buscar_campo,
    guardar_geojson,
    log,
)


def main() -> None:
    info = arcgis_service_info(URL_NUTS)
    capas = {c["id"]: c["name"] for c in info.get("layers", [])}
    log(f"Capas disponibles: {capas}")

    # Localizar la capa de NUTS-3 por su nombre
    capa_n3 = next(
        (lid for lid, nom in capas.items() if "NUTS3" in nom.upper()), None
    )
    if capa_n3 is None:
        capa_n3 = min(capas)
        log(f"No hallé capa 'NUTS3' por nombre; uso la capa {capa_n3}")

    campos = arcgis_layer_fields(URL_NUTS, capa_n3)
    nombres = {f["name"] for f in campos}
    log(f"Campos de la capa {capa_n3}: {sorted(nombres)}")
    (RAW / "nuts_campos.json").write_text(
        json.dumps(campos, indent=2), encoding="utf-8"
    )

    campo_pais = buscar_campo(nombres, "cntr_code", "cntr_id", "country")
    campo_id = buscar_campo(nombres, "nuts_id", "nuts_code")
    campo_nivel = buscar_campo(nombres, "levl_code", "level")

    if campo_pais is None or campo_id is None:
        raise SystemExit(f"Campos no identificados. Disponibles: {sorted(nombres)}")

    where = f"{campo_pais}='{PAIS}'"
    if campo_nivel:
        where += f" AND {campo_nivel}=3"

    log(f"Descargando NUTS-3 con where = {where}")
    feats = list(arcgis_query_paged(URL_NUTS, capa_n3, where=where))

    if not feats:
        log("0 registros con filtro de país. Reintento filtrando por nuts_id …")
        feats = list(
            arcgis_query_paged(
                URL_NUTS, capa_n3, where=f"{campo_id} LIKE '{PAIS}%'"
            )
        )

    # Normalizamos el nombre del identificador a NUTS_ID para que los
    # scripts de R no dependan de la capitalización del servicio.
    for f in feats:
        p = f["properties"]
        p["NUTS_ID"] = p.get(campo_id)
        p["NUTS_NAME"] = p.get(buscar_campo(nombres, "nuts_name", "name_latn"))

    guardar_geojson(feats, RAW / "nuts3_es.geojson")

    codigos = sorted(f["properties"]["NUTS_ID"] for f in feats)
    log(f"NUTS-3 descargadas: {len(codigos)}")
    log(f"Primeras: {codigos[:8]}")
    log(f"Últimas:  {codigos[-8:]}")

    # España tiene 59 NUTS-3 (50 provincias + 2 ciudades autónomas +
    # las 7 islas canarias contadas aparte según la versión).
    if not 55 <= len(codigos) <= 62:
        log(f"AVISO: se esperaban ~59 regiones, hay {len(codigos)}. Revisa el filtro.")

    # NUTS-2 por disolución del código
    n2 = sorted({c[:4] for c in codigos})
    log(f"NUTS-2 derivables por agregación: {len(n2)} -> {n2}")
    (RAW / "nuts2_desde_nuts3.json").write_text(
        json.dumps({"nuts2": n2}, indent=2), encoding="utf-8"
    )
    log(
        "NUTS-2 no se descarga: el servicio solo publica niveles 3 y 0. "
        "Se construirá en R disolviendo NUTS-3 por substr(NUTS_ID, 1, 4)."
    )


if __name__ == "__main__":
    main()
