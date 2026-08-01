"""
02 — Metadatos de estaciones de medición de NO2 en España.

Fuente : EEA, AirQualityDownloadServiceEUMonitoringStations (MapServer/0)
Salida : data/raw/estaciones_es.geojson  (EPSG:3035)

POR QUÉ ESTE PASO ES EL MÁS IMPORTANTE DEL PROYECTO
---------------------------------------------------
Este servicio trae, por estación, dos atributos que deciden la validez de
todo el análisis geoestadístico:

  * TIPO DE ESTACIÓN   (background / traffic / industrial)
  * TIPO DE ÁREA       (rural / suburban / urban)

Una estación de tráfico mide a escala de bocacalle: el NO2 cae un 60-80 %
en los primeros 50 m desde el borde de la calzada. Ese valor NO es una
realización de un campo aleatorio suave a escala regional. Si la metes en
el variograma obtendrás un nugget que devora la varianza estructural y una
superficie kriged sin sentido.

Estrategia adoptada:
  - Campo geoestadístico (Objetivo 1) -> SOLO estaciones background.
  - Estaciones de tráfico e industriales -> conjunto de validación externa
    y análisis de residuos (¿el modelo subestima donde hay tráfico? sí,
    y cuantificar ese sesgo es un resultado publicable).

Es lo mismo que hace el EEA para sus mapas de exposición europea.
"""

import json

from common import (
    URL_AQ_STATIONS,
    RAW,
    PAIS,
    CONTAMINANTE,
    arcgis_layer_fields,
    arcgis_query_paged,
    guardar_geojson,
    log,
)

LAYER = 0


def main() -> None:
    campos = arcgis_layer_fields(URL_AQ_STATIONS, LAYER)
    nombres = {f["name"] for f in campos}

    log("Campos de la capa AQ Stations:")
    for f in campos:
        log(f"    {f['name']:<40} {f['type']:<22} {f['alias']}")
    (RAW / "estaciones_campos.json").write_text(
        json.dumps(campos, indent=2), encoding="utf-8"
    )

    campo_pais = next(
        (
            c
            for c in (
                "CountryOrTerritory",
                "Countrycode",
                "CountryCode",
                "ISO2",
                "country",
            )
            if c in nombres
        ),
        None,
    )
    campo_poll = next(
        (
            c
            for c in ("AirPollutant", "Pollutant", "AirPollutantCode", "pollutant")
            if c in nombres
        ),
        None,
    )

    if campo_pais is None:
        raise SystemExit(
            "No identifiqué el campo de país. Mira data/raw/estaciones_campos.json "
            "y fija `campo_pais` manualmente."
        )

    # El campo de país puede contener el código ISO2 o el nombre completo.
    clausulas = [f"({campo_pais}='{PAIS}' OR {campo_pais}='Spain')"]
    if campo_poll:
        clausulas.append(f"{campo_poll} LIKE '%{CONTAMINANTE}%'")
    where = " AND ".join(clausulas)

    log(f"Consultando con where = {where}")
    feats = list(arcgis_query_paged(URL_AQ_STATIONS, LAYER, where=where))

    if not feats:
        log("0 registros. Reintentando sin filtro de contaminante …")
        feats = list(
            arcgis_query_paged(URL_AQ_STATIONS, LAYER, where=clausulas[0])
        )

    guardar_geojson(feats, RAW / "estaciones_es.geojson")

    # Resumen de tipologías: comprueba aquí que existen las tres clases.
    if feats:
        props = feats[0]["properties"].keys()
        log(f"Atributos por estación: {sorted(props)}")
        for clave in props:
            if "Type" in clave or "Area" in clave or "Classification" in clave:
                valores = {f["properties"].get(clave) for f in feats}
                log(f"  {clave}: {sorted(v for v in valores if v)}")


if __name__ == "__main__":
    main()
