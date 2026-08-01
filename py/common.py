"""
Utilidades compartidas para la capa de descarga.

Convenciones del proyecto
-------------------------
CRS de trabajo : EPSG:3035 (ETRS89 / LAEA Europe). Metros, equiareal.
                 Obligatorio para variogramas, kriging y procesos puntuales:
                 EPSG:4326 daría distancias en grados, sin sentido físico.
CRS de servicio: los ArcGIS REST del EEA sirven en EPSG:3857. Pedimos
                 siempre outSR=3035 para no reproyectar después.
Año de estudio : 2023 (último año con datos E1a verificados completos).
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Iterator

import requests

# --------------------------------------------------------------------------- #
# Configuración global
# --------------------------------------------------------------------------- #

PAIS = "ES"
ANIO = 2023
CONTAMINANTE = "NO2"
CRS_TRABAJO = 3035

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"

for _d in (RAW, INTERIM, PROCESSED):
    _d.mkdir(parents=True, exist_ok=True)

# Endpoints verificados (julio 2026)
URL_AQ_API = "https://eeadmz1-downloads-api-appservice.azurewebsites.net"
URL_AQ_STATIONS = (
    "https://air.discomap.eea.europa.eu/arcgis/rest/services/AirQuality/"
    "AirQualityDownloadServiceEUMonitoringStations/MapServer"
)
URL_NUTS = (
    "https://air.discomap.eea.europa.eu/arcgis/rest/services/EPRTR/"
    "NUTS_RG_01M_2021_3035_EU27/MapServer"
)
URL_IED = (
    "https://air.discomap.eea.europa.eu/arcgis/rest/services/Air/"
    "IED_SiteMap/MapServer"
)
URL_CLC = (
    "https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/"
    "CLC2018_WM/MapServer"
)
URL_EUROSTAT = (
    "https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data"
)
URL_OPENEO = "https://openeo.dataspace.copernicus.eu"

TIMEOUT = 120
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "proyecto-no2-espana/1.0 (academico)"})


# --------------------------------------------------------------------------- #
# Cliente ArcGIS REST
# --------------------------------------------------------------------------- #

def arcgis_service_info(base_url: str) -> dict[str, Any]:
    """Metadatos del MapServer: lista de capas con su id y nombre."""
    r = SESSION.get(base_url, params={"f": "pjson"}, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def arcgis_layer_fields(base_url: str, layer_id: int) -> list[dict[str, str]]:
    """Campos de una capa. Úsalo SIEMPRE antes de escribir un `where`:
    los nombres de campo del EEA cambian entre versiones del servicio."""
    r = SESSION.get(
        f"{base_url}/{layer_id}", params={"f": "pjson"}, timeout=TIMEOUT
    )
    r.raise_for_status()
    info = r.json()
    return [
        {"name": f["name"], "type": f["type"], "alias": f.get("alias", "")}
        for f in info.get("fields", [])
    ]


def buscar_campo(nombres: set[str], *candidatos: str) -> str | None:
    """
    Busca un campo por nombre, sin distinguir mayúsculas.

    Los servicios del EEA no son consistentes: el servicio NUTS expone
    'cntr_code' en minúsculas mientras que el de estaciones usa
    'CountryCode' en CamelCase. Comparar literalmente falla en uno de los dos.
    """
    indice = {n.lower(): n for n in nombres}
    for c in candidatos:
        if c.lower() in indice:
            return indice[c.lower()]
    return None


def arcgis_query_paged(
    base_url: str,
    layer_id: int,
    where: str = "1=1",
    out_fields: str = "*",
    out_sr: int = CRS_TRABAJO,
    page_size: int = 1000,
    geometry: bool = True,
) -> Iterator[dict[str, Any]]:
    """
    Consulta paginada de una capa ArcGIS REST.

    Los servicios del EEA tienen MaxRecordCount = 2000. Sin paginación
    truncarías silenciosamente el resultado: es el error clásico que
    convierte 3 500 instalaciones industriales en 2 000 y arruina
    el análisis de intensidad del proceso puntual.

    Devuelve features GeoJSON de una en una.
    """
    offset = 0
    while True:
        params = {
            "where": where,
            "outFields": out_fields,
            "outSR": out_sr,
            "returnGeometry": str(geometry).lower(),
            "resultOffset": offset,
            "resultRecordCount": page_size,
            "f": "geojson",
        }
        r = SESSION.get(
            f"{base_url}/{layer_id}/query", params=params, timeout=TIMEOUT
        )
        r.raise_for_status()
        payload = r.json()

        if "error" in payload:
            raise RuntimeError(f"ArcGIS devolvió error: {payload['error']}")

        feats = payload.get("features", [])
        if not feats:
            break

        yield from feats

        # `exceededTransferLimit` es la señal fiable de que hay más páginas.
        if not payload.get("exceededTransferLimit") and len(feats) < page_size:
            break
        offset += len(feats)
        time.sleep(0.3)  # cortesía con el servidor del EEA


def guardar_geojson(features: list[dict], destino: Path, crs: int = CRS_TRABAJO) -> None:
    """Escribe un FeatureCollection con el CRS declarado explícitamente."""
    fc = {
        "type": "FeatureCollection",
        "crs": {
            "type": "name",
            "properties": {"name": f"urn:ogc:def:crs:EPSG::{crs}"},
        },
        "features": features,
    }
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(json.dumps(fc), encoding="utf-8")
    print(f"  -> {destino.relative_to(ROOT)}  ({len(features)} registros)")


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
