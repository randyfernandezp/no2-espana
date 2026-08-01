"""
02b — Registro de estaciones con su CLASIFICACIÓN, vía DiscoData.

Fuente : EEA DiscoData, API SQL pública (https://discodata.eea.europa.eu/sql)
         Base 'Air Quality Data Flows', versión v1r1
Salida : data/raw/estaciones_es_meta.csv
         data/raw/discodata_catalogo.json

HISTORIAL DE ESTE SCRIPT (para el apartado de trazabilidad del informe)
----------------------------------------------------------------------
1. El fichero de metadatos de cmshare resultó ser un PDF, no datos.
2. DiscoData prohíbe consultar sys.databases e INFORMATION_SCHEMA
   ("system tables are not allowed"), así que no se puede descubrir el
   esquema por programación.
3. El explorador web de DiscoData reveló la base y sus tablas. Cuidado:
   si el navegador traduce la página, el nombre aparece castellanizado
   ("FlujosDatosCalidad del Aire"); el identificador real está en inglés.

Tablas visibles en la base: AirQualityStatistics, AssessmentRegimeMethods,
AssessmentRegimes, AttainmentMethods, Attainments, Measurements, Measures,
Models, Plans, Scenarios, SourceApportionment, Zones.

Las dos con metadatos de punto de muestreo son AirQualityStatistics
(estadísticos anuales por estación y contaminante, que es la fuente del
servicio ArcGIS del EEA) y Measurements (flujo D: métodos de evaluación
y descripción de los puntos de muestreo).

LO QUE BUSCAMOS
---------------
  AirQualityStationType : background | traffic | industrial
  AirQualityStationArea : urban | suburban | rural | rural-remote

Sin ellas no se separa el campo regional de fondo del recargo urbano por
tráfico, y el Objetivo 1 pierde su fundamento metodológico.

BONUS
-----
AirQualityStatistics contiene además la media anual ya calculada y validada
por el EEA. Sirve como CONTRASTE INDEPENDIENTE de las medias que calcula
03_series_no2.py a partir de los parquet diarios: si ambas coinciden, el
tratamiento de banderas de validez y el umbral de captura son correctos.
"""

from __future__ import annotations

import json
import time
from urllib.parse import quote

import pandas as pd
import requests

from common import RAW, PAIS, ANIO, log

URL_SQL = "https://discodata.eea.europa.eu/sql"

# El nombre mostrado en el explorador puede venir traducido por el navegador.
# Se prueban las variantes plausibles del identificador real.
BASES = [
    "AirQualityDataFlows",
    "Air Quality Data Flows",
    "AirQuality_DataFlows",
    "AQDataFlows",
    "Airquality_Dissem",
]
VERSIONES = ["v1r1", "latest", "v1"]
TABLAS = ["AirQualityStatistics", "Measurements", "AssessmentRegimes"]


def consultar(sql: str, pagina: int = 1, por_pagina: int = 1000) -> list[dict]:
    url = f"{URL_SQL}?query={quote(sql)}&p={pagina}&nrOfHits={por_pagina}"
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    payload = r.json()
    if "errors" in payload:
        raise RuntimeError(payload["errors"][0].get("error", payload["errors"]))
    return payload.get("results", [])


def sondear() -> list[dict]:
    """Prueba combinaciones base x version x tabla hasta que una responda."""
    hallazgos = []
    log("Sondeando combinaciones base / versión / tabla …")
    for base in BASES:
        for ver in VERSIONES:
            for tab in TABLAS:
                fqn = f"[{base}].[{ver}].[{tab}]"
                try:
                    filas = consultar(f"SELECT TOP 1 * FROM {fqn}", por_pagina=1)
                except Exception as exc:  # noqa: BLE001
                    msg = str(exc)[:70]
                    if "Invalid object name" not in msg:
                        log(f"  [!!] {fqn}: {msg}")
                    continue
                if not filas:
                    continue
                cols = list(filas[0].keys())
                bajas = {c.lower() for c in cols}
                tiene = any("stationtype" in c for c in bajas) and any(
                    "stationarea" in c for c in bajas
                )
                log(f"  [ok] {fqn}: {len(cols)} columnas"
                    f"{'  <- CON CLASIFICACIÓN' if tiene else ''}")
                log(f"       {cols}")
                hallazgos.append({"fqn": fqn, "columnas": cols, "clasifica": tiene})
                time.sleep(0.3)
            # Si ya encontramos algo en esta base/versión, no seguir probando otras
            if any(h["clasifica"] for h in hallazgos):
                return hallazgos
    return hallazgos


def descargar(fqn: str, cols: list[str]) -> pd.DataFrame:
    idx = {c.lower(): c for c in cols}
    col_pais = next(
        (idx[k] for k in ("countrycode", "country_code", "country") if k in idx), None
    )
    col_poll = next(
        (idx[k] for k in ("airpollutant", "pollutant", "airpollutantcode") if k in idx),
        None,
    )
    col_anio = next(
        (idx[k] for k in ("reportingyear", "year", "statisticsyear") if k in idx), None
    )

    partes = []
    if col_pais:
        partes.append(f"{col_pais} = '{PAIS}'")
    if col_poll:
        partes.append(f"{col_poll} LIKE '%NO2%'")
    if col_anio:
        partes.append(f"{col_anio} = {ANIO}")
    donde = " AND ".join(partes) if partes else "1=1"

    sql = f"SELECT * FROM {fqn} WHERE {donde}"
    log(f"Descargando:\n  {sql}")

    todas: list[dict] = []
    for pagina in range(1, 100):
        filas = consultar(sql, pagina=pagina, por_pagina=1000)
        if not filas:
            break
        todas.extend(filas)
        log(f"  página {pagina}: +{len(filas)} (total {len(todas)})")
        if len(filas) < 1000:
            break
        time.sleep(0.3)

    df = pd.DataFrame(todas)
    if df.empty and col_anio:
        log(f"Sin filas para {ANIO}. Reintento sin filtro de año …")
        partes = [p for p in partes if not p.startswith(str(col_anio))]
        sql = f"SELECT * FROM {fqn} WHERE {' AND '.join(partes) or '1=1'}"
        df = pd.DataFrame(consultar(sql, por_pagina=1000))
    return df


def main() -> None:
    hallazgos = sondear()
    (RAW / "discodata_catalogo.json").write_text(
        json.dumps(hallazgos, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    if not hallazgos:
        raise SystemExit(
            "Ninguna combinación respondió. En el explorador web de DiscoData, "
            "haz clic sobre una tabla: el nombre completo aparece escrito en el "
            "cuadro de consulta. Pásame esa línea literal."
        )

    buenos = [h for h in hallazgos if h["clasifica"]] or hallazgos
    elegido = buenos[0]
    log(f"Tabla elegida: {elegido['fqn']}")

    df = descargar(elegido["fqn"], elegido["columnas"])
    if df.empty:
        raise SystemExit("La tabla responde pero no hay filas para España.")

    destino = RAW / "estaciones_es_meta.csv"
    df.to_csv(destino, index=False, encoding="utf-8")
    log(f"-> {destino.name} ({len(df)} filas, {len(df.columns)} columnas)")

    col_tipo = next((c for c in df.columns if "stationtype" in c.lower()), None)
    col_area = next((c for c in df.columns if "stationarea" in c.lower()), None)
    col_eoi = next((c for c in df.columns if "eoicode" in c.lower()), df.columns[0])

    if col_tipo and col_area:
        print("\n" + "=" * 64)
        print("TABLA CRUZADA: tipo de estación x tipo de área")
        print("=" * 64)
        print(pd.crosstab(df[col_tipo], df[col_area], margins=True))
        fondo = df[df[col_tipo].astype(str).str.lower().str.contains("background")]
        print(f"\nEstaciones DE FONDO (base del variograma): {fondo[col_eoi].nunique()}")
        print(f"Estaciones totales con NO2: {df[col_eoi].nunique()}")
        print(
            "\nSi las de fondo bajan de ~60, habrá que relajar el criterio de\n"
            "inclusión antes de calcular el variograma."
        )
    else:
        log(f"Sin columnas de clasificación. Columnas: {list(df.columns)}")


if __name__ == "__main__":
    main()
