"""
04b — Masas de NOx emitidas, reconstruidas por unión en DiscoData.

Fuente : EEA DiscoData, base [IED].[v1r2]
Salida : data/raw/eprtr_emisiones.csv

EL ESQUEMA REAL (verificado con 04d_explorar.py)
------------------------------------------------
E-PRTR está normalizado y ninguna tabla tiene por sí sola lo que hace falta.
La cadena es:

  PollutantRelease            pollutant='NOX', mediumCode='AIR',
                              totalPollutantQuantityKg, facilityReportId
        |  facilityReportId -> Id
  ProductionFacilityReport    countryCode, reportingYear, localId, namespace
        |  (localId, namespace) -> (localId, namespace)
  ProductionFacility          coordenadas, siteId, NUTS3, actividad
        |  siteId -> id
  ProductionSite              localId, namespace del SITIO

El identificador INSPIRE se compone como  namespace + '/' + localId, que es
el formato que usa la capa IED_SiteMap de la que salió nuestro patrón de
puntos ("ES.CAED/12345.SITE").

DOS CAUTELAS
------------
1. ProductionFacility guarda una fila por instalación Y AÑO. Unir sin
   restringir el año duplicaría las emisiones tantas veces como años haya
   declarado la instalación.
2. Un sitio puede albergar varias instalaciones. Las masas se agregan por
   sitio antes de unir con el patrón de puntos, que está a nivel de sitio.

Se guardan además las coordenadas, para poder unir espacialmente si el
cruce por identificador fallara.
"""

from __future__ import annotations

import time
from urllib.parse import quote

import pandas as pd
import requests

from common import RAW, PAIS, ANIO, log

URL_SQL = "https://discodata.eea.europa.eu/sql"
B = "[IED].[v1r2]"
VENTANA_INI, VENTANA_FIN = ANIO - 4, ANIO   # 2019-2023, como en 04c


def consultar(sql: str, pagina: int = 1, por_pagina: int = 1000) -> list[dict]:
    url = f"{URL_SQL}?query={quote(sql)}&p={pagina}&nrOfHits={por_pagina}"
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    j = r.json()
    if "errors" in j:
        raise RuntimeError(j["errors"][0].get("error", j["errors"]))
    return j.get("results", [])


def descargar(sql: str, etiqueta: str) -> pd.DataFrame:
    todas: list[dict] = []
    for pagina in range(1, 150):
        filas = consultar(sql, pagina=pagina)
        if not filas:
            break
        todas.extend(filas)
        if pagina % 5 == 0:
            log(f"    {etiqueta}: {len(todas):,}")
        if len(filas) < 1000:
            break
        time.sleep(0.25)
    log(f"  {etiqueta}: {len(todas):,} filas")
    return pd.DataFrame(todas)


SQL = f"""
SELECT
    pfr.countryCode        AS pais,
    pfr.reportingYear      AS anio,
    pr.pollutant           AS contaminante,
    pr.mediumCode          AS medio,
    pr.totalPollutantQuantityKg AS cantidad_kg,
    pf.facilityName        AS instalacion,
    pf.localId             AS fac_localId,
    pf.namespace           AS fac_namespace,
    pf.EPRTRAnnexIMainActivity AS actividad,
    pf.NUTS3               AS nuts3,
    pf.x_4258              AS lon,
    pf.y_4258              AS lat,
    ps.localId             AS site_localId,
    ps.namespace           AS site_namespace,
    ps.siteName            AS sitio
FROM {B}.[PollutantRelease] pr
INNER JOIN {B}.[ProductionFacilityReport] pfr
        ON pr.facilityReportId = pfr.Id
INNER JOIN {B}.[ProductionFacility] pf
        ON pf.localId = pfr.localId
       AND pf.namespace = pfr.namespace
       AND pf.reportingYear = pfr.reportingYear
LEFT JOIN {B}.[ProductionSite] ps
        ON ps.id = pf.siteId
WHERE pr.pollutant = 'NOX'
  AND pr.mediumCode = 'AIR'
  AND pfr.countryCode = '{PAIS}'
  AND pfr.reportingYear BETWEEN {VENTANA_INI} AND {VENTANA_FIN}
""".strip()


def main() -> None:
    log("Consulta con unión en cuatro tablas:")
    print(SQL)
    df = descargar(SQL, "emisiones NOx")

    if df.empty:
        log("Sin filas. Probando sin restricción de año …")
        df = descargar(
            SQL.replace(
                f"  AND pfr.reportingYear BETWEEN {VENTANA_INI} AND {VENTANA_FIN}", ""
            ),
            "emisiones NOx (todos los años)",
        )
    if df.empty:
        raise SystemExit("La unión no devolvió filas. Revisa el esquema.")

    df["anio"] = pd.to_numeric(df["anio"], errors="coerce")
    df["cantidad_kg"] = pd.to_numeric(df["cantidad_kg"], errors="coerce")
    df = df[df["cantidad_kg"] > 0]

    log("\nDeclaraciones de NOx a aire por año:")
    for a, n in df.groupby("anio").size().sort_index().items():
        log(f"    {int(a)}: {n:>5} instalaciones")

    # ------------------------------------------------------------------ #
    # Identificador INSPIRE del sitio
    # ------------------------------------------------------------------ #
    tiene_sitio = df["site_localId"].notna()
    log(f"\nFilas con sitio asignado: {tiene_sitio.sum():,} de {len(df):,}")

    df["InspireSiteId"] = (
        df["site_namespace"].fillna(df["fac_namespace"]).astype(str)
        + "/"
        + df["site_localId"].fillna(df["fac_localId"]).astype(str)
    )
    log("Ejemplos de identificador construido:")
    for v in df["InspireSiteId"].dropna().unique()[:4]:
        log(f"    {v}")

    # ------------------------------------------------------------------ #
    # Agregación
    # ------------------------------------------------------------------ #
    # Primero se suman las instalaciones dentro de cada sitio y año, porque
    # un sitio puede albergar varias. Después se promedia entre años: sumar
    # los años contaría la misma capacidad emisora varias veces, cuando lo
    # que interesa es la magnitud anual característica.
    por_anio = (
        df.groupby(["InspireSiteId", "anio"], dropna=True)["cantidad_kg"]
        .sum()
        .reset_index()
    )
    agg = (
        por_anio.groupby("InspireSiteId")
        .agg(nox_kg=("cantidad_kg", "mean"),
             nox_kg_max=("cantidad_kg", "max"),
             n_anios=("anio", "nunique"),
             ultimo_anio=("anio", "max"))
        .reset_index()
    )
    agg["nox_t"] = agg["nox_kg"] / 1000.0

    # Coordenadas y metadatos, por si el cruce por identificador falla y hay
    # que recurrir a una unión espacial.
    meta = (
        df.sort_values("anio", ascending=False)
        .drop_duplicates("InspireSiteId")[
            ["InspireSiteId", "sitio", "instalacion", "actividad",
             "nuts3", "lon", "lat"]
        ]
    )
    agg = agg.merge(meta, on="InspireSiteId", how="left")

    destino = RAW / "eprtr_emisiones.csv"
    agg.to_csv(destino, index=False)
    log(f"\n-> {destino.name}: {len(agg)} sitios con NOx declarado")

    print("\nDistribución de la marca (t NOx/año):")
    print(agg["nox_t"].describe(percentiles=[0.5, 0.75, 0.9, 0.95, 0.99]))

    tot = agg["nox_t"].sum()
    log(f"\nNOx industrial total declarado: {tot:,.0f} t/año")
    log(f"Las 10 mayores concentran el "
        f"{100 * agg.nlargest(10, 'nox_t')['nox_t'].sum() / tot:.1f} %")

    print("\nLos 12 mayores emisores:")
    print(
        agg.nlargest(12, "nox_t")[["sitio", "actividad", "nuts3", "nox_t"]]
        .to_string(index=False, max_colwidth=38)
    )

    ratio = agg["nox_t"].quantile(0.99) / agg["nox_t"].median()
    log(f"\nRazón percentil 99 / mediana: {ratio:.0f}")
    if ratio > 50:
        log("Asimetría extrema: usa log(nox_t) como marca en el análisis de")
        log("correlación, o unas pocas térmicas dominarán la superficie.")


if __name__ == "__main__":
    main()
