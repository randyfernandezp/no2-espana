"""
04c — Depuración del patrón de puntos industrial.

Entrada: data/raw/ied_es.geojson        (103 299 registros brutos)
Salida : data/interim/ied_es_limpio.geojson
         data/interim/ied_diagnostico.csv

POR QUÉ HACE FALTA ESTE PASO
----------------------------
La descarga bruta del servicio IED_SiteMap trae 103 299 registros para
España. Ese número NO es el número de instalaciones industriales, y usarlo
tal cual invalidaría por completo el análisis de procesos puntuales.

Problema 1 — DUPLICACIÓN POR AÑO DE DECLARACIÓN
   El campo Site_reporting_year indica que el servicio publica una fila por
   sitio Y por año. Un proceso puntual espacial requiere un punto por
   entidad: con duplicados, la intensidad estimada lambda(s) se multiplica
   por el número de años declarados y la función K queda dominada por pares
   a distancia CERO (el mismo sitio consigo mismo en otro año). El resultado
   sería una agregación espuria masiva.
   -> Se conserva el registro del año más reciente por InspireSiteId.

Problema 2 — SECTORES IRRELEVANTES PARA EL NO2
   El desglose por eprtr_sectors del volcado bruto es revelador:
       INTENSIVE LIVESTOCK   49 936   <- 48 % del total
       MINERALS              10 537
       WASTE AND WASTEWATER  10 130
       METALS                 9 697
       CHEMICALS              7 112
       FOOD AND BEVERAGE      6 935
       ENERGY                 3 015
   Casi la mitad son explotaciones ganaderas intensivas (30 453 de porcino,
   7 673 de avicultura). Su emisión característica es amoníaco y metano,
   NO óxidos de nitrógeno de combustión. Incluirlas equivale a modelar la
   distribución de la ganadería española y llamarla "fuentes de NO2".

   El NOx procede de procesos de COMBUSTIÓN a alta temperatura. Los sectores
   relevantes son, por orden de importancia:
       ENERGY    (centrales térmicas, refinerías, cogeneración)
       MINERALS  (cemento, cal, vidrio, cerámica: hornos a >1400 C)
       METALS    (siderurgia, coquerías, sinterización)
       CHEMICALS (producción de ácido nítrico, craqueo)
       WASTE     (incineradoras)

ESTRATEGIA DE FILTRADO
----------------------
Se aplica en dos niveles y se documenta el efecto de cada uno:

  Nivel A (preferente): filtrar por el campo `pollutants`, que lista los
      contaminantes efectivamente declarados por el sitio. Es el criterio
      empírico: el propio operador declara si emite NOx.

  Nivel B (respaldo): si `pollutants` no resulta utilizable, filtrar por
      lista blanca de sectores de combustión.

Se guardan AMBAS versiones. Comparar el patrón obtenido con cada criterio
es en sí un análisis de sensibilidad que conviene reportar: si las
conclusiones del proceso puntual cambian según el criterio, hay que decirlo.

NOTA SOBRE EL SESGO RESULTANTE
------------------------------
Tras el filtrado, el patrón representa "grandes instalaciones de combustión
sujetas a la Directiva de Emisiones Industriales", no "toda fuente de NOx".
Quedan fuera el tráfico rodado (la fuente dominante de NO2 en España), la
calefacción residencial y las pymes bajo umbral. Esto NO es un defecto del
dato: es la definición de la población objetivo, y debe enunciarse así en
el informe.
"""

from __future__ import annotations

import json
from collections import Counter

import pandas as pd

from common import RAW, INTERIM, ANIO, log

# Sectores cuya actividad principal implica combustión a alta temperatura
SECTORES_NOX = {
    "ENERGY",
    "MINERALS",
    "METALS",
    "CHEMICALS",
    "WASTE AND WASTEWATER",
    "PAPER AND WOOD",
}

# Actividades del Anexo I especialmente intensivas en NOx
ACTIVIDADES_NOX = ["1(a)", "1(b)", "1(c)", "1(d)", "2(", "3(", "4(", "5(a)", "5(b)"]


def cargar() -> tuple[list[dict], dict]:
    ruta = RAW / "ied_es.geojson"
    if not ruta.exists():
        raise SystemExit(f"Falta {ruta}. Ejecuta antes 04_industrias.py")
    fc = json.loads(ruta.read_text(encoding="utf-8"))
    return fc["features"], fc.get("crs", {})


def main() -> None:
    feats, crs = cargar()
    log(f"Registros brutos: {len(feats):,}")

    df = pd.DataFrame([f["properties"] for f in feats])
    df["_idx"] = range(len(df))

    diagnostico = [("bruto", len(df))]

    # ---------------------------------------------------------------- #
    # Paso 1 — determinar la declaración de NOx SOBRE TODOS LOS AÑOS
    # ---------------------------------------------------------------- #
    # ORDEN IMPORTANTE: primero se evalúa qué declara cada sitio a lo largo
    # de toda su historia, y solo después se deduplica.
    #
    # Hacerlo al revés (deduplicar quedándose con el año más reciente y
    # luego mirar `pollutants`) introduce un sesgo de completitud: los
    # sitios cuyo último registro es 2024 se juzgan por un año que aún
    # está parcialmente reportado, mientras que los que terminan en 2023
    # se juzgan por un año cerrado. Dos sitios idénticos recibirían
    # tratamiento distinto según cuándo dejaron de declarar.
    PATRON_NOX = r"nitrogen ox|\bNOX\b|\bNO2\b"

    if "pollutants" in df.columns:
        df["_declara_nox"] = df["pollutants"].astype(str).str.contains(
            PATRON_NOX, case=False, na=False, regex=True
        )

        if "Site_reporting_year" in df.columns:
            log("\n  Sitios que declaran NOx por año de reporte:")
            porano = (
                df[df["_declara_nox"]]
                .groupby("Site_reporting_year")["InspireSiteId"]
                .nunique()
                .sort_index()
            )
            for a, n in porano.items():
                aviso = ""
                if n < 0.6 * porano.max():
                    aviso = "  <- cobertura baja, año probablemente incompleto"
                log(f"     {int(a)}: {n:>5}{aviso}")

        # Un sitio es fuente de NOx si lo declaró en CUALQUIER año.
        # Justificación: la localización industrial es muy estable; que una
        # planta quede un año por debajo del umbral de 100 t no la convierte
        # en no-fuente. Lo que modelamos es la ubicación de la capacidad
        # emisora, no la emisión de un ejercicio concreto.
        # --- Coherencia temporal con el año de las concentraciones ------ #
        # La serie de declarantes cae de forma monótona a lo largo de 18 años
        # (339 en 2007 -> 203 en 2023). Ese descenso NO es un artefacto de
        # reporte incompleto: refleja el cierre de centrales de carbón y la
        # instalación de sistemas de reducción catalítica en cementeras y
        # siderurgia. Solo el último año puede estar a medio consolidar.
        #
        # Consecuencia metodológica: el conjunto "declaró NOx alguna vez"
        # mezcla fuentes activas en 2007 con las activas en 2023. Una térmica
        # cerrada en 2012 no contribuye al NO2 medido en 2023, y meterla en
        # el patrón de puntos introduce ruido en la estimación de intensidad
        # y sesga a la baja cualquier asociación con las concentraciones.
        #
        # Se generan tres conjuntos con criterios temporales distintos:
        #   ANIO           -> coherencia estricta con el año del contaminante
        #   ventana 5 años -> fuentes persistentes, más robusto a que una
        #                     instalación quede un año bajo el umbral
        #   histórico      -> máxima cobertura, para sensibilidad
        VENTANA = 5
        if "Site_reporting_year" in df.columns:
            anio_col = pd.to_numeric(df["Site_reporting_year"], errors="coerce")

            nox_anio = set(
                df.loc[df["_declara_nox"] & (anio_col == ANIO), "InspireSiteId"]
                .dropna().unique()
            )
            nox_ventana = set(
                df.loc[
                    df["_declara_nox"] & (anio_col >= ANIO - VENTANA + 1)
                    & (anio_col <= ANIO),
                    "InspireSiteId",
                ].dropna().unique()
            )
            log(f"\n  Declaran NOx en {ANIO}: {len(nox_anio):,}")
            log(f"  Declaran NOx en {ANIO - VENTANA + 1}-{ANIO}: {len(nox_ventana):,}")
            globals()["_CONJUNTOS_TEMP"] = {
                "anio": nox_anio,
                "ventana": nox_ventana,
            }

        nox_alguna_vez = set(
            df.loc[df["_declara_nox"], "InspireSiteId"].dropna().unique()
        )
        log(f"\n  Sitios que declararon NOx en algún año: {len(nox_alguna_vez):,}")

    # ---------------------------------------------------------------- #
    # Paso 2 — deduplicar por sitio (año más reciente)
    # ---------------------------------------------------------------- #
    if "Site_reporting_year" in df.columns and "InspireSiteId" in df.columns:
        anios = sorted(int(a) for a in df["Site_reporting_year"].dropna().unique())
        log(f"\nAños presentes: {anios[0]}-{anios[-1]} ({len(anios)} años)")
        log(f"Sitios únicos: {df['InspireSiteId'].nunique():,}")

        antes = len(df)
        df = (
            df.sort_values("Site_reporting_year", ascending=False)
            .drop_duplicates(subset="InspireSiteId", keep="first")
            .copy()
        )
        log(f"Deduplicado: {antes:,} -> {len(df):,} registros")
        diagnostico.append(("tras deduplicar por sitio", len(df)))

    # ---------------------------------------------------------------- #
    # Paso 3A — patrón por declaración de NOx (criterio empírico)
    # ---------------------------------------------------------------- #
    df_a = None
    if "pollutants" in df.columns:
        muestra = df["pollutants"].dropna().head(2).tolist()
        log(f"\nEjemplos de 'pollutants': {muestra}")
        df_a = df[df["InspireSiteId"].isin(nox_alguna_vez)].copy()
        log(f"Nivel A — declara NOx en algún año: {len(df_a):,}")
        diagnostico.append(("nivel A: declara NOx", len(df_a)))

    # ---------------------------------------------------------------- #
    # Paso 3B — patrón por sector de combustión (criterio de respaldo)
    # ---------------------------------------------------------------- #
    df_b = None
    if "eprtr_sectors" in df.columns:
        def relevante(v) -> bool:
            return isinstance(v, str) and any(s in v.upper() for s in SECTORES_NOX)

        df_b = df[df["eprtr_sectors"].apply(relevante)].copy()
        log(f"Nivel B — sectores de combustión: {len(df_b):,}")
        diagnostico.append(("nivel B: sector de combustión", len(df_b)))

        log("\n  Composición del nivel B:")
        for sec, n in Counter(
            s.strip() for v in df_b["eprtr_sectors"].dropna() for s in str(v).split(",")
        ).most_common(10):
            log(f"    {n:>6}  {sec}")

    # Intersección: sector de combustión Y declara NOx. Es el núcleo duro.
    if df_a is not None and df_b is not None:
        inter = set(df_a["InspireSiteId"]) & set(df_b["InspireSiteId"])
        log(f"\n  Intersección A ∩ B: {len(inter):,} sitios")
        solo_a = set(df_a["InspireSiteId"]) - set(df_b["InspireSiteId"])
        if solo_a:
            log(f"  Declaran NOx pero fuera de sectores de combustión: {len(solo_a):,}")
            fuera = df_a[df_a["InspireSiteId"].isin(solo_a)]
            for sec, n in Counter(fuera["eprtr_sectors"].dropna()).most_common(5):
                log(f"      {n:>4}  {sec}")

    # ---------------------------------------------------------------- #
    # Paso 3 — elegir el patrón principal y guardar
    # ---------------------------------------------------------------- #
    # Preferimos el criterio empírico si dio un número razonable.
    if df_a is not None and 50 <= len(df_a) <= 20000:
        principal, etiqueta = df_a, "nivel A (declara NOx)"
    elif df_b is not None:
        principal, etiqueta = df_b, "nivel B (sector de combustión)"
    else:
        principal, etiqueta = df, "sin filtrar"

    log(f"\nPatrón principal: {etiqueta} -> {len(principal):,} instalaciones")

    def guardar(sub: pd.DataFrame, nombre: str) -> None:
        if sub is None or sub.empty:
            return
        sel = [feats[i] for i in sub["_idx"]]
        salida = {
            "type": "FeatureCollection",
            "crs": crs,
            "features": sel,
        }
        ruta = INTERIM / nombre
        ruta.write_text(json.dumps(salida), encoding="utf-8")
        log(f"  -> {ruta.name} ({len(sel):,} puntos)")

    temporales = globals().get("_CONJUNTOS_TEMP", {})
    if temporales and "InspireSiteId" in df.columns:
        for clave, conjunto in temporales.items():
            sub = df[df["InspireSiteId"].isin(conjunto)]
            guardar(sub, f"ied_es_nox_{clave}.geojson")
            diagnostico.append((f"nivel A restringido a {clave}", len(sub)))

    guardar(principal, "ied_es_limpio.geojson")
    guardar(df_a, "ied_es_nivelA_declara_nox.geojson")
    guardar(df_b, "ied_es_nivelB_sector.geojson")

    pd.DataFrame(diagnostico, columns=["paso", "n_registros"]).to_csv(
        INTERIM / "ied_diagnostico.csv", index=False
    )

    log(
        "\nGuardadas las dos versiones. Comparar el patrón que produce cada\n"
        "criterio es un análisis de sensibilidad reportable: si las\n"
        "conclusiones cambian según el filtro, hay que decirlo en el informe."
    )


if __name__ == "__main__":
    main()
