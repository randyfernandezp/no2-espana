"""
08 — Verificación cruzada de las medias anuales.

Entradas: data/interim/no2_anual_estacion.csv   (calculado por 03)
          data/raw/estaciones_es_meta.csv       (oficial EEA, vía DiscoData)
Salida  : data/interim/verificacion_medias.csv
          data/interim/estaciones_modelado.csv

POR QUÉ ESTE PASO
-----------------
El script 03 reconstruye las medias anuales a partir de 4,6 millones de
observaciones horarias, aplicando por el camino tres decisiones: el filtro
de validez, la agregación horaria->diaria y el umbral de captura del 75 %.
Cada una puede introducir sesgo si se implementa mal, y el resultado
seguiría pareciendo razonable.

Por suerte no hace falta confiar: la tabla AirQualityStatistics del EEA
contiene la media anual YA CALCULADA Y VALIDADA por el propio organismo
para cada punto de muestreo. Comparar ambas series es una verificación
independiente en el sentido fuerte: si coinciden, las tres decisiones son
correctas; si divergen sistemáticamente, el patrón de la discrepancia
señala cuál falló.

QUÉ SE COMPRUEBA ADEMÁS
-----------------------
1. UNIDADES. La columna Unit debe ser homogénea. Un subconjunto en mg/m3
   en vez de ug/m3 introduciría un factor 1000 en unas pocas estaciones,
   difícil de detectar a ojo en un resumen agregado.
2. TIPO DE AGREGACIÓN. La columna AggType distingue observaciones horarias
   de diarias. Si el fichero mezcla ambas para un mismo punto, la media
   estaría ponderando mal.
3. VERIFICACIÓN. La columna Verification indica si el dato pasó el control
   nacional. Pedimos dataset=2 (E1a), luego debería ser mayoritariamente
   verificado; si no, la petición no trajo lo que creíamos.

SALIDA PRINCIPAL
----------------
estaciones_modelado.csv: la tabla que consumirá R/10_prepare.R, con la
media anual, la clasificación de estación (fondo/tráfico/industrial), el
tipo de área y las coordenadas, todo unido y verificado.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import RAW, INTERIM, ANIO, log


def normalizar_id(s: pd.Series) -> pd.Series:
    """
    Homogeneiza el identificador de punto de muestreo.

    En los parquet aparece como 'ES/SPO-ES0118A_8_100' o similar; en las
    tablas de metadatos, a veces sin el prefijo de país. Se recorta el
    prefijo y se pasa a mayúsculas para que el cruce no falle por formato.
    """
    return (
        s.astype(str)
        .str.replace(r"^[A-Z]{2}/", "", regex=True)
        .str.strip()
        .str.upper()
    )


def revisar_parquet() -> None:
    """Comprueba unidades, tipo de agregación y estado de verificación."""
    ruta = INTERIM / "no2_diario.parquet"
    if not ruta.exists():
        log("No encuentro no2_diario.parquet; omito la revisión del crudo.")
        return

    df = pd.read_parquet(ruta)
    log(f"\nRevisión del fichero crudo ({len(df):,} registros)")

    for col in ("Unit", "AggType", "Verification", "Validity", "Pollutant"):
        if col in df.columns:
            conteo = df[col].value_counts(dropna=False).head(8)
            log(f"  {col}:")
            for v, n in conteo.items():
                pct = 100 * n / len(df)
                log(f"      {str(v)[:55]:<55} {n:>10,}  ({pct:5.1f} %)")

    # Alerta específica sobre unidades mezcladas
    if "Unit" in df.columns and df["Unit"].nunique() > 1:
        log(
            "  ATENCIÓN: hay más de una unidad en el fichero. Antes de\n"
            "  promediar hay que convertir todo a ug/m3."
        )


def main() -> None:
    # ------------------------------------------------------------------ #
    # 1. Series calculadas por nosotros
    # ------------------------------------------------------------------ #
    calc = pd.read_csv(INTERIM / "no2_anual_estacion.csv")
    calc["join_id"] = normalizar_id(calc["SamplingPointId"])
    log(f"Medias calculadas por 03: {len(calc)} puntos de muestreo")

    revisar_parquet()

    # ------------------------------------------------------------------ #
    # 2. Estadísticos oficiales del EEA
    # ------------------------------------------------------------------ #
    meta = pd.read_csv(RAW / "estaciones_es_meta.csv", low_memory=False)
    log(f"\nTabla AirQualityStatistics: {len(meta):,} filas")

    def hallar(*claves, obligatorio=True):
        for k in claves:
            for c in meta.columns:
                if k.lower() == c.lower():
                    return c
        for k in claves:
            for c in meta.columns:
                if k.lower() in c.lower():
                    return c
        if obligatorio:
            raise SystemExit(f"No hallo {claves}. Columnas: {list(meta.columns)}")
        return None

    c_sp = hallar("SamplingPointId", "Samplingpoint")
    c_tipo = hallar("AirQualityStationType")
    c_area = hallar("AirQualityStationArea")
    c_eoi = hallar("AirQualityStationEoICode", obligatorio=False)
    c_lon = hallar("Longitude", obligatorio=False)
    c_lat = hallar("Latitude", obligatorio=False)
    c_alt = hallar("Altitude", obligatorio=False)
    # Nombres reales del esquema AirQualityStatistics del EEA:
    #   AirPollutionLevel      -> el valor del estadístico
    #   DataAggregationProcess -> qué estadístico es (media anual, P50, etc.)
    #   YearOfStatistics       -> año de referencia
    c_stat = hallar("DataAggregationProcess", obligatorio=False)
    c_val = hallar("AirPollutionLevel", obligatorio=False)
    c_unidad = hallar("UnitOfAirpollutionLevel", obligatorio=False)
    c_anio = hallar("YearOfStatistics", "ReportingYear", obligatorio=False)
    c_captura = hallar("DataCapture", obligatorio=False)
    c_outlier = hallar("potentialOutlier", obligatorio=False)
    c_pob = hallar("CityPopulation", obligatorio=False)

    if c_unidad:
        log(f"  Unidades del estadístico: {meta[c_unidad].value_counts().to_dict()}")
    if c_outlier:
        n_out = meta[c_outlier].astype(str).str.upper().isin(["Y", "TRUE", "1"]).sum()
        log(f"  Filas marcadas como posible atípico por el EEA: {n_out:,}")

    log(f"  identificador: {c_sp} | tipo: {c_tipo} | área: {c_area}")
    log(f"  estadístico: {c_stat} | valor: {c_val}")

    if c_stat:
        log(f"\n  Estadísticos disponibles en la tabla:")
        for v, n in meta[c_stat].value_counts().head(10).items():
            log(f"      {str(v)[:50]:<50} {n:>8,}")

    meta["join_id"] = normalizar_id(meta[c_sp])

    # ------------------------------------------------------------------ #
    # 3. Extraer la media anual oficial
    # ------------------------------------------------------------------ #
    oficial = None
    if c_stat and c_val:
        m = meta.copy()
        if c_anio:
            m = m[pd.to_numeric(m[c_anio], errors="coerce") == ANIO]
        # La media anual aparece etiquetada como 'Annual mean' o 'P50'/'Mean'
        # El vocabulario EIONET codifica la media anual como 'P1Y'; algunas
        # versiones la escriben como 'Annual mean'. Se aceptan ambas.
        mascara = m[c_stat].astype(str).str.contains(
            r"P1Y|annual\s*mean|^mean$|yearly", case=False, na=False, regex=True
        )
        m = m[mascara]
        if len(m):
            oficial = (
                m.groupby("join_id")[c_val]
                .mean()
                .rename("no2_oficial")
                .reset_index()
            )
            log(f"\nMedias anuales oficiales para {ANIO}: {len(oficial)} puntos")

    # ------------------------------------------------------------------ #
    # 4. Comparación
    # ------------------------------------------------------------------ #
    if oficial is not None and len(oficial):
        comp = calc.merge(oficial, on="join_id", how="inner")
        log(f"Puntos cruzados: {len(comp)}")

        if len(comp):
            comp["dif"] = comp["no2_media"] - comp["no2_oficial"]
            comp["dif_pct"] = 100 * comp["dif"] / comp["no2_oficial"]

            r = np.corrcoef(comp["no2_media"], comp["no2_oficial"])[0, 1]
            print("\n" + "=" * 62)
            print("VERIFICACIÓN: media calculada vs media oficial del EEA")
            print("=" * 62)
            print(f"  n                  = {len(comp)}")
            print(f"  correlación        = {r:.5f}")
            print(f"  diferencia media   = {comp['dif'].mean():+.3f} ug/m3")
            print(f"  diferencia mediana = {comp['dif'].median():+.3f} ug/m3")
            print(f"  |dif| máxima       = {comp['dif'].abs().max():.3f} ug/m3")
            print(f"  RMSE               = {np.sqrt((comp['dif']**2).mean()):.3f}")
            print(f"\n  Nuestro máximo     = {comp['no2_media'].max():.2f} ug/m3")
            print(f"  Máximo oficial     = {comp['no2_oficial'].max():.2f} ug/m3")

            if r > 0.99 and comp["dif"].abs().median() < 0.5:
                print("\n  RESULTADO: las series coinciden. El procesamiento del")
                print("  script 03 queda validado contra la fuente oficial.")
            else:
                print("\n  RESULTADO: hay discrepancia. Revisa los 10 casos con")
                print("  mayor desviación en verificacion_medias.csv.")

            comp.sort_values("dif", key=abs, ascending=False).to_csv(
                INTERIM / "verificacion_medias.csv", index=False
            )
    else:
        log(
            "\nNo pude extraer la media anual oficial de la tabla. "
            "Mira los estadísticos listados arriba y dime cuál corresponde "
            "a la media anual."
        )

    # ------------------------------------------------------------------ #
    # 5. Tabla de modelado: media + clasificación + coordenadas
    # ------------------------------------------------------------------ #
    cols = ["join_id", c_tipo, c_area]
    for c in (c_eoi, c_lon, c_lat, c_alt, c_pob, c_outlier):
        if c:
            cols.append(c)

    clasif = meta[cols].drop_duplicates(subset="join_id")
    modelado = calc.merge(clasif, on="join_id", how="left")

    modelado = modelado.rename(
        columns={
            c_tipo: "tipo_estacion",
            c_area: "area_estacion",
            c_lon: "lon",
            c_lat: "lat",
            c_alt: "altitud_declarada",
            c_pob: "poblacion_ciudad",
            c_outlier: "atipico_eea",
        }
    )
    modelado["es_fondo"] = (
        modelado["tipo_estacion"].astype(str).str.lower().str.contains("background")
    )

    sin_clasif = modelado["tipo_estacion"].isna().sum()
    if sin_clasif:
        log(f"\nAVISO: {sin_clasif} puntos sin clasificación tras el cruce.")

    # Contraste adicional: nuestra captura frente a la oficial del EEA.
    # Divergencias grandes indicarían que el filtro de validez difiere del
    # criterio que aplica el propio organismo.
    if c_captura and c_anio:
        of = meta[pd.to_numeric(meta[c_anio], errors="coerce") == ANIO]
        of = of.groupby("join_id")[c_captura].mean().rename("captura_oficial")
        cap = modelado.merge(of, on="join_id", how="inner")
        if len(cap):
            cap["dif_cap"] = 100 * cap["captura"] - cap["captura_oficial"]
            log(
                f"\nCaptura calculada vs oficial (n={len(cap)}): "
                f"diferencia mediana {cap['dif_cap'].median():+.2f} puntos "
                f"porcentuales, máxima {cap['dif_cap'].abs().max():.1f}"
            )

    modelado.to_csv(INTERIM / "estaciones_modelado.csv", index=False)
    log(f"\n-> data/interim/estaciones_modelado.csv ({len(modelado)} filas)")

    print("\n" + "=" * 62)
    print("NO2 MEDIO ANUAL POR TIPO DE ESTACIÓN")
    print("=" * 62)
    print(
        modelado.groupby("tipo_estacion")["no2_media"]
        .agg(["count", "mean", "median", "max"])
        .round(2)
    )
    print("\nPor tipo de área:")
    print(
        modelado.groupby("area_estacion")["no2_media"]
        .agg(["count", "mean", "median", "max"])
        .round(2)
    )
    print(f"\nEstaciones de fondo (base del variograma): {modelado['es_fondo'].sum()}")
    print(
        "\nPatrón esperado si los datos son correctos:\n"
        "  tráfico > industrial > fondo, y urbano > suburbano > rural.\n"
        "Si no se cumple, hay un problema en el cruce o en la clasificación."
    )


if __name__ == "__main__":
    main()
