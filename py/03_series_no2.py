"""
03 — Series diarias de NO2 en España, 2023 (datos verificados E1a).

Fuente : EEA Air Quality Download Service, API Parquet
         https://eeadmz1-downloads-api-appservice.azurewebsites.net
Salida : data/raw/eea_no2_es_2023.zip
         data/interim/no2_diario.parquet
         data/interim/no2_anual_estacion.csv

POR QUÉ ESTA VERSIÓN DESCUBRE LOS PARÁMETROS
--------------------------------------------
El primer intento devolvió un ZIP de 0 MB sin ningún parquet dentro. El
servidor aceptó la petición y generó el archivo, luego la sintaxis del
cuerpo era válida; simplemente ningún registro casaba con el filtro.

La causa habitual es el identificador del contaminante: la API no espera
la cadena "NO2" sino el URI del vocabulario EIONET, del tipo
    http://dd.eionet.europa.eu/vocabulary/aq/pollutant/8
Como no conviene codificarlo a ciegas, el script pregunta primero a los
endpoints de descubrimiento (/Country, /Pollutant, /Property) y construye
el cuerpo con los valores que el propio servicio declara válidos.

Si aun así el ZIP viene vacío, prueba variantes del cuerpo de forma
sistemática y te dice cuál funcionó, para que quede documentado.

CRITERIO DE CAPTURA
-------------------
La Directiva 2008/50/CE exige cobertura temporal mínima del 75 % para que
una media anual sea legalmente válida. Una estación con cuatro meses de
datos tendría una media sesgada por la estacionalidad: el NO2 invernal en
España es del orden de 1,8 veces el estival, por inversión térmica y menor
fotólisis.
"""

from __future__ import annotations

import io
import json
import time
import zipfile
from datetime import datetime

import pandas as pd
import requests

from common import URL_AQ_API, RAW, INTERIM, PAIS, ANIO, CONTAMINANTE, SESSION, log

EMAIL = "randy.fernandez.p@uni.pe"

CAPTURA_MINIMA = 0.75
DIAS_MINIMOS = int(365 * CAPTURA_MINIMA)


# --------------------------------------------------------------------------- #
# Descubrimiento de parámetros válidos
# --------------------------------------------------------------------------- #
def explorar_endpoint(ruta: str) -> list | dict | None:
    try:
        r = SESSION.get(f"{URL_AQ_API}/{ruta}", timeout=120)
        r.raise_for_status()
        return r.json()
    except Exception as exc:  # noqa: BLE001
        log(f"  /{ruta}: {str(exc)[:70]}")
        return None


def identificador_no2() -> list[str]:
    """Devuelve los identificadores plausibles de NO2, en orden de preferencia."""
    log("Descubriendo el identificador de NO2 …")
    candidatos: list[str] = []

    for ruta in ("Pollutant", "Property", "Vocabulary/Pollutant"):
        datos = explorar_endpoint(ruta)
        if not isinstance(datos, list) or not datos:
            continue
        log(f"  /{ruta} devolvió {len(datos)} entradas")
        log(f"    ejemplo: {json.dumps(datos[0], ensure_ascii=False)[:200]}")

        for item in datos:
            if not isinstance(item, dict):
                continue
            texto = json.dumps(item, ensure_ascii=False).lower()
            # Buscar NO2 evitando NOX y otros óxidos
            if '"no2"' in texto or "nitrogen dioxide" in texto:
                for clave in ("id", "uri", "notation", "url", "value", "code"):
                    v = item.get(clave)
                    if isinstance(v, str) and v and v not in candidatos:
                        candidatos.append(v)
                log(f"    NO2 encontrado: {json.dumps(item, ensure_ascii=False)[:220]}")
        if candidatos:
            break

    # Respaldos conocidos del vocabulario EIONET
    for extra in (
        "http://dd.eionet.europa.eu/vocabulary/aq/pollutant/8",
        "8",
        "NO2",
    ):
        if extra not in candidatos:
            candidatos.append(extra)

    log(f"  Candidatos a probar: {candidatos}")
    return candidatos


# --------------------------------------------------------------------------- #
# Petición y descarga
# --------------------------------------------------------------------------- #
def cuerpo(poll: str, agregacion: str | None) -> dict:
    b = {
        "countries": [PAIS],
        "cities": [],
        "pollutants": [poll],
        "dataset": 2,          # E1a: datos verificados por los países
        "dateTimeStart": f"{ANIO}-01-01T00:00:00Z",
        "dateTimeEnd": f"{ANIO}-12-31T23:59:59Z",
        "email": EMAIL,
        "source": "API",
    }
    if agregacion:
        b["aggregationType"] = agregacion
    return b


def pedir(body: dict) -> str | None:
    try:
        r = SESSION.post(f"{URL_AQ_API}/ParquetFile/async", json=body, timeout=180)
        r.raise_for_status()
        return r.text.strip().strip('"')
    except Exception as exc:  # noqa: BLE001
        log(f"    POST falló: {str(exc)[:90]}")
        return None


def descargar(url: str, destino, max_minutos: int = 45) -> int:
    """Espera a que el servidor genere el ZIP y lo guarda. Devuelve su tamaño."""
    t0 = datetime.now()
    intento = 0
    while (datetime.now() - t0).total_seconds() < max_minutos * 60:
        intento += 1
        resp = requests.get(url, timeout=900)
        if resp.status_code == 404:
            if intento % 3 == 1:
                log(f"    generándose … (intento {intento})")
            time.sleep(20)
            continue
        resp.raise_for_status()
        destino.write_bytes(resp.content)
        return len(resp.content)
    raise TimeoutError("El EEA no generó el archivo en el plazo previsto.")


def contar_parquets(zip_path) -> int:
    try:
        with zipfile.ZipFile(zip_path) as z:
            return sum(1 for n in z.namelist() if n.lower().endswith(".parquet"))
    except zipfile.BadZipFile:
        return 0


def obtener_zip():
    zip_path = RAW / f"eea_{CONTAMINANTE.lower()}_{PAIS.lower()}_{ANIO}.zip"
    if zip_path.exists() and contar_parquets(zip_path) > 0:
        log(f"{zip_path.name} ya existe y contiene datos.")
        return zip_path

    for poll in identificador_no2():
        for agregacion in ("day", None, "hour"):
            etiqueta = f"pollutant={poll[:60]!r}, agg={agregacion!r}"
            log(f"Probando {etiqueta}")
            url = pedir(cuerpo(poll, agregacion))
            if not url:
                continue
            tam = descargar(url, zip_path)
            n = contar_parquets(zip_path)
            log(f"    ZIP de {tam / 1e6:.1f} MB con {n} parquet(s)")
            if n > 0:
                log(f"COMBINACIÓN VÁLIDA -> {etiqueta}")
                (RAW / "aq_api_parametros.json").write_text(
                    json.dumps(cuerpo(poll, agregacion), indent=2), encoding="utf-8"
                )
                return zip_path

    raise SystemExit(
        "Ninguna combinación devolvió datos.\n"
        "Comprueba manualmente en https://eeadmz1-downloads-webapp.azurewebsites.net/ "
        "qué años tienen datos E1a para España y NO2, y ajusta ANIO en common.py."
    )


# --------------------------------------------------------------------------- #
# Consolidación
# --------------------------------------------------------------------------- #
def consolidar(zip_path) -> pd.DataFrame:
    trozos = []
    with zipfile.ZipFile(zip_path) as z:
        nombres = [n for n in z.namelist() if n.lower().endswith(".parquet")]
        log(f"{len(nombres)} archivos parquet en el ZIP")
        for i, n in enumerate(nombres, 1):
            with z.open(n) as fh:
                trozos.append(pd.read_parquet(io.BytesIO(fh.read())))
            if i % 100 == 0:
                log(f"  leídos {i}/{len(nombres)}")
    df = pd.concat(trozos, ignore_index=True)
    log(f"Registros totales: {len(df):,}")
    log(f"Columnas: {list(df.columns)}")
    return df


def medias_anuales(df: pd.DataFrame) -> pd.DataFrame:
    def col(*nombres, obligatorio=True):
        idx = {c.lower(): c for c in df.columns}
        for n in nombres:
            if n.lower() in idx:
                return idx[n.lower()]
        for c in df.columns:
            if any(n.lower() in c.lower() for n in nombres):
                return c
        if obligatorio:
            raise SystemExit(f"No hallo columna {nombres}. Hay: {list(df.columns)}")
        return None

    c_val = col("Value")
    c_sp = col("Samplingpoint", "SamplingPointId")
    c_ini = col("Start")
    c_flag = col("Validity", obligatorio=False)

    d = df.copy()
    d[c_ini] = pd.to_datetime(d[c_ini], errors="coerce", utc=True)
    d = d[d[c_ini].dt.year == ANIO]
    d[c_val] = pd.to_numeric(d[c_val], errors="coerce")

    # Validity >= 1 marca observación válida en el vocabulario del EEA;
    # los negativos (-1, -99) señalan dato inválido o no verificado.
    if c_flag:
        d = d[pd.to_numeric(d[c_flag], errors="coerce") >= 1]
    d = d[d[c_val] >= 0]

    # Si el servidor devolvió datos horarios, agregamos a diario primero.
    d["fecha"] = d[c_ini].dt.date
    diario = d.groupby([c_sp, "fecha"], as_index=False)[c_val].mean()
    log(f"Observaciones diarias: {len(diario):,}")

    diario["mes"] = pd.to_datetime(diario["fecha"]).dt.month
    g = (
        diario.groupby(c_sp)
        .agg(
            no2_media=(c_val, "mean"),
            no2_mediana=(c_val, "median"),
            no2_p95=(c_val, lambda s: s.quantile(0.95)),
            n_dias=(c_val, "size"),
        )
        .reset_index()
        .rename(columns={c_sp: "SamplingPointId"})
    )

    inv = diario[diario["mes"].isin([12, 1, 2])].groupby(c_sp)[c_val].mean()
    ver = diario[diario["mes"].isin([6, 7, 8])].groupby(c_sp)[c_val].mean()
    g["no2_invierno"] = g["SamplingPointId"].map(inv)
    g["no2_verano"] = g["SamplingPointId"].map(ver)
    g["razon_inv_ver"] = g["no2_invierno"] / g["no2_verano"]
    g["captura"] = g["n_dias"] / 365

    antes = len(g)
    g = g[g["n_dias"] >= DIAS_MINIMOS].copy()
    log(f"Puntos de muestreo: {antes} -> {len(g)} con captura >= {CAPTURA_MINIMA:.0%}")
    return g


def main() -> None:
    zip_path = obtener_zip()
    df = consolidar(zip_path)
    df.to_parquet(INTERIM / "no2_diario.parquet", index=False)

    anual = medias_anuales(df)
    anual.to_csv(INTERIM / "no2_anual_estacion.csv", index=False)
    log(f"-> data/interim/no2_anual_estacion.csv ({len(anual)} estaciones)")
    print(anual["no2_media"].describe())
    print(
        "\nContraste: la media española de NO2 de fondo suele situarse entre\n"
        "8 y 15 ug/m3, y las estaciones de tráfico urbano entre 30 y 50.\n"
        "Si tus cifras se alejan mucho, revisa las unidades del campo Value."
    )


if __name__ == "__main__":
    main()
