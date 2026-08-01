"""
05 — Densidad poblacional por NUTS-3 (Eurostat demo_r_d3dens).

Fuente : https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/demo_r_d3dens
Salida : data/raw/densidad_nuts3_es.csv

La API devuelve JSON-stat 2.0: un array plano `value` indexado por el
producto cartesiano de las dimensiones, con los índices en `dimension`.
No es un formato tabular; hay que reconstruir la tabla. Se hace a mano
abajo para no depender del paquete `pyjstat`.

Se descargan varios años para poder usar la TASA DE CAMBIO de densidad
como covariable: una provincia que se despuebla y otra que se densifica
pueden tener hoy la misma densidad pero trayectorias de emisión opuestas.
"""

from __future__ import annotations

import itertools

import pandas as pd

from common import URL_EUROSTAT, RAW, SESSION, TIMEOUT, log

DATASET = "demo_r_d3dens"
ANIOS = ["2015", "2019", "2021", "2022", "2023"]


def descargar_jsonstat(anios: list[str]) -> dict:
    params = [("format", "JSON"), ("lang", "EN")]
    params += [("time", a) for a in anios]
    r = SESSION.get(f"{URL_EUROSTAT}/{DATASET}", params=params, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def jsonstat_a_dataframe(js: dict) -> pd.DataFrame:
    """Reconstruye la tabla larga a partir del array `value` disperso."""
    dims = js["id"]                      # p.ej. ['freq','unit','geo','time']
    tamanos = js["size"]
    etiquetas = {}
    for d in dims:
        idx = js["dimension"][d]["category"]["index"]
        if isinstance(idx, dict):
            # dict {codigo: posicion} -> lista ordenada por posición
            etiquetas[d] = [k for k, _ in sorted(idx.items(), key=lambda kv: kv[1])]
        else:
            etiquetas[d] = list(idx)

    valores = js["value"]                # dict {indice_plano: valor}
    combinaciones = list(itertools.product(*[etiquetas[d] for d in dims]))

    filas = []
    for pos, combo in enumerate(combinaciones):
        v = valores.get(str(pos), valores.get(pos))
        if v is None:
            continue
        fila = dict(zip(dims, combo))
        fila["value"] = v
        filas.append(fila)

    df = pd.DataFrame(filas)
    assert len(combinaciones) == pd.Series(tamanos).prod(), "Desajuste de dimensiones"
    return df


def main() -> None:
    js = descargar_jsonstat(ANIOS)
    df = jsonstat_a_dataframe(js)
    log(f"Registros recibidos: {len(df):,}  columnas: {list(df.columns)}")

    # NUTS-3 español: código de 5 caracteres que empieza por ES.
    # (ES30 es NUTS-2 Comunidad de Madrid; ES300 es la provincia.)
    es = df[df["geo"].str.startswith("ES") & (df["geo"].str.len() == 5)].copy()
    es = es.rename(columns={"geo": "NUTS_ID", "time": "anio", "value": "densidad"})
    es["anio"] = es["anio"].astype(int)

    ancho = es.pivot_table(
        index="NUTS_ID", columns="anio", values="densidad"
    ).add_prefix("dens_")
    ancho = ancho.reset_index()

    if "dens_2015" in ancho and "dens_2023" in ancho:
        ancho["dens_cambio_pct"] = (
            100 * (ancho["dens_2023"] / ancho["dens_2015"] - 1)
        )

    ancho.to_csv(RAW / "densidad_nuts3_es.csv", index=False)
    log(f"-> data/raw/densidad_nuts3_es.csv ({len(ancho)} regiones NUTS-3)")
    print(ancho.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
