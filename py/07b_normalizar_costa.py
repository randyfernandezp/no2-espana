"""
07b — Normalización costera de las fracciones focales.

Entrada: data/interim/focal/frac_*.tif   (generados por 07_corine.py)
Salida : data/interim/focal_norm/frac_*.tif
         data/interim/focal_norm/frac_tierra_<r>m.tif
         data/interim/clc_estaciones.csv  (reescrito con ambas versiones)

EL PROBLEMA QUE CORRIGE
-----------------------
El ráster CORINE codifica el mar y el exterior de la cobertura como -128.
La media focal calculada en 07_corine.py divide entre TODOS los píxeles del
disco, mar incluido. Para una estación costera eso significa que su fracción
artificial se calcula sobre un disco que es medio agua:

    Estación en Barcelona, radio 5 km
        superficie artificial : 40 km2
        disco completo        : 78,5 km2  (de los cuales ~35 km2 son mar)
        fracción cruda        : 0,51
        fracción sobre tierra : 0,92

La versión cruda no mide "cuán urbanizado está el entorno" sino "cuánta
ciudad hay en el disco", que confunde dos cosas distintas. En España el
efecto no es menor: Barcelona, Valencia, Bilbao, Málaga, Gijón, A Coruña,
Cartagena y Tarragona son costeras y están entre las de mayor NO2. El sesgo
va justo en la dirección que más daño hace, deprimiendo la covariable
precisamente donde la concentración es alta, lo que atenúa el coeficiente
estimado en el kriging con deriva externa.

LA CORRECCIÓN
-------------
Basta una convolución adicional por radio, la de la máscara de tierra:

    fraccion_normalizada = fraccion_cruda / fraccion_tierra

Ambas comparten denominador (el disco completo), así que el cociente da
directamente la proporción sobre superficie terrestre. No hace falta
recalcular las 35 convoluciones ya hechas.

QUÉ VERSIÓN USAR
----------------
Se conservan las dos, y la decisión es sustantiva, no técnica:

  - NORMALIZADA para covariables de composición del entorno (¿qué proporción
    del suelo circundante es urbano o industrial?). Es la que entra en el
    kriging y en el modelo de intensidad.

  - CRUDA cuando lo relevante sea la cantidad absoluta de fuente en el
    entorno, no su proporción. Media ciudad y medio mar emite la mitad que
    una ciudad entera.

Además, frac_tierra es en sí una covariable útil: es una medida continua de
"proximidad al mar" que captura la ventilación por brisa marina, un factor
real de dispersión del NO2 en las ciudades costeras españolas.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio

from common import INTERIM, CRS_TRABAJO, log

DIR_FOCAL = INTERIM / "focal"
DIR_NORM = INTERIM / "focal_norm"
CLC_ES = INTERIM / "clc_es_100m.tif"
RES = 100
RADIOS_M = [500, 1000, 2000, 5000, 10000]

# Umbral por debajo del cual el cociente es inestable: si menos del 5 % del
# disco es tierra, la fracción normalizada se dispara por división entre casi
# cero. Esos píxeles (mar abierto) se marcan como sin dato.
MIN_TIERRA = 0.05


def kernel_disco(radio_px: int) -> np.ndarray:
    y, x = np.ogrid[-radio_px : radio_px + 1, -radio_px : radio_px + 1]
    disco = (x**2 + y**2) <= radio_px**2
    return disco.astype(np.float32) / disco.sum()


def convolucionar(mascara: np.ndarray, radio_m: int) -> np.ndarray:
    from scipy.signal import oaconvolve

    k = kernel_disco(max(1, int(round(radio_m / RES))))
    return oaconvolve(mascara.astype(np.float32), k, mode="same").astype("float32")


def main() -> None:
    if not DIR_FOCAL.exists():
        raise SystemExit(f"No existe {DIR_FOCAL}. Ejecuta antes 07_corine.py")
    DIR_NORM.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Máscara de tierra
    # ------------------------------------------------------------------ #
    with rasterio.open(CLC_ES) as src:
        datos = src.read(1)
        perfil = src.profile.copy()
        perfil.update(dtype="float32", count=1, compress="lzw", nodata=np.nan)

    # CORINE NO codifica el mar como nodata: lo clasifica como una clase más.
    # En la nomenclatura de nivel 3 son
    #     521 lagunas costeras, 522 estuarios, 523 mar y océano
    # que en la codificación por índice corresponden a 42, 43 y 44.
    # El valor -128 significa "fuera de la cobertura CORINE", no mar.
    #
    # Las aguas continentales (511 cursos de agua, 512 masas de agua;
    # índices 40 y 41) SÍ se cuentan como superficie terrestre: forman parte
    # del paisaje continental y de su hidrología, aunque no sean fuente de
    # emisión. Excluirlas convertiría la covariable en algo distinto, y
    # penalizaría injustamente a estaciones junto a embalses o rías.
    if datos.max() >= 111:
        MARINAS = [521, 522, 523]
        valido = (datos >= 111) & (datos <= 523)
    else:
        MARINAS = [42, 43, 44]
        valido = (datos >= 1) & (datos <= 44)

    tierra = valido & ~np.isin(datos, MARINAS)

    n_mar = int(np.isin(datos, MARINAS).sum())
    n_fuera = int((~valido).sum())
    log(f"Composición del recorte ({datos.size:,} píxeles):")
    log(f"    tierra firme        {100 * tierra.mean():5.1f} %")
    log(f"    mar y aguas costeras {100 * n_mar / datos.size:5.1f} %")
    log(f"    fuera de cobertura   {100 * n_fuera / datos.size:5.1f} %")
    log(f"  Valores fuera de cobertura: "
        f"{np.unique(datos[~valido])[:8].tolist()}")

    # ------------------------------------------------------------------ #
    # 2. Fracción de tierra por radio
    # ------------------------------------------------------------------ #
    fracs_tierra = {}
    for r in RADIOS_M:
        destino = DIR_NORM / f"frac_tierra_{r}m.tif"
        if destino.exists():
            with rasterio.open(destino) as src:
                fracs_tierra[r] = src.read(1)
            log(f"  {destino.name} ya existe")
            continue
        log(f"  calculando fracción de tierra a {r} m …")
        ft = convolucionar(tierra, r)
        fracs_tierra[r] = ft
        with rasterio.open(destino, "w", **perfil) as dst:
            dst.write(ft, 1)
        log(f"  -> {destino.name}")

    # ------------------------------------------------------------------ #
    # 3. Normalizar cada fracción existente
    # ------------------------------------------------------------------ #
    log("\nNormalizando las fracciones de clase …")
    for ruta in sorted(DIR_FOCAL.glob("frac_*.tif")):
        partes = ruta.stem.split("_")
        try:
            radio = int(partes[-1].replace("m", ""))
        except ValueError:
            continue
        if radio not in fracs_tierra:
            continue

        destino = DIR_NORM / ruta.name
        if destino.exists():
            continue

        with rasterio.open(ruta) as src:
            cruda = src.read(1)

        ft = fracs_tierra[radio]
        with np.errstate(divide="ignore", invalid="ignore"):
            norm = np.where(ft >= MIN_TIERRA, cruda / ft, np.nan).astype("float32")
        # El cociente puede exceder 1 por ruido numérico en los bordes
        norm = np.clip(norm, 0.0, 1.0)

        with rasterio.open(destino, "w", **perfil) as dst:
            dst.write(norm, 1)
        log(f"  -> {destino.name}")

    # ------------------------------------------------------------------ #
    # 4. Remuestrear en las estaciones, con ambas versiones
    # ------------------------------------------------------------------ #
    est = pd.read_csv(INTERIM / "estaciones_modelado.csv")
    g = gpd.GeoDataFrame(
        est,
        geometry=gpd.points_from_xy(est["lon"], est["lat"]),
        crs=4326,
    ).to_crs(CRS_TRABAJO)
    coords = [(p.x, p.y) for p in g.geometry]

    salida = est[["join_id"]].copy()
    salida["x_3035"] = [c[0] for c in coords]
    salida["y_3035"] = [c[1] for c in coords]

    for ruta in sorted(DIR_FOCAL.glob("frac_*.tif")):
        nombre = ruta.stem.replace("frac_", "")
        with rasterio.open(ruta) as src:
            salida[f"{nombre}_crudo"] = [v[0] for v in src.sample(coords)]
        norm = DIR_NORM / ruta.name
        if norm.exists():
            with rasterio.open(norm) as src:
                salida[nombre] = [v[0] for v in src.sample(coords)]

    for r in RADIOS_M:
        ruta = DIR_NORM / f"frac_tierra_{r}m.tif"
        if ruta.exists():
            with rasterio.open(ruta) as src:
                salida[f"tierra_{r}m"] = [v[0] for v in src.sample(coords)]

    salida.to_csv(INTERIM / "clc_estaciones.csv", index=False)
    log(f"\n-> clc_estaciones.csv ({len(salida)} filas, {len(salida.columns)} col.)")

    # ------------------------------------------------------------------ #
    # 5. Diagnóstico del efecto de la corrección
    # ------------------------------------------------------------------ #
    print("\n" + "=" * 66)
    print("EFECTO DE LA NORMALIZACIÓN COSTERA")
    print("=" * 66)
    for r in RADIOS_M:
        c, n, t = f"artificial_{r}m_crudo", f"artificial_{r}m", f"tierra_{r}m"
        if c in salida and n in salida:
            costeras = salida[salida[t] < 0.9]
            print(
                f"  {r:>6} m | cruda {salida[c].mean():.4f} -> "
                f"normalizada {salida[n].mean():.4f} | "
                f"estaciones con mar en el disco: {len(costeras)}"
            )

    if "tierra_5000m" in salida:
        cost = salida.nsmallest(8, "tierra_5000m")[
            ["join_id", "tierra_5000m", "artificial_5000m_crudo", "artificial_5000m"]
        ]
        print("\n  Las ocho estaciones más costeras (radio 5 km):")
        print(cost.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(
            "\n  Comprueba que la columna normalizada sube claramente respecto\n"
            "  a la cruda: ahí estaba el sesgo que se acaba de corregir."
        )


if __name__ == "__main__":
    main()
