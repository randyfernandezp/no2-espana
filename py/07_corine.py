"""
07 — Cobertura del suelo CORINE 2018 -> covariables de uso del suelo.

Fuente : Copernicus Land Monitoring Service, CLC2018 ráster 100 m (V2020_20u1)
Entrada: cualquier fichero U2018_CLC2018*.tif dentro de data/raw/ (se busca solo)
Salida : data/interim/clc_es_100m.tif           ráster recortado a España
         data/interim/focal/frac_<clase>_<r>m.tif  fracciones por radio
         data/interim/clc_estaciones.csv        muestreo en las estaciones
         data/interim/clc_fracciones_nuts3.csv  agregado por región

CAMBIO DE ESTRATEGIA RESPECTO A LA VERSIÓN ANTERIOR
---------------------------------------------------
La versión inicial recortaba un buffer alrededor de cada estación y contaba
píxeles: 510 estaciones x 5 radios = 2 550 operaciones de máscara sobre un
ráster de 201 MB. Lento, y además solo servía para las estaciones.

Ahora se calcula UNA VEZ, para cada radio, un ráster continuo con la
fracción de cada clase en el entorno de cada píxel (media focal por
convolución). Ese ráster se muestrea después donde haga falta.

Dos ventajas decisivas:
  1. Es de uno a dos órdenes de magnitud más rápido: la convolución por
     solapamiento y suma es O(n log n) sea cual sea el radio.
  2. Resuelve las covariables de la REJILLA DE PREDICCIÓN, no solo de las
     estaciones. El kriging con deriva externa exige conocer la covariable
     en todos los puntos donde se predice, no solo donde se observa; sin
     estos rásteres el KDE es inviable.

POR QUÉ VARIOS RADIOS ANIDADOS
------------------------------
El uso del suelo entra en el modelo a dos escalas distintas:
  - En el kriging representa el entorno de emisión de la estación. La
    literatura de Land Use Regression para NO2 sitúa el radio informativo
    entre 300 m y 5 km, variable según el tipo de estación.
  - En el modelo de intensidad del proceso puntual representa
    disponibilidad de suelo industrial, a escala de decenas de km.
Calcular varios radios y dejar que la validación cruzada elija es más
honesto que fijar 1 km por convención. El radio ganador ES un resultado:
estima la escala espacial a la que opera el proceso.

AGRUPACIÓN DE CLASES
--------------------
Los 44 códigos de nivel 3 se colapsan en seis grupos. Usar los 44 como
variables ficticias daría un modelo con más parámetros que estaciones.
"""

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.windows import from_bounds
from shapely.geometry import mapping

from common import RAW, INTERIM, CRS_TRABAJO, log

RADIOS_M = [500, 1000, 2000, 5000, 10000]
RES = 100  # resolución nativa de CLC en metros

CLC_ES = INTERIM / "clc_es_100m.tif"
DIR_FOCAL = INTERIM / "focal"

# Codificación por índice de clase (1-44), que es la del ráster oficial.
# La tabla .tif.vat.dbf que acompaña al fichero confirma esta correspondencia.
GRUPOS_INDICE = {
    "urbano_denso": [1],                        # tejido urbano continuo
    "urbano_disperso": [2],                     # tejido urbano discontinuo
    "industrial": [3],                          # zonas industriales/comerciales
    "transporte": [4, 5, 6],                    # viario, ferrocarril, puertos, aeropuertos
    "otro_artificial": [7, 8, 9, 10, 11],       # extracción, vertederos, obras, verde urbano
    "agricola": list(range(12, 23)),
    "natural": list(range(23, 45)),
}
# Codificación alternativa por código de nivel 3 (111-523)
GRUPOS_CODIGO = {
    "urbano_denso": [111],
    "urbano_disperso": [112],
    "industrial": [121],
    "transporte": [122, 123, 124],
    "otro_artificial": [131, 132, 133, 141, 142],
    "agricola": list(range(211, 245)),
    "natural": list(range(311, 336)) + list(range(411, 424))
    + [511, 512, 521, 522, 523],
}
ARTIFICIALES = ["urbano_denso", "urbano_disperso", "industrial", "transporte",
                "otro_artificial"]


# --------------------------------------------------------------------------- #
# Localización del ráster y recorte a España
# --------------------------------------------------------------------------- #
def localizar_clc() -> Path:
    patrones = ["**/U2018_CLC2018*.tif", "**/*CLC2018*.tif", "**/clc*100m.tif"]
    for pat in patrones:
        for p in sorted(RAW.glob(pat)):
            # Descartar los auxiliares (.tif.ovr, .tif.aux.xml, etc.)
            if p.name.count(".") == 1 and p.stat().st_size > 10_000_000:
                log(f"Ráster CLC localizado: {p.name} "
                    f"({p.stat().st_size / 1e6:.0f} MB)")
                return p
    raise SystemExit(
        f"No encuentro el ráster CLC en {RAW}.\n"
        "Copia el fichero U2018_CLC2018_V2020_20u1.tif (el de ~200 MB) "
        "desde la carpeta DATA del ZIP a data/raw/.\n"
        "Los ficheros auxiliares (.tfw, .ovr, .vat.dbf) no hacen falta."
    )


def recortar_espana(origen: Path) -> Path:
    """Recorta a la extensión peninsular + Baleares, con margen para los buffers."""
    if CLC_ES.exists():
        log(f"{CLC_ES.name} ya existe, omito el recorte.")
        return CLC_ES

    nuts = gpd.read_file(RAW / "nuts3_es.geojson").to_crs(CRS_TRABAJO)
    fuera = nuts["NUTS_ID"].str.startswith(("ES70", "ES63", "ES64"))
    nuts = nuts[~fuera]

    margen = max(RADIOS_M) + 5000
    minx, miny, maxx, maxy = nuts.total_bounds
    caja = (minx - margen, miny - margen, maxx + margen, maxy + margen)
    log(f"Recortando a {[round(v) for v in caja]} (margen {margen / 1000:.0f} km)")

    with rasterio.open(origen) as src:
        log(f"  CRS del origen: {src.crs} | tamaño {src.width}x{src.height}")
        if src.crs.to_epsg() != CRS_TRABAJO:
            log(f"  AVISO: el ráster no está en EPSG:{CRS_TRABAJO}.")
        ventana = from_bounds(*caja, transform=src.transform)
        datos = src.read(1, window=ventana)
        perfil = src.profile.copy()
        perfil.update(
            height=datos.shape[0],
            width=datos.shape[1],
            transform=src.window_transform(ventana),
            compress="lzw",
        )

    CLC_ES.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(CLC_ES, "w", **perfil) as dst:
        dst.write(datos, 1)
    log(f"-> {CLC_ES.name}  {datos.shape[1]}x{datos.shape[0]} px "
        f"({CLC_ES.stat().st_size / 1e6:.0f} MB)")
    return CLC_ES


def detectar_codificacion(datos: np.ndarray) -> dict[str, list[int]]:
    valores = np.unique(datos)
    log(f"Valores presentes en el ráster: {valores[:20].tolist()}"
        f"{' …' if len(valores) > 20 else ''}")
    if valores.max() >= 111:
        log("Codificación por código de nivel 3 (111-523).")
        return GRUPOS_CODIGO
    log("Codificación por índice de clase (1-44), la del producto oficial.")
    return GRUPOS_INDICE


# --------------------------------------------------------------------------- #
# Fracciones focales por convolución
# --------------------------------------------------------------------------- #
def kernel_disco(radio_px: int) -> np.ndarray:
    """Kernel circular normalizado. Circular, no cuadrado: la distancia a la
    que un uso del suelo influye no depende de la orientación."""
    y, x = np.ogrid[-radio_px : radio_px + 1, -radio_px : radio_px + 1]
    disco = (x**2 + y**2) <= radio_px**2
    return disco.astype(np.float32) / disco.sum()


def fraccion_focal(mascara: np.ndarray, radio_m: int) -> np.ndarray:
    radio_px = max(1, int(round(radio_m / RES)))
    k = kernel_disco(radio_px)
    try:
        from scipy.signal import oaconvolve

        return oaconvolve(mascara.astype(np.float32), k, mode="same").astype("float32")
    except ImportError:
        # Respaldo sin scipy: ventana cuadrada mediante imagen integral.
        # Aproximación aceptable; la diferencia frente al disco es pequeña
        # a estas escalas, pero conviene dejar constancia de que se usó.
        log("    (scipy no disponible: uso ventana cuadrada como aproximación)")
        pad = radio_px
        m = np.pad(mascara.astype(np.float32), pad, mode="edge")
        integral = m.cumsum(0).cumsum(1)
        lado = 2 * radio_px + 1
        h, w = mascara.shape
        a = integral[0:h, 0:w]
        b = integral[0:h, lado : lado + w]
        c = integral[lado : lado + h, 0:w]
        d = integral[lado : lado + h, lado : lado + w]
        return ((d - b - c + a) / (lado * lado)).astype("float32")


def generar_focales(ruta: Path) -> dict[tuple[str, int], Path]:
    DIR_FOCAL.mkdir(parents=True, exist_ok=True)
    generados: dict[tuple[str, int], Path] = {}

    with rasterio.open(ruta) as src:
        datos = src.read(1)
        perfil = src.profile.copy()
        perfil.update(dtype="float32", count=1, compress="lzw", nodata=np.nan)

    grupos = detectar_codificacion(datos)

    # Máscara de superficie artificial total, además de cada grupo por separado
    mascaras = {n: np.isin(datos, c) for n, c in grupos.items()}
    mascaras["artificial"] = np.logical_or.reduce(
        [mascaras[n] for n in ARTIFICIALES]
    )

    log("\nComposición del territorio recortado:")
    total = datos.size
    for nombre, m in sorted(mascaras.items(), key=lambda kv: -kv[1].sum()):
        log(f"    {nombre:<18} {100 * m.sum() / total:5.2f} %")

    clases_utiles = ["artificial", "industrial", "urbano_denso", "urbano_disperso",
                     "transporte", "agricola", "natural"]

    for clase in clases_utiles:
        if clase not in mascaras:
            continue
        for r in RADIOS_M:
            destino = DIR_FOCAL / f"frac_{clase}_{r}m.tif"
            generados[(clase, r)] = destino
            if destino.exists():
                continue
            frac = fraccion_focal(mascaras[clase], r)
            with rasterio.open(destino, "w", **perfil) as dst:
                dst.write(frac, 1)
            log(f"  -> {destino.name}")

    return generados


# --------------------------------------------------------------------------- #
# Muestreo en puntos y agregación por región
# --------------------------------------------------------------------------- #
def muestrear_estaciones(focales: dict[tuple[str, int], Path]) -> None:
    est = pd.read_csv(INTERIM / "estaciones_modelado.csv")
    if not {"lon", "lat"} <= set(est.columns):
        raise SystemExit("estaciones_modelado.csv no tiene lon/lat. Ejecuta 08 antes.")

    g = gpd.GeoDataFrame(
        est,
        geometry=gpd.points_from_xy(est["lon"], est["lat"]),
        crs=4326,
    ).to_crs(CRS_TRABAJO)
    coords = [(p.x, p.y) for p in g.geometry]
    log(f"\nMuestreando {len(coords)} estaciones en los rásteres focales …")

    salida = est[["join_id"]].copy()
    salida["x_3035"] = [c[0] for c in coords]
    salida["y_3035"] = [c[1] for c in coords]

    for (clase, r), ruta in sorted(focales.items()):
        with rasterio.open(ruta) as src:
            vals = [v[0] for v in src.sample(coords)]
        salida[f"{clase}_{r}m"] = vals

    salida.to_csv(INTERIM / "clc_estaciones.csv", index=False)
    log(f"-> clc_estaciones.csv ({len(salida)} filas, {len(salida.columns)} columnas)")

    log("\n  Fracción artificial media por radio (control de coherencia):")
    for r in RADIOS_M:
        col = f"artificial_{r}m"
        if col in salida:
            log(f"    {r:>6} m: {salida[col].mean():.4f}")
    log("  Debe decrecer al ampliar el radio: las estaciones están en zonas")
    log("  más artificiales que su entorno amplio. Si crece, hay un error.")


def agregar_nuts3(ruta_clc: Path) -> None:
    nuts = gpd.read_file(RAW / "nuts3_es.geojson").to_crs(CRS_TRABAJO)
    log(f"\nAgregando por {len(nuts)} regiones NUTS-3 …")

    with rasterio.open(ruta_clc) as src:
        muestra = src.read(1, window=((0, 500), (0, 500)))
        grupos = detectar_codificacion(muestra)

        filas = []
        for i, (nid, geom) in enumerate(zip(nuts["NUTS_ID"], nuts.geometry), 1):
            try:
                recorte, _ = mask(src, [mapping(geom)], crop=True, filled=False)
            except ValueError:
                continue
            pix = recorte.compressed()
            if pix.size == 0:
                continue
            fila = {"NUTS_ID": nid, "n_pixeles": int(pix.size)}
            for nombre, codigos in grupos.items():
                fila[f"frac_{nombre}"] = float(np.isin(pix, codigos).mean())
            fila["frac_artificial"] = sum(fila[f"frac_{n}"] for n in ARTIFICIALES)
            # Superficie absoluta en km2: una provincia grande con 2 % de suelo
            # industrial tiene más suelo industrial que una pequeña con 5 %,
            # y para el conteo de fuentes lo relevante es lo absoluto.
            fila["area_industrial_km2"] = fila["frac_industrial"] * pix.size * 0.01
            fila["area_artificial_km2"] = fila["frac_artificial"] * pix.size * 0.01
            filas.append(fila)
            if i % 20 == 0:
                log(f"    {i}/{len(nuts)}")

    df = pd.DataFrame(filas)
    df.to_csv(INTERIM / "clc_fracciones_nuts3.csv", index=False)
    log(f"-> clc_fracciones_nuts3.csv ({len(df)} regiones)")
    log("\n  Provincias con mayor superficie industrial (km2):")
    for _, r in df.nlargest(8, "area_industrial_km2").iterrows():
        log(f"    {r['NUTS_ID']}: {r['area_industrial_km2']:8.1f} km2  "
            f"({100 * r['frac_industrial']:.2f} % del territorio)")


def main() -> None:
    origen = localizar_clc()
    recorte = recortar_espana(origen)
    focales = generar_focales(recorte)
    (INTERIM / "clc_focales.json").write_text(
        json.dumps({f"{c}_{r}": str(p) for (c, r), p in focales.items()}, indent=2),
        encoding="utf-8",
    )
    muestrear_estaciones(focales)
    agregar_nuts3(recorte)
    log(
        "\nLos rásteres de data/interim/focal/ sirven tanto para las estaciones "
        "como para la rejilla de predicción del kriging con deriva externa."
    )


if __name__ == "__main__":
    main()
