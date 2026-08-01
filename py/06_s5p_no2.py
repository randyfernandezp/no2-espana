"""
06 — Columna troposférica de NO2 (Sentinel-5P / TROPOMI) sobre España, 2023.

Fuente : Copernicus Data Space Ecosystem, API openEO
Salida : data/raw/s5p_no2_2023_mensual/*.tif  (12 composiciones mensuales)
         data/raw/s5p_no2_2023_media.tif      (media anual)

REQUIERE CUENTA GRATUITA
------------------------
Regístrate en https://dataspace.copernicus.eu/ antes de ejecutar.
La primera ejecución abre el navegador para autenticación OIDC; el token
queda cacheado en ~/.local/share/openeo-python-client/.

Instalación:  pip install openeo rioxarray

LIMITACIÓN CONOCIDA DEL SERVICIO
--------------------------------
En la implementación actual de openEO y Sentinel Hub, la colección
Sentinel-5P admite UNA SOLA BANDA por petición. Pedir NO2 y CO a la vez
devuelve error. Aquí solo pedimos NO2, así que no afecta.

POR QUÉ COMPOSICIONES MENSUALES Y NO UNA SOLA MEDIA ANUAL
---------------------------------------------------------
Tres razones:
  1. Robustez operativa: una petición anual sobre toda España suele agotar
     el tiempo de procesamiento. Doce peticiones mensuales siempre pasan.
  2. La media anual de la columna pondera implícitamente por número de
     observaciones válidas, que no es uniforme (más nubes en invierno en
     la cornisa cantábrica). Con mensuales puedes promediar los doce meses
     con peso igual y evitar ese sesgo estacional-geográfico.
  3. Te deja hacer kriging estacional después sin volver a descargar.

QUÉ MIDE Y QUÉ NO
-----------------
TROPOMI mide la COLUMNA troposférica integrada verticalmente (mol/m2), no
la concentración en superficie (ug/m3). La relación entre ambas depende de
la altura de la capa de mezcla, que varía con la estación y la orografía.
Por eso la columna entra en el modelo como COVARIABLE en kriging con deriva
externa, no como sustituto de la medición en superficie. El coeficiente
estimado es, de hecho, un resultado interpretable: mide el factor de
conversión columna-superficie efectivo para España.

Además el paso local es de ~3,5 x 5,5 km en nadir: no resuelve gradientes
urbanos. Su valor añadido está en las zonas SIN estaciones, que es
exactamente donde el kriging más lo necesita.
"""

from __future__ import annotations

import calendar
from pathlib import Path

import openeo

from common import RAW, ANIO, URL_OPENEO, log

# Extensión peninsular + Baleares en EPSG:4326
EXTENSION = {"west": -9.8, "south": 35.8, "east": 4.6, "north": 44.0}

COLECCION = "SENTINEL_5P_L2"
BANDA = "NO2"
RESOLUCION_GRADOS = 0.045  # ~5 km, cercano a la resolución nativa de TROPOMI

DESTINO = RAW / f"s5p_no2_{ANIO}_mensual"


def conectar() -> openeo.Connection:
    log(f"Conectando a {URL_OPENEO} …")
    con = openeo.connect(URL_OPENEO)
    con.authenticate_oidc()  # abre el navegador la primera vez
    log("Autenticado.")
    return con


def descargar_mes(con: openeo.Connection, mes: int) -> Path:
    ultimo = calendar.monthrange(ANIO, mes)[1]
    inicio = f"{ANIO}-{mes:02d}-01"
    fin = f"{ANIO}-{mes:02d}-{ultimo:02d}"
    salida = DESTINO / f"s5p_no2_{ANIO}{mes:02d}.tif"

    if salida.exists():
        log(f"  {salida.name} ya existe, omito.")
        return salida

    cubo = con.load_collection(
        COLECCION,
        spatial_extent=EXTENSION,
        temporal_extent=[inicio, fin],
        bands=[BANDA],
    )
    # Media temporal del mes -> imagen 2D
    cubo = cubo.reduce_dimension(dimension="t", reducer="mean")
    cubo = cubo.resample_spatial(resolution=RESOLUCION_GRADOS, projection=4326)

    log(f"  procesando {inicio} .. {fin} …")
    cubo.download(str(salida), format="GTiff")
    log(f"  -> {salida.name}")
    return salida


def media_anual(rutas: list[Path]) -> None:
    """
    Promedio con peso igual por mes.

    Por qué peso igual y no un promedio de todos los píxeles del año:
    el número de observaciones válidas de TROPOMI no es uniforme en el
    tiempo ni en el espacio. En invierno hay más nubosidad, y la cornisa
    cantábrica pierde muchas más pasadas que el sureste peninsular. Un
    promedio directo ponderaría implícitamente por número de pasadas y
    daría un mapa sesgado hacia los meses y regiones despejados, que son
    justo los de menor NO2 (más fotólisis, capa de mezcla más alta).
    Promediando primero dentro de cada mes y luego entre meses con peso
    igual, ese sesgo estacional-geográfico desaparece.
    """
    import numpy as np
    import rioxarray as rxr

    capas, nombres = [], []
    for p in sorted(rutas):
        if p.exists():
            capas.append(rxr.open_rasterio(p, masked=True).squeeze())
            nombres.append(p.stem[-6:])
    if not capas:
        log("No hay capas mensuales para promediar.")
        return

    apilado = np.stack([c.values for c in capas])

    # Diagnóstico de cobertura: qué fracción del dominio tiene dato cada mes.
    # Conviene documentarlo: si algún mes baja del 50 %, la media anual se
    # apoya en menos meses de los que parece en esas zonas.
    log("Cobertura de píxeles válidos por mes:")
    for nombre, capa in zip(nombres, apilado):
        valido = np.isfinite(capa).mean()
        aviso = "  <- cobertura baja" if valido < 0.5 else ""
        log(f"    {nombre}: {valido:6.1%}{aviso}")

    media = np.nanmean(apilado, axis=0)
    n_meses = np.isfinite(apilado).sum(axis=0)

    log(f"Meses válidos por píxel: mediana {np.median(n_meses):.0f}, "
        f"mínimo {n_meses.min():.0f}")
    finitos = media[np.isfinite(media)]
    if finitos.size:
        log(f"Columna troposférica de NO2 (mol/m2): "
            f"min {finitos.min():.3e}, mediana {np.median(finitos):.3e}, "
            f"max {finitos.max():.3e}")

    plantilla = capas[0].copy(data=media)

    # rioxarray falla si '_FillValue' está a la vez en attrs y en encoding.
    # Se limpia de ambos sitios antes de fijar el nodata definitivo.
    for objetivo in (plantilla.attrs, plantilla.encoding):
        for clave in ("_FillValue", "missing_value", "scale_factor", "add_offset"):
            objetivo.pop(clave, None)

    plantilla.rio.write_nodata(np.nan, encoded=False, inplace=True)
    destino = RAW / f"s5p_no2_{ANIO}_media.tif"
    plantilla.rio.to_raster(destino)
    log(f"-> {destino.name}  (media con peso igual de {len(capas)} meses)")

    # Mapa auxiliar del número de meses válidos: sirve para ponderar la
    # confianza de la covariable en el kriging, o para excluir zonas con
    # cobertura insuficiente.
    conteo = capas[0].copy(data=n_meses.astype("float32"))
    for objetivo in (conteo.attrs, conteo.encoding):
        for clave in ("_FillValue", "missing_value", "scale_factor", "add_offset"):
            objetivo.pop(clave, None)
    conteo.rio.write_nodata(np.nan, encoded=False, inplace=True)
    conteo.rio.to_raster(RAW / f"s5p_no2_{ANIO}_nmeses.tif")
    log(f"-> s5p_no2_{ANIO}_nmeses.tif  (meses válidos por píxel)")


def main() -> None:
    DESTINO.mkdir(parents=True, exist_ok=True)
    con = conectar()
    rutas = []
    for mes in range(1, 13):
        try:
            rutas.append(descargar_mes(con, mes))
        except Exception as exc:  # noqa: BLE001
            log(f"  FALLO en el mes {mes}: {exc}")
            log("  (si es error de cuota, espera y reejecuta: es incremental)")
    media_anual(rutas)


if __name__ == "__main__":
    main()
