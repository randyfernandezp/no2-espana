# =============================================================================
# 00 — Instalación y verificación de paquetes
#
# Ejecutar UNA VEZ antes que cualquier otro script de R.
# Desde RStudio: abrir este fichero y pulsar Source (Ctrl+Shift+S).
# Desde PowerShell: & "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" R/00_paquetes.R
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# -----------------------------------------------------------------------------
# Paquetes por bloque, con su función en el proyecto
# -----------------------------------------------------------------------------
PAQUETES <- list(
  # Manipulación de datos y utilidades
  base = c("dplyr", "readr", "tidyr", "stringr", "purrr", "tibble"),

  # Datos espaciales: sf para vectorial, terra para ráster.
  # Ambos compilan contra GDAL/PROJ/GEOS, igual que rasterio en Python.
  # En Windows CRAN sirve binarios ya compilados, así que no debería haber
  # problemas de compilación.
  espacial = c("sf", "terra", "units"),

  # OBJETIVO 1 — geoestadística
  #   gstat  : variogramas empíricos y direccionales, ajuste, kriging, CV
  #   automap: ajuste automático de variograma, útil como contraste
  #   fields : utilidades de campos espaciales
  geoestadistica = c("gstat", "automap", "fields"),

  # OBJETIVO 3 — procesos puntuales
  #   spatstat es un metapaquete; se instalan sus componentes
  puntual = c("spatstat", "spatstat.geom", "spatstat.explore",
              "spatstat.model", "spatstat.random"),

  # OBJETIVO 2 — datos de área
  #   spdep      : matrices de vecindad, Moran, LISA
  #   spatialreg : modelos SAR, SEM, SDM
  areal = c("spdep", "spatialreg"),

  # Gráficos y mapas
  graficos = c("ggplot2", "tmap", "patchwork", "scales", "viridis",
               "classInt", "RColorBrewer")
)

# -----------------------------------------------------------------------------
# Instalación
# -----------------------------------------------------------------------------
instalar_si_falta <- function(pkgs, etiqueta) {
  cat("\n===", etiqueta, "===\n")
  faltan <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(faltan) == 0) {
    cat("  todos presentes\n")
    return(invisible(TRUE))
  }
  cat("  instalando:", paste(faltan, collapse = ", "), "\n")
  install.packages(faltan, dependencies = TRUE)
}

for (nm in names(PAQUETES)) instalar_si_falta(PAQUETES[[nm]], nm)

# -----------------------------------------------------------------------------
# INLA: no está en CRAN, tiene repositorio propio
# -----------------------------------------------------------------------------
# Se usa para dos cosas distintas y ambas centrales:
#   - Objetivo 1: modelo geoestadístico bayesiano vía SPDE, que da la
#     distribución posterior completa del campo. Sin ella no se puede
#     calcular P(NO2 > 20 ug/m3) ni propagar la incertidumbre al agregar
#     por región.
#   - Objetivo 2: modelo areal BYM2, que separa la variación espacialmente
#     estructurada del ruido independiente.
# Si la instalación falla, el proyecto sigue siendo viable con gstat y
# spatialreg; se pierde la cuantificación bayesiana de incertidumbre.
if (!requireNamespace("INLA", quietly = TRUE)) {
  cat("\n=== INLA (repositorio externo) ===\n")
  tryCatch({
    install.packages(
      "INLA",
      repos = c(getOption("repos"),
                INLA = "https://inla.r-inla-download.org/R/stable"),
      dependencies = TRUE
    )
  }, error = function(e) {
    cat("  falló:", conditionMessage(e), "\n")
    cat("  No es bloqueante. Puedes continuar sin INLA y añadirla después.\n")
  })
} else {
  cat("\n=== INLA === presente\n")
}

# -----------------------------------------------------------------------------
# Verificación final
# -----------------------------------------------------------------------------
cat("\n\n", strrep("=", 62), "\n", sep = "")
cat("VERIFICACIÓN\n")
cat(strrep("=", 62), "\n", sep = "")

todos <- c(unlist(PAQUETES, use.names = FALSE), "INLA")
estado <- vapply(todos, requireNamespace, logical(1), quietly = TRUE)

for (p in todos[estado]) cat(sprintf("  [ok] %-22s %s\n", p,
                                     as.character(packageVersion(p))))
if (any(!estado)) {
  cat("\n  FALTAN:\n")
  for (p in todos[!estado]) cat("   -", p, "\n")
}

# Comprobación de que sf y terra encuentran sus librerías geoespaciales:
# es donde suelen aparecer los problemas en Windows.
cat("\n--- Librerías geoespaciales ---\n")
if (requireNamespace("sf", quietly = TRUE)) {
  print(sf::sf_extSoftVersion())
  # Prueba real de reproyección al CRS de trabajo del proyecto
  pt <- sf::st_sfc(sf::st_point(c(-3.70, 40.42)), crs = 4326)
  pr <- sf::st_transform(pt, 3035)
  cat(sprintf("\nPrueba EPSG:4326 -> 3035 (Madrid): %.0f, %.0f\n",
              sf::st_coordinates(pr)[1], sf::st_coordinates(pr)[2]))
  cat("Debe dar aproximadamente 3162000, 2029000\n")
}
