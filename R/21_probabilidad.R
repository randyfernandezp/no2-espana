# =============================================================================
# 21 — Probabilidad de superar los valores límite, y puente al Objetivo 2
#
# Convierte la superficie kriged en el producto de política pública del
# proyecto: dónde España incumpliría la norma que le será exigible en 2030.
#
# ENTRADA : data/processed/superficie_kde.rds, modelo_geo.rds, nuts3.rds
# SALIDAS : output/figuras/10_prob_20.png, 11_prob_40.png, 12_exposicion.png
#           output/tablas/21_exposicion_nuts3.csv
#           data/processed/exposicion_nuts3.rds
#
# VALORES DE REFERENCIA
#   40 ug/m3  Directiva 2008/50/CE, vigente hoy
#   20 ug/m3  Directiva (UE) 2024/2881, exigible desde 2030
#   10 ug/m3  Guía OMS 2021
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(ggplot2)
  library(terra); library(tidyr)
})

PROC <- file.path("data", "processed")
INT  <- file.path("data", "interim")
FIG  <- file.path("output", "figuras")
TAB  <- file.path("output", "tablas")

msg <- function(...) cat(sprintf(...), "\n")
sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n",
                       strrep("=", 68), "\n", sep = "")

kde   <- readRDS(file.path(PROC, "superficie_kde.rds"))
mod   <- readRDS(file.path(PROC, "modelo_geo.rds"))
nuts3 <- readRDS(file.path(PROC, "nuts3.rds"))

UMBRALES <- c(OMS = 10, Directiva_2030 = 20, Directiva_vigente = 40)

# =============================================================================
# 1. PROBABILIDAD DE SUPERACIÓN
# =============================================================================
sec("1. PROBABILIDAD DE SUPERAR CADA UMBRAL")

# El kriging entrega, en cada celda, una predicción y su varianza. Bajo el
# supuesto de normalidad EN LA ESCALA DEL MODELO, eso define una distribución
# predictiva completa, no solo un valor puntual. La probabilidad de superar
# un umbral se obtiene directamente de ella.
#
# Esto es lo que un mapa de valores puntuales no puede dar: distingue entre
# "predigo 19 ug/m3 con mucha certeza" y "predigo 19 pero podría ser 26".
# Para gestión pública, la segunda situación exige actuar y la primera no.
#
# Se hace en la escala del modelo (logarítmica) porque es donde el supuesto
# de normalidad es defendible; el umbral se transforma, no la predicción.
prob_superar <- function(umbral, pred, varianza, log_escala) {
  u <- if (log_escala) log(umbral) else umbral
  1 - pnorm((u - pred) / sqrt(varianza))
}

for (nm in names(UMBRALES)) {
  kde[[paste0("p_", UMBRALES[nm])]] <-
    prob_superar(UMBRALES[nm], kde$var1.pred, kde$var1.var, mod$usar_log)
}

for (nm in names(UMBRALES)) {
  u <- UMBRALES[nm]
  p <- kde[[paste0("p_", u)]]
  msg("%-18s (%2d ug/m3): superficie con P>0,5 = %5.2f %% | con P>0,9 = %5.2f %%",
      nm, u, 100 * mean(p > 0.5), 100 * mean(p > 0.9))
}

# =============================================================================
# 2. MAPAS
# =============================================================================
mapa_prob <- function(campo, titulo, subtitulo) {
  ggplot() +
    geom_sf(data = kde, aes(colour = .data[[campo]]), size = 0.32, shape = 15) +
    geom_sf(data = st_boundary(nuts3), fill = NA, colour = "white",
            linewidth = 0.1) +
    scale_colour_viridis_c(
      option = "inferno", limits = c(0, 1), name = "P(superar)",
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    labs(title = titulo, subtitle = subtitulo) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

ggsave(file.path(FIG, "10_prob_20.png"),
       mapa_prob("p_20", "Probabilidad de superar 20 ug/m3",
                 "Valor límite anual de NO2 exigible desde 2030 (Directiva UE 2024/2881)"),
       width = 7.5, height = 6, dpi = 150)

ggsave(file.path(FIG, "11_prob_40.png"),
       mapa_prob("p_40", "Probabilidad de superar 40 ug/m3",
                 "Valor límite vigente (Directiva 2008/50/CE)"),
       width = 7.5, height = 6, dpi = 150)

ggsave(file.path(FIG, "12_prob_10.png"),
       mapa_prob("p_10", "Probabilidad de superar 10 ug/m3",
                 "Guía de calidad del aire de la OMS (2021)"),
       width = 7.5, height = 6, dpi = 150)

# =============================================================================
# 3. PESO POBLACIONAL DASIMÉTRICO
# =============================================================================
sec("3. REDISTRIBUCIÓN DASIMÉTRICA DE LA POBLACIÓN")

# Eurostat da la población por provincia, pero dentro de cada provincia la
# gente no se reparte de forma uniforme: se concentra en núcleos urbanos.
# Promediar la superficie de NO2 con peso igual por celda respondería a
# "¿cuánto NO2 hay en el territorio?", que no es la pregunta de salud
# pública. La pregunta relevante es "¿a cuánto NO2 está expuesta la gente?".
#
# El método dasimétrico reparte la población conocida de cada provincia
# proporcionalmente a la superficie urbanizada de cada celda. Es una
# aproximación estándar en estudios de exposición y mucho mejor que asumir
# uniformidad, porque las zonas urbanas son a la vez las más pobladas y las
# más contaminadas: ignorar esa coincidencia subestima la exposición real.

# krige() devuelve únicamente predicción, varianza y geometría: no arrastra
# las covariables de la rejilla. Hay que recuperarlas de mod$rejilla, que se
# guardó en 20 y conserva el mismo orden de filas.
if (!"artificial_1000m" %in% names(kde)) {
  rej <- mod$rejilla
  stopifnot(nrow(rej) == nrow(kde))
  cols <- setdiff(names(st_drop_geometry(rej)), names(st_drop_geometry(kde)))
  kde <- bind_cols(kde, st_drop_geometry(rej)[, cols, drop = FALSE])
  msg("Covariables recuperadas de la rejilla: %d columnas", length(cols))
}

peso_urb <- NULL
for (v in c("urbano_denso_1000m", "artificial_1000m")) {
  if (v %in% names(kde) && sum(!is.na(kde[[v]])) > 0.9 * nrow(kde)) {
    peso_urb <- kde[[v]]
    msg("Proxy de población: %s", v)
    break
  }
}
if (is.null(peso_urb)) {
  warning("Sin proxy urbano; se usará peso uniforme.")
  peso_urb <- rep(1, nrow(kde))
}
# Suelo sin urbanizar recibe peso residual, no cero: hay población dispersa.
kde$peso_pob <- pmax(peso_urb, 0.01)

# =============================================================================
# 4. AGREGACIÓN POR PROVINCIA
# =============================================================================
sec("4. EXPOSICIÓN POR PROVINCIA (NUTS-3)")

# Agregar la superficie predicha, y no la media de las estaciones, tiene dos
# ventajas decisivas para el Objetivo 2:
#   - Cubre las provincias sin ninguna estación, que de otro modo quedarían
#     como datos ausentes en el modelo areal.
#   - Pondera por territorio o por población, no por densidad de red. La
#     media de estaciones de Madrid está dominada por sus muchas estaciones
#     urbanas; la de Soria depende de una sola.
asig <- st_join(kde, nuts3["NUTS_ID"], join = st_within)

expos <- asig |>
  st_drop_geometry() |>
  filter(!is.na(NUTS_ID)) |>
  group_by(NUTS_ID) |>
  summarise(
    n_celdas       = n(),
    no2_areal      = mean(no2, na.rm = TRUE),
    no2_poblacion  = weighted.mean(no2, peso_pob, na.rm = TRUE),
    no2_p90        = quantile(no2, 0.90, na.rm = TRUE),
    no2_max        = max(no2, na.rm = TRUE),
    # Superficie y población que superan cada umbral
    sup_sobre_20   = mean(p_20 > 0.5, na.rm = TRUE),
    sup_sobre_40   = mean(p_40 > 0.5, na.rm = TRUE),
    pob_sobre_20   = weighted.mean(p_20 > 0.5, peso_pob, na.rm = TRUE),
    pob_sobre_10   = weighted.mean(p_10 > 0.5, peso_pob, na.rm = TRUE),
    # Probabilidad media de superación, que integra la incertidumbre en
    # lugar de dicotomizarla con un corte en 0,5
    prob_media_20  = weighted.mean(p_20, peso_pob, na.rm = TRUE),
    .groups = "drop"
  )

nuts3 <- nuts3 |> left_join(expos, by = "NUTS_ID")

msg("Provincias con exposición estimada: %d de %d",
    sum(!is.na(nuts3$no2_poblacion)), nrow(nuts3))

# --- Contraste: media areal frente a media ponderada por población ----------
# La brecha entre ambas mide cuánto se subestima la exposición al ignorar
# que la gente vive donde más contaminación hay. Debe ser positiva y mayor
# en las provincias más urbanizadas.
nuts3$brecha_exposicion <- nuts3$no2_poblacion - nuts3$no2_areal
msg("Brecha poblacional-areal: mediana %+.2f | máxima %+.2f ug/m3",
    median(nuts3$brecha_exposicion, na.rm = TRUE),
    max(nuts3$brecha_exposicion, na.rm = TRUE))

cat("\nDiez provincias con mayor exposición poblacional:\n")
print(
  nuts3 |> st_drop_geometry() |>
    select(NUTS_ID, NUTS_NAME, no2_areal, no2_poblacion, brecha_exposicion,
           pob_sobre_20, prob_media_20) |>
    arrange(desc(no2_poblacion)) |> head(10) |>
    mutate(across(where(is.numeric), ~ round(.x, 3))) |>
    as.data.frame()
)

cat("\nProvincias donde más del 10 % de la población superaría 20 ug/m3:\n")
riesgo <- nuts3 |> st_drop_geometry() |>
  filter(pob_sobre_20 > 0.10) |>
  select(NUTS_ID, NUTS_NAME, no2_poblacion, pob_sobre_20, sup_sobre_20) |>
  arrange(desc(pob_sobre_20))
if (nrow(riesgo) == 0) {
  msg("Ninguna. El umbral de 2030 no se superaría de forma generalizada,")
  msg("aunque conviene mirar prob_media_20: puede haber mucha población")
  msg("con probabilidad intermedia, que la dicotomía en 0,5 oculta.")
} else {
  print(riesgo |> mutate(across(where(is.numeric), ~ round(.x, 3))) |>
          as.data.frame())
}

# =============================================================================
# 5. MAPAS DE EXPOSICIÓN
# =============================================================================
mapa_nuts <- function(campo, titulo, etiqueta, paleta = "C") {
  ggplot(nuts3) +
    geom_sf(aes(fill = .data[[campo]]), colour = "white", linewidth = 0.15) +
    scale_fill_viridis_c(option = paleta, name = etiqueta, na.value = "grey85") +
    labs(title = titulo) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

ggsave(file.path(FIG, "13_exposicion_poblacional.png"),
       mapa_nuts("no2_poblacion",
                 "NO2 medio ponderado por población",
                 expression(mu*g/m^3)),
       width = 7, height = 6, dpi = 150)

ggsave(file.path(FIG, "14_poblacion_sobre_20.png"),
       mapa_nuts("pob_sobre_20",
                 "Fracción de población que superaría 20 ug/m3",
                 "fracción", "B"),
       width = 7, height = 6, dpi = 150)

ggsave(file.path(FIG, "15_brecha_exposicion.png"),
       mapa_nuts("brecha_exposicion",
                 "Brecha entre exposición poblacional y media territorial",
                 expression(mu*g/m^3), "D"),
       width = 7, height = 6, dpi = 150)

# =============================================================================
# 6. GUARDAR
# =============================================================================
write_csv(st_drop_geometry(nuts3), file.path(TAB, "21_exposicion_nuts3.csv"))
saveRDS(nuts3, file.path(PROC, "nuts3_exposicion.rds"))
saveRDS(kde, file.path(PROC, "superficie_kde_prob.rds"))

sec("RESUMEN")
msg("Superficie nacional con P(NO2>20) > 0,5 : %.2f %%",
    100 * mean(kde$p_20 > 0.5))
msg("Población nacional con P(NO2>20) > 0,5 : %.2f %%",
    100 * weighted.mean(kde$p_20 > 0.5, kde$peso_pob))
msg("Población nacional con P(NO2>10) > 0,5 : %.2f %%",
    100 * weighted.mean(kde$p_10 > 0.5, kde$peso_pob))
# =============================================================================
# 7. CONTRASTE CON LAS MEDICIONES REALES
# =============================================================================
sec("7. QUÉ MIDE ESTA SUPERFICIE Y QUÉ NO")

# La superficie interpola el campo de FONDO regional, porque se ajustó solo
# con estaciones de fondo (ver la justificación en 10_prepare.R). No
# representa los máximos de exposición junto a vías de tráfico intenso, que
# ocurren a escala de decenas de metros y quedan por debajo de la resolución
# de 2 km de la rejilla.
#
# Interpretar el mapa como "exposición máxima" sería un error grave: las
# cifras de superación de fondo son bajas, pero la población que vive junto a
# ejes viarios está expuesta al fondo MÁS el recargo local, que el script 20
# estimó entre 2,9 y 3,8 ug/m3 en estaciones de tráfico.
#
# El contraste directo con las mediciones sitúa el mapa en su lugar.
est <- readRDS(file.path(PROC, "estaciones.rds"))
cat("\nEstaciones que superan cada umbral (medición directa, 2023):\n")
print(
  est |> st_drop_geometry() |>
    group_by(tipo) |>
    summarise(
      n = n(),
      media = round(mean(no2_media), 2),
      maxima = round(max(no2_media), 2),
      `sobre_10` = round(mean(no2_media > 10), 3),
      `sobre_20` = round(mean(no2_media > 20), 3),
      `sobre_40` = round(mean(no2_media > 40), 3),
      .groups = "drop"
    ) |> as.data.frame()
)

msg("\nSuperficie de fondo con P(NO2>20)>0,5 : %.2f %%", 100 * mean(kde$p_20 > 0.5))
msg("Estaciones de tráfico que miden >20   : %.1f %%",
    100 * mean(est$no2_media[est$tipo == "traffic"] > 20, na.rm = TRUE))
cat("La distancia entre ambas cifras es la magnitud del fenómeno de\n")
cat("microescala que la geoestadística regional no puede capturar, y es\n")
cat("justo lo que el proceso puntual y el modelo areal deben complementar.\n")

cat("\nLa comparación entre estas tres cifras es el mensaje de política\n")
cat("pública: España cumple hoy con holgura, tiene margen ajustado frente\n")
cat("al umbral de 2030, y está lejos de la recomendación de la OMS.\n")
cat("\nSiguiente: R/40_procesos_puntuales.R (Objetivo 3)\n")
