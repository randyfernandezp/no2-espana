# =============================================================================
# 20 — OBJETIVO 1: variación espacial continua
#
#   A. Exploratorio, transformación y desagrupamiento
#   B. Tendencia de primer orden
#   C. Variograma empírico: omnidireccional y direccional (anisotropía)
#   D. Ajuste de modelos (exponencial, esférico, Matérn)
#   E. Rejilla de predicción con covariables
#   F. Kriging ordinario (KO) y con deriva externa (KDE)
#   G. Validación cruzada dejando uno fuera + validación externa
#
# ENTRADA : data/processed/estaciones.rds, nuts3.rds, dominio.rds
# SALIDAS : output/figuras/*.png, output/tablas/*.csv
#           data/processed/superficie_ko.rds, superficie_kde.rds
#           data/processed/modelo_geo.rds
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(gstat); library(dplyr)
  library(readr); library(tidyr); library(ggplot2); library(stringr)
})

CRS_T <- 3035
PROC  <- file.path("data", "processed")
INT   <- file.path("data", "interim")
RAW   <- file.path("data", "raw")
FIG   <- file.path("output", "figuras")
TAB   <- file.path("output", "tablas")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf(...), "\n")
sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n",
                       strrep("=", 68), "\n", sep = "")

est     <- readRDS(file.path(PROC, "estaciones.rds"))
nuts3   <- readRDS(file.path(PROC, "nuts3.rds"))
dominio <- readRDS(file.path(PROC, "dominio.rds"))

fondo   <- est |> filter(es_fondo, !is.na(no2_media))
externa <- est |> filter(!es_fondo, !is.na(no2_media))
msg("Ajuste: %d estaciones de fondo | Validación externa: %d",
    nrow(fondo), nrow(externa))

# =============================================================================
# A. EXPLORATORIO
# =============================================================================
sec("A. ANÁLISIS EXPLORATORIO")

print(summary(fondo$no2_media))
asim <- function(x) mean((x - mean(x))^3) / sd(x)^3
msg("Asimetría: bruto %.3f | log %.3f",
    asim(fondo$no2_media), asim(log(fondo$no2_media)))

sw_b <- shapiro.test(fondo$no2_media)
sw_l <- shapiro.test(log(fondo$no2_media))
msg("Shapiro-Wilk: bruto W=%.4f p=%.3g | log W=%.4f p=%.3g",
    sw_b$statistic, sw_b$p.value, sw_l$statistic, sw_l$p.value)

# El kriging no exige normalidad para ser el mejor predictor lineal
# insesgado, pero SÍ la asumen los intervalos de predicción, y de ellos
# depende el mapa de probabilidad de superar el valor límite. Por eso la
# decisión se toma con el contraste, no por costumbre.
USAR_LOG <- sw_l$statistic > sw_b$statistic
fondo$z  <- if (USAR_LOG) log(fondo$no2_media) else fondo$no2_media
msg("Variable de trabajo: %s", if (USAR_LOG) "log(NO2)" else "NO2")

# --- Desagrupamiento por celdas ---------------------------------------------
# Las estaciones no son una muestra aleatoria del territorio: se instalan
# donde hay población y donde se sospecha un problema. Madrid concentra
# decenas en pocos km2 mientras que Soria tiene una en 10 000. Sin ponderar,
# la media muestral estima la media de la RED, no la del TERRITORIO, y el
# variograma queda dominado por pares urbanos a corta distancia.
CELDA <- 25000
celdas <- st_make_grid(dominio, cellsize = CELDA)
idx <- st_intersects(fondo, celdas) |> as.numeric()
fondo <- fondo |>
  group_by(celda = idx) |>
  mutate(peso = 1 / n()) |>
  ungroup()

msg("Celdas de %d km ocupadas: %d", CELDA / 1000, n_distinct(fondo$celda))
msg("Media aritmética %.3f | desagrupada %.3f ug/m3",
    mean(fondo$no2_media), weighted.mean(fondo$no2_media, fondo$peso))

p <- ggplot(st_drop_geometry(fondo), aes(no2_media)) +
  geom_histogram(bins = 30, fill = "grey35") +
  labs(x = expression(NO[2]~"medio anual"~(mu*g/m^3)), y = "estaciones",
       title = "Distribución en estaciones de fondo") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG, "01_histograma.png"), p, width = 6, height = 4, dpi = 150)

# =============================================================================
# B. TENDENCIA DE PRIMER ORDEN
# =============================================================================
sec("B. TENDENCIA Y COVARIABLES")

fondo$xk <- fondo$x / 1000; fondo$yk <- fondo$y / 1000  # km, mejor condicionado

tend <- lm(z ~ xk + yk + I(xk^2) + I(yk^2) + I(xk * yk), data = fondo)
msg("Tendencia cuadrática en coordenadas: R2 = %.4f",
    summary(tend)$r.squared)

CANDIDATAS <- c("s5p_no2", "altitud", "artificial_1000m", "artificial_5000m",
                "industrial_5000m", "urbano_denso_1000m", "transporte_1000m",
                "tierra_5000m", "poblacion_ciudad")
disp <- intersect(CANDIDATAS, names(fondo))

cat("\nCorrelación de Spearman con", if (USAR_LOG) "log(NO2)" else "NO2", ":\n")
correl <- sapply(disp, function(v) {
  x <- st_drop_geometry(fondo)[[v]]
  if (all(is.na(x))) NA else cor(x, fondo$z, method = "spearman",
                                 use = "complete.obs")
})
print(round(sort(correl, decreasing = TRUE), 3))
write_csv(tibble(covariable = names(correl), rho = as.numeric(correl)),
          file.path(TAB, "20_correlaciones.csv"))

# Modelo lineal con las covariables disponibles en TODA la rejilla.
# La altitud se excluye del KDE salvo que exista un modelo digital: el
# kriging con deriva externa exige conocer la covariable en cada punto de
# predicción, no solo donde se observa.
f_lm <- as.formula(paste("z ~", paste(setdiff(disp, "altitud"), collapse = " + ")))
mod_lm <- lm(f_lm, data = fondo)
cat("\nModelo lineal de deriva:\n")
print(summary(mod_lm))

# =============================================================================
# C. VARIOGRAMA EMPÍRICO
# =============================================================================
sec("C. VARIOGRAMA EMPÍRICO")

bb <- st_bbox(dominio)
diag_dom <- sqrt((bb$xmax - bb$xmin)^2 + (bb$ymax - bb$ymin)^2)
# El cutoff debe cubrir la escala a la que ocurre la estructura, no la del
# dominio. Con un cutoff de 475 km y bines de 19 km, la estructura real
# (decenas de km) cae en los dos o tres primeros puntos y el ajuste queda
# determinado por bines de meseta plana que no aportan información.
# Se toma un cutoff mucho menor y bines finos, lo que da ~30 puntos
# concentrados donde el variograma cambia.
CUTOFF <- min(diag_dom / 10, 200000)
ANCHO  <- CUTOFF / 30
msg("Diagonal del dominio %.0f km | cutoff %.0f km | ancho %.1f km",
    diag_dom / 1000, CUTOFF / 1000, ANCHO / 1000)

# Rango PRÁCTICO: distancia a la que el variograma alcanza el 95 % de su
# meseta. A diferencia del parámetro `range`, es comparable entre familias:
# el `range` de un Matérn con kappa alto y el de un esférico no significan
# lo mismo, y compararlos directamente induce a error.
rango_practico <- function(m, maxd = 600000) {
  d <- seq(1, maxd, length.out = 4000)
  g <- gstat::variogramLine(m, dist_vector = d)$gamma
  meseta <- sum(m$psill)
  i <- which(g >= 0.95 * meseta)[1]
  if (is.na(i)) NA_real_ else d[i]
}

vg <- variogram(z ~ 1, data = fondo, cutoff = CUTOFF, width = ANCHO)
png(file.path(FIG, "02_variograma_omni.png"), 900, 600, res = 110)
print(plot(vg, main = "Variograma empírico omnidireccional"))
dev.off()

# Variograma direccional. En España la anisotropía no es un tecnicismo: los
# vientos dominantes y la orientación de los valles del Ebro (NO-SE) y del
# Guadalquivir (NE-SO) generan continuidad espacial preferente. Si aparece,
# es un resultado interpretable, no un problema a corregir.
vg_dir <- variogram(z ~ 1, data = fondo, cutoff = CUTOFF,
                    width = ANCHO * 1.5, alpha = c(0, 45, 90, 135),
                    tol.hor = 22.5)
png(file.path(FIG, "03_variograma_direccional.png"), 1100, 800, res = 110)
print(plot(vg_dir, main = "Variogramas direccionales (0/45/90/135 grados)"))
dev.off()

rangos_dir <- vg_dir |>
  group_by(dir.hor) |>
  summarise(gamma_max = max(gamma), n_pares = sum(np), .groups = "drop")
print(rangos_dir)
msg("Si gamma_max difiere mucho entre direcciones, hay anisotropía zonal.")

# =============================================================================
# D. AJUSTE DE MODELOS
# =============================================================================
sec("D. AJUSTE DEL VARIOGRAMA")

s0 <- var(fondo$z); r0 <- CUTOFF / 4
candidatos <- list(
  Exp = vgm(0.75 * s0, "Exp", r0, 0.25 * s0),
  Sph = vgm(0.75 * s0, "Sph", r0, 0.25 * s0),
  Gau = vgm(0.75 * s0, "Gau", r0, 0.25 * s0),
  Mat = vgm(0.75 * s0, "Mat", r0, 0.25 * s0, kappa = 1)
)

ajustes <- lapply(names(candidatos), function(nm) {
  fit <- tryCatch(
    # kappa se restringe a [0.3, 2.5]. Sin acotar, el ajuste se pega al tope
    # de 5 que explora gstat, señal de inestabilidad: un Matérn con kappa 5
    # es casi gaussiano, implica un campo diferenciable infinitas veces y
    # produce matrices de covarianza mal condicionadas en el kriging.
    fit.variogram(vg, candidatos[[nm]],
                  fit.kappa = if (nm == "Mat") seq(0.3, 2.5, 0.1) else FALSE),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(fit) || any(fit$psill < 0)) return(NULL)
  list(nombre = nm, m = fit, sce = attr(fit, "SSErr"))
})
ajustes <- Filter(Negate(is.null), ajustes)
if (length(ajustes) == 0) stop("Ningún modelo de variograma convergió.")

tabla <- bind_rows(lapply(ajustes, function(a) {
  m <- a$m; k <- nrow(m)
  tibble(modelo = a$nombre,
         nugget = m$psill[1],
         sill_parcial = m$psill[k],
         sill_total = sum(m$psill),
         # Pepita relativa: por debajo de 0,25 hay estructura espacial
         # fuerte; por encima de 0,75 el campo es casi ruido y el kriging
         # aporta poco sobre la media.
         pepita_rel = m$psill[1] / sum(m$psill),
         rango_par_km = m$range[k] / 1000,
         rango_practico_km = rango_practico(m) / 1000,
         kappa = if (!is.null(m$kappa)) m$kappa[k] else NA_real_,
         SCE = a$sce)
})) |> arrange(SCE)
print(as.data.frame(tabla))
write_csv(tabla, file.path(TAB, "20_ajuste_variogramas.csv"))

mejor <- ajustes[[which.min(sapply(ajustes, `[[`, "sce"))]]
vg_fit <- mejor$m
msg("Seleccionado: %s | rango práctico %.0f km | pepita relativa %.1f%%",
    mejor$nombre, rango_practico(vg_fit) / 1000,
    100 * vg_fit$psill[1] / sum(vg_fit$psill))
if (vg_fit$psill[1] < 1e-8)
  msg("AVISO: pepita nula. Implica ausencia de error de medida y de")
  msg("variación por debajo de la menor distancia entre estaciones,")
  msg("supuesto poco creíble en calidad del aire. Revisar el primer bin.")

png(file.path(FIG, "04_variograma_ajustado.png"), 900, 600, res = 110)
print(plot(vg, vg_fit, main = sprintf("Ajuste %s", mejor$nombre)))
dev.off()

# =============================================================================
# E. REJILLA DE PREDICCIÓN CON COVARIABLES
# =============================================================================
sec("E. REJILLA DE PREDICCIÓN")

RES <- 2000
rej <- st_make_grid(dominio, cellsize = RES, what = "centers") |>
  st_as_sf() |> st_filter(dominio)
names(rej)[1] <- "geometry"; st_geometry(rej) <- "geometry"
xy <- st_coordinates(rej)
rej$x <- xy[, 1]; rej$y <- xy[, 2]
rej$xk <- rej$x / 1000; rej$yk <- rej$y / 1000
msg("Rejilla: %d celdas de %d m", nrow(rej), RES)

# --- Covariables sobre la rejilla -------------------------------------------
# Aquí se paga el trabajo de 07_corine.py: los rásteres focales permiten
# evaluar la composición del entorno en cualquier punto, no solo en las
# estaciones. Sin ellos el KDE sería inviable.
rej_v <- vect(rej)

ruta_s5p <- file.path(RAW, "s5p_no2_2023_media.tif")
if (file.exists(ruta_s5p)) {
  s5p <- rast(ruta_s5p)
  rej$s5p_no2 <- terra::extract(s5p, vect(st_transform(rej, crs(s5p))))[, 2] * 1e5
}

for (v in setdiff(disp, c("s5p_no2", "altitud", "poblacion_ciudad"))) {
  f <- file.path(INT, "focal_norm", paste0("frac_", v, ".tif"))
  if (!file.exists(f)) f <- file.path(INT, "focal", paste0("frac_", v, ".tif"))
  if (file.exists(f)) {
    rej[[v]] <- terra::extract(rast(f), rej_v)[, 2]
    msg("  %s: %d/%d válidos en la rejilla",
        v, sum(!is.na(rej[[v]])), nrow(rej))
  }
}

covs_kde <- intersect(names(rej), setdiff(disp, c("altitud", "poblacion_ciudad")))
covs_kde <- covs_kde[sapply(covs_kde, function(v) sum(is.na(rej[[v]])) < 0.02 * nrow(rej))]
msg("Covariables utilizables en KDE: %s", paste(covs_kde, collapse = ", "))

completas <- complete.cases(st_drop_geometry(rej)[, covs_kde, drop = FALSE])
rej_ok <- rej[completas, ]
msg("Celdas con todas las covariables: %d de %d", nrow(rej_ok), nrow(rej))

# =============================================================================
# F. KRIGING
# =============================================================================
sec("F. KRIGING")

msg("Kriging ordinario …")
ko <- krige(z ~ 1, fondo, rej_ok, model = vg_fit, debug.level = 0)

f_kde <- as.formula(paste("z ~", paste(covs_kde, collapse = " + ")))
msg("Kriging con deriva externa: %s", deparse(f_kde))

# El variograma del KDE debe estimarse sobre los RESIDUOS de la deriva, no
# sobre la variable bruta. Que el rango residual sea menor que el bruto es
# la confirmación de que las covariables han absorbido la estructura de gran
# escala, que es justo su cometido.
vg_res <- variogram(f_kde, data = fondo, cutoff = CUTOFF, width = ANCHO)
s_res <- var(residuals(lm(f_kde, data = fondo)))
vg_res_fit <- fit.variogram(vg_res,
                            vgm(0.7 * s_res, mejor$nombre, r0 * 0.6, 0.3 * s_res))
msg("Rango práctico: residual %.0f km  vs  bruto %.0f km",
    rango_practico(vg_res_fit) / 1000, rango_practico(vg_fit) / 1000)
msg("Varianza: bruta %.4f -> residual %.4f (%.0f %% absorbido por la deriva)",
    var(fondo$z), s_res, 100 * (1 - s_res / var(fondo$z)))

png(file.path(FIG, "05_variograma_residual.png"), 900, 600, res = 110)
print(plot(vg_res, vg_res_fit, main = "Variograma de los residuos de la deriva"))
dev.off()

kde <- krige(f_kde, fondo, rej_ok, model = vg_res_fit, debug.level = 0)

# --- Retrotransformación -----------------------------------------------------
# Con logaritmos, exp(mu) estima la MEDIANA, no la media. El corrector
# lognormal exp(mu + sigma2/2) devuelve la media, que es la magnitud que
# compara con el valor límite normativo.
retro <- function(pred, varianza) {
  if (USAR_LOG) exp(pred + varianza / 2) else pred
}
ko$no2 <- retro(ko$var1.pred, ko$var1.var)
kde$no2 <- retro(kde$var1.pred, kde$var1.var)

# Diagnóstico de extrapolación. El KDE aplica el modelo de deriva a celdas
# cuyas covariables pueden caer fuera del rango observado en las estaciones.
# La red española se concentra en zonas urbanas, mientras que la mayor parte
# del territorio es rural con fracción artificial cercana a cero: buena parte
# de la rejilla queda por debajo del mínimo calibrado, y ahí la predicción es
# extrapolación, no interpolación.
cat("\nExtrapolación de covariables fuera del rango calibrado:\n")
for (v in covs_kde) {
  rg <- range(st_drop_geometry(fondo)[[v]], na.rm = TRUE)
  x <- rej_ok[[v]]
  fuera <- mean(x < rg[1] | x > rg[2], na.rm = TRUE)
  msg("  %-20s estaciones [%.3f, %.3f] | %.1f %% de la rejilla fuera",
      v, rg[1], rg[2], 100 * fuera)
}

msg("KO  -> media %.2f | rango %.2f - %.2f ug/m3",
    mean(ko$no2), min(ko$no2), max(ko$no2))
msg("KDE -> media %.2f | rango %.2f - %.2f ug/m3",
    mean(kde$no2), min(kde$no2), max(kde$no2))

# =============================================================================
# G. VALIDACIÓN
# =============================================================================
sec("G. VALIDACIÓN CRUZADA")

# Se calculan por separado en dos escalas, porque mezclarlas no tiene
# sentido: la varianza de kriging vive en la escala del MODELO (logarítmica
# si se transformó), mientras que RMSE y sesgo solo son interpretables en
# la escala ORIGINAL, en ug/m3.
metricas_orig <- function(obs, pred) {
  e <- pred - obs
  c(n = length(obs), RMSE = sqrt(mean(e^2)), MAE = mean(abs(e)),
    sesgo = mean(e), R2 = 1 - sum(e^2) / sum((obs - mean(obs))^2))
}

# En la escala del modelo se evalúa la CALIBRACIÓN de la incertidumbre.
# La cobertura empírica del intervalo nominal al 95 % debería rondar 0,95 y
# la desviación de los residuos estandarizados debería rondar 1. Si la
# desviación es mucho mayor que 1, el modelo es sobreconfiado y el mapa de
# probabilidad de superar el valor límite resultaría engañoso.
metricas_modelo <- function(obs_z, pred_z, var_z) {
  e <- pred_z - obs_z
  zt <- e / sqrt(var_z)
  c(RMSE_log = sqrt(mean(e^2)),
    cobertura95 = mean(abs(zt) < 1.96),
    ESD = sd(zt))
}

cv_ko <- krige.cv(z ~ 1, fondo, model = vg_fit, nfold = nrow(fondo))
cv_kde <- krige.cv(f_kde, fondo, model = vg_res_fit, nfold = nrow(fondo))

obs_o <- fondo$no2_media
obs_z <- fondo$z

res_cv <- bind_rows(
  tibble(modelo = "Kriging ordinario",
         !!!as.list(metricas_orig(obs_o, retro(cv_ko$var1.pred, cv_ko$var1.var))),
         !!!as.list(metricas_modelo(obs_z, cv_ko$var1.pred, cv_ko$var1.var))),
  tibble(modelo = "Kriging con deriva externa",
         !!!as.list(metricas_orig(obs_o, retro(cv_kde$var1.pred, cv_kde$var1.var))),
         !!!as.list(metricas_modelo(obs_z, cv_kde$var1.pred, cv_kde$var1.var)))
)
print(as.data.frame(res_cv))
write_csv(res_cv, file.path(TAB, "20_validacion_cruzada.csv"))

mejora <- 100 * (1 - res_cv$RMSE[2] / res_cv$RMSE[1])
msg("Las covariables reducen el RMSE un %.1f %%", mejora)

# --- Validación externa: el recargo local ------------------------------------
# El modelo se ajustó SIN estaciones de tráfico e industriales. El sesgo
# medido en ellas no es un fallo: cuantifica cuánto excede la concentración
# local al fondo regional, separando "fondo" de "recargo". Es un resultado
# sustantivo para el Objetivo 5.
if (nrow(externa) > 0) {
  pe <- krige(f_kde, fondo, externa, model = vg_res_fit, debug.level = 0)
  externa$pred_fondo <- retro(pe$var1.pred, pe$var1.var)
  externa$recargo <- externa$no2_media - externa$pred_fondo

  cat("\nRecargo local sobre el fondo regional (ug/m3):\n")
  print(externa |> st_drop_geometry() |>
          group_by(tipo, area) |>
          summarise(n = n(), mediana = median(recargo),
                    media = mean(recargo), .groups = "drop") |>
          mutate(across(where(is.numeric), ~ round(.x, 2))) |>
          as.data.frame())

  write_csv(externa |> st_drop_geometry() |>
              select(join_id, tipo, area, no2_media, pred_fondo, recargo),
            file.path(TAB, "20_recargo_local.csv"))
}

# =============================================================================
# MAPAS Y GUARDADO
# =============================================================================
mapa <- function(sfobj, campo, titulo, paleta = "C") {
  ggplot() +
    geom_sf(data = sfobj, aes(colour = .data[[campo]]), size = 0.35,
            shape = 15) +
    geom_sf(data = st_boundary(nuts3), fill = NA, colour = "white",
            linewidth = 0.12) +
    scale_colour_viridis_c(option = paleta, name = expression(mu*g/m^3)) +
    labs(title = titulo) +
    theme_void(base_size = 11)
}

ggsave(file.path(FIG, "06_mapa_ko.png"),
       mapa(ko, "no2", "Kriging ordinario"), width = 7, height = 6, dpi = 150)
ggsave(file.path(FIG, "07_mapa_kde.png"),
       mapa(kde, "no2", "Kriging con deriva externa"),
       width = 7, height = 6, dpi = 150)

kde$ee <- sqrt(kde$var1.var)
ggsave(file.path(FIG, "08_mapa_error.png"),
       mapa(kde, "ee", "Error estándar de predicción (escala del modelo)", "A"),
       width = 7, height = 6, dpi = 150)

saveRDS(ko,  file.path(PROC, "superficie_ko.rds"))
saveRDS(kde, file.path(PROC, "superficie_kde.rds"))
saveRDS(list(vg_fit = vg_fit, vg_res_fit = vg_res_fit, usar_log = USAR_LOG,
             formula_kde = f_kde, covs = covs_kde, rejilla = rej_ok,
             cutoff = CUTOFF, modelo = mejor$nombre),
        file.path(PROC, "modelo_geo.rds"))

sec("OBJETIVO 1 COMPLETADO")
cat("Figuras en output/figuras/ y tablas en output/tablas/\n")
cat("Siguiente: R/30_spde_inla.R para la version bayesiana y el mapa de\n")
cat("probabilidad de superar los 20 ug/m3 exigibles desde 2030.\n")
