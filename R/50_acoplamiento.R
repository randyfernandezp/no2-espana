# =============================================================================
# 50 — OBJETIVO 4: acoplamiento entre el proceso puntual y el campo continuo
#
# Pregunta: ¿la presión de emisión industrial explica la concentración de NO2
# más allá de lo que ya explican el uso del suelo y la señal satelital?
# Y si la explica, ¿a qué ESCALA ESPACIAL opera esa influencia?
#
# CONTRIBUCIÓN ORIGINAL
# ---------------------
# Se construye la superficie de intensidad del proceso puntual (Objetivo 3)
# con distintos anchos de banda y se introduce cada una como covariable en
# el modelo geoestadístico (Objetivo 1). El ancho que minimiza el error de
# validación cruzada ESTIMA la escala espacial a la que la emisión industrial
# influye sobre la concentración de fondo.
#
# Ese número tiene interpretación física directa: es la distancia
# característica de transporte y dilución del NOx emitido antes de dejar de
# ser distinguible del fondo. Y puede contrastarse con la escala del
# conglomerado estimada por el modelo de Thomas: si coinciden, la unidad de
# influencia es el distrito industrial; si la de influencia es mucho mayor,
# el mecanismo relevante es el transporte atmosférico regional.
#
# ENTRADA : data/processed/pp_industrias.rds, estaciones.rds, modelo_geo.rds
# SALIDAS : output/figuras/3x_*.png, output/tablas/50_*.csv
#           data/processed/acoplamiento.rds
# =============================================================================

suppressPackageStartupMessages({
  library(spatstat.geom); library(spatstat.explore)
  library(sf); library(terra); library(gstat); library(dplyr)
  library(readr); library(ggplot2); library(tidyr)
})

PROC <- file.path("data", "processed")
FIG  <- file.path("output", "figuras")
TAB  <- file.path("output", "tablas")

msg <- function(...) cat(sprintf(...), "\n")
sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n",
                       strrep("=", 68), "\n", sep = "")

pp   <- readRDS(file.path(PROC, "pp_industrias.rds"))
est  <- readRDS(file.path(PROC, "estaciones.rds"))
mod  <- readRDS(file.path(PROC, "modelo_geo.rds"))
vent <- readRDS(file.path(PROC, "ventana.rds"))

fondo <- est |> filter(es_fondo, !is.na(no2_media))
fondo$z <- if (mod$usar_log) log(fondo$no2_media) else fondo$no2_media
msg("Estaciones de fondo: %d | Instalaciones: %d", nrow(fondo), npoints(pp))

# =============================================================================
# 1. SUPERFICIES DE PRESIÓN DE EMISIÓN A DISTINTAS ESCALAS
# =============================================================================
sec("1. BARRIDO DE ANCHOS DE BANDA")

ANCHOS <- c(5, 10, 15, 25, 40, 60, 100) * 1000   # metros
covs_base <- mod$covs

# Conversión de imagen de spatstat a ráster de terra.
# spatstat ordena las filas de sur a norte y terra de norte a sur, de modo
# que hay que invertir el orden o la superficie saldrá reflejada.
im_a_rast <- function(x) {
  m <- as.matrix(x)                 # filas: yrow ascendente
  m <- m[nrow(m):1, ]               # invertir a orientación de terra
  r <- rast(nrows = nrow(m), ncols = ncol(m),
            xmin = x$xrange[1], xmax = x$xrange[2],
            ymin = x$yrange[1], ymax = x$yrange[2],
            crs = "EPSG:3035")
  values(r) <- as.vector(t(m))
  r
}

# ¿Hay marcas de masa? Si las hay, la superficie pondera por toneladas de
# NOx; si no, todas las instalaciones pesan igual y el resultado es una
# versión degradada que confunde una térmica con una cerámica.
tiene_marcas <- is.marked(pp) &&
  (if (is.data.frame(marks(pp))) any(!is.na(marks(pp)$nox_t)) else
     any(!is.na(marks(pp))))

if (tiene_marcas) {
  m <- marks(pp)
  w <- if (is.data.frame(m)) m$nox_t else as.numeric(m)
  ok <- !is.na(w) & w > 0
  pp_w <- pp[ok]; pesos <- w[ok]
  msg("Ponderando por masa de NOx: %d instalaciones, %.0f t totales",
      sum(ok), sum(pesos))
} else {
  pp_w <- unmark(pp); pesos <- NULL
  msg("SIN marcas de masa: todas las instalaciones pesan igual.")
  msg("Descarga eprtr_emisiones.csv para ponderar por toneladas de NOx.")
}

# density.ppp deja sin valor los píxeles exteriores a la ventana. Tras
# simplificar la frontera a 1 km, algunas estaciones costeras caen en esos
# píxeles y quedan como NA. Un solo NA basta para que lm y krige.cv se
# ajusten sobre submuestras distintas y las comparaciones dejen de ser
# válidas, así que se rellenan por vecindad antes de seguir.
rellenar_na <- function(r, pts, intentos = 4) {
  v <- terra::extract(r, pts)[, 2]
  k <- 0
  while (any(is.na(v)) && k < intentos) {
    k <- k + 1
    r <- terra::focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE)
    v[is.na(v)] <- terra::extract(r, pts[is.na(v)])[, 2]
  }
  list(v = v, n_rellenados = k)
}

presiones <- list()
pts_fondo <- vect(fondo)

for (h in ANCHOS) {
  lam <- density(pp_w, sigma = h, weights = pesos, positive = TRUE,
                 edge = TRUE)
  r <- im_a_rast(lam)
  nombre <- sprintf("presion_%dkm", h / 1000)

  ext <- rellenar_na(r, pts_fondo)
  v <- ext$v
  n_na <- sum(is.na(v))

  # Transformación logarítmica con desplazamiento PROPORCIONAL, no un
  # epsilon arbitrario. Con sigma pequeño la mayor parte del territorio
  # tiene intensidad casi nula, y log(v + 1e-14) amontonaría el 80 % de las
  # estaciones en un mismo valor de suelo, creando un artefacto. Sumar el
  # 1 % de la intensidad media reparte esos casos de forma continua y hace
  # comparables entre sí las distintas escalas.
  c_desp <- 0.01 * mean(v, na.rm = TRUE)
  fondo[[nombre]] <- log(v + c_desp)
  presiones[[nombre]] <- r

  msg("  sigma = %3d km | NA %d | log-presión de %.2f a %.2f | en el suelo %.0f %%",
      h / 1000, n_na,
      min(fondo[[nombre]], na.rm = TRUE), max(fondo[[nombre]], na.rm = TRUE),
      100 * mean(v < c_desp, na.rm = TRUE))
}

# Comprobación imprescindible antes de comparar modelos: todos deben
# ajustarse sobre las mismas estaciones.
cols_p <- names(presiones)
n_inc <- sum(!complete.cases(st_drop_geometry(fondo)[, c(covs_base, cols_p)]))
if (n_inc > 0) {
  msg("Descartando %d estaciones con algún valor ausente", n_inc)
  fondo <- fondo[complete.cases(
    st_drop_geometry(fondo)[, c(covs_base, cols_p)]), ]
}
msg("Estaciones usadas en todas las comparaciones: %d", nrow(fondo))

# =============================================================================
# 2. SELECCIÓN DE ESCALA POR VALIDACIÓN CRUZADA
# =============================================================================
sec("2. ¿A QUÉ ESCALA INFLUYE LA INDUSTRIA?")

f_base <- as.formula(paste("z ~", paste(covs_base, collapse = " + ")))

vg_de <- function(f) {
  v <- variogram(f, data = fondo, cutoff = mod$cutoff, width = mod$cutoff / 30)
  s <- var(residuals(lm(f, data = fondo)))
  fit.variogram(v, vgm(0.7 * s, mod$modelo, mod$cutoff / 6, 0.3 * s))
}

evaluar <- function(f, etiqueta) {
  vgf <- tryCatch(vg_de(f), error = function(e) {
    msg("  [%s] fallo al ajustar el variograma: %s", etiqueta,
        conditionMessage(e)); NULL })
  if (is.null(vgf)) return(NULL)
  cv <- tryCatch(
    krige.cv(f, fondo, model = vgf, nfold = nrow(fondo), verbose = FALSE),
    error = function(e) {
      msg("  [%s] fallo en la validación cruzada: %s", etiqueta,
          conditionMessage(e)); NULL })
  if (is.null(cv)) return(NULL)
  pred <- if (mod$usar_log) exp(cv$var1.pred + cv$var1.var / 2) else cv$var1.pred
  e <- pred - fondo$no2_media
  tibble(modelo = etiqueta,
         RMSE = sqrt(mean(e^2)),
         MAE = mean(abs(e)),
         R2 = 1 - sum(e^2) / sum((fondo$no2_media - mean(fondo$no2_media))^2),
         AIC_lm = AIC(lm(f, data = fondo)))
}

resultados <- evaluar(f_base, "Sin presión industrial")

for (nombre in names(presiones)) {
  f <- update(f_base, paste("~ . +", nombre))
  r <- evaluar(f, nombre)
  if (!is.null(r)) resultados <- bind_rows(resultados, r)
}

resultados <- resultados |>
  mutate(mejora_RMSE = 100 * (resultados$RMSE[1] - RMSE) / resultados$RMSE[1],
         delta_AIC = AIC_lm - min(AIC_lm))
print(as.data.frame(resultados))
write_csv(resultados, file.path(TAB, "50_barrido_escalas.csv"))

mejor_i <- which.min(resultados$RMSE)
if (mejor_i == 1) {
  msg("\nNINGUNA escala mejora el modelo base.")
  msg("Lectura: una vez controlado el uso del suelo y la señal satelital,")
  msg("la presión industrial no aporta información sobre el NO2 de fondo.")
  msg("Es un resultado, no un fallo: concuerda con el recargo casi nulo")
  msg("de las estaciones industriales medido en el script 20.")
  H_OPT <- NA
} else {
  H_OPT <- ANCHOS[mejor_i - 1]
  msg("\nEscala óptima: sigma = %.0f km", H_OPT / 1000)
  msg("Mejora del RMSE: %.2f %%", resultados$mejora_RMSE[mejor_i])
}

p <- ggplot(resultados[-1, ] |>
              mutate(h_km = ANCHOS[seq_len(n())] / 1000),
            aes(h_km, RMSE)) +
  geom_hline(yintercept = resultados$RMSE[1], linetype = 2, colour = "grey40") +
  geom_line(colour = "steelblue") + geom_point(size = 2) +
  scale_x_log10(breaks = ANCHOS / 1000) +
  labs(x = "ancho de banda del núcleo (km, escala log)",
       y = "RMSE de validación cruzada",
       title = "Escala espacial de influencia industrial",
       subtitle = "La línea discontinua es el modelo sin presión industrial") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG, "30_escala_influencia.png"), p,
       width = 7, height = 4.5, dpi = 150)

# =============================================================================
# 3. MODELO EN LA ESCALA ÓPTIMA
# =============================================================================
sec("3. EFECTO ESTIMADO")

nombre_opt <- if (is.na(H_OPT)) names(presiones)[1] else
  sprintf("presion_%dkm", H_OPT / 1000)
f_opt <- update(f_base, paste("~ . +", nombre_opt))

lm_opt <- lm(f_opt, data = fondo)
lm_base <- lm(f_base, data = fondo)
cat("\nCoeficientes del modelo con presión industrial:\n")
print(summary(lm_opt)$coefficients)

cat("\nContraste F frente al modelo sin presión industrial:\n")
print(anova(lm_base, lm_opt))

msg("\nR2 ajustado: base %.4f -> con presión %.4f",
    summary(lm_base)$adj.r.squared, summary(lm_opt)$adj.r.squared)

b <- coef(lm_opt)[nombre_opt]
if (!is.na(b) && mod$usar_log) {
  # En escala log, un aumento de una unidad en log-presión (es decir,
  # multiplicar la presión por e) cambia el NO2 en un factor exp(b).
  msg("Multiplicar por e la presión industrial cambia el NO2 en %+.2f %%",
      100 * (exp(b) - 1))
}

# =============================================================================
# 4. CONTRASTE CON LA ESCALA DEL PROCESO PUNTUAL
# =============================================================================
sec("4. DIÁLOGO ENTRE LAS DOS RAMAS")

ruta_kppm <- file.path(PROC, "modelo_kppm.rds")
if (file.exists(ruta_kppm) && !is.na(H_OPT)) {
  kp <- readRDS(ruta_kppm)
  escala_cluster <- tryCatch(kp$clustpar[["scale"]], error = function(e) NA)
  if (!is.na(escala_cluster)) {
    msg("Escala del conglomerado industrial (Thomas) : %.1f km",
        escala_cluster / 1000)
    msg("Diámetro efectivo del distrito              : %.1f km",
        4 * escala_cluster / 1000)
    msg("Escala de influencia sobre la concentración : %.1f km", H_OPT / 1000)
    razon <- H_OPT / (4 * escala_cluster)
    msg("Razón entre ambas: %.1f", razon)
    if (razon > 2) {
      msg("La influencia se extiende mucho más allá del distrito: domina el")
      msg("transporte atmosférico regional sobre la proximidad inmediata.")
    } else if (razon < 0.5) {
      msg("La influencia es más local que el propio distrito: el efecto se")
      msg("agota antes de abarcar el conjunto de instalaciones agrupadas.")
    } else {
      msg("Ambas escalas coinciden: la unidad de influencia sobre el NO2 es")
      msg("el distrito industrial, no la instalación aislada.")
    }
  }
}

# =============================================================================
# 5. SUPERFICIE DE PRESIÓN PARA LA REJILLA
# =============================================================================
# Se guarda la superficie en la escala elegida para que el script 60 pueda
# agregarla por provincia y confrontarla con el modelo areal.
r_opt <- presiones[[nombre_opt]]
writeRaster(r_opt, file.path(PROC, "presion_industrial.tif"), overwrite = TRUE)

png(file.path(FIG, "31_presion_industrial.png"), 1000, 800, res = 110)
plot(log(r_opt + 1e-14),
     main = sprintf("log presión de emisión industrial (sigma = %s km)",
                    ifelse(is.na(H_OPT), "—", H_OPT / 1000)))
points(st_coordinates(fondo), pch = 20, cex = 0.3, col = "white")
dev.off()

saveRDS(list(resultados = resultados, h_opt = H_OPT,
             formula = f_opt, modelo_lm = lm_opt,
             tiene_marcas = tiene_marcas),
        file.path(PROC, "acoplamiento.rds"))

sec("OBJETIVO 4 COMPLETADO")
if (!tiene_marcas) {
  cat("AVISO: sin marcas de masa el análisis trata igual a todas las\n")
  cat("instalaciones. Con eprtr_emisiones.csv el resultado sería más nítido.\n")
}
cat("Siguiente: R/60_areal.R (Objetivo 2)\n")
