# =============================================================================
# 40 — OBJETIVO 3: procesos puntuales de fuentes industriales
#
#   A. Descriptiva del patrón e intensidad de primer orden
#   B. Contraste de aleatoriedad espacial completa
#   C. Funciones K y g, homogéneas e inhomogéneas, con envolventes
#   D. Modelo de intensidad de Poisson inhomogéneo (ppm) con covariables
#   E. Diagnóstico de residuos y agregación no explicada
#   F. Modelo de conglomerados si queda estructura
#   G. Análisis marcado por masa de NOx
#
# ENTRADA : data/processed/pp_industrias.rds, ventana.rds, industrias.rds
#           data/interim/focal_norm/*.tif
# SALIDAS : output/figuras/2x_*.png, output/tablas/40_*.csv
#           data/processed/modelo_ppm.rds, intensidad.rds
# =============================================================================

suppressPackageStartupMessages({
  library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
  library(spatstat.random); library(sf); library(terra); library(dplyr)
  library(readr)
})

PROC <- file.path("data", "processed")
INT  <- file.path("data", "interim")
FIG  <- file.path("output", "figuras")
TAB  <- file.path("output", "tablas")

msg <- function(...) cat(sprintf(...), "\n")
sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n",
                       strrep("=", 68), "\n", sep = "")

pp      <- readRDS(file.path(PROC, "pp_industrias.rds"))

# El patrón llega MARCADO con la masa de NOx (así lo deja 10_prepare.R desde
# que las emisiones se descargan correctamente). Las marcas son
# imprescindibles para la superficie de presión del bloque H y del Objetivo
# 4, pero estorban en todo el análisis de intensidad: ppm() rechaza patrones
# marcados no multitipo, y las envolventes comparan un patrón marcado contra
# simulaciones sin marcas.
#
# La distinción es conceptual, no técnica: la intensidad de primer orden
# modela DÓNDE hay instalaciones, no cuánto emiten. Quitar las marcas ahí es
# lo correcto, no un rodeo.
pp_u <- if (is.marked(pp)) unmark(pp) else pp
ventana <- readRDS(file.path(PROC, "ventana.rds"))
nuts3   <- readRDS(file.path(PROC, "nuts3.rds"))

set.seed(20230101)
NSIM <- 199   # 199 simulaciones dan una prueba exacta al 5 % con envolventes

# =============================================================================
# A. DESCRIPTIVA
# =============================================================================
sec("A. DESCRIPTIVA DEL PATRÓN")

area_km2 <- area.owin(ventana) / 1e6
msg("Instalaciones: %d", npoints(pp_u))
msg("Ventana: %.0f km2", area_km2)
msg("Intensidad media: %.4f instalaciones por 1000 km2",
    1000 * npoints(pp_u) / area_km2)
msg("Distancia media al vecino más próximo: %.1f km",
    mean(nndist(pp_u)) / 1000)
# Bajo aleatoriedad espacial completa, la distancia esperada al vecino más
# próximo es 1/(2*sqrt(lambda)). Compararla con la observada da una primera
# señal: menor que la esperada indica agregación; mayor, regularidad.
lambda <- npoints(pp_u) / area.owin(ventana)
msg("Esperada bajo ACE: %.1f km", (1 / (2 * sqrt(lambda))) / 1000)

png(file.path(FIG, "20_patron_industrial.png"), 900, 800, res = 110)
plot(pp_u, main = "Instalaciones industriales con emisión de NOx",
     pch = 20, cex = 0.5, cols = "firebrick")
dev.off()

# =============================================================================
# B. INTENSIDAD DE PRIMER ORDEN
# =============================================================================
sec("B. INTENSIDAD DE PRIMER ORDEN")

# El ancho de banda se elige por validación cruzada en lugar de fijarlo a
# ojo. bw.diggle minimiza el error cuadrático medio de la estimación y es
# adecuado cuando interesa detectar agregación; bw.ppl (verosimilitud de
# punto) tiende a dar anchos mayores y superficies más suaves. Se calculan
# ambos y se reporta la diferencia, porque la elección afecta a todo lo que
# venga después.
bw_d <- bw.diggle(pp_u)
bw_p <- bw.ppl(pp_u)
msg("Ancho de banda: Diggle %.1f km | verosimilitud %.1f km",
    bw_d / 1000, bw_p / 1000)

# La discrepancia entre ambos criterios es informativa: si difieren mucho,
# el patrón tiene estructura a dos escalas distintas (agrupamiento local
# fino dentro de una heterogeneidad regional amplia). Conviene declararlo
# y comprobar que las conclusiones no dependen de la elección.
if (max(bw_d, bw_p) / min(bw_d, bw_p) > 3)
  msg("AVISO: los dos criterios difieren en un factor %.1f. El patrón tiene",
      max(bw_d, bw_p) / min(bw_d, bw_p))
  msg("estructura a dos escalas; se repetirá el análisis con ambos.")

BW <- bw_p   # el mayor: más estable para su uso posterior como covariable
lam <- density(pp_u, sigma = BW, positive = TRUE)

png(file.path(FIG, "21_intensidad.png"), 1000, 800, res = 110)
plot(lam, main = sprintf("Intensidad estimada (sigma = %.0f km)", BW / 1000))
plot(pp_u, add = TRUE, pch = 20, cex = 0.28, cols = "white")
dev.off()

# Cociente entre máximo y media: cuantifica lo concentrado del patrón
msg("Intensidad: media %.3g | máxima %.3g | cociente %.1f",
    mean(lam), max(lam), max(lam) / mean(lam))

# =============================================================================
# C. CONTRASTE DE ALEATORIEDAD ESPACIAL COMPLETA
# =============================================================================
sec("C. ¿ES EL PATRÓN ALEATORIO?")

# Con 296 puntos en 64 cuadrantes hay celdas con conteo esperado muy bajo,
# y la aproximación chi-cuadrado deja de ser fiable. La versión Monte Carlo
# no depende de esa aproximación.
qt <- quadrat.test(pp_u, nx = 8, ny = 8, method = "MonteCarlo", nsim = 999)
print(qt)

png(file.path(FIG, "22_cuadrantes.png"), 900, 800, res = 110)
plot(qt, main = "Contraste de cuadrantes (8 x 8)")
dev.off()

# =============================================================================
# D. FUNCIONES DE SEGUNDO ORDEN
# =============================================================================
sec("D. FUNCIONES K Y g")

# La K homogénea contrasta contra un Poisson de intensidad CONSTANTE. Para
# la industria ese contraste está garantizado a rechazar, y su rechazo no
# informa: nadie espera que las fábricas se distribuyan uniformemente por
# España. Se calcula solo como referencia.
#
# La versión INHOMOGÉNEA es la que responde la pregunta interesante:
# condicionando en la intensidad estimada, ¿queda agregación residual? Si la
# respuesta es sí, hay un mecanismo de atracción entre instalaciones
# (distritos industriales, infraestructura compartida) que va más allá de
# que unas zonas sean más propicias que otras.
RMAX <- 100000   # 100 km: escala de interés para distritos industriales

# ATENCIÓN AL DISEÑO DE LA SIMULACIÓN.
# envelope() simula por defecto desde un Poisson HOMOGÉNEO. Comparar contra
# esa referencia una K inhomogénea calculada condicionando en la intensidad
# real mezcla dos cosas distintas: la K observada queda deflactada por el
# condicionamiento mientras que las simuladas no, y el contraste pierde
# sentido (típicamente no detecta nada aunque haya agregación evidente).
#
# La referencia correcta es un Poisson INHOMOGÉNEO con la misma función de
# intensidad estimada. Así la hipótesis nula es "los puntos son
# independientes DADA la intensidad", que es justo lo que se quiere
# contrastar: si se rechaza, hay interacción entre instalaciones y no solo
# heterogeneidad del territorio.
msg("Calculando envolventes de K inhomogénea (%d simulaciones)…", NSIM)
env_kinh <- envelope(pp_u, Kinhom, lambda = lam, nsim = NSIM,
                     simulate = expression(rpoispp(lam)),
                     correction = "border", rmax = RMAX, verbose = FALSE)

png(file.path(FIG, "23_Kinhom.png"), 900, 700, res = 110)
plot(env_kinh, main = "K inhomogénea con envolventes",
     xlab = "distancia (m)", legend = FALSE)
dev.off()

msg("Calculando envolventes de g inhomogénea…")
env_ginh <- envelope(pp_u, pcfinhom, lambda = lam, nsim = NSIM,
                     simulate = expression(rpoispp(lam)),
                     rmax = RMAX, verbose = FALSE)

png(file.path(FIG, "24_ginhom.png"), 900, 700, res = 110)
plot(env_ginh, main = "Función de correlación de pares inhomogénea",
     xlab = "distancia (m)", legend = FALSE, ylim = c(0, 5))
dev.off()

# Escala de agregación: distancia hasta la que la K observada supera la
# envolvente superior. Es una estimación directa del alcance del mecanismo
# de atracción entre instalaciones.
fuera <- with(env_kinh, obs > hi)
if (any(fuera)) {
  r_max_agr <- max(env_kinh$r[fuera])
  msg("Agregación significativa hasta %.0f km", r_max_agr / 1000)
} else {
  msg("No se detecta agregación fuera de las envolventes.")
}

# =============================================================================
# E. MODELO DE INTENSIDAD CON COVARIABLES
# =============================================================================
sec("E. MODELO DE POISSON INHOMOGÉNEO")

# Convertir rásteres a imágenes de spatstat.
# terra recorre las filas de norte a sur y spatstat de sur a norte, de modo
# que hay que voltear la matriz o el mapa saldrá invertido.
raster_a_im <- function(ruta, res_m = 5000) {
  r <- rast(ruta)
  f <- max(1, round(res_m / res(r)[1]))
  if (f > 1) r <- terra::aggregate(r, fact = f, fun = "mean", na.rm = TRUE)
  m <- as.matrix(r, wide = TRUE)
  m <- m[nrow(m):1, ]
  e <- as.vector(ext(r))
  im(m,
     xcol = seq(e[1], e[2], length.out = ncol(m)),
     yrow = seq(e[3], e[4], length.out = nrow(m)))
}

covs <- list()
mapa_covs <- c(
  suelo_industrial = "frac_industrial_5000m.tif",
  suelo_artificial = "frac_artificial_5000m.tif",
  transporte       = "frac_transporte_5000m.tif",
  tierra           = "frac_tierra_10000m.tif"
)
for (nm in names(mapa_covs)) {
  ruta <- file.path(INT, "focal_norm", mapa_covs[nm])
  if (!file.exists(ruta)) ruta <- file.path(INT, "focal", mapa_covs[nm])
  if (file.exists(ruta)) {
    covs[[nm]] <- raster_a_im(ruta)
    msg("  covariable %s: rango %.3f - %.3f", nm,
        min(covs[[nm]], na.rm = TRUE), max(covs[[nm]], na.rm = TRUE))
  }
}

# Densidad poblacional por provincia, rasterizada.
# Es la covariable que separa "la industria está donde hay gente" de "la
# industria está donde hay suelo industrial": sin ella no se puede
# distinguir urbanización de industrialización, que es justo la pregunta
# del Objetivo 5.
if ("dens_2023" %in% names(nuts3)) {
  plantilla <- rast(ext(vect(nuts3)), resolution = 5000, crs = "EPSG:3035")
  rp <- rasterize(vect(nuts3), plantilla, field = "dens_2023")
  tmp <- tempfile(fileext = ".tif"); writeRaster(rp, tmp, overwrite = TRUE)
  covs$dens_pob <- raster_a_im(tmp)
  msg("  covariable dens_pob añadida")
}

if (length(covs) == 0) stop("Sin covariables: revisa data/interim/focal_norm/")

f_ppm <- as.formula(paste("~", paste(names(covs), collapse = " + ")))
msg("Modelo: %s", deparse(f_ppm))

mod_ppm <- ppm(pp_u, trend = f_ppm, covariates = covs)
print(summary(mod_ppm))

# exp(beta) es el factor por UNIDAD de covariable, pero las fracciones de
# uso del suelo varían entre 0 y 0,6 como mucho: un factor "por unidad"
# describe un cambio que nunca ocurre y sale en cifras absurdas (1e8).
# Se reporta además el efecto sobre el rango intercuartílico observado de
# cada covariable, que sí es un contraste realista.
rango_iqr <- sapply(names(covs), function(nm) {
  v <- covs[[nm]][pp_u]           # valor de la covariable en los puntos
  as.numeric(diff(quantile(v, c(0.25, 0.75), na.rm = TRUE)))
})

coefs <- tibble(
  termino = names(coef(mod_ppm)),
  estimado = as.numeric(coef(mod_ppm)),
  ee = sqrt(diag(vcov(mod_ppm)))
) |>
  mutate(z = estimado / ee,
         p = 2 * pnorm(-abs(z)),
         factor_unidad = exp(estimado),
         iqr = rango_iqr[termino],
         # Factor multiplicativo sobre la intensidad al pasar del primer al
         # tercer cuartil de la covariable: la magnitud interpretable
         factor_iqr = ifelse(is.na(iqr), NA, exp(estimado * iqr)))
print(as.data.frame(coefs))
write_csv(coefs, file.path(TAB, "40_coeficientes_ppm.csv"))

# =============================================================================
# F. DIAGNÓSTICO DE RESIDUOS
# =============================================================================
sec("F. RESIDUOS DEL MODELO")

png(file.path(FIG, "25_residuos_ppm.png"), 1000, 900, res = 110)
diagnose.ppm(mod_ppm, which = "smooth",
             main = "Residuos suavizados del modelo de intensidad")
dev.off()

# La pregunta clave: ¿queda agregación DESPUÉS de condicionar en las
# covariables? Si la K inhomogénea calculada con la intensidad AJUSTADA
# sigue por encima de la envolvente, existe un mecanismo de atracción entre
# instalaciones que las covariables observadas no explican. La lectura
# sustantiva es la existencia de distritos industriales: aglomeración por
# infraestructura compartida, cadenas de suministro o herencia histórica.
lam_ppm <- predict(mod_ppm, type = "trend")
msg("Envolventes de K inhomogénea sobre la intensidad ajustada…")
env_res <- envelope(pp_u, Kinhom, lambda = lam_ppm, nsim = NSIM,
                    correction = "border", rmax = RMAX,
                    simulate = expression(rpoispp(lam_ppm)), verbose = FALSE)

png(file.path(FIG, "26_Kinhom_residual.png"), 900, 700, res = 110)
plot(env_res, main = "K inhomogénea condicionada en el modelo ajustado",
     xlab = "distancia (m)", legend = FALSE)
dev.off()

queda <- with(env_res, obs > hi)
if (any(queda)) {
  msg("QUEDA agregación no explicada hasta %.0f km",
      max(env_res$r[queda]) / 1000)
  msg("Interpretación: distritos industriales no capturados por el uso del")
  msg("suelo ni la población. Procede un modelo de conglomerados.")
} else {
  msg("Las covariables explican la agregación: un Poisson inhomogéneo basta.")
}

# =============================================================================
# G. MODELO DE CONGLOMERADOS
# =============================================================================
if (any(queda)) {
  sec("G. PROCESO DE CONGLOMERADOS DE THOMAS")
  # kappa  : densidad de centros de conglomerado (distritos por unidad de área)
  # scale  : dispersión de las instalaciones alrededor de su centro
  # mu     : número esperado de instalaciones por conglomerado
  mod_th <- tryCatch(
    kppm(pp_u, trend = f_ppm, covariates = covs, clusters = "Thomas"),
    error = function(e) { msg("kppm falló: %s", conditionMessage(e)); NULL }
  )
  if (!is.null(mod_th)) {
    print(mod_th)
    # En kppm, kappa y scale viven en $clustpar, pero el tamaño medio de
    # conglomerado está en $mu (no dentro de clustpar), y su nombre cambia
    # entre versiones de spatstat. Se extrae con tolerancia.
    cp <- mod_th$clustpar
    escala <- if ("scale" %in% names(cp)) cp[["scale"]] else NA_real_
    kappa  <- if ("kappa" %in% names(cp)) cp[["kappa"]] else NA_real_
    mu <- tryCatch(as.numeric(mod_th$mu), error = function(e) NA_real_)

    if (!is.na(escala))
      msg("Escala del conglomerado (desv. típica de dispersión): %.1f km",
          escala / 1000)
    # El diámetro efectivo del distrito industrial ronda 4 veces la escala,
    # que es donde se concentra el 95 % de la descendencia gaussiana.
    if (!is.na(escala))
      msg("Diámetro efectivo del distrito: %.1f km", 4 * escala / 1000)
    if (!is.na(kappa))
      msg("Conglomerados esperados en España: %.0f",
          kappa * area.owin(ventana))
    if (length(mu) == 1 && !is.na(mu))
      msg("Instalaciones esperadas por conglomerado: %.2f", mu)
    saveRDS(mod_th, file.path(PROC, "modelo_kppm.rds"))
  }
}

# =============================================================================
# H. ANÁLISIS MARCADO
# =============================================================================
sec("H. PATRÓN MARCADO POR MASA DE NOx")

if (is.marked(pp)) {
  m <- marks(pp)
  nox <- if (is.data.frame(m)) m$nox_t else as.numeric(m)
  nox <- nox[!is.na(nox)]
  if (length(nox) > 10) {
    print(summary(nox))
    msg("Asimetría de la marca: %.2f",
        mean((nox - mean(nox))^3) / sd(nox)^3)
    msg("Las 10 mayores concentran el %.1f %% del NOx total",
        100 * sum(sort(nox, decreasing = TRUE)[1:10]) / sum(nox))
    msg("Se trabaja con log(masa): sin transformar, unas pocas centrales")
    msg("térmicas dominarían por completo la superficie ponderada.")

    ppm_marcado <- pp[!is.na(if (is.data.frame(m)) m$nox_t else m)]
    marks(ppm_marcado) <- log(nox)

    # Función de correlación de marcas: ¿los grandes emisores se agrupan
    # entre sí? Valores por encima de 1 a distancias cortas indicarían que
    # las instalaciones próximas tienden a tener emisiones similares.
    km <- markcorr(ppm_marcado, r = seq(0, RMAX, length.out = 100))
    png(file.path(FIG, "27_markcorr.png"), 900, 700, res = 110)
    plot(km, main = "Correlación de marcas: log(masa de NOx)",
         xlab = "distancia (m)")
    dev.off()

    # Intensidad ponderada por masa: es la superficie de PRESIÓN DE EMISIÓN
    # que el Objetivo 4 usará como covariable del campo de concentración.
    lam_w <- density(ppm_marcado, sigma = BW, weights = exp(marks(ppm_marcado)),
                     positive = TRUE)
    png(file.path(FIG, "28_presion_emision.png"), 1000, 800, res = 110)
    plot(lam_w, main = "Presión de emisión (intensidad ponderada por masa)")
    plot(pp_u, add = TRUE, pch = 20, cex = 0.25, cols = "white")
    dev.off()
    saveRDS(lam_w, file.path(PROC, "presion_emision.rds"))
    msg("-> presion_emision.rds, insumo del Objetivo 4")
  }
} else {
  msg("El patrón no tiene marcas: falta data/raw/eprtr_emisiones.csv.")
  msg("Sin ellas, el Objetivo 4 tratará por igual una térmica de 20 000 t")
  msg("y una cerámica de 110 t. Descarga el Excel de industry.eea.europa.eu")
  # Sin marcas, la intensidad simple sirve como versión degradada
  saveRDS(lam, file.path(PROC, "presion_emision.rds"))
}

saveRDS(mod_ppm, file.path(PROC, "modelo_ppm.rds"))
saveRDS(lam, file.path(PROC, "intensidad.rds"))

sec("OBJETIVO 3 COMPLETADO")
cat("Siguiente: R/50_acoplamiento.R (Objetivo 4)\n")
