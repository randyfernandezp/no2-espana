# =============================================================================
# 60 — OBJETIVO 2: variación espacial discreta (datos de área)
#
#   A. Variable respuesta y covariables por provincia
#   B. Matriz de vecindad y autocorrelación global (Moran)
#   C. Agrupamientos locales (LISA)
#   D. Regresión OLS y diagnóstico de dependencia residual
#   E. Modelos espaciales SAR, SEM y SDM
#   F. Descomposición de varianza (LMG y bloques temáticos)
#   G. Chequeo de MAUP: escala (NUTS-2) y zonificación (k-medias)
#
# ENTRADA : data/processed/nuts3_exposicion.rds, nuts2.rds,
#           presion_industrial.tif, instalaciones filtradas (E-PRTR)
# SALIDAS : output/figuras/4x_*.png, output/tablas/60_*.csv
#           data/processed/modelos_areales.rds
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(spdep); library(spatialreg); library(dplyr)
  library(readr); library(ggplot2); library(terra); library(tidyr)
})
if (!requireNamespace("car", quietly = TRUE))
  stop("Falta el paquete 'car' (VIF). install.packages('car')")

PROC <- file.path("data", "processed")
FIG  <- file.path("output", "figuras")
TAB  <- file.path("output", "tablas")
for (p in c(FIG, TAB)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf(...), "\n")
sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n",
                       strrep("=", 68), "\n", sep = "")

n3 <- readRDS(file.path(PROC, "nuts3_exposicion.rds"))
msg("Provincias: %d", nrow(n3))

# 10_prepare.R deja duplicadas las columnas de identificación en distinta
# caja (nuts_id / NUTS_ID). Se comprueba y se elimina la redundante.
for (par in list(c("nuts_id", "NUTS_ID"), c("nuts_name", "NUTS_NAME"))) {
  if (all(par %in% names(n3))) {
    if (identical(n3[[par[1]]], n3[[par[2]]])) {
      n3[[par[1]]] <- NULL
    } else {
      stop("Las columnas ", par[1], " y ", par[2], " difieren. ",
           "Revisar los joins de 10_prepare.R antes de continuar.")
    }
  }
}

# =============================================================================
# A. VARIABLE RESPUESTA Y COVARIABLES
# =============================================================================
sec("A. CONSTRUCCIÓN DEL CONJUNTO")

# --- A.1 Emisiones industriales por provincia ---------------------------------
#
# CORRECCIÓN. La columna nox_total_t que trae el RDS venía en cero para las
# 50 provincias: aguas arriba se agregaba por `InspireSiteId` (la clave de
# identificación) en lugar de por la marca de emisión. Aquí se reconstruye
# por vía GEOMÉTRICA, sin ningún join por atributo: cada instalación cae en
# el polígono que la contiene y se suman sus toneladas.
#
# El cruce de las cuatro tablas de DiscoData ya lo resolvió
# 04c_filtrar_industrias.py; no se duplica esa lógica aquí.

ruta_ied <- file.path(PROC, "industrias.rds")
if (!file.exists(ruta_ied))
  stop("No se encuentra ", ruta_ied, ". Ejecuta primero R/10_prepare.R.")
ied <- readRDS(ruta_ied)
msg("Instalaciones leídas de %s: %d", basename(ruta_ied), nrow(ied))

stopifnot("nox_t" %in% names(ied))
ied <- st_transform(ied, st_crs(n3))

nox_prov <- st_join(ied, n3["NUTS_ID"], join = st_within) |>
  st_drop_geometry() |>
  filter(!is.na(NUTS_ID)) |>
  group_by(NUTS_ID) |>
  summarise(nox_total_t = sum(nox_t, na.rm = TRUE), .groups = "drop")

n3$nox_total_t <- NULL
n3 <- left_join(n3, nox_prov, by = "NUTS_ID") |>
  mutate(nox_total_t = coalesce(nox_total_t, 0))

# Guardia imprescindible: un cero legítimo (provincia sin grandes emisores) y
# un cero por join fallido son indistinguibles a simple vista. Sin este
# control el fallo vuelve a pasar en silencio.
perdidas <- sum(ied$nox_t, na.rm = TRUE) - sum(n3$nox_total_t)
msg("NOx total: %.1f t en instalaciones | %.1f t asignado a provincias",
    sum(ied$nox_t, na.rm = TRUE), sum(n3$nox_total_t))
msg("Provincias con cero emisiones declaradas: %d", sum(n3$nox_total_t == 0))
if (abs(perdidas) > 1)
  warning("Se pierden ", round(perdidas, 1), " t en el cruce espacial: ",
          "hay instalaciones fuera de todo polígono. Revisar geometrías.")

# --- A.2 Presión industrial ---------------------------------------------------
# Escala de 10 km identificada como óptima en el script 50.
ruta_pres <- file.path(PROC, "presion_industrial.tif")
if (file.exists(ruta_pres)) {
  pres <- rast(ruta_pres)
  n3$presion_ind <- terra::extract(pres, vect(n3), fun = mean,
                                   na.rm = TRUE)[, 2]
  # log1p en lugar de log(x + 1e-14): el offset minúsculo mandaba a -32 las
  # provincias con presión nula y creaba valores atípicos artificiales.
  n3$log_presion <- log1p(n3$presion_ind)
  msg("Presión industrial agregada por provincia")
}

# --- A.3 Variable respuesta ---------------------------------------------------
#
# CORRECCIÓN. La respuesta es la MEDIA AREAL de bloque, no la exposición
# ponderada por población.
#
# La exposición ponderada se construye con la misma distribución de población
# que después entra como predictor: correlaciona 0,948 con log_dens_pob por
# construcción, así que regresarla sobre densidad es circular. Sigue siendo
# la métrica de política pública correcta y se conserva para el informe y
# como análisis de sensibilidad, pero no puede ser la variable dependiente
# de un modelo cuyo objetivo es atribuir causas.

d <- n3 |>
  mutate(
    y            = no2_areal,
    y_pob        = no2_poblacion,
    log_dens_pob = log(dens_2023),
    log_nox      = log1p(nox_total_t),
    nox_por_km2  = nox_total_t / area_km2,
    dens_ind     = 1000 * n_instalaciones / area_km2,
    frac_urbano  = frac_urbano_denso + frac_urbano_disperso,
    log_area     = log(area_km2)
  ) |>
  filter(!is.na(y))

msg("Provincias con respuesta válida: %d", nrow(d))
print(summary(d$y))

VARS <- c("log_dens_pob", "frac_artificial", "frac_urbano", "frac_industrial",
          "dens_ind", "log_nox", "log_presion", "log_area")
VARS <- intersect(VARS, names(d))

# Detectar variables sin varianza ANTES de correlacionar. Spearman sobre un
# vector constante devuelve NA, ese NA viaja a la matriz y rompe el `if`.
D <- st_drop_geometry(d)[, VARS, drop = FALSE]
const <- vapply(D, function(x) {
  x <- x[is.finite(x)]
  length(x) < 3L || length(unique(x)) == 1L
}, logical(1))
if (any(const)) {
  warning("Variables sin varianza, excluidas del análisis: ",
          paste(names(const)[const], collapse = ", "),
          ". Revisar su construcción aguas arriba.")
  VARS <- VARS[!const]
  D <- D[, VARS, drop = FALSE]
}

cat("\nCorrelación con la media areal:\n")
cm <- sapply(VARS, function(v) cor(D[[v]], d$y, use = "pairwise.complete.obs",
                                   method = "spearman"))
print(round(sort(cm, decreasing = TRUE), 3))

cat("\nCorrelación con la exposición poblacional (referencia):\n")
cmp <- sapply(VARS, function(v) cor(D[[v]], d$y_pob,
                                    use = "pairwise.complete.obs",
                                    method = "spearman"))
print(round(sort(cmp, decreasing = TRUE), 3))

# Colinealidad entre predictores: con 50 unidades, dos covariables muy
# correlacionadas hacen inestables los coeficientes y hacen imposible
# atribuir el efecto a una u otra, que es justo lo que pide el Objetivo 5.
cat("\nCorrelaciones entre predictores (|rho| > 0,7):\n")
M <- cor(D, use = "pairwise.complete.obs", method = "spearman")
for (i in seq_len(nrow(M) - 1)) for (j in (i + 1):ncol(M)) {
  # isTRUE() para que un eventual NA no rompa el bucle
  if (isTRUE(abs(M[i, j]) > 0.7))
    msg("  %-18s <-> %-18s rho = %+.3f", rownames(M)[i], colnames(M)[j], M[i, j])
}
write_csv(as.data.frame(round(M, 3)) |> mutate(variable = rownames(M), .before = 1),
          file.path(TAB, "60_correlaciones.csv"))

# El área es el confusor latente del bloque: las provincias pequeñas son
# densas, urbanas, costeras e industriales a la vez, así que toda variable
# INTENSIVA (por unidad de superficie) hereda ese eje común. log_nox es un
# total, no una densidad, y por eso es la única medida industrial que no
# resulta ser un alias de urbanización.
cat("\nCorrelación de cada predictor con log_area (confusor de escala):\n")
print(round(sort(M[, "log_area"]), 3))

# =============================================================================
# B. VECINDAD Y AUTOCORRELACIÓN GLOBAL
# =============================================================================
sec("B. ESTRUCTURA ESPACIAL")

# CORRECCIÓN. Baleares entra como TRES unidades NUTS-3 (Mallorca, Menorca,
# Eivissa-Formentera). Enlazarlas a sus dos vecinos más próximos por
# centroide las conecta entre sí y deja un triángulo aislado: el grafo sigue
# teniendo varias componentes y los contrastes dejan de ser fiables.
#
# Se usa un grafo de k = 5 vecinos más próximos como estructura principal.
# Conecta el archipiélago con el levante peninsular, garantiza una sola
# componente y no depende de la precisión de las fronteras. La contigüidad
# tipo reina se conserva como análisis de sensibilidad.

cent <- st_coordinates(st_point_on_surface(st_geometry(d)))
nb5 <- knn2nb(knearneigh(cent, k = 5))
lw <- nb2listw(nb5, style = "W")
msg("Grafo principal: k = 5 vecinos más próximos | componentes: %d",
    n.comp.nb(make.sym.nb(nb5))$nc)

nb_q <- poly2nb(d, queen = TRUE)
msg("Contigüidad reina (sensibilidad): %d componentes, %d unidades aisladas",
    n.comp.nb(nb_q)$nc, sum(card(nb_q) == 0))
if (sum(card(nb_q) == 0) > 0)
  msg("  Aisladas: %s", paste(d$NUTS_NAME[card(nb_q) == 0], collapse = ", "))
lw_q <- nb2listw(nb_q, style = "W", zero.policy = TRUE)

mi <- moran.test(d$y, lw)
print(mi)
msg("Moran I = %.4f (esperado bajo independencia %.4f), p = %.4g",
    mi$estimate[1], mi$estimate[2], mi$p.value)

# ADVERTENCIA para el informe: y es un agregado de bloque de un campo
# continuo con rango práctico de 71 km sobre polígonos de ~113 km de
# diámetro equivalente. Agregar un campo suave genera autocorrelación por
# construcción, aunque el proceso areal no la tenga. El Moran observado no
# debe leerse como evidencia de un mecanismo de contagio entre provincias.

png(file.path(FIG, "40_moran_scatter.png"), 800, 700, res = 110)
moran.plot(d$y, lw, labels = d$NUTS_ID,
           xlab = "media areal de NO2 (ug/m3)",
           ylab = "media de los vecinos")
dev.off()

# =============================================================================
# C. AGRUPAMIENTOS LOCALES
# =============================================================================
sec("C. INDICADORES LOCALES (LISA)")

lisa <- localmoran(d$y, lw)
d$lisa_I <- lisa[, 1]
d$lisa_p <- p.adjust(lisa[, 5], method = "BH")   # 50 contrastes simultáneos

z <- scale(d$y)[, 1]
zl <- lag.listw(lw, z)
d$cuadrante <- case_when(
  d$lisa_p > 0.05           ~ "no significativo",
  z > 0  & zl > 0           ~ "alto-alto",
  z < 0  & zl < 0           ~ "bajo-bajo",
  z > 0  & zl < 0           ~ "alto-bajo",
  TRUE                      ~ "bajo-alto"
)
print(table(d$cuadrante))

cat("\nProvincias en agrupamiento alto-alto:\n")
print(d |> st_drop_geometry() |> filter(cuadrante == "alto-alto") |>
        select(NUTS_ID, NUTS_NAME, y, lisa_p) |>
        mutate(across(where(is.numeric), ~ round(.x, 3))) |> as.data.frame())

pal <- c("alto-alto" = "#b2182b", "bajo-bajo" = "#2166ac",
         "alto-bajo" = "#ef8a62", "bajo-alto" = "#67a9cf",
         "no significativo" = "grey88")
p <- ggplot(d) +
  geom_sf(aes(fill = cuadrante), colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(title = "Agrupamientos locales de NO2 de fondo",
       subtitle = "LISA sobre la media areal de bloque (p ajustada por BH)") +
  theme_void(base_size = 11)
ggsave(file.path(FIG, "41_lisa.png"), p, width = 7, height = 6, dpi = 150)

# =============================================================================
# D. REGRESIÓN Y DIAGNÓSTICO
# =============================================================================
sec("D. MODELO OLS Y DEPENDENCIA RESIDUAL")

# CORRECCIÓN de la especificación. Los candidatos anteriores (frac_artificial
# y log_presion) son alias de urbanización: frac_artificial va con
# log_dens_pob a rho = 0,955 y log_presion con dens_ind a 0,943. Con n = 50
# eso no da un modelo interpretable.
#
# Se elige un representante por mecanismo:
#   log_dens_pob -> urbanización
#   log_nox      -> industria (total de emisiones; la única medida
#                   industrial que no correlaciona con el área)
#   log_area     -> confusor de escala
#
# Llevar log_nox y log_area juntos en logaritmos resuelve además el
# desajuste dimensional entre una respuesta intensiva (ug/m3) y un predictor
# extensivo (toneladas): b2*log(NOx) + b3*log(area) equivale a modelar
# log(NOx / area^theta) con theta = -b3/b2. No se impone si la industria
# actúa como total o como densidad; se deja que los datos estimen el
# exponente. (Solo interpretable si b2 es significativo.)
CAND <- c("log_dens_pob", "log_nox", "log_area")
CAND <- intersect(CAND, names(d))
f <- as.formula(paste("y ~", paste(CAND, collapse = " + ")))
msg("Modelo: %s", deparse(f))

ols <- lm(f, data = d)
print(summary(ols))

cat("\nFactores de inflación de varianza:\n")
print(round(car::vif(ols), 3))

# Diagnóstico de observaciones influyentes. Con n = 50 el umbral habitual es
# 4/n = 0,08.
ck <- cooks.distance(ols)
cat("\nDistancias de Cook más altas (umbral 4/n = ", round(4/nrow(d), 3), "):\n",
    sep = "")
print(round(sort(ck, decreasing = TRUE)[1:5], 3))
infl <- which(ck > 4 / nrow(d))
if (length(infl) > 0) {
  msg("Influyentes: %s", paste(d$NUTS_NAME[infl], collapse = ", "))
  ols_sin <- lm(f, data = d[-infl, ])
  comp_infl <- tibble(
    termino = names(coef(ols)),
    completo = as.numeric(coef(ols)),
    sin_influyentes = as.numeric(coef(ols_sin)[names(coef(ols))]),
    p_completo = summary(ols)$coefficients[, 4],
    p_sin = summary(ols_sin)$coefficients[names(coef(ols)), 4]
  )
  cat("\nSensibilidad al borrado de influyentes:\n")
  print(as.data.frame(comp_infl))
  write_csv(comp_infl, file.path(TAB, "60_influyentes.csv"))
}

# Moran sobre los residuos: si sigue habiendo autocorrelación, el OLS viola
# el supuesto de independencia y sus errores estándar están subestimados.
mr <- lm.morantest(ols, lw)
msg("\nMoran de los residuos: I = %.4f, p = %.4g", mr$estimate[1], mr$p.value)

# Contrastes del multiplicador de Lagrange. La regla de Anselin-Florax: si el
# robusto del retardo es significativo y el del error no, procede SAR; si es
# al revés, SEM. lm.LMtests está deprecada desde spdep 1.3; lm.RStests es la
# función vigente (mismos contrastes, nombres RSerr / RSlag).
lm_tests <- lm.RStests(ols, lw, test = "all")
print(summary(lm_tests))

cat("\nMismos contrastes con contigüidad reina (sensibilidad):\n")
print(summary(lm.RStests(ols, lw_q, test = "all", zero.policy = TRUE)))

# =============================================================================
# E. MODELOS ESPACIALES
# =============================================================================
sec("E. MODELOS AUTORREGRESIVOS")

sar <- lagsarlm(f, data = d, listw = lw)      # retardo en la respuesta
sem <- errorsarlm(f, data = d, listw = lw)    # retardo en el error
sdm <- lagsarlm(f, data = d, listw = lw, type = "mixed")  # Durbin

cat("\n--- SAR (retardo espacial) ---\n"); print(summary(sar))
cat("\n--- SEM (error espacial) ---\n");   print(summary(sem))

comp <- tibble(
  modelo = c("OLS", "SAR", "SEM", "SDM"),
  AIC = c(AIC(ols), AIC(sar), AIC(sem), AIC(sdm)),
  logLik = c(as.numeric(logLik(ols)), as.numeric(logLik(sar)),
             as.numeric(logLik(sem)), as.numeric(logLik(sdm))),
  parametro = c(NA, sar$rho, sem$lambda, sdm$rho)
) |> mutate(delta_AIC = AIC - min(AIC)) |> arrange(AIC)
print(as.data.frame(comp))
write_csv(comp, file.path(TAB, "60_comparacion_modelos.csv"))

mejor <- comp$modelo[1]
msg("Modelo preferido por AIC: %s", mejor)

# ¿Queda dependencia sin modelar en el mejor ajuste?
mejor_obj <- switch(mejor, OLS = ols, SAR = sar, SEM = sem, SDM = sdm)
mt <- moran.test(residuals(mejor_obj), lw)
msg("Moran de los residuos de %s: I = %.4f, p = %.4g",
    mejor, mt$estimate[1], mt$p.value)

# En modelos con retardo en la respuesta, los coeficientes NO son efectos
# marginales: un cambio en una provincia se propaga a sus vecinas y regresa.
# Los impactos directos, indirectos y totales sí son interpretables.
if (mejor %in% c("SAR", "SDM")) {
  obj <- if (mejor == "SAR") sar else sdm
  W <- as(lw, "CsparseMatrix")
  imp <- impacts(obj, tr = trW(W, type = "mult"), R = 500)
  cat("\nImpactos directos, indirectos y totales:\n")
  print(summary(imp, zstats = TRUE, short = TRUE))
}

# =============================================================================
# F. DESCOMPOSICIÓN DE VARIANZA
# =============================================================================
sec("F. ¿QUÉ EXPLICA EL NO2 REGIONAL?")

r2 <- function(vars, dat = d) {
  if (length(vars) == 0) return(0)
  summary(lm(as.formula(paste("y ~", paste(vars, collapse = " + "))),
             data = dat))$r.squared
}

# --- F.1 Reparto entre predictores (LMG) --------------------------------------
# El R2 secuencial depende del orden de entrada. LMG promedia sobre todas las
# ordenaciones posibles y da una partición exacta que suma el R2 del modelo.
cat("\nDescomposición LMG del modelo ajustado:\n")
perms <- function(v) {
  if (length(v) == 1) return(matrix(v))
  do.call(rbind, lapply(seq_along(v),
                        function(i) cbind(v[i], perms(v[-i]))))
}
P <- perms(seq_along(CAND))
lmg <- setNames(numeric(length(CAND)), CAND)
for (i in seq_len(nrow(P))) for (k in seq_along(CAND)) {
  j <- P[i, k]
  lmg[j] <- lmg[j] +
    (r2(CAND[P[i, seq_len(k)]]) - r2(CAND[P[i, seq_len(k - 1)]])) / nrow(P)
}
desc_lmg <- tibble(predictor = CAND, R2 = as.numeric(lmg),
                   pct = 100 * as.numeric(lmg) / sum(lmg)) |> arrange(desc(R2))
print(as.data.frame(desc_lmg |> mutate(across(where(is.numeric), ~ round(.x, 4)))))
write_csv(desc_lmg, file.path(TAB, "60_lmg.csv"))

# --- F.2 Descomposición del modelo espacial -----------------------------------
# y = Xb + u, con u = lambda*W*u + eps. El residuo se separa en campo
# espacial estructurado e innovación. Como los componentes no son
# ortogonales, se reparte por covarianzas: los términos suman var(y) de
# forma exactamente aditiva.
{
  Xb  <- as.numeric(model.matrix(ols) %*% sem$coefficients)
  u   <- d$y - Xb
  eps <- as.numeric(u - sem$lambda * lag.listw(lw, u))
  sp  <- u - eps
  comp_var <- c(fijos    = cov(d$y, Xb),
                espacial = cov(d$y, sp),
                ruido    = cov(d$y, eps)) / var(d$y)
  cat("\nDescomposición de varianza del SEM (proporción de var(y)):\n")
  print(round(rbind(proporcion = comp_var, pct = 100 * comp_var), 3))
  msg("Suma de componentes: %.6f (debe ser 1)", sum(comp_var))
  write_csv(tibble(componente = names(comp_var), proporcion = comp_var),
            file.path(TAB, "60_descomposicion_sem.csv"))
}

# --- F.3 Bloques temáticos ----------------------------------------------------
# La pregunta del Objetivo 5 no es qué covariable es significativa, sino
# cuánta varianza aporta cada MECANISMO. Como los bloques están
# correlacionados entre sí, se separa la contribución única (lo que un bloque
# explica y ningún otro puede) de la compartida (lo que varios explican
# indistintamente). Atribuir la parte compartida a un solo bloque sería
# arbitrario, y es el error habitual al interpretar coeficientes.
BLOQUES <- list(
  urbanizacion = intersect(c("log_dens_pob", "frac_urbano"), names(d)),
  uso_suelo    = intersect(c("frac_artificial", "frac_industrial"), names(d)),
  industria    = intersect(c("log_presion", "dens_ind", "log_nox"), names(d)),
  geografia    = intersect(c("log_area"), names(d))
)
BLOQUES <- BLOQUES[lengths(BLOQUES) > 0]

todos <- unlist(BLOQUES, use.names = FALSE)
R2_total <- r2(todos)
msg("\nR2 con todos los bloques: %.4f", R2_total)

desc <- tibble(
  bloque = names(BLOQUES),
  R2_solo = sapply(BLOQUES, r2),
  R2_unico = sapply(names(BLOQUES), function(b)
    R2_total - r2(setdiff(todos, BLOQUES[[b]])))
) |>
  mutate(R2_compartido = R2_solo - R2_unico,
         pct_unico = 100 * R2_unico / R2_total) |>
  arrange(desc(R2_unico))
print(as.data.frame(desc))
write_csv(desc, file.path(TAB, "60_descomposicion_varianza.csv"))

msg("\nVarianza única total: %.4f de %.4f (%.1f %%)",
    sum(desc$R2_unico), R2_total, 100 * sum(desc$R2_unico) / R2_total)
msg("El resto es varianza compartida: no atribuible a un solo mecanismo.")

# --- F.4 Techo del efecto industrial ------------------------------------------
# Si log_nox no es significativa, el resultado informativo no es el p-valor
# sino la COTA: cuánto NO2 podría como máximo explicar la industria sin ser
# incompatible con los datos. Un techo cuantificado es una afirmación mucho
# más fuerte que una ausencia de significancia.
if ("log_nox" %in% CAND) {
  ic <- confint(ols)["log_nox", ]
  rango <- diff(range(d$log_nox))
  techo <- max(abs(ic)) * rango
  cat("\nCota del efecto industrial:\n")
  msg("  IC 95%% de log_nox: [%.4f, %.4f]", ic[1], ic[2])
  msg("  Efecto máximo compatible al recorrer todo el rango observado de")
  msg("  emisiones: %.3f ug/m3 (%.1f %% de la media areal nacional)",
      techo, 100 * techo / mean(d$y))
}

p <- ggplot(desc, aes(reorder(bloque, R2_unico))) +
  geom_col(aes(y = R2_unico, fill = "única"), width = 0.6) +
  geom_col(aes(y = R2_solo, fill = "total (única + compartida)"),
           width = 0.6, alpha = 0.28) +
  coord_flip() +
  scale_fill_manual(values = c("única" = "#2166ac",
                               "total (única + compartida)" = "#92c5de"),
                    name = NULL) +
  labs(x = NULL, y = expression(R^2),
       title = "Aportación de cada mecanismo al NO2 de fondo regional") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG, "42_descomposicion.png"), p, width = 7.5, height = 4.5,
       dpi = 150)

# =============================================================================
# G. PROBLEMA DE LA UNIDAD ESPACIAL MODIFICABLE
# =============================================================================
sec("G. CHEQUEO DE MAUP")

# El MAUP tiene dos componentes y hay que documentar los dos:
#   ESCALA       — mismo criterio de trazado, distinto número de unidades
#   ZONIFICACIÓN — mismo número de unidades, distinto trazado
# Reportar solo el primero es la omisión habitual.

# CORRECCIÓN de la agregación. La versión anterior promediaba las variables
# YA TRANSFORMADAS: la media ponderada de log_area entre provincias no es el
# logaritmo del área total, ni la de log_nox el logaritmo de las emisiones
# sumadas. Hay que agregar en la escala CRUDA y volver a transformar.
#
# La respuesta se agrega ponderando por ÁREA (no por población): la media
# areal de la unión de polígonos es exactamente la media de las partes
# ponderada por superficie, sin aproximación.

agregar <- function(dat, g) {
  tibble(g = g,
         y           = dat$y,
         poblacion   = dat$poblacion,
         area_km2    = dat$area_km2,
         nox_total_t = dat$nox_total_t) |>
    group_by(g) |>
    summarise(y           = weighted.mean(y, area_km2, na.rm = TRUE),
              poblacion   = sum(poblacion, na.rm = TRUE),
              area_km2    = sum(area_km2, na.rm = TRUE),
              nox_total_t = sum(nox_total_t, na.rm = TRUE),
              .groups = "drop") |>
    mutate(log_dens_pob = log(poblacion / area_km2),
           log_nox      = log1p(nox_total_t),
           log_area     = log(area_km2))
}

dd <- st_drop_geometry(d)

# --- G.1 Escala: NUTS-2 -------------------------------------------------------
d2 <- agregar(dd, substr(dd$NUTS_ID, 1, 4))
msg("Regiones NUTS-2 con datos: %d", nrow(d2))

ols2 <- lm(f, data = d2)
print(summary(ols2))
cat("\nVIF a NUTS-2:\n"); print(round(car::vif(ols2), 3))

ck2 <- cooks.distance(ols2)
cat("\nDistancias de Cook a NUTS-2 (umbral 4/n = ", round(4/nrow(d2), 3), "):\n",
    sep = "")
print(round(sort(ck2, decreasing = TRUE)[1:4], 3))
infl2 <- which(ck2 > 4 / nrow(d2))
if (length(infl2) > 0) {
  msg("Influyentes a NUTS-2: %s", paste(d2$g[infl2], collapse = ", "))
  ols2b <- lm(f, data = d2[-infl2, ])
  cat("\nAjuste a NUTS-2 sin las unidades influyentes:\n")
  print(summary(ols2b))
  msg("R2 ajustado: %.3f con influyentes -> %.3f sin ellas",
      summary(ols2)$adj.r.squared, summary(ols2b)$adj.r.squared)
  msg("Una caída fuerte indica que el ajuste aparente de la escala gruesa")
  msg("se apoyaba en pocas unidades, no en información adicional.")
}

# Coeficientes estandarizados: los b crudos NO son comparables entre escalas
# porque la dispersión de los predictores cambia al agregar.
zb <- function(m, dat) {
  v <- all.vars(formula(m))
  s <- lm(reformulate(v[-1], v[1]),
          data = as.data.frame(scale(as.data.frame(dat)[, v])))
  coef(s)[-1]
}
comp_esc <- tibble(
  termino = CAND,
  beta_z_NUTS3 = as.numeric(zb(ols, dd)[CAND]),
  beta_z_NUTS2 = as.numeric(zb(ols2, d2)[CAND])
) |> mutate(cambio_signo = sign(beta_z_NUTS3) != sign(beta_z_NUTS2))
cat("\nCoeficientes estandarizados por escala:\n")
print(as.data.frame(comp_esc))
write_csv(comp_esc, file.path(TAB, "60_maup_escala.csv"))

msg("\nR2: NUTS-3 %.3f -> NUTS-2 %.3f",
    summary(ols)$r.squared, summary(ols2)$r.squared)

# --- G.2 Zonificación: particiones alternativas -------------------------------
# Se fija el número de unidades en el de NUTS-2 y se varía SOLO el trazado,
# agrupando las provincias por k-medias sobre centroides. Cada réplica con
# nstart = 1 converge a un óptimo local distinto, lo que da un muestreo de
# zonificaciones plausibles a escala constante.
set.seed(1)
K <- nrow(d2)
R <- 1000

una_rep <- function(i) {
  g <- tryCatch(kmeans(cent, centers = K, nstart = 1, iter.max = 50)$cluster,
                error = function(e) NULL)
  if (is.null(g)) return(rep(NA_real_, 2 * length(CAND) + 1))
  m <- lm(f, data = agregar(dd, g))
  s <- summary(m)$coefficients
  c(setNames(s[CAND, 1], paste0("b_", CAND)),
    setNames(s[CAND, 4], paste0("p_", CAND)),
    r2 = summary(m)$adj.r.squared)
}
zon <- as.data.frame(t(vapply(seq_len(R), una_rep,
                              numeric(2 * length(CAND) + 1))))
zon <- zon[complete.cases(zon), ]
msg("Zonificaciones evaluadas: %d de %d", nrow(zon), R)

qs <- t(apply(zon[, paste0("b_", CAND), drop = FALSE], 2, quantile,
              c(.025, .25, .5, .75, .975)))
cat("\nDistribución de los coeficientes sobre zonificaciones alternativas:\n")
print(round(qs, 3))

tasa <- sapply(CAND, function(v) mean(zon[[paste0("p_", v)]] < 0.05))
cat("\nFracción de particiones en que cada predictor sale significativo:\n")
print(round(tasa, 3))
msg("Referencia: con efecto nulo se espera 0,05. Muy por encima indicaría")
msg("que la zonificación fabrica significancia; muy por debajo, que la")
msg("ausencia de efecto es estable frente al redibujado de fronteras.")

# ¿Cae el coeficiente observado en NUTS-2 dentro del rango de lo que produce
# una partición cualquiera del mismo tamaño?
zon_tab <- tibble(
  termino = CAND,
  beta_NUTS2 = as.numeric(coef(ols2)[CAND]),
  q025 = qs[, 1], mediana = qs[, 3], q975 = qs[, 5],
  pct_significativo = 100 * as.numeric(tasa[CAND])
) |> mutate(fuera_del_rango = beta_NUTS2 < q025 | beta_NUTS2 > q975)
cat("\nCoeficiente de NUTS-2 frente a la distribución por zonificación:\n")
print(as.data.frame(zon_tab))
write_csv(zon_tab, file.path(TAB, "60_maup_zonificacion.csv"))

if (any(zon_tab$fuera_del_rango)) {
  msg("\nATENCIÓN: %s queda fuera del intervalo del 95 %% generado por",
      paste(zon_tab$termino[zon_tab$fuera_del_rango], collapse = ", "))
  msg("particiones alternativas del mismo tamaño. El valor estimado sobre la")
  msg("división administrativa real es atípico y debe reportarse como no")
  msg("robusto, con el diagnóstico de influyentes como evidencia.")
}

p <- zon |>
  select(all_of(paste0("b_", CAND))) |>
  pivot_longer(everything(), names_to = "termino", values_to = "beta") |>
  mutate(termino = sub("^b_", "", termino)) |>
  ggplot(aes(beta, termino)) +
  geom_violin(fill = "grey85", colour = NA) +
  geom_point(data = zon_tab, aes(beta_NUTS2, termino),
             colour = "#b2182b", size = 2.5) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  labs(x = "coeficiente", y = NULL,
       title = "Efecto de zonificación sobre los coeficientes",
       subtitle = paste0(nrow(zon), " particiones de ", K,
                         " grupos; en rojo, la división NUTS-2 real")) +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG, "43_maup_zonificacion.png"), p, width = 7.5, height = 4,
       dpi = 150)

# =============================================================================
saveRDS(list(datos = d, lw = lw, lw_contig = lw_q, formula = f,
             ols = ols, sar = sar, sem = sem, sdm = sdm,
             comparacion = comp, lmg = desc_lmg, descomposicion = desc,
             maup_escala = comp_esc, maup_zonificacion = zon_tab,
             zonificaciones = zon),
        file.path(PROC, "modelos_areales.rds"))

sec("OBJETIVO 2 COMPLETADO")
cat("Siguiente: R/70_sintesis.R (Objetivo 5)\n")
