# =============================================================================
# 47_nuts2_especificacion.R
# Recalcula el modelo areal a NUTS-2 con la MISMA especificación que NUTS-3,
# y produce la tabla de sensibilidad de especificación.
#
# La agregación se elige por VALIDACIÓN: NUTS-3 ya trae la respuesta que usó
# el pipeline (`no2_areal`), así que se prueba cuál de las tres formas de
# promediar la superficie la reproduce, y esa misma se aplica a NUTS-2.
#
# Genera: output/tablas/nuts2_especificacion.csv
# Uso:    source("R/47_nuts2_especificacion.R")
# =============================================================================

library(sf)

RUTA <- "data/processed"
leer <- function(n) {
  f <- file.path(RUTA, paste0(n, ".rds"))
  if (!file.exists(f)) stop("No existe ", f)
  readRDS(f)
}

nuts2 <- leer("nuts2")
nuts3 <- leer("nuts3_exposicion")
sup   <- leer("superficie_kde")

ID2 <- "NUTS2_ID"; ID3 <- "NUTS_ID"
cat("NUTS-2:", nrow(nuts2), "| NUTS-3:", nrow(nuts3), "\n\n")

# --- Media de bloque, tres convenciones -------------------------------------
col_sup <- intersect(c("pred", "var1.pred", "no2", "z"), names(sup))
if (!length(col_sup))
  stop("No identifico la columna predicha. Columnas: ",
       paste(names(sup), collapse = ", "))
col_sup <- col_sup[1]

rango <- range(sup[[col_sup]], na.rm = TRUE)
ES_LOG <- max(rango) < 6
cat("Superficie: columna", col_sup, "| rango", round(rango, 3),
    "| escala", if (ES_LOG) "logarítmica" else "cruda", "\n\n")

medias <- function(poligonos, id) {
  pol <- poligonos[, id, drop = FALSE]
  j <- st_join(st_transform(sup[, col_sup], st_crs(pol)), pol)
  v <- j[[col_sup]]; g <- j[[id]]
  crudo <- if (ES_LOG) exp(v) else v
  a <- tapply(v,     g, mean, na.rm = TRUE)   # media tal cual
  b <- tapply(crudo, g, mean, na.rm = TRUE)   # media en escala cruda
  data.frame(id = names(b),
             yA = as.numeric(a),              # A: directa
             yB = log(as.numeric(b)),         # B: cruda y volver a transformar
             yC = as.numeric(b),              # C: cruda sin transformar
             stringsAsFactors = FALSE)
}

# --- Validación: ¿cuál reproduce el `no2_areal` del pipeline? ----------------
m3 <- medias(nuts3, ID3)
k3 <- match(nuts3[[ID3]], m3$id)
ref <- as.numeric(nuts3$no2_areal)

cat("=== ¿Qué agregación reproduce `no2_areal`? ===\n")
comp <- do.call(rbind, lapply(c("yA", "yB", "yC"), function(v) {
  x <- m3[[v]][k3]
  data.frame(agregacion = v,
             correlacion = round(cor(x, ref, use = "complete.obs"), 4),
             dif_media   = round(mean(x - ref, na.rm = TRUE), 4),
             dif_max_abs = round(max(abs(x - ref), na.rm = TRUE), 4))
}))
print(comp, row.names = FALSE)

AGREGACION <- comp$agregacion[which.min(comp$dif_max_abs)]
cat("\nAgregación elegida:", AGREGACION,
    "(la de menor discrepancia máxima)\n\n")

# --- Referencia: el modelo del paper a NUTS-3 -------------------------------
d3 <- st_drop_geometry(nuts3)
d3$y        <- ref
d3$log_nox  <- log1p(as.numeric(d3$nox_total_t))
d3$log_area <- log(as.numeric(d3$area_km2))
d3$log_dens <- as.numeric(d3$log_dens_pob)

f3 <- lm(y ~ log_dens + log_nox + log_area, data = d3)
cat("=== NUTS-3, especificación del paper (n =", nrow(d3), ") ===\n")
print(round(summary(f3)$coefficients, 5))
cat("R2 ajustado:", round(summary(f3)$adj.r.squared, 4), "\n\n")

# --- El mismo modelo a NUTS-2 -----------------------------------------------
m2 <- medias(nuts2, ID2)
d  <- st_drop_geometry(nuts2)
d$y        <- m2[[AGREGACION]][match(d[[ID2]], m2$id)]
d$log_nox  <- log1p(as.numeric(d$nox_total_t))
d$log_area <- log(as.numeric(d$area_km2))
d$log_dens <- as.numeric(d$log_dens_pob)
d$log_conteo <- log1p(as.numeric(d$n_instalaciones))
d$madrid   <- substr(d[[ID2]], 1, 4) == "ES30"

base <- subset(d, !is.na(y))
cat("=== NUTS-2, misma especificación (n =", nrow(base), ") ===\n")
m <- lm(y ~ log_dens + log_nox + log_area, data = base)
print(round(summary(m)$coefficients, 5))
cat("R2 ajustado:", round(summary(m)$adj.r.squared, 4), "\n")

# --- Diagnóstico de influyentes ---------------------------------------------
cook   <- cooks.distance(m)
umbral <- 4 / nrow(base)
infl   <- sort(cook, decreasing = TRUE)[1:3]
cat("\nUmbral de Cook 4/n =", round(umbral, 4), "\n")
print(data.frame(codigo       = base[[ID2]][as.integer(names(infl))],
                 cook         = round(infl, 3),
                 veces_umbral = round(infl / umbral, 1)),
      row.names = FALSE)

sin_mad <- subset(base, !madrid)
m_sm <- lm(y ~ log_dens + log_nox + log_area, data = sin_mad)
cat("\n=== NUTS-2 sin Madrid ===\n")
print(round(summary(m_sm)$coefficients["log_nox", ], 5))
cat("R2 ajustado:", round(summary(m_sm)$adj.r.squared, 4), "\n")

# --- Sensibilidad de especificación -----------------------------------------
fila <- function(etiqueta, datos, formula, termino) {
  fit <- lm(formula, data = datos)
  cf  <- summary(fit)$coefficients
  if (!termino %in% rownames(cf))
    stop("Término '", termino, "' ausente en '", etiqueta, "'. Hay: ",
         paste(rownames(cf), collapse = ", "))
  data.frame(especificacion = etiqueta,
             n = nrow(datos), termino = termino,
             coef = round(cf[termino, 1], 4),
             p    = signif(cf[termino, 4], 3),
             r2_aj = round(summary(fit)$adj.r.squared, 3))
}

sens <- rbind(
  fila("NUTS-3, referencia", d3, y ~ log_dens + log_nox + log_area, "log_nox"),
  fila("NUTS-2, especificación del paper",
       base, y ~ log_dens + log_nox + log_area, "log_nox"),
  fila("NUTS-2, sin el término de área",
       base, y ~ log_dens + log_nox, "log_nox"),
  fila("NUTS-2, conteo de instalaciones (log)",
       base, y ~ log_dens + log_conteo + log_area, "log_conteo"),
  fila("NUTS-2, conteo sin transformar",
       base, y ~ log_dens + n_instalaciones + log_area, "n_instalaciones"),
  fila("NUTS-2, excluyendo Madrid",
       sin_mad, y ~ log_dens + log_nox + log_area, "log_nox")
)

cat("\n=== Sensibilidad de especificación ===\n")
print(sens, row.names = FALSE)

dir.create("output/tablas", recursive = TRUE, showWarnings = FALSE)
write.csv(sens, "output/tablas/nuts2_especificacion.csv", row.names = FALSE)
cat("\nTabla escrita en output/tablas/nuts2_especificacion.csv\n")
