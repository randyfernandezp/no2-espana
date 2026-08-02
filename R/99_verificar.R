# =============================================================================
# 99 — Verificación de integridad del proyecto
#
# No comprueba que el análisis sea correcto: comprueba que los datos que
# alimentan el análisis son los que se cree que son. Está escrito contra la
# clase de fallo que ya se produjo una vez — ceros silenciosos y ficheros
# derivados obsoletos — no contra errores de sintaxis.
#
# Uso:  source("R/99_verificar.R")
# Sale con un recuento de fallos. Cero fallos no demuestra que todo esté
# bien; un fallo sí demuestra que algo está mal.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(spatstat.geom)
})

PROC <- file.path("data", "processed")
RAW  <- file.path("data", "raw")
INT  <- file.path("data", "interim")

FALLOS <- 0L
AVISOS <- 0L

chk <- function(cond, texto, detalle = "") {
  ok <- isTRUE(cond)
  if (!ok) FALLOS <<- FALLOS + 1L
  cat(sprintf("[%s] %s%s\n", if (ok) " OK " else "FALLO", texto,
              if (!ok && nzchar(detalle)) paste0("\n         ", detalle) else ""))
  invisible(ok)
}
avi <- function(cond, texto) {
  if (!isTRUE(cond)) {
    AVISOS <<- AVISOS + 1L
    cat(sprintf("[aviso] %s\n", texto))
  }
  invisible(TRUE)
}
sec <- function(t) cat("\n", strrep("-", 68), "\n", t, "\n",
                       strrep("-", 68), "\n", sep = "")

# =============================================================================
# 1. FRESCURA: ningún derivado más antiguo que sus fuentes
# =============================================================================
# Este es el bloque que habría detectado el fallo original. nuts3.rds se
# generó ANTES de que existiera eprtr_emisiones.csv, así que arrastraba
# nox_total_t = 0, y nuts3_exposicion.rds heredó el defecto. Un fichero
# derivado más viejo que su fuente es sospechoso siempre.
sec("1. FRESCURA DE LOS FICHEROS DERIVADOS")

edad <- function(p) if (file.exists(p)) file.mtime(p) else NA

CADENA <- list(
  list(hijo = file.path(PROC, "nuts3.rds"),
       padres = c(file.path(RAW, "eprtr_emisiones.csv"),
                  file.path(RAW, "nuts3_es.geojson"),
                  file.path(RAW, "densidad_nuts3_es.csv"),
                  file.path(INT, "clc_fracciones_nuts3.csv"))),
  list(hijo = file.path(PROC, "industrias.rds"),
       padres = c(file.path(RAW, "eprtr_emisiones.csv"),
                  file.path(INT, "ied_es_nox_ventana.geojson"))),
  list(hijo = file.path(PROC, "pp_industrias.rds"),
       padres = file.path(PROC, "industrias.rds")),
  list(hijo = file.path(PROC, "nuts3_exposicion.rds"),
       padres = file.path(PROC, "nuts3.rds")),
  list(hijo = file.path(PROC, "presion_industrial.tif"),
       padres = file.path(PROC, "pp_industrias.rds")),
  list(hijo = file.path(PROC, "modelos_areales.rds"),
       padres = c(file.path(PROC, "nuts3_exposicion.rds"),
                  file.path(PROC, "presion_industrial.tif")))
)

for (e in CADENA) {
  th <- edad(e$hijo)
  if (is.na(th)) { avi(FALSE, paste("No existe:", e$hijo)); next }
  for (p in e$padres) {
    tp <- edad(p)
    if (is.na(tp)) { avi(FALSE, paste("Fuente ausente:", p)); next }
    chk(th >= tp,
        sprintf("%s más reciente que %s", basename(e$hijo), basename(p)),
        sprintf("derivado %s | fuente %s -> REGENERAR", th, tp))
  }
}

# =============================================================================
# 2. INVARIANTES DE CONSERVACIÓN
# =============================================================================
# Cantidades que deben cuadrar por construcción. Si no cuadran, hay un cruce
# roto en alguna parte, aunque todo se haya ejecutado sin errores.
sec("2. INVARIANTES")

nuts3 <- readRDS(file.path(PROC, "nuts3.rds"))
ied   <- readRDS(file.path(PROC, "industrias.rds"))
pp    <- readRDS(file.path(PROC, "pp_industrias.rds"))
vent  <- readRDS(file.path(PROC, "ventana.rds"))

chk(!all(is.na(nuts3$nox_total_t)),
    "nox_total_t no es NA en todas las provincias")
chk(!all(nuts3$nox_total_t == 0, na.rm = TRUE),
    "nox_total_t no es cero en todas las provincias",
    "Este es EXACTAMENTE el fallo original. Revisar 10_prepare.R")

if (!all(is.na(ied$nox_t))) {
  dif <- abs(sum(ied$nox_t, na.rm = TRUE) - sum(nuts3$nox_total_t, na.rm = TRUE))
  chk(dif < 1,
      sprintf("NOx conservado en el cruce (dif = %.2f t)", dif),
      "Hay instalaciones fuera de todo polígono o asignadas dos veces")
}

chk(npoints(pp) == nrow(ied),
    sprintf("ppp e ied coinciden (%d vs %d)", npoints(pp), nrow(ied)),
    "La deduplicación se aplicó a uno y no al otro")

chk(sum(nuts3$n_instalaciones) == nrow(ied),
    sprintf("Instalaciones asignadas: %d de %d",
            sum(nuts3$n_instalaciones), nrow(ied)))

area_dom <- area(vent) / 1e6
dif_area <- abs(sum(nuts3$area_km2) - area_dom) / area_dom
chk(dif_area < 0.02,
    sprintf("Área: suma provincial %.0f km2 vs ventana %.0f km2 (%.1f %%)",
            sum(nuts3$area_km2), area_dom, 100 * dif_area),
    "La ventana se simplifica a 1 km; una desviación mayor del 2 % es anómala")

chk(all(nuts3$poblacion > 0),
    "Todas las provincias tienen población positiva")
chk(!any(duplicated(nuts3$NUTS_ID)),
    "No hay NUTS_ID duplicados")
chk(sum(grepl("^nuts_id$|^nuts_name$", names(nuts3))) == 0,
    "Sin columnas de identificación duplicadas en minúscula")

# =============================================================================
# 3. RANGOS PLAUSIBLES
# =============================================================================
# Un valor dentro de rango no prueba nada; uno fuera de rango prueba que hay
# un fallo. Los límites son deliberadamente generosos.
sec("3. RANGOS")

est <- readRDS(file.path(PROC, "estaciones.rds"))

chk(all(est$no2_media > 0 & est$no2_media < 120, na.rm = TRUE),
    "NO2 en estaciones dentro de 0-120 ug/m3")
chk(min(ied$nox_t, na.rm = TRUE) >= 100,
    sprintf("Marca mínima de NOx = %.1f t (umbral E-PRTR: 100)",
            min(ied$nox_t, na.rm = TRUE)),
    "Por debajo de 100 t no debería haber declaraciones")
chk(st_crs(nuts3)$epsg == 3035 && st_crs(ied)$epsg == 3035,
    "CRS 3035 en geometrías")
chk(!any(grepl("^ES70|^ES63|^ES64", nuts3$NUTS_ID)),
    "Canarias, Ceuta y Melilla excluidas del dominio")

# Provincias con cero emisiones: unas pocas es normal, cincuenta es el bug.
n0 <- sum(nuts3$nox_total_t == 0, na.rm = TRUE)
chk(n0 < 10, sprintf("Provincias sin emisores declarados: %d", n0))

# El ranking industrial debe ser reconocible. Si no lo es, el cruce falló
# aunque los totales cuadren.
top <- nuts3 |> st_drop_geometry() |>
  arrange(desc(nox_total_t)) |> head(8)
cat("\nTop 8 en NOx industrial (debe incluir Asturias, A Coruña, Tarragona,\n",
    "Cádiz, Huelva — los focos documentados):\n", sep = "")
print(data.frame(provincia = top$NUTS_NAME, t = round(top$nox_total_t)))
esperados <- c("Asturias", "Coru", "Tarragona", "Huelva")
n_esp <- sum(sapply(esperados, function(p) any(grepl(p, top$NUTS_NAME))))
chk(n_esp >= 3,
    sprintf("Focos documentados en el top 8: %d de %d", n_esp, length(esperados)),
    "Si el ranking es irreconocible, el join asignó mal")

# =============================================================================
# 4. REGRESIÓN CONTRA LOS RESULTADOS DOCUMENTADOS
# =============================================================================
# Los valores esperados son los VIGENTES, posteriores a la deduplicación por
# InspireSiteId, y son los que cita el artículo. Se leen de los objetos
# ajustados, no se transcriben a mano: si un modelo se reajusta y cambia, este
# bloque lo detecta.
#
# Si algún extractor no reconoce la estructura del objeto, avisa y sigue en
# vez de romper: un verificador que no se puede ejecutar no verifica nada.
sec("4. COMPARACIÓN CON LOS VALORES DOCUMENTADOS")

`%||%` <- function(a, b) if (is.null(a)) b else a

leer_si <- function(n) {
  f <- file.path(PROC, paste0(n, ".rds"))
  if (file.exists(f)) readRDS(f) else NULL
}

# --- sigma del proceso agregado -> diámetro del distrito = 4 sigma ----------
extraer_sigma_km <- function(k) {
  if (is.null(k)) return(NA_real_)
  for (via in list(function(x) as.numeric(x$clustpar[["scale"]]),
                   function(x) as.numeric(x$modelpar[["sigma"]]),
                   function(x) as.numeric(x$par[["scale"]]))) {
    v <- try(via(k), silent = TRUE)
    if (!inherits(v, "try-error") && length(v) == 1 && is.finite(v))
      return(if (v > 1000) v / 1000 else v)   # metros -> km
  }
  NA_real_
}

# --- rango práctico del variograma: rho(r*) = 0,05 --------------------------
extraer_rango_km <- function(g, cual = "vg_fit") {
  if (is.null(g)) return(NA_real_)
  if (!is.data.frame(g) && !is.null(g[[cual]]) && is.data.frame(g[[cual]]))
    g <- g[[cual]]
  vm <- g
  if (!is.data.frame(vm))
    for (nm in c("vg_fit", "vg_res_fit", "variograma", "vgm", "fit"))
      if (!is.null(g[[nm]]) && is.data.frame(g[[nm]])) { vm <- g[[nm]]; break }
  if (!is.data.frame(vm) || !all(c("model", "range") %in% names(vm)))
    return(NA_real_)
  i <- which(as.character(vm$model) != "Nug")
  if (!length(i)) return(NA_real_)
  i <- i[which.max(vm$psill[i])]
  r <- as.numeric(vm$range[i])
  mult <- switch(as.character(vm$model[i]),
                 Exp = 3, Gau = sqrt(3), Sph = 1, Mat = 4, 3)
  r <- r * mult
  if (r > 1000) r / 1000 else r
}

# --- escala óptima del barrido de anchos de banda ---------------------------
extraer_escala_km <- function(a) {
  if (is.null(a)) return(NA_real_)
  if (!is.null(a$h_opt)) {
    v <- as.numeric(a$h_opt)[1]
    if (is.finite(v)) return(if (v > 1000) v / 1000 else v)
  }
  df <- if (is.data.frame(a)) a else
        Filter(is.data.frame, a)[[1]] %||% NULL
  if (is.null(df) || !is.data.frame(df)) return(NA_real_)
  ch <- grep("^(h|banda|bw|ancho)", names(df), ignore.case = TRUE, value = TRUE)
  cr <- grep("rmse", names(df), ignore.case = TRUE, value = TRUE)
  if (!length(ch) || !length(cr)) return(NA_real_)
  h <- as.numeric(df[[ch[1]]]); r <- as.numeric(df[[cr[1]]])
  ok <- is.finite(h) & is.finite(r)
  if (!any(ok)) return(NA_real_)
  h[ok][which.min(r[ok])]
}
geo       <- leer_si("modelo_geo")
sigma_km  <- extraer_sigma_km(leer_si("modelo_kppm"))
rango_km  <- extraer_rango_km(geo, "vg_fit")
rango_res <- extraer_rango_km(geo, "vg_res_fit")
escala_km <- extraer_escala_km(leer_si("acoplamiento"))

REF <- tribble(
  ~magnitud,                             ~esperado, ~tol, ~observado,
  "instalaciones en el patrón",              297,     5,   as.numeric(npoints(pp)),
  "rango práctico del variograma (km)",       71,     8,   rango_km,
  "escala óptima de acoplamiento (km)",       10,     3,   escala_km,
  "diámetro del distrito, 4 sigma (km)",    10.7,     3,   4 * sigma_km,
  "rango residual tras la deriva (km)",       24,      6,   rango_res
)
REF$desvio <- round(REF$observado - REF$esperado, 2)
print(as.data.frame(REF))

for (i in seq_len(nrow(REF))) {
  if (is.na(REF$observado[i])) {
    avi(FALSE, paste0("No pude leer '", REF$magnitud[i],
                      "' del objeto ajustado; revisar el extractor."))
    next
  }
  chk(abs(REF$observado[i] - REF$esperado[i]) <= REF$tol[i],
      sprintf("%s: %.2f (esperado %.2f +/- %.2f)", REF$magnitud[i],
              REF$observado[i], REF$esperado[i], REF$tol[i]),
      "El valor documentado y el ajustado se han separado: actualizar uno u otro")
}

# --- El hallazgo central: la razón entre las dos estimaciones ---------------
if (all(is.finite(c(escala_km, sigma_km)))) {
  razon <- escala_km / (4 * sigma_km)
  cat(sprintf("\nRazón escala/diámetro = %.3f\n", razon))
  chk(abs(razon - 1) < 0.25,
      sprintf("Convergencia de escalas dentro del 25 %% (razón %.3f)", razon),
      "La convergencia es el hallazgo central: si se aleja, hay que decirlo")
} else {
  avi(FALSE, "No se pudo calcular la razón escala/diámetro.")
}

# =============================================================================
sec("RESUMEN")
msg <- sprintf("%d fallos, %d avisos", FALLOS, AVISOS)
cat(msg, "\n")
if (FALLOS == 0) {
  cat("Sin fallos. Esto NO demuestra que el análisis sea correcto:\n")
  cat("demuestra que los datos son internamente consistentes y están\n")
  cat("actualizados. La validez de los modelos se juzga aparte.\n")
} else {
  cat("Hay fallos. No usar ningún resultado hasta resolverlos.\n")
}
