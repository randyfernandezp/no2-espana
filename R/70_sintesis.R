# =============================================================================
# 70 — OBJETIVO 5: síntesis de las tres ramas
#
# Este script NO estima nada nuevo. Lee lo que las tres ramas dejaron en
# disco y construye el argumento conjunto:
#
#   A. Insumos y verificación de procedencia
#   B. Convergencia de escalas (el hallazgo central)
#   C. Curva escala-señal: el efecto industrial frente a la agregación
#   D. ¿Qué explica el NO2 regional?
#   E. Magnitudes en unidades físicas, no en significancia
#   F. Qué aporta cada rama que las otras no pueden dar
#   G. Figura de síntesis
#
# ENTRADA : data/processed/*.rds (productos de 20_, 21_, 40_, 50_, 60_)
# SALIDAS : output/tablas/70_*.csv, output/figuras/50_*.png
#
# Las dos preguntas del objetivo son:
#   ¿El NO2 regional se explica por densidad industrial, urbanización, uso
#    del suelo, o una combinación?
#   ¿Qué evidencia aporta cada rama que las otras dos no pueden dar?
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(readr); library(ggplot2); library(tidyr)
})

PROC <- file.path("data", "processed")
FIG  <- file.path("output", "figuras")
TAB  <- file.path("output", "tablas")
for (p in c(FIG, TAB)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf(...), "\n")
sec <- function(t) cat("\n", strrep("=", 68), "\n", t, "\n",
                       strrep("=", 68), "\n", sep = "")

# =============================================================================
# A. INSUMOS
# =============================================================================
sec("A. INSUMOS Y PROCEDENCIA")

leer <- function(f) {
  p <- file.path(PROC, f)
  if (!file.exists(p)) { msg("AUSENTE: %s", f); return(NULL) }
  readRDS(p)
}

areal <- leer("modelos_areales.rds")   # de 60_areal.R
kppm_ <- leer("modelo_kppm.rds")       # de 40_procesos_puntuales.R
acopl <- leer("acoplamiento.rds")      # de 50_acoplamiento.R
n3    <- leer("nuts3_exposicion.rds")  # de 21_probabilidad.R

if (is.null(areal)) stop("Falta modelos_areales.rds. Ejecuta R/60_areal.R.")

# Procedencia: la síntesis solo vale si los insumos son de la MISMA cadena.
# Un producto más antiguo que otro significa que se están mezclando
# ejecuciones distintas, que es exactamente el fallo que ya se produjo una vez.
arch <- c("nuts3.rds", "industrias.rds", "pp_industrias.rds",
          "nuts3_exposicion.rds", "modelo_kppm.rds", "acoplamiento.rds",
          "modelos_areales.rds")
arch <- arch[file.exists(file.path(PROC, arch))]
tiempos <- file.mtime(file.path(PROC, arch))
print(data.frame(fichero = arch, modificado = format(tiempos)))
span <- as.numeric(difftime(max(tiempos), min(tiempos), units = "hours"))
if (span > 6)
  warning("Los insumos abarcan ", round(span, 1), " horas. Verifica que ",
          "proceden de la misma ejecución (R/99_verificar.R).")

# --- Extracción tolerante --------------------------------------------------
# Los objetos de 40_ y 50_ cambian de forma entre versiones de spatstat. Se
# intenta extraer; si no se consigue, se usa el valor documentado y se marca
# como tal, en lugar de fallar o de fingir que se ha leído.
DOC <- list(escala_acoplamiento_km = 10.0,
            diametro_distrito_km   = 10.7,
            F_acoplamiento         = 13.86,
            p_acoplamiento         = 2.5e-4,
            R2_base                = 0.7540,
            R2_presion             = 0.7686,
            efecto_e_pct           = 6.35,
            rango_variograma_km    = 71.0,
            mejora_rmse_pct        = 0.41)

origen <- c()
val <- function(nombre, expr) {
  v <- tryCatch(suppressWarnings(eval(expr)), error = function(e) NULL)
  if (is.null(v) || !is.finite(v)) {
    origen[[nombre]] <<- "documentado"
    DOC[[nombre]]
  } else {
    origen[[nombre]] <<- "extraído"
    as.numeric(v)
  }
}

diam_km <- val("diametro_distrito_km",
               quote(4 * kppm_$clustpar[["scale"]] / 1000))
esc_km  <- val("escala_acoplamiento_km",
               quote(acopl$h_opt / 1000))   # h_opt viene en metros (EPSG:3035)

msg("\nDiámetro del distrito industrial : %.1f km (%s)",
    diam_km, origen[["diametro_distrito_km"]])
msg("Escala de acoplamiento           : %.1f km (%s)",
    esc_km, origen[["escala_acoplamiento_km"]])

# =============================================================================
# B. CONVERGENCIA DE ESCALAS
# =============================================================================
sec("B. EL HALLAZGO CENTRAL: DOS RAMAS, UNA ESCALA")

# Las dos estimaciones son INDEPENDIENTES en el sentido fuerte: usan datos
# distintos y responden a preguntas distintas.
#
#   Objetivo 3 — solo coordenadas de instalaciones. Ninguna concentración
#     entra en el cálculo. Pregunta: ¿a qué distancia se agrupan las
#     instalaciones entre sí?
#   Objetivo 4 — instalaciones y estaciones. Pregunta: ¿a qué distancia deja
#     de notarse la emisión sobre el campo de NO2?
#
# No hay ninguna razón aritmética para que coincidan. Que lo hagan implica
# que la unidad de influencia sobre la concentración es el DISTRITO
# industrial, no la instalación aislada.

razon <- esc_km / diam_km
conv <- tibble(
  magnitud = c("Diámetro efectivo del distrito industrial",
               "Escala óptima de influencia sobre la concentración"),
  rama = c("Objetivo 3 — proceso puntual (Thomas)",
           "Objetivo 4 — acoplamiento (barrido de sigma)"),
  datos_usados = c("solo coordenadas de instalaciones",
                   "instalaciones + estaciones de fondo"),
  km = c(diam_km, esc_km)
)
print(as.data.frame(conv))
msg("\nRazón escala/diámetro: %.2f", razon)
write_csv(conv |> mutate(razon = razon), file.path(TAB, "70_convergencia.csv"))

if (abs(log(razon)) < log(1.35)) {
  msg("Las dos escalas coinciden dentro de un 35 %%.")
  msg("Lectura: la unidad de influencia sobre el NO2 es el distrito")
  msg("industrial, no la instalación aislada. Ninguna de las dos ramas")
  msg("puede establecer esto por sí sola.")
} else {
  msg("Las escalas difieren; la convergencia no se sostiene y hay que")
  msg("reportarlo como resultado negativo, no omitirlo.")
}

# Robustez: la convergencia se verificó dos veces sobre patrones distintos.
rob <- tibble(
  version = c("antes de depurar el pipeline", "tras regenerar la cadena"),
  instalaciones = c("296-299", "297"),
  diametro_km = c(9.5, diam_km),
  escala_km = c(10.0, esc_km)
) |> mutate(razon = escala_km / diametro_km)
cat("\nRobustez de la convergencia:\n")
print(as.data.frame(rob))
msg("Que sobreviva a un cambio en el patrón de puntos es evidencia de que")
msg("no es un artefacto de un ajuste concreto.")

# =============================================================================
# C. CURVA ESCALA-SEÑAL
# =============================================================================
sec("C. EL EFECTO INDUSTRIAL DESAPARECE AL AGREGAR")

# El mismo mecanismo físico medido a tres resoluciones espaciales. La
# industria es detectable a 10 km y se disuelve al promediar sobre unidades
# administrativas. No es que el efecto no exista: es que la unidad de
# análisis lo destruye.

d <- areal$datos
area_media <- mean(d$area_km2)
diam_equiv <- 2 * sqrt(area_media / pi)          # diámetro del disco equivalente
cobertura  <- 100 * (pi * (esc_km / 2)^2) / area_media

p_nox_n3 <- summary(areal$ols)$coefficients["log_nox", 4]
p_nox_sem <- summary(areal$sem)$Coef["log_nox", 4]
zon <- areal$maup_zonificacion
tasa_nox <- zon$pct_significativo[zon$termino == "log_nox"] / 100

escala <- tibble(
  soporte = c("Continuo (superficie a 2 km)",
              "NUTS-3 (50 provincias)",
              "NUTS-2 (16 comunidades)"),
  extension_km = c(esc_km, diam_equiv, diam_equiv * sqrt(nrow(d) / 16)),
  industria_significativa = c("sí", "no", "no (robusto)"),
  evidencia = c(sprintf("F = %.2f, p = %.1e", DOC$F_acoplamiento,
                        DOC$p_acoplamiento),
                sprintf("p = %.2f (OLS), %.2f (SEM)", p_nox_n3, p_nox_sem),
                sprintf("significativa en %.1f %% de 1000 zonificaciones",
                        100 * tasa_nox))
)
print(as.data.frame(escala))
write_csv(escala, file.path(TAB, "70_curva_escala.csv"))

msg("\nLa provincia media mide %.0f km2: un disco equivalente de %.0f km de",
    area_media, diam_equiv)
msg("diámetro. Un radio de influencia de %.0f km cubre el %.1f %% de ella.",
    esc_km, cobertura)
msg("Agregar a NUTS-3 promedia la señal industrial contra dos órdenes de")
msg("magnitud de territorio donde no ocurre nada.")
msg("")
msg("La urbanización sobrevive a la agregación porque OPERA a escala")
msg("provincial. La industria no. La diferencia no está en la intensidad")
msg("del mecanismo sino en su alcance espacial.")

# =============================================================================
# D. ¿QUÉ EXPLICA EL NO2 REGIONAL?
# =============================================================================
sec("D. RESPUESTA A LA PRIMERA PREGUNTA DEL OBJETIVO")

cat("\nReparto entre predictores (LMG, suma exacta el R2):\n")
print(as.data.frame(areal$lmg |> mutate(across(where(is.numeric), ~round(.x, 4)))))

cat("\nBloques temáticos (única frente a compartida):\n")
print(as.data.frame(areal$descomposicion |>
                      mutate(across(where(is.numeric), ~round(.x, 4)))))

unica  <- sum(areal$descomposicion$R2_unico)
# La tabla ya trae pct_unico = 100 * R2_unico / R2_total, así que el R2 total
# se recupera de cualquier fila sin volver a ajustar el modelo.
R2_tot <- areal$descomposicion$R2_unico[1] /
          (areal$descomposicion$pct_unico[1] / 100)
msg("\nR2 con todos los bloques: %.4f", R2_tot)
msg("Varianza única: %.4f (%.1f %% del R2). El resto es compartida.",
    unica, 100 * unica / R2_tot)
msg("")
msg("ATENCIÓN a la aparente contradicción entre las dos descomposiciones:")
msg("  LMG (3 predictores)  -> urbanización 84,5 %% del R2")
msg("  Bloques (todas)      -> uso de suelo aporta la mayor varianza ÚNICA")
msg("No se contradicen: miden cosas distintas. La LMG reparte el R2 de un")
msg("modelo donde la urbanización es el ÚNICO representante antrópico. La de")
msg("bloques mide qué aporta cada mecanismo que ningún otro puede aportar, y")
msg("ahí frac_artificial (rho = 0,955 con log_dens_pob) capta la señal")
msg("antrópica con más resolución que la densidad sola.")
msg("La lectura conjunta es que urbanización y uso de suelo son la MISMA")
msg("dimensión medida de dos formas, no dos mecanismos competidores.")
msg("")
msg("RESPUESTA: una combinación dominada por el eje antrópico")
msg("(urbanización/uso de suelo, indistinguibles entre sí a esta escala),")
msg("con la industria aportando menos del 1 %% de varianza única y una")
msg("fracción mayoritaria NO atribuible a ningún mecanismo por separado.")
msg("El dato honesto no es 'gana la urbanización' sino que los mecanismos")
msg("son inseparables a escala provincial: la industria se ubica donde hay")
msg("ciudades, puertos y autopistas, que son las mismas variables que")
msg("generan NO2 de tráfico.")
msg("")
msg("Esa endogeneidad no es un defecto del análisis: es una propiedad del")
msg("territorio, y el objetivo 3 la cuantifica (ver bloque F).")

# =============================================================================
# E. MAGNITUDES FÍSICAS
# =============================================================================
sec("E. CUÁNTO, NO SOLO SI")

# Un p-valor no dice cuánto. Estas tres cifras sí, y en las mismas unidades.
ic <- confint(areal$ols)["log_nox", ]
techo <- max(abs(ic)) * diff(range(d$log_nox))
media_nac <- mean(d$y)

mag <- tibble(
  concepto = c("Cota superior del efecto industrial areal",
               "Recargo local medido en estaciones industriales",
               "Recargo local medido en estaciones de tráfico",
               "Media areal nacional de NO2 de fondo"),
  ug_m3 = c(techo, 1.0, 3.5, media_nac),
  fuente = c("IC 95 % de log_nox sobre todo el rango de emisiones",
             "objetivo 1, validación externa (0,5-1,5)",
             "objetivo 1, validación externa (2,9-4,8)",
             "objetivo 1, agregación de bloque")
)
print(as.data.frame(mag |> mutate(ug_m3 = round(ug_m3, 3))))
write_csv(mag, file.path(TAB, "70_magnitudes.csv"))

msg("\nEl techo industrial (%.2f ug/m3, %.1f %% de la media nacional) es una",
    techo, 100 * techo / media_nac)
msg("afirmación más fuerte que una no-significancia: no es que ignoremos si")
msg("la industria influye, es que sabemos que no puede influir más que eso.")
msg("")
msg("Y el contraste con el tráfico es el mensaje de política pública: el")
msg("recargo de tráfico triplica al industrial y ocurre por debajo de la")
msg("resolución de la superficie, allí donde vive la gente.")

# =============================================================================
# F. QUÉ APORTA CADA RAMA
# =============================================================================
sec("F. RESPUESTA A LA SEGUNDA PREGUNTA DEL OBJETIVO")

aporta <- tribble(
  ~rama, ~evidencia_unica, ~por_que_las_otras_no_pueden,

  "Geoestadística (obj. 1)",
  "Concentración en puntos no muestreados, con incertidumbre propagada; y la separación entre campo de fondo y recargo microescala (tráfico 2,9-4,8 ug/m3 frente a industrial 0,5-1,5).",
  "El modelo areal parte de un agregado y no puede recuperar el interior de la unidad. El proceso puntual no observa concentraciones en absoluto.",

  "Procesos puntuales (obj. 3)",
  "Dónde se ubican las fuentes y por qué: suelo industrial x3,94 sobre la intensidad, densidad de población NO significativa. Y la agregación residual hasta 38 km, que revela distritos industriales tras condicionar en covariables.",
  "La geoestadística trata las fuentes como covariable dada, no como resultado. El modelo areal solo ve conteos por provincia y no puede distinguir 20 instalaciones agrupadas de 20 dispersas.",

  "Modelo areal (obj. 2)",
  "Exposición por unidad de gestión y su brecha frente a la media territorial (Madrid 7,52 -> 15,01). La cota del efecto industrial a escala administrativa. Y la dependencia espacial residual (lambda = 0,446).",
  "Ninguna de las otras dos produce cifras sobre la unidad en que se legisla ni pondera por población.",

  "Acoplamiento (obj. 4)",
  "La escala física de influencia industrial, 10 km, estimada por barrido y no impuesta. Es lo que permite comparar con el diámetro del distrito.",
  "El objetivo 3 estima agrupamiento de fuentes sin mirar concentraciones; el objetivo 1 usa covariables a escala fija; el objetivo 2 no tiene resolución para ver 10 km."
)
print(as.data.frame(aporta))
write_csv(aporta, file.path(TAB, "70_aportacion_ramas.csv"))

msg("\nY lo que NINGUNA rama da por separado: la coincidencia entre la")
msg("escala de agrupamiento de las fuentes (%.1f km) y la escala a la que",
    diam_km)
msg("influyen sobre la concentración (%.1f km). Esa comparación exige dos",
    esc_km)
msg("ramas metodológicamente independientes, y es el argumento que justifica")
msg("integrar las tres en lugar de elegir una.")

# =============================================================================
# G. FIGURA DE SÍNTESIS
# =============================================================================
sec("G. FIGURA")

# Panel único: el efecto industrial en función del soporte espacial, con las
# dos escalas convergentes marcadas.
gd <- tibble(
  soporte = factor(c("continuo\n(2 km)", "NUTS-3\n(50 unidades)",
                     "NUTS-2\n(16 unidades)"),
                   levels = c("continuo\n(2 km)", "NUTS-3\n(50 unidades)",
                              "NUTS-2\n(16 unidades)")),
  extension = c(esc_km, diam_equiv, diam_equiv * sqrt(nrow(d) / 16)),
  detectable = c("sí", "no", "no")
)

p <- ggplot(gd, aes(soporte, extension, fill = detectable)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = diam_km, linetype = 2, colour = "#b2182b") +
  annotate("text", x = 2, y = diam_km, vjust = -0.8, size = 3.4,
           colour = "#b2182b",
           label = sprintf("diámetro del distrito industrial: %.1f km", diam_km)) +
  scale_y_log10() +
  scale_fill_manual(values = c("sí" = "#2166ac", "no" = "grey75"),
                    name = "¿efecto industrial\ndetectable?") +
  labs(x = NULL, y = "extensión característica (km, escala log)",
       title = "La señal industrial se disuelve al agregar",
       subtitle = paste0("El mecanismo opera a ~", round(esc_km),
                         " km; las unidades administrativas promedian sobre",
                         " un orden de magnitud más")) +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG, "50_curva_escala.png"), p, width = 7.5, height = 5,
       dpi = 150)

# Descomposición por bloques
p2 <- areal$descomposicion |>
  select(bloque, unica = R2_unico, compartida = R2_compartido) |>
  pivot_longer(-bloque, names_to = "tipo", values_to = "R2") |>
  ggplot(aes(reorder(bloque, R2), R2, fill = tipo)) +
  geom_col(width = 0.62) + coord_flip() +
  scale_fill_manual(values = c("unica" = "#2166ac", "compartida" = "#bdd7e7"),
                    name = NULL) +
  labs(x = NULL, y = expression(R^2),
       title = "Ningún mecanismo explica el NO2 regional por sí solo",
       subtitle = "La mayor parte de la varianza es compartida entre bloques") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG, "51_sintesis_bloques.png"), p2, width = 7.5, height = 4.2,
       dpi = 150)

msg("Figuras: 50_curva_escala.png, 51_sintesis_bloques.png")

# =============================================================================
sec("OBJETIVO 5 COMPLETADO")
cat("Tablas en output/tablas/70_*.csv\n")
cat("Siguiente: redacción del informe.\n\n")
cat("Valores usados y su procedencia:\n")
print(data.frame(magnitud = names(origen), origen = unlist(origen)))
