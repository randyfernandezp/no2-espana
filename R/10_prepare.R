# =============================================================================
# 10 — Preparación del conjunto de modelado
#
# Une todo lo descargado por los scripts de py/ en objetos listos para las
# tres ramas del análisis.
#
# ENTRADAS
#   data/interim/estaciones_modelado.csv    NO2 anual + clasificación + coords
#   data/interim/clc_estaciones.csv         fracciones de uso del suelo
#   data/interim/ied_es_nox_ventana.geojson patrón industrial 2019-2023
#   data/interim/clc_fracciones_nuts3.csv   uso del suelo por provincia
#   data/raw/nuts3_es.geojson               geometría NUTS-3
#   data/raw/s5p_no2_2023_media.tif         columna troposférica satelital
#   data/raw/densidad_nuts3_es.csv          densidad poblacional
#   data/raw/eprtr_emisiones.csv            masas de NOx (opcional)
#
# SALIDAS (data/processed/)
#   estaciones.rds     sf de puntos con NO2 y todas las covariables
#   nuts3.rds          sf de polígonos con covariables agregadas
#   nuts2.rds          idem por disolución, para el chequeo de MAUP
#   industrias.rds     sf de puntos industriales
#   pp_industrias.rds  objeto ppp de spatstat, con marcas si las hay
#   ventana.rds        owin del dominio de estudio
#
# CRS de trabajo: EPSG:3035 (ETRS89 / LAEA Europe), métrico y equiareal.
# Obligatorio: en EPSG:4326 los rangos del variograma saldrían en grados y
# las áreas del proceso puntual estarían distorsionadas por la latitud.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(dplyr); library(readr)
  library(stringr); library(tidyr); library(spatstat.geom)
})

CRS_T <- 3035
ANIO  <- 2023

RAW  <- file.path("data", "raw")
INT  <- file.path("data", "interim")
PROC <- file.path("data", "processed")
dir.create(PROC, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf(...), "\n")

# =============================================================================
# 1. DOMINIO DE ESTUDIO
# =============================================================================
# Se excluyen Canarias (ES70*), Ceuta (ES630) y Melilla (ES640).
#
# No es una decisión estética: las tres ramas presuponen un dominio conexo.
#   - El variograma estima una covarianza sobre distancias. Pares de
#     estaciones separados por 1 800 km de océano no informan sobre la
#     continuidad espacial del campo peninsular.
#   - El proceso puntual necesita una ventana con área bien definida;
#     incluir islas remotas obligaría a contar el Atlántico intermedio
#     como territorio observado.
#   - El modelo areal usa contigüidad, y las islas no tienen vecinos:
#     producirían regiones aisladas en la matriz de pesos W.
# Canarias podría analizarse aparte como dominio propio.

nuts3_all <- st_read(file.path(RAW, "nuts3_es.geojson"), quiet = TRUE) |>
  st_transform(CRS_T) |>
  st_make_valid()

EXCLUIR <- "^ES70|^ES63|^ES64"
nuts3 <- nuts3_all |> filter(!str_detect(NUTS_ID, EXCLUIR))
msg("NUTS-3: %d totales -> %d en el dominio peninsular + Baleares",
    nrow(nuts3_all), nrow(nuts3))

dominio <- st_union(nuts3) |> st_make_valid()
msg("Área del dominio: %.0f km2", as.numeric(st_area(dominio)) / 1e6)

# =============================================================================
# 2. ESTACIONES
# =============================================================================
est_tab <- read_csv(file.path(INT, "estaciones_modelado.csv"),
                    show_col_types = FALSE)
msg("Estaciones en la tabla: %d", nrow(est_tab))

if (!all(c("lon", "lat") %in% names(est_tab)))
  stop("Faltan lon/lat en estaciones_modelado.csv. Reejecuta py/08_verificar.py")

est <- est_tab |>
  filter(!is.na(lon), !is.na(lat), !is.na(no2_media)) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(CRS_T)

antes <- nrow(est)
est <- est |> st_filter(dominio)
msg("Dentro del dominio: %d de %d", nrow(est), antes)

# --- Clasificación: la decisión metodológica central -------------------------
# Las estaciones de tráfico miden a microescala. El NO2 cae entre un 60 y un
# 80 % en los primeros 50 m desde el borde de la calzada, de modo que su
# lectura NO es una realización del campo regional suave que queremos
# interpolar. Si entran en el variograma, el efecto pepita absorbe la
# varianza estructural y la superficie kriged pierde sentido.
#
# Estrategia: solo estaciones de fondo alimentan el modelo geoestadístico.
# Las de tráfico e industriales quedan como validación externa, y el sesgo
# medido en ellas CUANTIFICA el recargo local, que es un resultado en sí.
est <- est |>
  mutate(
    tipo = str_to_lower(as.character(tipo_estacion)),
    area = str_to_lower(as.character(area_estacion)),
    es_fondo = str_detect(tipo, "background"),
    x = st_coordinates(geometry)[, 1],
    y = st_coordinates(geometry)[, 2]
  )

print(table(est$tipo, est$area, useNA = "ifany"))

# =============================================================================
# 3. COVARIABLES EN LAS ESTACIONES
# =============================================================================

# --- 3a. Uso del suelo (fracciones normalizadas por superficie terrestre) ----
clc <- read_csv(file.path(INT, "clc_estaciones.csv"), show_col_types = FALSE) |>
  select(-any_of(c("x_3035", "y_3035")))
est <- est |> left_join(clc, by = "join_id")
msg("Covariables CLC añadidas: %d", ncol(clc) - 1)

# --- 3b. Columna troposférica de NO2 (Sentinel-5P) ---------------------------
# Se extrae transformando los PUNTOS al CRS del ráster, no al revés.
# Reproyectar el ráster remuestrea e introduce suavizado artificial;
# transformar puntos es exacto.
ruta_s5p <- file.path(RAW, sprintf("s5p_no2_%d_media.tif", ANIO))
if (file.exists(ruta_s5p)) {
  s5p <- rast(ruta_s5p)
  pts_s5p <- est |> st_transform(crs(s5p)) |> vect()
  # Escalado a 1e-5 mol/m2: evita coeficientes en notación científica
  # ilegible en las tablas del informe.
  est$s5p_no2 <- terra::extract(s5p, pts_s5p)[, 2] * 1e5
  msg("S5P: %d válidos | rango %.2f - %.2f (1e-5 mol/m2)",
      sum(!is.na(est$s5p_no2)),
      min(est$s5p_no2, na.rm = TRUE), max(est$s5p_no2, na.rm = TRUE))

  r_n <- file.path(RAW, sprintf("s5p_no2_%d_nmeses.tif", ANIO))
  if (file.exists(r_n))
    est$s5p_nmeses <- terra::extract(rast(r_n), pts_s5p)[, 2]
} else {
  warning("Falta el ráster S5P. El KDE quedará sin covariable satelital.")
  est$s5p_no2 <- NA_real_
}

# --- 3c. Altitud -------------------------------------------------------------
# Se usa la declarada por el operador, no un DEM. El operador conoce la cota
# exacta del emplazamiento; un DEM a 90 m promedia una celda entera, lo que
# en terreno abrupto (buena parte de España) yerra por decenas de metros.
if ("altitud_declarada" %in% names(est)) {
  est$altitud <- suppressWarnings(as.numeric(est$altitud_declarada))
  msg("Altitud: %d válidas | rango %.0f - %.0f m",
      sum(!is.na(est$altitud)),
      min(est$altitud, na.rm = TRUE), max(est$altitud, na.rm = TRUE))
}

# =============================================================================
# 4. PATRÓN DE PUNTOS INDUSTRIAL
# =============================================================================
# Variante temporal 2019-2023 (ver py/04c_filtrar_industrias.py): coherencia
# con el año de las concentraciones sin ser tan estricto como para perder
# instalaciones que quedaron un año por debajo del umbral de 100 t/año.
ruta_ied <- file.path(INT, "ied_es_nox_ventana.geojson")
if (!file.exists(ruta_ied)) ruta_ied <- file.path(INT, "ied_es_limpio.geojson")

ied <- st_read(ruta_ied, quiet = TRUE) |>
  st_transform(CRS_T) |>
  st_filter(dominio)
msg("Instalaciones en el dominio: %d (fuente: %s)",
    nrow(ied), basename(ruta_ied))

# --- Marcas: masa de NOx emitida --------------------------------------------
# Sin marcas, una térmica de 20 000 t y una cerámica de 110 t pesan igual en
# la estimación de intensidad. La marca es imprescindible para la superficie
# de presión de emisión del Objetivo 4.
ruta_emis <- file.path(RAW, "eprtr_emisiones.csv")
if (file.exists(ruta_emis)) {
  emis <- read_csv(ruta_emis, show_col_types = FALSE)
  col_geo <- names(ied)[str_detect(names(ied),
                                   regex("InspireSiteId", ignore_case = TRUE))][1]
  col_emis <- names(emis)[str_detect(names(emis),
                                     regex("InspireSiteId|facilityInspire",
                                           ignore_case = TRUE))][1]
  # La columna de masa NO se infiere. El CSV trae nox_kg, nox_kg_max y nox_t:
  # tres columnas que casan con cualquier patrón laxo y difieren en un factor
  # 1000. Adivinar aquí produce un error silencioso de tres órdenes de
  # magnitud, que es exactamente lo que pasó al detectarla por expresión
  # regular (eligió nox_kg y el total nacional salió en 1,2e8).
  col_masa <- "nox_t"
  if (!col_masa %in% names(emis))
    stop("Falta la columna 'nox_t' en ", basename(ruta_emis),
         ". Candidatas presentes: ",
         paste(grep("nox", names(emis), value = TRUE, ignore.case = TRUE),
               collapse = ", "),
         ". No se infiere automáticamente: difieren en un factor 1000.")

  if (!is.na(col_geo) && !is.na(col_emis)) {
    # CORRECCIÓN. Antes se unía emis directamente. El volcado del E-PRTR trae
    # una fila por declaración, de modo que un sitio con varias corrientes o
    # varios años aparece repetido: el left_join DUPLICA filas de ied, y esos
    # duplicados inflan tanto el conteo por provincia como la intensidad del
    # proceso puntual. Hay que colapsar a una fila por sitio ANTES de unir.
    emis_sitio <- emis |>
      select(.clave = all_of(col_emis), .masa = all_of(col_masa)) |>
      filter(!is.na(.clave)) |>
      group_by(.clave) |>
      summarise(nox_t = sum(.masa, na.rm = TRUE), .groups = "drop")
    msg("Emisiones: %d filas -> %d sitios únicos", nrow(emis), nrow(emis_sitio))

    n_antes <- nrow(ied)
    ied <- ied |>
      left_join(emis_sitio, by = setNames(".clave", col_geo))
    if (nrow(ied) != n_antes)
      stop("El join de emisiones cambió el número de instalaciones de ",
           n_antes, " a ", nrow(ied), ". Hay claves duplicadas sin colapsar.")

    n_marca <- sum(!is.na(ied$nox_t))
    # Guardia de unidades. España declara del orden de 1e5 t/año de NOx
    # industrial al E-PRTR. Un total cuatro órdenes por encima solo puede ser
    # kilogramos leídos como toneladas.
    total_t <- sum(ied$nox_t, na.rm = TRUE)
    if (total_t > 5e5)
      stop("Total de NOx = ", format(round(total_t), big.mark = " "),
           " t. Implausible para España (esperado ~1e5). ",
           "Probable confusión de unidades: revisar col_masa.")
    msg("Marcas de NOx: %d de %d instalaciones (%.0f %%) | total %.0f t",
        n_marca, nrow(ied), 100 * n_marca / nrow(ied),
        sum(ied$nox_t, na.rm = TRUE))
    if (n_marca < 0.5 * nrow(ied)) {
      msg("AVISO: menos de la mitad casan. Revisa que ambos identificadores")
      msg("sean del mismo nivel jerárquico (sitio frente a instalación).")
    }
  } else {
    msg("No identifico las columnas de unión; patrón sin marcas.")
    ied$nox_t <- NA_real_
  }
} else {
  warning("Sin data/raw/eprtr_emisiones.csv: patrón SIN marcas de masa.")
  ied$nox_t <- NA_real_
}

# CORRECCIÓN CRÍTICA. sum(x, na.rm = TRUE) sobre un vector íntegramente NA
# devuelve 0, no NA. Cuando este script se ejecutaba antes de descargar
# eprtr_emisiones.csv, las 50 provincias recibían nox_total_t = 0 sin un solo
# aviso, y el replace_na(.x, 0) de más abajo remataba el disimulo. Ese cero
# viajaba a nuts3.rds y de ahí a todo lo aguas abajo, donde es
# indistinguible de una provincia sin grandes emisores.
#
# La bandera evita que el conjunto se guarde en silencio con marcas ausentes.
HAY_MARCAS <- !all(is.na(ied$nox_t))
if (!HAY_MARCAS)
  warning("SIN MARCAS DE MASA. nox_total_t se guardará como NA, no como 0. ",
          "Ejecuta py/04b_emisiones.py y vuelve a lanzar este script antes ",
          "de usar cualquier resultado industrial.")

# CORRECCIÓN. Antes se deduplicaba el ppp pero no ied, de modo que
# npoints(pp) y nrow(ied) divergían y los conteos por provincia del modelo
# areal no coincidían con el patrón puntual. La deduplicación se hace ahora
# sobre ied, y el ppp se construye ya limpio.
#
# El criterio es la coordenada, no la coordenada más la marca: dos registros
# en el mismo punto son el mismo emplazamiento aunque declaren masas
# distintas, y sus masas deben sumarse, no competir.
xy <- st_coordinates(ied)
clave_xy <- paste(round(xy[, 1], 1), round(xy[, 2], 1))
n_dup <- sum(duplicated(clave_xy))
if (n_dup > 0) {
  msg("Colapsando %d instalaciones en coordenadas duplicadas", n_dup)
  ied <- ied |>
    mutate(.clave_xy = clave_xy) |>
    group_by(.clave_xy) |>
    mutate(nox_t = if (all(is.na(nox_t))) NA_real_ else sum(nox_t, na.rm = TRUE)) |>
    slice(1) |>
    ungroup() |>
    select(-.clave_xy)
}
msg("Instalaciones tras deduplicar: %d", nrow(ied))

# --- Ventana de observación --------------------------------------------------
# Debe ser la frontera real, no el rectángulo envolvente: la corrección de
# borde de la K de Ripley depende de ella, y con bbox se asignaría intensidad
# al mar y a territorio francés y portugués.
# Se simplifica a 1 km para que spatstat no se ahogue con miles de vértices;
# a escala nacional el error es despreciable.
#
# Va DESPUÉS de la deduplicación porque necesita comprobar qué instalaciones
# quedan dentro.
ventana_sf <- st_simplify(dominio, dTolerance = 1000) |> st_make_valid()

# La simplificación desplaza la costa hasta 1 km hacia dentro y puede dejar
# fuera instalaciones portuarias que sí pertenecen al dominio. spatstat las
# rechaza en silencio ("point rejected as lying outside the window") y el
# patrón queda con menos puntos que ied. Se dilata la ventana lo mínimo para
# recuperarlas: el área cambia menos del 1 % y la corrección de borde de la K
# de Ripley no se ve afectada a escala nacional.
fuera <- lengths(st_within(ied, ventana_sf)) == 0
if (any(fuera)) {
  msg("%d instalaciones fuera de la ventana simplificada; dilatando 1 km",
      sum(fuera))
  ventana_sf <- st_buffer(ventana_sf, 1000) |> st_make_valid()
  fuera <- lengths(st_within(ied, ventana_sf)) == 0
  if (any(fuera)) {
    nom <- if ("siteName" %in% names(ied)) ied$siteName[fuera] else which(fuera)
    msg("Siguen fuera y se descartan: %s", paste(nom, collapse = ", "))
    ied <- ied[!fuera, ]
  }
}
ventana <- as.owin(ventana_sf)
msg("Ventana: %.0f km2 | dominio real: %.0f km2",
    area(ventana) / 1e6, as.numeric(st_area(dominio)) / 1e6)

marcas <- if (!HAY_MARCAS) NULL else data.frame(nox_t = ied$nox_t)
pp <- ppp(
  x = st_coordinates(ied)[, 1],
  y = st_coordinates(ied)[, 2],
  window = ventana,
  marks = marcas
)
stopifnot(npoints(pp) == nrow(ied))
msg("Patrón: %d instalaciones | intensidad %.3f por 1000 km2",
    npoints(pp), 1000 * npoints(pp) / (area(ventana) / 1e6))

# =============================================================================
# 5. AGREGADOS POR REGIÓN
# =============================================================================
dens <- read_csv(file.path(RAW, "densidad_nuts3_es.csv"), show_col_types = FALSE)
clc3 <- read_csv(file.path(INT, "clc_fracciones_nuts3.csv"), show_col_types = FALSE)

# CORRECCIÓN. st_join sin argumento usa st_intersects: un punto que caiga
# exactamente sobre una frontera se asigna a DOS provincias y se cuenta dos
# veces. st_within es unívoco para puntos.
ind_nuts <- ied |>
  st_join(nuts3["NUTS_ID"], join = st_within) |>
  st_drop_geometry() |>
  filter(!is.na(NUTS_ID)) |>
  group_by(NUTS_ID) |>
  summarise(n_instalaciones = n(),
            # Si no hay marcas, NA. Nunca 0: un cero por dato ausente es
            # indistinguible de una provincia sin grandes emisores.
            nox_total_t = if (!HAY_MARCAS) NA_real_ else sum(nox_t, na.rm = TRUE),
            .groups = "drop")

# Verificación de conservación: toda instalación del dominio debe caer en
# alguna provincia y toda tonelada debe quedar asignada.
n_asignadas <- sum(ind_nuts$n_instalaciones)
if (n_asignadas != nrow(ied))
  warning("Solo ", n_asignadas, " de ", nrow(ied), " instalaciones caen ",
          "dentro de una provincia. Revisar geometrías costeras.")
if (HAY_MARCAS) {
  perdidas <- sum(ied$nox_t, na.rm = TRUE) - sum(ind_nuts$nox_total_t)
  msg("NOx: %.1f t en instalaciones | %.1f t asignadas a provincias",
      sum(ied$nox_t, na.rm = TRUE), sum(ind_nuts$nox_total_t))
  if (abs(perdidas) > 1)
    warning("Se pierden ", round(perdidas, 1), " t en el cruce espacial.")
}

# Media por estaciones: PROVISIONAL.
# El Objetivo 2 NO debe usarla. Promediar estaciones pondera por densidad de
# red, no por territorio: Madrid tiene 24 estaciones en 8 000 km2 y Soria
# una en 10 000. La media correcta se obtiene agregando la superficie kriged
# (script 60), lo que además propaga la incertidumbre.
est_nuts <- est |>
  st_join(nuts3["NUTS_ID"], join = st_within) |>
  st_drop_geometry() |>
  filter(!is.na(NUTS_ID)) |>
  group_by(NUTS_ID) |>
  summarise(n_estaciones = n(),
            no2_estaciones = mean(no2_media, na.rm = TRUE),
            # which() y no es_fondo a secas: si es_fondo tiene algún NA, el
            # subíndice lógico devuelve filas NA y contamina la media.
            no2_fondo = mean(no2_media[which(es_fondo)], na.rm = TRUE),
            .groups = "drop")

nuts3 <- nuts3 |>
  left_join(dens, by = "NUTS_ID") |>
  left_join(clc3, by = "NUTS_ID") |>
  left_join(ind_nuts, by = "NUTS_ID") |>
  left_join(est_nuts, by = "NUTS_ID") |>
  mutate(
    # n_instalaciones y n_estaciones sí van a 0: un conteo ausente ES un cero.
    # nox_total_t NO: solo se rellena con 0 cuando hay marcas y esa provincia
    # simplemente no tiene emisores por encima del umbral E-PRTR.
    across(c(n_instalaciones, n_estaciones), ~ replace_na(.x, 0)),
    nox_total_t  = if (HAY_MARCAS) replace_na(nox_total_t, 0) else NA_real_,
    area_km2     = as.numeric(st_area(geometry)) / 1e6,
    dens_instal  = 1000 * n_instalaciones / area_km2,
    nox_por_km2  = nox_total_t / area_km2,
    log_dens_pob = log(dens_2023),
    poblacion    = dens_2023 * area_km2
  )

# Los CSV de densidad y de fracciones CLC traen sus propias columnas de
# identificación en minúsculas, que quedan duplicadas tras los joins. Se
# comprueba que son idénticas y se elimina la redundante, para que los
# scripts posteriores no unan por la versión equivocada.
for (par in list(c("nuts_id", "NUTS_ID"), c("nuts_name", "NUTS_NAME"))) {
  if (all(par %in% names(nuts3))) {
    if (identical(nuts3[[par[1]]], nuts3[[par[2]]])) {
      nuts3[[par[1]]] <- NULL
      msg("Columna duplicada eliminada: %s", par[1])
    } else {
      stop("Las columnas ", par[1], " y ", par[2], " difieren. ",
           "Revisar las claves de densidad_nuts3_es.csv y clc_fracciones_nuts3.csv")
    }
  }
}

msg("NUTS-3 con covariables: %d x %d", nrow(nuts3), ncol(nuts3))
msg("Provincias sin ninguna estación: %d", sum(nuts3$n_estaciones == 0))

# --- NUTS-2 por disolución ---------------------------------------------------
# El servicio del EEA solo publica los niveles 0 y 3, pero la jerarquía NUTS
# está codificada en el identificador: ES300 (provincia de Madrid) pertenece
# a ES30 (Comunidad de Madrid). La agregación por los cuatro primeros
# caracteres es exacta, no aproximada.
# Sirve para el chequeo de MAUP: si los coeficientes del modelo areal cambian
# de signo o magnitud al pasar de provincia a comunidad, la conclusión depende
# de la escala de agregación y hay que declararlo.
# Se separa en dos pasos deliberadamente:
#   (a) atributos sobre la tabla sin geometría
#   (b) disolución geométrica
# Motivo: dentro de un mismo summarise(), dplyr permite que una expresión
# vea las columnas creadas por las anteriores. Si se calcula primero
# area_km2 = sum(area_km2), las medias ponderadas posteriores reciben como
# pesos un escalar de longitud 1 en lugar del vector original, y fallan.
# Separarlo evita ese acoplamiento y además es más rápido, porque no
# arrastra la geometría en los cálculos de atributos.

attr3 <- st_drop_geometry(nuts3) |>
  mutate(NUTS2_ID = str_sub(NUTS_ID, 1, 4))

nuts2_attr <- attr3 |>
  group_by(NUTS2_ID) |>
  summarise(
    # Las medias ponderadas van PRIMERO, con los pesos aún vectoriales
    frac_artificial = weighted.mean(frac_artificial, area_km2, na.rm = TRUE),
    frac_industrial = weighted.mean(frac_industrial, area_km2, na.rm = TRUE),
    # y después los totales
    area_km2        = sum(area_km2, na.rm = TRUE),
    poblacion       = sum(poblacion, na.rm = TRUE),
    n_instalaciones = sum(n_instalaciones, na.rm = TRUE),
    nox_total_t     = sum(nox_total_t, na.rm = TRUE),
    n_estaciones    = sum(n_estaciones, na.rm = TRUE),
    n_provincias    = n(),
    .groups = "drop"
  ) |>
  mutate(
    dens_pob     = poblacion / area_km2,
    log_dens_pob = log(dens_pob),
    dens_instal  = 1000 * n_instalaciones / area_km2,
    nox_por_km2  = nox_total_t / area_km2
  )

nuts2_geom <- nuts3 |>
  mutate(NUTS2_ID = str_sub(NUTS_ID, 1, 4)) |>
  group_by(NUTS2_ID) |>
  summarise(.groups = "drop")

nuts2 <- nuts2_geom |> left_join(nuts2_attr, by = "NUTS2_ID")

msg("NUTS-2 derivadas: %d", nrow(nuts2))

# =============================================================================
# 6. GUARDAR
# =============================================================================
# Guardia final: no dejar en disco un conjunto con las marcas industriales
# vacías, porque nuts3.rds alimenta a todos los scripts posteriores y el
# defecto sería invisible aguas abajo.
if (!HAY_MARCAS)
  warning("Se guarda nuts3.rds SIN masas de NOx (nox_total_t = NA). ",
          "Cualquier resultado industrial obtenido a partir de este fichero ",
          "carece de validez.")
if (HAY_MARCAS && all(nuts3$nox_total_t == 0, na.rm = TRUE))
  stop("nox_total_t es cero en TODAS las provincias pese a haber marcas. ",
       "El cruce espacial ha fallado; no se guarda nada.")

saveRDS(est,     file.path(PROC, "estaciones.rds"))
saveRDS(nuts3,   file.path(PROC, "nuts3.rds"))
saveRDS(nuts2,   file.path(PROC, "nuts2.rds"))
saveRDS(ied,     file.path(PROC, "industrias.rds"))
saveRDS(pp,      file.path(PROC, "pp_industrias.rds"))
saveRDS(ventana, file.path(PROC, "ventana.rds"))
saveRDS(dominio, file.path(PROC, "dominio.rds"))

msg("IMPORTANTE: nuts3_exposicion.rds es un derivado de nuts3.rds. Tras")
msg("regenerar este fichero hay que volver a ejecutar 20_ y 21_ antes de 60_.")

# =============================================================================
# 7. RESUMEN
# =============================================================================
cat("\n", strrep("=", 66), "\n", sep = "")
cat("CONJUNTO DE MODELADO\n")
cat(strrep("=", 66), "\n", sep = "")

msg("Estaciones totales      : %d", nrow(est))
msg("  de fondo (variograma) : %d", sum(est$es_fondo, na.rm = TRUE))
msg("  de validación externa : %d", sum(!est$es_fondo, na.rm = TRUE))
msg("Instalaciones (proceso) : %d", npoints(pp))
msg("Regiones NUTS-3         : %d", nrow(nuts3))
msg("Regiones NUTS-2         : %d", nrow(nuts2))

cat("\nNO2 medio anual (ug/m3) por tipo y área:\n")
print(
  est |> st_drop_geometry() |>
    group_by(tipo, area) |>
    summarise(n = n(), media = mean(no2_media), max = max(no2_media),
              .groups = "drop") |>
    mutate(across(where(is.numeric), ~ round(.x, 2))) |>
    as.data.frame()
)

cat("\nCobertura de covariables en estaciones de fondo:\n")
covs <- c("s5p_no2", "altitud", "artificial_1000m", "artificial_5000m",
          "industrial_5000m", "tierra_5000m", "poblacion_ciudad")
for (v in intersect(covs, names(est))) {
  x <- st_drop_geometry(est)[est$es_fondo, v, drop = TRUE]
  msg("  %-20s %3d/%3d validos | mediana %8.3f",
      v, sum(!is.na(x)), length(x), median(x, na.rm = TRUE))
}

cat("\nSiguiente paso: R/20_variograma_kriging.R\n")
