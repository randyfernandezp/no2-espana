# =============================================================================
# 46_tres_escalas.R
# Panel comparativo de las tres escalas de agregación sobre España:
#   NUTS-2 (16 comunidades) | NUTS-3 (50 provincias) | distritos industriales
#
# El tercer panel NO dibuja los 520 padres del proceso de Thomas: esos son
# latentes y el modelo estima su número y dispersión, no sus posiciones.
# Dibuja la unión de discos del diámetro estimado (10,7 km) en torno a las
# instalaciones observadas, que es la huella espacial del mecanismo.
#
# Genera: output/figuras/45_tres_escalas.png
# Uso:    source("R/46_tres_escalas.R")
# =============================================================================

library(sf)

DIAMETRO_KM <- 10.7      # 4*sigma del proceso agregado ajustado
N_DISTRITOS <- 10        # cuántos distritos se rotulan (los de mayor emisión)
RUTA        <- "data/processed"

leer <- function(nombre) {
  f <- file.path(RUTA, paste0(nombre, ".rds"))
  if (!file.exists(f)) stop("No existe ", f)
  readRDS(f)
}

# --- Nombres a partir del código NUTS ----------------------------------------
# Las capas del pipeline guardan el código (NUTS2_ID / NUTS3_ID), no el nombre.
# Los códigos NUTS son estándar y estables, así que la tabla va aquí fijada.

NOMBRES_NUTS2 <- c(
  ES11 = "Galicia",            ES12 = "Asturias",
  ES13 = "Cantabria",          ES21 = "País Vasco",
  ES22 = "Navarra",            ES23 = "La Rioja",
  ES24 = "Aragón",             ES30 = "Madrid",
  ES41 = "Castilla y León",    ES42 = "Castilla-La Mancha",
  ES43 = "Extremadura",        ES51 = "Cataluña",
  ES52 = "C. Valenciana",      ES53 = "Baleares",
  ES61 = "Andalucía",          ES62 = "Murcia",
  ES63 = "Ceuta",              ES64 = "Melilla",
  ES70 = "Canarias")

NOMBRES_NUTS3 <- c(
  ES111 = "A Coruña",   ES112 = "Lugo",        ES113 = "Ourense",
  ES114 = "Pontevedra", ES120 = "Asturias",    ES130 = "Cantabria",
  ES211 = "Álava",      ES212 = "Gipuzkoa",    ES213 = "Bizkaia",
  ES220 = "Navarra",    ES230 = "La Rioja",    ES241 = "Huesca",
  ES242 = "Teruel",     ES243 = "Zaragoza",    ES300 = "Madrid",
  ES411 = "Ávila",      ES412 = "Burgos",      ES413 = "León",
  ES414 = "Palencia",   ES415 = "Salamanca",   ES416 = "Segovia",
  ES417 = "Soria",      ES418 = "Valladolid",  ES419 = "Zamora",
  ES421 = "Albacete",   ES422 = "Ciudad Real", ES423 = "Cuenca",
  ES424 = "Guadalajara",ES425 = "Toledo",      ES431 = "Badajoz",
  ES432 = "Cáceres",    ES511 = "Barcelona",   ES512 = "Girona",
  ES513 = "Lleida",     ES514 = "Tarragona",   ES521 = "Alicante",
  ES522 = "Castellón",  ES523 = "Valencia",    ES531 = "Eivissa-Form.",
  ES532 = "Mallorca",   ES533 = "Menorca",     ES611 = "Almería",
  ES612 = "Cádiz",      ES613 = "Córdoba",     ES614 = "Granada",
  ES615 = "Huelva",     ES616 = "Jaén",        ES617 = "Málaga",
  ES618 = "Sevilla",    ES620 = "Murcia",      ES630 = "Ceuta",
  ES640 = "Melilla")

col_codigo <- function(x) {
  cand <- c("NUTS2_ID", "NUTS3_ID", "NUTS_ID", "nuts2_id", "nuts3_id", "id")
  hit <- cand[cand %in% names(x)]
  if (!length(hit))
    stop("No encuentro columna de código NUTS. Disponibles: ",
         paste(names(x), collapse = ", "))
  hit[1]
}

# Traduce el código a nombre; si algún código no está en la tabla, deja el
# código tal cual y avisa, para que el fallo sea visible y no silencioso.
etiquetar <- function(x, tabla) {
  cod <- as.character(x[[col_codigo(x)]])
  nom <- unname(tabla[cod])
  faltan <- is.na(nom)
  if (any(faltan)) {
    warning("Códigos sin nombre en la tabla: ",
            paste(unique(cod[faltan]), collapse = ", "))
    nom[faltan] <- cod[faltan]
  }
  nom
}

nuts2 <- leer("nuts2")
nuts3 <- leer("nuts3")
ind   <- leer("industrias")

if (!inherits(ind, "sf")) {
  cols <- names(ind)
  xy <- cols[tolower(cols) %in% c("x", "y", "lon", "lat", "este", "norte")]
  if (length(xy) < 2)
    stop("`industrias.rds` no es sf. Columnas: ", paste(cols, collapse = ", "))
  ind <- st_as_sf(ind, coords = xy[1:2], crs = st_crs(nuts3))
}
ind <- st_transform(ind, st_crs(nuts3))

nom2 <- etiquetar(nuts2, NOMBRES_NUTS2)
nom3 <- etiquetar(nuts3, NOMBRES_NUTS3)
cat("NUTS-2:", nrow(nuts2), "| NUTS-3:", nrow(nuts3),
    "| instalaciones:", nrow(ind), "\n")

# --- Distritos: unión de discos del diámetro estimado ------------------------
radio_m <- DIAMETRO_KM * 1000 / 2
union_d <- st_union(st_buffer(st_geometry(ind), radio_m))
dominio <- st_union(st_geometry(nuts3))
piezas  <- st_cast(st_sfc(st_intersection(union_d, dominio),
                          crs = st_crs(nuts3)), "POLYGON")

area_dom <- as.numeric(sum(st_area(dominio))) / 1e6
area_pz  <- as.numeric(st_area(piezas)) / 1e6
cat("Distritos conexos:", length(piezas),
    "| superficie:", round(sum(area_pz)), "km2",
    sprintf("(%.1f%% del dominio)\n", 100 * sum(area_pz) / area_dom))

# --- Caracterizar cada distrito: emisión total e instalación dominante -------
# Se ordena por NOx acumulado, no por superficie: la superficie mide cuántas
# plantas se agrupan, no cuánto emiten. Los focos que cita el texto (Huelva,
# Tarragona, Puertollano) son pocos emisores muy grandes, no muchos pequeños.

dentro <- st_within(st_geometry(ind), piezas)          # a qué pieza va cada una
idx    <- vapply(dentro, function(z) if (length(z)) z[1] else NA_integer_, 1L)

nox <- if ("nox_t" %in% names(ind)) as.numeric(ind$nox_t) else rep(NA_real_,
                                                                   nrow(ind))
nox[is.na(nox)] <- 0

emision <- tapply(nox, factor(idx, levels = seq_along(piezas)), sum)
emision[is.na(emision)] <- 0
emision <- as.numeric(emision)

# Nombre del distrito = instalación de mayor emisión que contiene
nombre_pieza <- function(k) {
  q <- which(idx == k)
  if (!length(q)) return(paste0("distrito ", k))
  campo <- if ("siteName" %in% names(ind)) "siteName" else
           if ("cities"   %in% names(ind)) "cities"   else NULL
  if (is.null(campo)) return(paste0("distrito ", k))
  nm <- as.character(ind[[campo]][q[which.max(nox[q])]])
  if (is.na(nm) || !nzchar(nm)) paste0("distrito ", k) else limpiar(nm)
}

# Los nombres del registro vienen en mayúsculas y con forma jurídica: se
# recortan a algo legible sobre el mapa.
limpiar <- function(x) {
  x <- toupper(x)
  x <- gsub("\\b(S\\.?A\\.?U?|S\\.?L\\.?U?|C\\.?I\\.?)\\b", "", x)
  x <- gsub("\\b(ESPAÑA|IBERIA|PLANTA?|FACTORÍA)\\b.*$", "", x)
  x <- gsub("[[:space:],.-]+$", "", x)
  x <- gsub("\\s+", " ", trimws(x))
  x <- paste0(substr(x, 1, 1), tolower(substring(x, 2)))
  substr(x, 1, 22)
}

ord   <- order(emision, decreasing = TRUE)[seq_len(min(N_DISTRITOS,
                                                       length(piezas)))]
cen_d <- suppressWarnings(st_centroid(piezas[ord]))
etiq  <- vapply(ord, nombre_pieza, character(1))

cat("\nDistritos por emisión acumulada:\n")
print(data.frame(distrito = etiq,
                 nox_t    = round(emision[ord]),
                 km2      = round(area_pz[ord]),
                 plantas  = vapply(ord, function(k) sum(idx == k, na.rm = TRUE),
                                   1L)))

media_km2 <- function(x) round(mean(as.numeric(st_area(st_geometry(x)))) / 1e6)
centros   <- function(x) suppressWarnings(st_coordinates(
                            st_point_on_surface(st_geometry(x))))

# --- Figura ------------------------------------------------------------------
dir.create("output/figuras", recursive = TRUE, showWarnings = FALSE)
png("output/figuras/45_tres_escalas.png",
    width = 3900, height = 1250, res = 200)
op <- par(mfrow = c(1, 3), mar = c(0.5, 0.5, 3.2, 0.5))

encabezado <- function(t, s) {
  title(main = t, cex.main = 1.25, font.main = 2, line = 1.8)
  mtext(s, side = 3, line = 0.4, cex = 0.78, col = "grey35")
}

# Panel 1 — NUTS-2
plot(st_geometry(nuts3), border = NA, col = NA)
plot(st_geometry(nuts2), col = "grey93", border = "grey35", lwd = 0.7,
     add = TRUE)
text(centros(nuts2), labels = nom2,
     cex = 0.62, col = "grey15")
encabezado("NUTS-2: comunidades autónomas",
           sprintf("%d unidades  ·  %s km² de media",
                   nrow(nuts2), media_km2(nuts2)))

# Panel 2 — NUTS-3
plot(st_geometry(nuts3), border = NA, col = NA)
plot(st_geometry(nuts3), col = "grey93", border = "grey35", lwd = 0.5,
     add = TRUE)
text(centros(nuts3), labels = nom3,
     cex = 0.46, col = "grey15")
encabezado("NUTS-3: provincias",
           sprintf("%d unidades  ·  %s km² de media",
                   nrow(nuts3), media_km2(nuts3)))

# Panel 3 — distritos industriales
plot(st_geometry(nuts3), border = NA, col = NA)
plot(st_geometry(nuts3), border = "grey82", lwd = 0.3, add = TRUE)
plot(piezas, col = "#A0192D", border = NA, add = TRUE)
points(st_coordinates(ind), pch = 19, cex = 0.18, col = "black")
xy   <- st_coordinates(cen_d)
lado <- ifelse(xy[, 1] > mean(range(xy[, 1])), 2, 4)   # derecha -> texto a izq.
text(xy, labels = etiq, cex = 0.62, col = "#7A1222",
     pos = lado, offset = 0.55, font = 2)
encabezado("Distritos industriales",
           sprintf("%d conexos  ·  %.1f km de diámetro  ·  %.1f%% del territorio",
                   length(piezas), DIAMETRO_KM,
                   100 * sum(area_pz) / area_dom))

par(op)
dev.off()

cat("\nFigura escrita en output/figuras/45_tres_escalas.png\n")
