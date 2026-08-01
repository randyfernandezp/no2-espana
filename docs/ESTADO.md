# Estado del proyecto — NO₂ en España

Documento de continuidad. Si retomas el trabajo en otra sesión, empieza por aquí.

**Ruta local:** `C:\Users\randy.fernandez\Documents\no2-espana`
**Entorno Python:** `.venv` (Python 3.12) — activar con `.venv\Scripts\Activate.ps1`
**R:** 4.5.1 + RStudio. Directorio de trabajo: la raíz del proyecto.

---

## Qué está hecho

| Fase | Estado |
|---|---|
| Descarga de datos (7 fuentes) | Completa y verificada |
| Objetivo 1 — geoestadística | Completo |
| Objetivo 3 — procesos puntuales | Completo |
| Objetivo 4 — acoplamiento | Completo |
| Objetivo 2 — modelo areal | Completo (SEM + descomposición + MAUP escala y zonificación) |
| Objetivo 5 — síntesis | Pendiente de escribir |
| Integridad del pipeline | `99_verificar.R`: 0 fallos, 0 avisos (01/08/2026) |

---

## Decisiones de diseño (no reabrir sin motivo)

**País y contaminante:** España, NO₂, año 2023. España se eligió porque sus focos industriales (Huelva, Tarragona, Puertollano, As Pontes, Avilés, Algeciras) están espacialmente desacoplados de los grandes focos urbanos, lo que permite separar el efecto industrial del de urbanización. En Italia serían colineales.

**CRS de trabajo:** EPSG:3035 en todo. Métrico y equiareal.

**Dominio:** peninsular + Baleares. Se excluyen Canarias, Ceuta y Melilla porque las tres ramas presuponen un dominio conexo.

**Estaciones de fondo para el variograma.** Las de tráfico miden a microescala (el NO₂ cae 60-80 % en los primeros 50 m desde la calzada) y no son realizaciones del campo regional. Se usan como validación externa; el sesgo medido ahí cuantifica el recargo local.

**Desagrupamiento** por celdas de 25 km: la media pasa de 11,36 a 8,58 µg/m³.

**Patrón de puntos:** variante temporal 2019-2023, filtrado a sitios que declaran NOx. Deduplicado por `InspireSiteId` (el volcado bruto traía 103 299 registros para 9 168 sitios, con 18 años de declaraciones).

**Ventana del proceso puntual:** frontera real simplificada a 1 km y dilatada 1 km para no perder instalaciones portuarias. 504 327 km² frente a 498 522 del dominio real (+1,2 %), que es el sesgo a la baja de la intensidad estimada.

**Normalización costera de CORINE:** el mar es la clase 523 (índice 44), no `nodata`. Sin corregir, Cádiz aparecía con fracción artificial 0,15 cuando es 0,72.

---

## Resultados principales

**Objetivo 1.** Variograma gaussiano, rango práctico 71 km, pepita 9,6 %. Las covariables absorben el 76 % de la varianza (rango residual 24 km). KDE frente a KO: R² 0,776 vs 0,481, RMSE −34,3 %, sesgo 0,16 vs 0,93. Cobertura del intervalo al 95 %: 0,91.

**Exposición.** Madrid 7,52 µg/m³ areal → 15,01 ponderado por población. Barcelona 6,64 → 14,21. Tres provincias con más del 10 % de población sobre 20 µg/m³: Madrid (37,3 %), Barcelona (30,6 %), Bizkaia (22,7 %). El 43,3 % de las estaciones de tráfico ya superan hoy los 20 µg/m³, frente al 0,16 % de la superficie de fondo.

**Objetivo 3.** 297 instalaciones (patrón deduplicado). Test de cuadrantes Monte Carlo X² = 531,64, p = 0,002. Ancho de banda: Diggle 3,7 km vs verosimilitud 21,8 km — difieren en un factor 5,9, el patrón tiene estructura a dos escalas y se analiza con ambos. Efectos por rango intercuartílico sobre la intensidad: suelo industrial ×3,94, superficie artificial ×1,32, transporte ×0,748, **densidad de población ×0,999 (p = 0,997, no significativa)**. Queda agregación no explicada hasta 38 km. Modelo de Thomas: escala 2,7 km, **diámetro efectivo del distrito 10,7 km**, 520 distritos, psib = 0,915, phi = 10,74. La marca de NOx es muy asimétrica (3,30): las 10 mayores concentran el 17,2 % del total, de ahí que la superficie de presión use log(masa).

**Objetivo 4.** Barrido de siete anchos de banda (5 a 100 km). **Escala óptima 10 km**, con curva de RMSE en U limpia (peor a 5 y a 15). F = 13,86, p = 2,5 × 10⁻⁴. R² ajustado 0,754 → 0,769. Multiplicar por *e* la presión industrial sube el NO₂ un 6,35 %. **La mejora del RMSE es solo del 0,41 %**: significativa pero de poco poder predictivo, coherente con que Canarias concentre el 24 % del NOx industrial fuera del dominio.

**Objetivo 2.** Respuesta: `no2_areal` (media de bloque), **no** `no2_poblacion` — la exposición ponderada correlaciona 0,948 con `log_dens_pob` por construcción (la ponderación usa la misma distribución de población que el predictor), así que regresarla sobre densidad es circular. Queda como métrica de política pública, no como variable dependiente.

Modelo: `no2_areal ~ log_dens_pob + log_nox + log_area`, n = 50. Colinealidad resuelta: el bloque antrópico (`log_dens_pob`, `frac_artificial`, `frac_urbano`, `frac_industrial`) va junto a ρ 0,77–0,97, y `dens_ind` ↔ `log_presion` = 0,943. El confusor latente es el área: todas las variables intensivas van contra `log_area` (−0,54 a −0,69) porque las provincias pequeñas son densas, urbanas, costeras e industriales a la vez. `log_nox` es un total, no una densidad, y da −0,046 con el área: es la única medida industrial que no es alias de urbanización. VIF < 1,8.

OLS: `log_dens_pob` β = 0,575 (p = 2,1 × 10⁻⁹), `log_nox` p = 0,63, `log_area` p = 0,18. R² aj. 0,601. Moran sobre residuos I = 0,185 (p = 0,0085). Con contigüidad hay 4 subgrafos (península + Mallorca + Menorca + Eivissa-Formentera), así que el diagnóstico se hace con k = 5 vecinos: `RSlag` p = 0,92 frente a `RSerr` p = 0,027, y al robustecer `adjRSerr` 0,0015 vs `adjRSlag` 0,023. Firma limpia de Anselin–Florax para **error espacial**.

SEM (`errorsarlm`, lw5): lambda = 0,446 (LR p = 0,031; Wald p = 0,0090). Moran de residuos p = 0,54 — dependencia absorbida. AIC 69,92 vs 72,57 del OLS. `log_dens_pob` **sube** a 0,601 (z = 9,04): el campo espacial residual enmascaraba el efecto de urbanización, no lo inflaba. `log_nox` sigue sin significancia (p = 0,68).

**Decisión: se reporta SEM, con SDM como sensibilidad.** El AIC prefiere por los pelos el Durbin (69,67 vs 69,92 del SEM, 72,57 del OLS, 74,56 del SAR), pero ΔAIC = 0,25 es indistinguible y tres argumentos apuntan al error espacial: los contrastes robustos separan limpiamente (`adjRSerr` 0,0015 vs `adjRSlag` 0,023), el rho del SAR no es significativo (p = 0,91), y **ninguno de los impactos directos, indirectos o totales del SDM alcanza significancia** (p entre 0,11 y 0,74). Elegir SDM por 0,25 de AIC contra el criterio de los contrastes sería dejarse llevar por el ruido.

**Descomposición LMG** (promedio sobre las 3! ordenaciones, suma exacta el R²): urbanización **84,5 %**, área 10,8 %, industria **4,9 %**. En R²: 0,528 / 0,066 / 0,031.

**Descomposición del SEM** por covarianzas (aditiva exacta, suma 1,000000): efectos fijos **68,3 %**, campo espacial **0,2 %**, ruido **31,5 %**. Lambda es significativa pero el campo estructurado apenas aporta varianza — la dependencia existe y conviene modelarla para no subestimar los errores estándar, pero no es un mecanismo sustantivo.

**Descomposición por bloques temáticos** (R² total 0,9397 con todo dentro): único de uso de suelo 0,1556, urbanización 0,0224, industria 0,0050, geografía 0,0018. **Varianza única total: 19,7 %**; el 80,3 % restante es compartida y no atribuible a un solo mecanismo. Ese reparto es en sí un resultado del objetivo 5.

**Cota del efecto industrial — el número más fuerte del objetivo 2.** IC 95 % de `log_nox`: [−0,1006, +0,0617]. Al recorrer todo el rango observado de emisiones, el efecto máximo compatible con los datos es **0,956 µg/m³, el 22,3 % de la media areal nacional**. Un techo cuantificado dice mucho más que un p-valor: no es que no sepamos si la industria influye, es que sabemos que no puede influir más que eso.

**La industria no aparece en cuatro especificaciones**: OLS, con `log_presion` en lugar de `log_nox` (p = 0,36, R² aj. idéntico 0,601), sin Madrid (p = 0,90), y bajo SEM (p = 0,68).

**MAUP — escala.** A NUTS-2 (16 unidades) `log_nox` sale significativa y negativa (β = −0,562, p = 0,0062), pero es íntegramente Madrid: Cook = 4,64, diecinueve veces el umbral 4/n y diecisiete veces la siguiente. Sin ES30 el coeficiente cae a −0,063 (p = 0,68) y el R² aj. de 0,775 a 0,597, casi idéntico al 0,601 de NUTS-3. La agregación infla el ajuste aparente sin aportar información. **No hay inversión de signo con la escala: la señal industrial se disuelve, no se invierte.**

**MAUP — zonificación.** 1000 particiones de las 50 provincias en 16 grupos por k-medias sobre centroides. β_nox: IC 95 % [−0,146, +0,003], mediana −0,053 — el −0,562 de NUTS-2 queda **fuera del rango completo por un factor de casi cuatro**. Significancia de `log_nox` en el 2,5 % de las particiones (por debajo del nominal); de `log_dens_pob`, en el 100 % (β entre 0,290 y 0,633). Madrid **nunca** queda aislada en las 1000 réplicas: k-medias agrupa ~3 provincias por centro y Madrid está rodeada en el centro de la meseta. NUTS-2 la aísla por razones históricas y políticas, no geométricas — de ahí el artefacto.

**El hallazgo central:** la escala de influencia industrial sobre la concentración (10,0 km) coincide con el diámetro del distrito industrial estimado por el proceso puntual (10,7 km). **Razón 0,9.** Dos ramas independientes convergen en la misma escala. La unidad de influencia sobre el NO₂ es el distrito industrial, no la instalación aislada.

La convergencia se ha verificado dos veces sobre datos distintos: antes de depurar el pipeline daba 10 / 9,5 = 1,05, y tras regenerar toda la cadena con el patrón deduplicado da 10 / 10,7 = 0,90. Que sobreviva a un cambio en el patrón de puntos es evidencia de que no es un artefacto de un ajuste concreto.

**Recargo local sobre el fondo predicho:** tráfico 2,9-3,8 µg/m³, industrial 0,5-1,2. En entorno rural las estaciones industriales miden un 52 % más que las de fondo (5,95 vs 3,91); en entorno urbano, menos. La industria eleva el NO₂ donde el fondo es bajo, pero queda sepultada por el tráfico en las ciudades.

---

## Limitaciones a declarar en el informe

- **Umbral E-PRTR.** Solo se declara por encima de 100 t/año de NOx. La marca mínima observada es 102 t. El patrón representa grandes emisores, no toda fuente industrial.
- **Canarias concentra el 24 % del NOx industrial español** en siete centrales diésel, fuera del dominio. Por eso el efecto industrial peninsular es modesto.
- **CORINE clase 121** mezcla industria con comercio y logística, por lo que está casi colineal con urbanización. La señal industrial limpia es el patrón E-PRTR.
- **La superficie modela NO₂ de fondo**, no exposición máxima. Los picos de tráfico ocurren por debajo de la resolución de 2 km.
- **Pesos dasimétricos** sin escalar por población provincial: las medias dentro de cada provincia son correctas, el agregado nacional está sesgado.
- **Envolventes con lambda fijada** de los propios datos: contraste ligeramente anticonservador.
- **Curva escala–señal (material central del objetivo 5).** La industria es significativa a 10 km en el continuo (F = 13,86, p = 2,5 × 10⁻⁴), ausente a NUTS-3 (~113 km de diámetro equivalente, n = 50) y ausente a NUTS-2 (n = 15 sin Madrid). La provincia media peninsular ronda 10 000 km²; un radio de 10 km cubre el 3 % de ella. Agregar a NUTS-3 promedia la señal industrial contra dos órdenes de magnitud de territorio donde no pasa nada. La urbanización sobrevive a la agregación porque opera a escala provincial; la industria no.
- **Reportar m2 completo (16 unidades), con advertencia explícita.** Excluir la capital de un análisis de política pública no se sostiene, pero publicar β_nox = −0,56 significativa sin marcarla como no robusta equivaldría a afirmar que la industria limpia el aire. m2b va como evidencia de robustez.
- **No hay autocorrelación global en la respuesta, pero sí en los residuos.** Moran de `no2_areal` con k = 5: I = 0,013, p = 0,32 — no significativa, y LISA no detecta ni un agrupamiento alto-alto. Sobre los residuos del OLS, en cambio, I = 0,187 con p = 0,0019. Son cosas distintas y conviene no confundirlas en el informe: la dependencia aparece **después** de condicionar en densidad y área, no antes. Justifica el SEM. Y refuerza el punto anterior: agregar un campo continuo de rango 71 km sobre polígonos de ~113 km podría haber generado autocorrelación por construcción, y sin embargo la respuesta no la muestra.
- **Barcelona quinta en NOx industrial (5 666 t).** El desacoplamiento entre focos industriales y urbanos que justificó elegir España es cierto en el continuo, pero se pierde al agregar a NUTS-3: la provincia se traga el Baix Llobregat. Es MAUP puro y va en la discusión.
- **VIF a NUTS-2 sin Madrid sube a ~3** (con ella era ~1,7–2,1). ES30 era el punto que rompía la asociación densidad–área–NOx; sin ella las comunidades vuelven al patrón español general. Aceptable pero conviene declararlo.
- El modelo gaussiano ganó por poco (SCE 1,86 frente a 2,27 del esférico); conviene repetir con el esférico como sensibilidad.

---

## Orden de ejecución

```powershell
# Python (con .venv activado)
python py/00_smoke_test.py
python py/01_nuts.py ; python py/02_estaciones.py ; python py/02b_metadatos.py
python py/03_series_no2.py ; python py/04_industrias.py ; python py/04b_emisiones.py
python py/04c_filtrar_industrias.py ; python py/05_eurostat.py
python py/06_s5p_no2.py ; python py/07_corine.py ; python py/07b_normalizar_costa.py
python py/08_verificar.py
```

```r
source("R/10_prepare.R")
source("R/20_variograma_kriging.R")
source("R/21_probabilidad.R")
source("R/40_procesos_puntuales.R")
source("R/50_acoplamiento.R")
source("R/60_areal.R")
source("R/99_verificar.R")  # integridad: debe dar 0 fallos
```

---

## Lo que falta

1. Escribir `70_sintesis.R` (objetivo 5): tabla de qué aporta cada rama que las otras no pueden dar. El argumento central ya está montado — ver la curva escala–señal en las limitaciones.
2. Revisar la columna `NA 0` del barrido de anchos de banda en `50_acoplamiento.R`: aparece en las siete filas, no impidió nada, pero es un estadístico que no se está calculando.
3. Opcional: `30_spde_inla.R` para la versión bayesiana con SPDE.
4. Opcional: repetir el variograma con modelo esférico como sensibilidad (el gaussiano ganó por poco).
5. Redactar el informe.

---

## Bugs detectados y corregidos

Cadena regenerada por completo el 01/08/2026 y verificada con `R/99_verificar.R`: **0 fallos, 0 avisos**.

### El fallo original NO era de código, era de orden de ejecución

Diagnóstico definitivo, confirmado por las fechas de modificación:

| Fichero | Modificado |
|---|---|
| `data/raw/eprtr_emisiones.csv` | 31/07 19:19:48 |
| `data/processed/nuts3.rds` | 31/07 19:24:24 |
| `data/processed/nuts3_exposicion.rds` | **31/07 18:27:45** |

`nuts3_exposicion.rds` se generó **casi una hora antes que su propia fuente**, cuando `nuts3.rds` aún era la versión anterior a la descarga de emisiones. `nuts3.rds` tenía las 122 392,5 t correctas; el derivado tenía 0. `10_prepare.R` nunca escribió el cero.

Se descartaron por el camino dos hipótesis equivocadas que llegaron a quedar escritas aquí: que `col_emis <- "InspireSiteId"` estaba mal (se usa correctamente como clave del join, `by = setNames(col_emis, col_geo)`) y que la agregación se hacía por atributo en vez de geométricamente. Ninguna era el problema.

**Moraleja operativa:** un fallo de proceso se previene regenerando derivados y comprobando frescura, no parcheando funciones. De ahí el bloque 1 de `99_verificar.R`, que compara `file.mtime()` de cada derivado contra sus fuentes.

### Defectos latentes corregidos

Ninguno causó el fallo observado, pero todos podían causarlo:

1. **`sum(x, na.rm = TRUE)` devuelve 0 sobre un vector todo-NA**, no NA. Si `eprtr_emisiones.csv` faltara, las 50 provincias recibirían cero sin aviso. Corregido con la bandera `HAY_MARCAS`, que propaga NA y aborta si todas las provincias salen a cero pese a haber marcas.
2. **Detección de la columna de masa por expresión regular.** El CSV trae `nox_kg`, `nox_kg_max` y `nox_t`; el patrón laxo eligió `nox_kg` y el total nacional salió en 1,2 × 10⁸ t. Ahora se exige el nombre exacto `nox_t` y se falla con la lista de candidatas. Añadido un **guardia de plausibilidad**: si el total supera 5 × 10⁵ t (España declara ~1 × 10⁵), el script para.
3. **El join de emisiones podía duplicar filas de `ied`.** Corregido colapsando `emis` por sitio antes de unir (316 filas → 316 sitios únicos, así que no había duplicados reales), con `stopifnot` sobre `nrow(ied)`.
4. **`pp` se deduplicaba pero `ied` no**, así que `npoints(pp)` y `nrow(ied)` divergían — el origen de la discrepancia 296-299 vs 300. Ahora se deduplica `ied` por coordenada sumando masas (3 colapsadas → **297**) y el `ppp` se construye limpio.
5. **La ventana simplificada a 1 km dejaba fuera una instalación portuaria**, que spatstat rechazaba en silencio. Se dilata 1 km si detecta puntos fuera. La ventana pasa de 498 522 a 504 327 km², un **1,2 % más** — afecta a la intensidad del proceso puntual en esa proporción.
6. **`ppm()` rechaza patrones marcados no multitipo.** Al empezar a llegar marcas de verdad, `40_procesos_puntuales.R` falló. Corregido con `pp_u <- unmark(pp)` en todo el análisis de intensidad (16 usos) y `pp` marcado conservado solo en el bloque H, que construye la superficie de presión. La distinción es conceptual: la intensidad de primer orden modela **dónde** hay instalaciones, no cuánto emiten.
7. **`st_join` sin `join = st_within`** usaba `st_intersects`: un punto sobre una frontera se asignaba a dos provincias. Corregido en instalaciones y estaciones.
8. **`nuts_id`/`NUTS_ID` y `nuts_name`/`NUTS_NAME` duplicadas** en distinta caja, de las claves de `densidad_nuts3_es.csv` y `clc_fracciones_nuts3.csv`. Se comprueba `identical()` y se elimina la redundante.
9. **El bloque de correlaciones se rompía con variables constantes.** `cor(..., "spearman")` sobre un vector constante devuelve `NA`, y `if (abs(M[i,j]) > 0.7)` recibe `NA`. Corregido con `isTRUE()` y filtro previo.
10. **`if` sin llaves** en el aviso de cobertura de marcas: el segundo `msg` se ejecutaba siempre.
11. **`log(presion_ind + 1e-14)`** mandaba a −32 las provincias con presión nula. Cambiado a `log1p`.
12. **Baleares entra como TRES unidades NUTS-3.** El parche que las enlazaba a sus dos vecinos más próximos las conectaba entre sí y dejaba el archipiélago aislado igual (4 componentes). Sustituido por grafo **k = 5, una sola componente**; contigüidad reina conservada como sensibilidad.
13. **La agregación a NUTS-2 promediaba variables ya transformadas.** La media ponderada de `log_area` no es el log del área total. Se agrega en escala cruda y se re-transforma.
14. **`lm.LMtests` deprecada** desde spdep 1.3 → `lm.RStests`.
15. **`mean(no2_media[es_fondo])`**: si `es_fondo` tiene algún NA el subíndice lógico devuelve filas NA. Cambiado a `which(es_fondo)`.
16. **`Sys.glob` buscaba `ied*.rds`/`instalaciones*.rds`** pero el fichero se llama `industrias.rds`. Sustituido por ruta directa.

### Qué cambió al regenerar

El objetivo 1 reproduce **exactamente** (no dependía de las marcas). Los objetivos 3 y 4 se mueven ~10 %, coherente con haber colapsado 3 instalaciones duplicadas:

| | Antes | Ahora |
|---|---|---|
| Instalaciones | 296-299 | 297 |
| Escala de Thomas | 2,4 km | 2,7 km |
| Diámetro del distrito | 9,5 km | 10,7 km |
| Conglomerados | 551 | 520 |
| psib | 0,927 | 0,915 |
| Suelo industrial (IQR) | ×3,58 | ×3,94 |
| F acoplamiento | 14,86 | 13,86 |
| Efecto de *e* sobre NO₂ | 6,6 % | 6,35 % |
| **Razón escala/diámetro** | **1,05** | **0,90** |

El objetivo 2 reproduce sin cambios: lambda 0,446 idéntica, `log_nox` no significativa en todas las especificaciones, zonificación con `log_dens_pob` al 100 % y `log_nox` al 2,5 %.

**Nombres reales** (para no volver a inventarlos): respuesta `no2_areal` (también `no2_estaciones`, `no2_fondo`, `no2_poblacion`, `no2_p90`, `no2_max`, `prob_media_20`, `brecha_exposicion`); densidad `dens_2023` (verificado: `log_dens_pob == log(dens_2023)`); población `poblacion`; masa industrial `nox_t` (**no** `nox_kg`); fichero de puntos `data/processed/industrias.rds`.

---

## Fuentes y trazabilidad

| Variable | Servicio exacto | Registro |
|---|---|---|
| NUTS 2021 | EEA `EPRTR/NUTS_RG_01M_2021_3035_EU27/MapServer` (solo niveles 0 y 3; NUTS-2 por disolución) | no |
| NO₂ por estación | EEA `ParquetFile/async`, dataset 2 (E1a), pollutant URI `.../pollutant/8`, sin `aggregationType` | no |
| Clasificación de estaciones | DiscoData `[AirQualityDataFlows].[v1r1].[AirQualityStatistics]` | no |
| Instalaciones IED | ArcGIS `Air/IED_SiteMap/MapServer` (MaxRecordCount 1000, paginar) | no |
| Emisiones NOx | DiscoData `[IED].[v1r2]`, unión de PollutantRelease + ProductionFacilityReport + ProductionFacility + ProductionSite. Código del contaminante: `NOX` | no |
| S5P NO₂ | CDSE openEO, colección `SENTINEL_5P_L2`, una banda por petición | **sí** |
| CORINE CLC2018 | `land.copernicus.eu`, ráster 100 m V2020_20u1 | **sí** |
| Densidad poblacional | Eurostat `demo_r_d3dens` (JSON-stat) | no |

Licencias: EEA y Eurostat bajo CC-BY 4.0; Sentinel bajo el aviso legal de Copernicus.
