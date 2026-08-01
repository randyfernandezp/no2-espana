# Esqueleto del informe — NO₂ en España

Documento de trabajo para la redacción. No es el informe: es el argumento que
el informe tiene que sostener, con la asignación de qué figura y qué tabla
respalda cada afirmación.

**Principio rector.** La estructura NO sigue el orden de los objetivos. Un
informe organizado como «objetivo 1, objetivo 2, objetivo 3…» reproduce el
orden en que se hizo el análisis, no el orden en que se entiende, y entierra
el hallazgo central en la cuarta sección. Las ramas entran donde hacen falta
como evidencia, y el mapeo a los cinco objetivos se declara al principio para
que el lector (o el evaluador) pueda localizarlos.

**Extensión orientativa:** 25-35 páginas más anexos. Ocho figuras en el
cuerpo; el resto va a anexo.

---

## Tabla de correspondencia con los objetivos

Se coloca al final de la introducción, en una caja. Permite evaluar el
cumplimiento formal sin distorsionar la narración.

| Objetivo del enunciado | Dónde se resuelve |
|---|---|
| 1 — variación continua, variograma e interpolación | §3 y anexo A |
| 2 — agregación areal y relación con covariables | §5 |
| 3 — patrón de puntos e intensidad modelada | §4 |
| 4 — acoplamiento intensidad ↔ concentración | §4.3 |
| 5 — discusión integrada | §6 y §7 |

---

## §1. La pregunta

Arranque con la tensión real, no con un párrafo de contexto sobre la
contaminación del aire.

**Afirmación que sostiene la sección:** el NO₂ es el único contaminante donde
la pregunta «¿industria o urbanización?» está genuinamente abierta, y España
el único país europeo grande donde ambos efectos no son colineales por
construcción.

- Por qué NO₂ y no PM2.5 u O₃: es el único con producto satelital utilizable
  (S5P/TROPOMI) **y** con firma industrial y de tráfico simultáneas.
- Por qué España: los focos industriales (Huelva, Tarragona, Puertollano, As
  Pontes, Avilés, Algeciras) están espacialmente desacoplados de los grandes
  focos urbanos. En Italia o Polonia serían colineales.
- El coste de esa elección, dicho ya aquí y no escondido en limitaciones:
  59 NUTS-3 es poco para un modelo areal, y hay vacíos de muestreo.

**Sin figuras.** Máximo dos páginas.

---

## §2. Datos y decisiones

Sección corta y honesta. Las decisiones metodológicas que condicionan todo lo
demás van aquí, no dispersas.

**Afirmación:** tres decisiones determinan los resultados y son defendibles.

1. **Solo estaciones de fondo alimentan el modelo geoestadístico.** Las de
   tráfico miden a microescala (el NO₂ cae 60-80 % en los primeros 50 m) y no
   son realizaciones del campo regional. Se reservan como validación externa,
   y el sesgo medido en ellas *cuantifica el recargo local* — que deja de ser
   un problema para convertirse en un resultado.
2. **Dominio peninsular + Baleares.** Las tres ramas presuponen un dominio
   conexo. Canarias se excluye y eso tiene una consecuencia cuantificada que
   hay que declarar aquí: concentra el 24 % del NOx industrial español.
3. **Desagrupamiento por celdas de 25 km**: la media pasa de 11,36 a
   8,58 µg/m³. Sin corregir, la red urbana sesga toda la superficie.

| Elemento | Recurso |
|---|---|
| Tabla de fuentes y trazabilidad | `ESTADO.md` → sección «Fuentes» |
| Mapa de estaciones por tipo | figura de `20_` (confirmar nombre) |

---

## §3. El campo de fondo

**Afirmación:** existe un campo regional suave, de rango 71 km, y las
covariables satelitales y de uso de suelo explican tres cuartas partes de su
varianza.

- Variograma gaussiano, rango práctico **71 km**, pepita **9,6 %**.
- Las covariables absorben el **76 %** de la varianza (rango residual 24 km).
- KDE frente a KO: R² **0,776 vs 0,481**, RMSE **−34,3 %**, sesgo 0,16 vs
  0,93, cobertura del intervalo al 95 %: **0,91**.
- **El resultado que importa para lo que viene:** el recargo local medido en
  las estaciones excluidas. Tráfico 2,9-4,8 µg/m³, industrial 0,5-1,5. En
  entorno rural las industriales miden un 52 % más que las de fondo (5,95 vs
  3,91); en urbano, menos. *La industria eleva el NO₂ donde el fondo es bajo,
  pero queda sepultada por el tráfico en las ciudades.*

Esa última frase es la primera aparición del argumento central. Conviene que
el lector la retenga.

| Contenido | Recurso |
|---|---|
| Variograma ajustado | figura de `20_` |
| Superficie KDE | figura de `20_` |
| Validación cruzada KO vs KDE | tabla de `20_` |
| Recargo por tipo de estación | tabla de `20_` |

**Advertencia a incluir:** el propio ajuste señala variación por debajo de la
menor distancia entre estaciones. Es un supuesto poco creíble en calidad del
aire y hay que decirlo, no ocultarlo.

---

## §4. Dónde están las fuentes y hasta dónde llegan

El corazón metodológico. Dos subsecciones que preparan la tercera.

### §4.1 El patrón industrial no es aleatorio

**Afirmación:** las instalaciones se agrupan, y no solo porque se agrupe la
población.

- 297 instalaciones. Test de cuadrantes Monte Carlo X² = 531,64, **p = 0,002**.
- Intensidad modelada: suelo industrial **×3,94** por rango intercuartílico,
  superficie artificial ×1,32, transporte ×0,748, y **densidad de población
  ×0,999, p = 0,997 — no significativa**.
- Ese último dato es más interesante de lo que parece y merece un párrafo
  propio: *dónde vive la gente no predice dónde se pone la industria*, una vez
  se controla por uso del suelo. Es la primera evidencia directa contra la
  hipótesis de colinealidad total.

### §4.2 Quedan distritos

**Afirmación:** tras condicionar en todas las covariables observadas, sigue
habiendo agrupamiento — hay un mecanismo no capturado.

- Agregación residual hasta **38 km**.
- Modelo de Thomas: escala 2,7 km, **diámetro efectivo del distrito 10,7 km**,
  520 conglomerados, psib 0,915.
- Interpretación: infraestructura compartida, cadenas de suministro, herencia
  histórica. No es ruido, es estructura.

### §4.3 A qué distancia se nota

**Afirmación:** la influencia industrial sobre la concentración tiene una
escala física estimable, y son 10 km.

- Barrido de siete anchos de banda (5 a 100 km), curva de RMSE en U limpia.
- **Escala óptima 10 km.** F = 13,86, p = 2,5 × 10⁻⁴. R² ajustado
  0,754 → 0,769. Multiplicar por *e* la presión industrial sube el NO₂
  un 6,35 %.
- **Y la mejora del RMSE es solo del 0,41 %.** Decirlo aquí, no esconderlo:
  significativa pero de escaso poder predictivo. Coherente con que el 24 % del
  NOx industrial español esté fuera del dominio.

| Contenido | Recurso |
|---|---|
| Intensidad KDE con instalaciones | figura de `40_` |
| K inhomogénea residual | `26_Kinhom_residual.png` |
| Superficie de presión de emisión | `28_presion_emision.png` |
| Coeficientes del ppm | `40_coeficientes_ppm.csv` |
| Barrido de anchos de banda | tabla de `50_` |

---

## §5. La misma pregunta a escala administrativa

**Afirmación:** a escala provincial la industria desaparece, y no por casualidad.

Estructura de la sección como una investigación, no como un informe de
resultados: se plantea el modelo, se obtiene un no-resultado, y se somete ese
no-resultado a cuatro intentos de refutación.

1. **Especificación.** Respuesta `no2_areal`, no `no2_poblacion` — la
   exposición ponderada correlaciona 0,948 con la densidad *por construcción*,
   y regresarla sobre densidad sería circular. Explicar esto bien: es un error
   frecuente y detectarlo es un punto a favor.
2. **Colinealidad y el confusor de escala.** El bloque antrópico va junto a
   ρ 0,77-0,97. El área es el eje común: las provincias pequeñas son densas,
   urbanas, costeras e industriales a la vez. `log_nox` da −0,046 con el
   área — es la única medida industrial que no es alias de urbanización.
3. **El no-resultado.** `log_nox` p = 0,63. Y los cuatro intentos de tumbarlo:
   con `log_presion` en lugar de emisiones (p = 0,36), sin Madrid (p = 0,90),
   bajo SEM (p = 0,68), y en 1000 zonificaciones alternativas (significativa
   en el 2,5 %, por debajo del nominal).
4. **Dependencia espacial.** Moran de la respuesta **no** significativo
   (I = 0,013, p = 0,32); sobre residuos sí (I = 0,187, p = 0,0019). La
   dependencia aparece *después* de condicionar. SEM con lambda = 0,446.
   Reportar SEM y no SDM pese al AIC: ΔAIC = 0,25 es indistinguible, los
   contrastes robustos separan limpiamente, el rho del SAR no es significativo
   y ningún impacto del Durbin lo es.
5. **MAUP.** Escala y zonificación, los dos componentes. A NUTS-2 aparece un
   β negativo significativo que es íntegramente Madrid (Cook = 4,64) y que
   queda fuera del intervalo del 95 % generado por 1000 particiones del mismo
   tamaño. Madrid **nunca** queda aislada en esas réplicas: NUTS-2 la aísla por
   razones históricas, no geométricas.

| Contenido | Recurso |
|---|---|
| LISA | `41_lisa.png` (aunque no haya agrupamientos: el resultado negativo es informativo) |
| Comparación de modelos | `60_comparacion_modelos.csv` |
| Sensibilidad a influyentes | `60_influyentes.csv` |
| MAUP zonificación | `43_maup_zonificacion.png`, `60_maup_zonificacion.csv` |

---

## §6. Por qué desaparece

La sección que justifica haber hecho tres análisis en lugar de uno.

**Afirmación 1 — no es falta de potencia, es la unidad de análisis.**

La provincia media mide 9 970 km²: un disco equivalente de 113 km de diámetro.
Un radio de influencia de 10 km cubre el **0,8 %** de ella. Agregar a NUTS-3
promedia la señal industrial contra dos órdenes de magnitud de territorio donde
no ocurre nada. La urbanización sobrevive a la agregación porque *opera* a
escala provincial; la industria no. **La diferencia no está en la intensidad
del mecanismo sino en su alcance espacial.**

**Afirmación 2 — la escala de 10 km no es arbitraria.**

Dos estimaciones independientes convergen:

| Magnitud | Rama | Datos usados | Valor |
|---|---|---|---|
| Diámetro del distrito industrial | Proceso puntual (Thomas) | solo coordenadas de instalaciones | 10,7 km |
| Escala de influencia sobre la concentración | Acoplamiento (barrido) | instalaciones + estaciones | 10,0 km |

**Razón 0,93.** No hay ninguna razón aritmética para que coincidan: usan datos
distintos y responden a preguntas distintas. Que coincidan implica que **la
unidad de influencia sobre el NO₂ es el distrito industrial, no la instalación
aislada**.

Y sobrevivió a una depuración del pipeline que cambió el patrón de puntos:
antes daba 10 / 9,5 = 1,05, después 10 / 10,7 = 0,93. Mencionarlo — la
robustez frente a un cambio no planeado vale más que cualquier bootstrap.

**Afirmación 3 — los mecanismos son inseparables, y eso es un hallazgo.**

Descomposición de varianza: **el 80,3 % es varianza compartida**, no
atribuible a ningún mecanismo por separado. Y la aparente contradicción entre
las dos descomposiciones (LMG da urbanización 84,5 %; bloques da uso de suelo
como aporte único mayor) **hay que explicarla explícitamente**: miden cosas
distintas, y urbanización y uso de suelo son la misma dimensión medida de dos
formas. No dejar las dos tablas juntas sin comentario.

La endogeneidad no es un defecto del análisis: es una propiedad del territorio.
Y el objetivo 3 la cuantifica al mostrar que la densidad de población **no**
predice dónde están las instalaciones.

| Contenido | Recurso |
|---|---|
| Curva escala-señal | `50_curva_escala.png` ← **figura principal del informe** |
| Convergencia de escalas | `70_convergencia.csv` |
| Descomposición por bloques | `51_sintesis_bloques.png`, `60_descomposicion_varianza.csv` |
| Qué aporta cada rama | `70_aportacion_ramas.csv` ← **tabla que responde el objetivo 5** |

---

## §7. Qué significa para la gestión

Cerrar en unidades físicas y en la unidad en que se legisla.

**Afirmación:** la política de calidad del aire en España es un problema de
tráfico urbano, no de industria, y el análisis lo cuantifica.

| Magnitud | Valor |
|---|---|
| Cota superior del efecto industrial areal | **0,96 µg/m³** (22,3 % de la media nacional) |
| Recargo local en estaciones industriales | ~1,0 µg/m³ |
| Recargo local en estaciones de tráfico | ~3,5 µg/m³ |
| Media areal nacional de fondo | 4,28 µg/m³ |

El techo industrial es una afirmación **más fuerte** que una no-significancia:
no es que ignoremos si la industria influye, es que sabemos que no puede
influir más que eso. Insistir en esta distinción — es lo que separa un
resultado nulo bien hecho de uno mal hecho.

- **Brecha de exposición:** Madrid 7,52 → 15,01 µg/m³ al ponderar por
  población. Barcelona 6,64 → 14,21. La media territorial no es la métrica de
  política pública.
- **Provincias con más del 10 % de población sobre 20 µg/m³:** Madrid 37,3 %,
  Barcelona 30,6 %, Bizkaia 22,7 %.
- **El contraste que cierra el informe:** el 43,3 % de las estaciones de
  tráfico ya superan hoy los 20 µg/m³, frente al 0,16 % de la superficie de
  fondo. España cumple hoy con holgura y tiene margen frente al umbral de
  2030, pero está lejos de la recomendación de la OMS, y el problema vive
  por debajo de la resolución a la que este trabajo puede verlo.

---

## §8. Limitaciones

No como trámite. Las que de verdad condicionan la lectura, en orden de
importancia:

1. **Umbral E-PRTR de 100 t/año.** La marca mínima observada es 102 t. El
   patrón representa grandes emisores, no toda fuente industrial.
2. **Canarias fuera del dominio con el 24 % del NOx industrial español.** Por
   eso el efecto peninsular es modesto — es una limitación del diseño, no un
   hallazgo sobre la industria española en conjunto.
3. **La superficie modela NO₂ de fondo, no exposición máxima.** Los picos de
   tráfico ocurren por debajo de la resolución de 2 km. El propio análisis
   señala su propio punto ciego.
4. **CORINE clase 121** mezcla industria con comercio y logística; casi
   colineal con urbanización. La señal industrial limpia es el patrón E-PRTR.
5. **Pesos dasimétricos** sin escalar por población provincial: las medias
   dentro de cada provincia son correctas, el agregado nacional está sesgado.
6. **Envolventes con lambda fijada** de los propios datos: contraste
   ligeramente anticonservador.
7. **La ventana del proceso puntual** se dilató 1 km para no perder
   instalaciones portuarias: +1,2 % de área, sesgo a la baja de la intensidad
   en esa proporción.
8. **59 NUTS-3, 50 en el dominio.** Es poco para un modelo areal, y el
   diagnóstico de influyentes lo confirma.
9. **El modelo gaussiano ganó por poco** (SCE 1,86 vs 2,27 del esférico).

---

## Anexos

- **A.** Detalle geoestadístico: variogramas direccionales, comparación de
  modelos, validación cruzada completa.
- **B.** Diagnósticos del proceso puntual: residuos, correlación de marcas.
- **C.** Reproducibilidad: orden de ejecución, `99_verificar.R`, y la nota
  sobre el fallo de frescura de derivados. Incluirlo es un punto a favor, no
  en contra: demuestra control sobre el pipeline.
- **D.** Tablas completas.

---

## Pendiente antes de redactar

- Confirmar los nombres de las figuras de `20_`, `21_` y `50_`:
  `list.files("output/figuras")`.
- Decidir formato: Quarto (`.qmd`) permite regenerar el informe con los
  números en vivo desde los `.rds` y evita transcribir cifras a mano, que es
  donde se cuelan los errores.
