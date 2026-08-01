# NO₂ en España: escala de la influencia industrial

Análisis espacial del dióxido de nitrógeno de fondo en la España peninsular y
balear (2023), combinando geoestadística, procesos puntuales y modelos areales.

**Resultado principal.** Dos estimaciones metodológicamente independientes
sitúan la escala de la influencia industrial en torno a 10 km: un proceso de
conglomerados de Thomas ajustado *solo a coordenadas* de instalaciones estima
un diámetro de distrito de 10,7 km, y un barrido de anchos de banda que acopla
emisión con concentraciones observadas sitúa la escala óptima en 10,0 km
(razón 0,93). A escala provincial ese efecto es indetectable, y el trabajo
cuantifica por qué.

---

## Requisitos

| Componente | Versión probada |
|---|---|
| R | 4.5.1 |
| Python | 3.12 |
| TeX Live | 2023 o posterior (solo para el paper y las diapositivas) |

Paquetes de R: `sf`, `terra`, `gstat`, `spatstat`, `spdep`, `spatialreg`,
`dplyr`, `readr`, `ggplot2`, `tidyr`, `stringr`, `car`.

```r
install.packages(c("sf","terra","gstat","spatstat","spdep","spatialreg",
                   "dplyr","readr","ggplot2","tidyr","stringr","car"))
```

Paquetes de Python: ver `requirements.txt`.

```bash
python -m venv .venv
# Windows:  .venv\Scripts\Activate.ps1
# Linux/macOS:  source .venv/bin/activate
pip install -r requirements.txt
```

---

## Ejecución

El orden importa. Los guiones de R consumen productos de los de Python, y
entre sí forman una cadena de dependencias estricta.

### 1. Descarga de datos

```bash
python py/00_smoke_test.py          # comprueba conectividad con los servicios
python py/01_nuts.py                # geometrías NUTS 2021
python py/02_estaciones.py          # estaciones de medición
python py/02b_metadatos.py          # clasificación de estaciones
python py/03_series_no2.py          # series de NO2 (dataset E1a)
python py/04_industrias.py          # instalaciones IED
python py/04b_emisiones.py          # emisiones de NOx (E-PRTR)
python py/04c_filtrar_industrias.py # cruce y filtrado por umbral
python py/05_eurostat.py            # densidad de población
python py/06_s5p_no2.py             # Sentinel-5P (requiere registro CDSE)
python py/07_corine.py              # CORINE CLC2018 (requiere registro)
python py/07b_normalizar_costa.py   # corrección de la clase mar
python py/08_verificar.py           # verificación de las descargas
```

**Dos fuentes requieren registro previo:** Copernicus Data Space Ecosystem
(Sentinel-5P) y `land.copernicus.eu` (CORINE). Las credenciales se leen de
variables de entorno; ver `py/README_datos.md`.

### 2. Análisis

```r
source("R/10_prepare.R")             # preparación y control de calidad
source("R/20_variograma_kriging.R")  # objetivo 1: campo continuo
source("R/21_probabilidad.R")        # probabilidad de superar umbrales
source("R/40_procesos_puntuales.R")  # objetivo 3: patrón industrial
source("R/50_acoplamiento.R")        # objetivo 4: escala de influencia
source("R/60_areal.R")               # objetivo 2: modelo areal y MAUP
source("R/70_sintesis.R")            # objetivo 5: síntesis
source("R/99_verificar.R")           # integridad: debe dar 0 fallos
```

### 3. Documentos

```bash
cd paper  && latexmk -pdf paper.tex
cd slides && latexmk -pdf slides.tex
```

---

## Verificación de integridad

`R/99_verificar.R` no comprueba que el análisis sea correcto: comprueba que
los datos que lo alimentan son los que se cree que son. Ejecutarlo tras
cualquier cambio en la cadena.

Verifica tres cosas:

1. **Frescura.** Ningún fichero derivado puede ser más antiguo que sus
   fuentes. Este bloque existe porque el fallo real que se produjo durante el
   desarrollo no fue de código sino de orden de ejecución: un `.rds` derivado
   se generó antes de que se descargaran las emisiones y arrastró ceros
   durante semanas.
2. **Invariantes de conservación.** El NOx asignado a provincias debe igualar
   al de las instalaciones; el número de puntos del patrón debe igualar las
   filas del objeto espacial; la suma de áreas provinciales debe aproximar el
   área del dominio.
3. **Rangos plausibles.** Incluido un control sobre el ranking industrial: si
   Asturias, A Coruña, Tarragona y Huelva no aparecen entre los ocho primeros,
   el cruce espacial ha fallado aunque los totales cuadren.

**Cero fallos no demuestra que el análisis sea correcto. Un fallo sí demuestra
que algo está mal.**

---

## Estructura

```
no2-espana/
├── py/                 descarga y preprocesado
├── R/                  análisis
│   ├── 10_prepare.R
│   ├── 20_variograma_kriging.R
│   ├── 21_probabilidad.R
│   ├── 40_procesos_puntuales.R
│   ├── 50_acoplamiento.R
│   ├── 60_areal.R
│   ├── 70_sintesis.R
│   └── 99_verificar.R
├── data/
│   ├── raw/            descargas sin modificar   (no versionado)
│   ├── interim/        productos intermedios     (no versionado)
│   └── processed/      conjuntos de análisis     (no versionado)
├── output/
│   ├── figuras/        .png
│   └── tablas/         .csv
├── paper/              paper.tex, refs.bib
├── slides/             slides.tex
└── docs/               ESTADO.md (continuidad del proyecto)
```

Los directorios de datos **no se versionan**: se regeneran ejecutando los
guiones de Python. Descargar todo lleva del orden de una hora, dominada por
Sentinel-5P.

---

## Decisiones de diseño

Documentadas en detalle en `docs/ESTADO.md`. Las tres que más condicionan los
resultados:

**Solo estaciones de fondo alimentan el modelo geoestadístico** (212 de 462).
Las de tráfico miden a microescala —el NO₂ cae 60-80 % en los primeros 50 m
desde la calzada— y no son realizaciones del campo regional. Se reservan como
validación externa, y el sesgo medido en ellas cuantifica el recargo local.

**La respuesta del modelo areal es la media de bloque, no la exposición
ponderada por población.** Esta última se construye con la misma distribución
de población que después entra como predictor (ρ = 0,948), de modo que
regresarla sobre densidad sería circular. Se conserva como métrica de política
pública, no como variable dependiente.

**El dominio excluye Canarias, Ceuta y Melilla** porque las tres ramas
presuponen conexidad. Canarias concentra el 24 % del NOx industrial español:
el efecto peninsular estimado es modesto por construcción del dominio, no
necesariamente por debilidad del mecanismo.

---

## Fuentes y licencias

| Variable | Fuente | Licencia |
|---|---|---|
| Geometrías NUTS 2021 | EEA | CC-BY 4.0 |
| NO₂ por estación | EEA (dataset E1a) | CC-BY 4.0 |
| Instalaciones IED y emisiones | EEA DiscoData | CC-BY 4.0 |
| NO₂ troposférico | Copernicus Sentinel-5P | Aviso legal Copernicus |
| Cobertura del suelo | CORINE CLC2018 | Aviso legal Copernicus |
| Densidad de población | Eurostat | CC-BY 4.0 |

El código de este repositorio se publica bajo licencia MIT (ver `LICENSE`).

---

## Cita

Si este trabajo le resulta útil:

```bibtex
@misc{no2espana2026,
  author = {Fernández, Randy},
  title  = {La escala importa: convergencia entre el agrupamiento de fuentes
            industriales y su influencia sobre el NO2 de fondo en España},
  year   = {2026},
  url    = {https://github.com/randyfernandezp/no2-espana}
}
```
