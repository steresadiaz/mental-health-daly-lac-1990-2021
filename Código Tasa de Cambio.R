##Código con Tasa de Cambio 
# ============================================================
#  SALUD MENTAL EN AMÉRICA LATINA Y EL CARIBE — PIPELINE R
#  Sofí - Tesis Anáhuac  
# ============================================================

# -------- 0. PAQUETES ---------------------------------------------------------
to_install <- c(
  "tidyverse", "data.table", "janitor", "readr", "readxl", "here", "glue",
  "broom", "ggcorrplot", "car", "DescTools", "rstatix", "lme4", "lmerTest",
  "performance", "parameters", "effectsize", "MASS", "robustbase", "sandwich",
  "modelsummary", "segmented", "splines", "rsample", "boot", "lavaan", "semPlot",
  "sf", "rnaturalearth", "rnaturalearthdata", "spdep", "tmap", "viridis", "patchwork", 
  "countrycode", "ragg", "wbstats", "dplyr", "stringr", "moments", "pheatmap", 
  "grid", "ggrepel", "boot", "lme4"
)

# Instalar los paquetes si no están instalados
new_pk <- setdiff(to_install, rownames(installed.packages()))
if (length(new_pk)) install.packages(new_pk, Ncpus = parallel::detectCores())

# Cargar los paquetes
lapply(to_install, library, character.only = TRUE)


set.seed(12345)
options(dplyr.summarise.inform = FALSE)
theme_set(theme_minimal(base_size = 12))
REF_YEAR <- 2000
BREAKS   <- c(2004)


# -------- 1. RUTAS / PARÁMETROS ----------------------------------------------
dir.create(here("output"), showWarnings = FALSE)

# AJUSTA ESTAS RUTAS A TUS ARCHIVOS
PATH_GBD_WIDE  <- here("data","gbd_latam.csv")          # country, iso3, year, cause, metric, value
PATH_CAUSAS    <- here("data","dalys_causas_latam.csv") # country, iso3, year, cause, dalys
PATH_DET       <- here("data","determinantes_wb_z.csv")   # country, iso3, year, pobreza, gini, violencia, desempleo, educacion, gasto_salud, etc.

# Lista de países soberanos en LAC (WB region = "Latin America & Caribbean")
LAC_ISO3 <- wbstats::wb_countries() %>%
  dplyr::filter(region == "Latin America & Caribbean", !is.na(iso3c)) %>%
  dplyr::pull(iso3c) %>% unique()


# -------- 2. UTILIDADES -------------------------------------------------------
cv <- function(x) sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE)

year_splines <- function(year, breaks = BREAKS) {
  tibble(
    year_c = year - REF_YEAR,
    s1 = pmax(0, year - breaks[1]),
    s2 = pmax(0, year - breaks[length(breaks)])
  )
}

winsorize_vec <- function(x, probs = c(0.01, 0.99)) DescTools::Winsorize(x, probs = probs, na.rm = TRUE)

# -------- 3. CARGA Y LIMPIEZA GBD (LAC) --------------------------------------
message("Cargando GBD LAC…")
stopifnot(file.exists(PATH_GBD_WIDE))
gbd <- fread(PATH_GBD_WIDE) |> clean_names()

gbd <- fread(PATH_GBD_WIDE) |>
  clean_names() |>
  rename(
    country = location,
    value   = val
  )

# Estructura esperada:
# country (chr), iso3 (chr), year (int), cause (chr),
# metric (DALYs/YLLs/YLDs/Prevalence/Incidence/Deaths), value (dbl)
req_cols <- c("country","iso3","year","cause","metric","value")
if (!all(req_cols %in% names(gbd))) stop("gbd_latam.csv debe tener: country, iso3, year, cause, metric, value")


# ---- Completar columna iso3 en tu base GBD a partir de gbd$country ---------
# Por si no existe aún el vector de ALC del WB (nota: no incluye PRI/VIR/BMU)
if (!exists("LAC_ISO3")) {
  LAC_ISO3 <- wb_countries() %>%
    dplyr::filter(region == "Latin America & Caribbean", !is.na(iso3c)) %>%
    dplyr::pull(iso3c) %>% unique()
}

# Diccionario para las etiquetas que mostraste (y equivalentes)
custom_names <- c(
  "Venezuela (Bolivarian Republic of)" = "VEN",
  "Bolivia (Plurinational State of)"   = "BOL",
  "Saint Kitts and Nevis"              = "KNA",
  "Saint Lucia"                        = "LCA",
  "Saint Vincent and the Grenadines"   = "VCT",
  "Dominican Republic"                 = "DOM",
  "Trinidad and Tobago"                = "TTO",
  "Antigua and Barbuda"                = "ATG",
  "Bahamas"                            = "BHS",
  "Barbados"                           = "BRB",
  "Belize"                             = "BLZ",
  "Grenada"                            = "GRD",
  "Suriname"                           = "SUR",
  "Guyana"                             = "GUY",
  "Haiti"                              = "HTI",
  "Jamaica"                            = "JAM",
  "Costa Rica"                         = "CRI",
  "El Salvador"                        = "SLV",
  "Guatemala"                          = "GTM",
  "Honduras"                           = "HND",
  "Mexico"                             = "MEX",
  "Nicaragua"                          = "NIC",
  "Panama"                             = "PAN",
  "Paraguay"                           = "PRY",
  "Peru"                               = "PER",
  "Ecuador"                            = "ECU",
  "Colombia"                           = "COL",
  "Chile"                              = "CHL",
  "Brazil"                             = "BRA",
  "Cuba"                               = "CUB",
  "Dominica"                           = "DMA",
  # Territorios/complementarios en tu lista:
  "Puerto Rico"                        = "PRI",
  "United States Virgin Islands"       = "VIR",
  "Bermuda"                            = "BMU"
)

# Asegura columna iso3
if (!"iso3" %in% names(gbd)) gbd$iso3 <- NA_character_
gbd$iso3 <- as.character(gbd$iso3)

gbd <- gbd %>%
  dplyr::mutate(
    country_raw   = country,
    country_clean = stringr::str_squish(stringr::str_remove(country_raw, "\\s*\\([^)]*\\)\\s*$")),
    iso3 = dplyr::coalesce(
      dplyr::na_if(iso3, ""),
      countrycode::countrycode(country_raw,   "country.name.en", "iso3c", warn = FALSE),
      countrycode::countrycode(country_clean, "country.name.en", "iso3c", warn = FALSE),
      dplyr::recode(country_raw,   !!!custom_names, .default = NA_character_),
      dplyr::recode(country_clean, !!!custom_names, .default = NA_character_)
    )
  )



# QA: nombres que aún no mapearon
faltan <- gbd %>% filter(is.na(iso3)) %>% distinct(country_raw = country) %>% arrange(country_raw)
if (nrow(faltan)) {
  message("⚠️ Países/territorios sin ISO3 tras el mapeo (añádelos al diccionario):")
  print(faltan, n = Inf)
} else {
  message("✅ Todas las filas tienen ISO3.")
}

# (Opcional) quedarte solo con ALC + territorios que tú sí necesitas
extra_iso3 <- c("PRI","VIR","BMU")  # los de tu lista que no siempre están en ALC del WB
keep_iso3  <- unique(c(LAC_ISO3, extra_iso3))
gbd    <- gbd %>% filter(iso3 %in% keep_iso3)
# Filtra a LAC
gbd <- gbd |> filter(iso3 %in% LAC_ISO3)



message("Validando DALYs = YLLs + YLDs (solo casos completos y DALYs>0)…")

# 1) Identificadores = todas las columnas salvo metric y value
id_cols <- setdiff(names(gbd), c("metric","value"))

# 2) Reetiquetar métricas largas → cortas
gbd_eval <- gbd %>%
  mutate(metric_short = case_when(
    str_starts(metric, "DALYs") ~ "DALYs",
    str_starts(metric, "YLDs")  ~ "YLDs",
    str_starts(metric, "YLLs")  ~ "YLLs",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(metric_short)) %>%
  # 3) Pivot sin rellenar con ceros (queremos NAs para poder filtrar)
  pivot_wider(
    id_cols    = all_of(id_cols),
    names_from = metric_short,
    values_from = value,
    values_fill = NA_real_,                           # ← importante
    values_fn   = ~ if (all(is.na(.x))) NA_real_ else sum(.x, na.rm = TRUE)
  ) %>%
  # 4) Solo casos completos y con DALYs > 0
  tidyr::drop_na(DALYs, YLLs, YLDs) %>%
  filter(DALYs > 0) %>%
  # 5) Chequeo
  mutate(
    diff_abs = abs(DALYs - (YLLs + YLDs)),
    diff_rel = diff_abs / DALYs
  )

consistencia <- sum(gbd_eval$diff_abs, na.rm = TRUE)
message(sprintf("Diferencia absoluta total (casos completos, DALYs>0): %.6g (ideal≈0)", consistencia))

# Top discrepancias
gbd_eval %>%
  arrange(desc(diff_abs)) %>%
  dplyr::select(any_of(c("country","iso3","year","cause","DALYs","YLLs","YLDs","diff_abs","diff_rel"))) %>%
  head(10)

# Tolerancia: relativa al tamaño de DALYs
gbd_eval <- gbd_eval %>%
  mutate(
    diff_abs = abs(DALYs - (YLLs + YLDs)),
    diff_rel = diff_abs / pmax(DALYs, 1e-12),
    ok = diff_rel <= 1e-8 | diff_abs <= 1e-6   # tolerancias seguras
  )

# 1) Afirmar que todo pasa
stopifnot(all(gbd_eval$ok))


# -------- 4. MÉTRICAS PRINCIPALES --------------------------------------------
unique(gbd$metric)
gbd <- gbd %>% dplyr::mutate(metric = stringr::str_replace(metric, "^(DALYs|YLDs|YLLs).*", "\\1"))

message("Construyendo tablas por indicador LAC…")
dalys_total <- gbd |> filter(metric == "DALYs") |>
  group_by(country, iso3, year) |>
  summarise(DALYs = sum(value, na.rm=TRUE), .groups="drop")
colnames(dalys_total)


# 1. Filtrar los años de interés (2004 y 2021)
dalys_filtrado <- dalys_total %>%
  filter(year %in% c(2004, 2021))

# 2. Pivotar (hacer ancho) el dataframe para tener los DALYs de 2004 y 2021 en columnas separadas
dalys_ancho <- dalys_filtrado %>%
  dplyr::select(country, year, DALYs) %>%
  pivot_wider(
    names_from = year,
    values_from = DALYs,
    names_prefix = "DALYs_"
  )

# 3. Calcular la tasa de cambio porcentual
tasa_cambio_2004_2021 <- dalys_ancho %>%
  # Asegurarse de que no haya divisiones por cero (aunque es poco probable para DALYs)
  filter(!is.na(DALYs_2004) & DALYs_2004 != 0) %>%
  mutate(
    # Diferencia entre los años
    cambio_absoluto = DALYs_2021 - DALYs_2004,
    # Tasa de cambio porcentual
    tasa_cambio_porcentaje = (cambio_absoluto / DALYs_2004) * 100
  ) %>%
  # Seleccionar las columnas relevantes para la salida
  dplyr::select(country, DALYs_2004, DALYs_2021, tasa_cambio_porcentaje)

# 4. Mostrar el resultado
print(tasa_cambio_2004_2021)

# --- PASO PREVIO (RE-CÁLCULO DE LA TASA DE CAMBIO 2004-2021) ---
# Asumo que 'dalys_total' está cargado.
dalys_filtrado <- dalys_total %>%
  filter(year %in% c(2004, 2021))

dalys_ancho <- dalys_filtrado %>%
  dplyr::select(country, iso3, year, DALYs) %>%
  pivot_wider(
    names_from = year,
    values_from = DALYs,
    names_prefix = "DALYs_"
  )

# Crear el nuevo dataframe 'roc_total' con la Tasa de Cambio como la nueva variable de interés
roc_total <- dalys_ancho %>%
  # El DALYs_2004 se convierte en el denominador. Filtramos NA y 0s.
  filter(!is.na(`DALYs_2004`), `DALYs_2004` != 0) %>%
  mutate(
    # Variable clave: Tasa de Cambio Porcentual (RoC) 2004 a 2021
    RoC_2004_2021 = ((`DALYs_2021` - `DALYs_2004`) / `DALYs_2004`) * 100
  ) %>%
  # Seleccionamos las columnas relevantes para este nuevo análisis
  dplyr::select(country, iso3, RoC_2004_2021)

# El nuevo 'dataset' tiene una fila por país. La variable clave es RoC_2004_2021.
print(summary(roc_total$RoC_2004_2021))
# -----------------------------------------------------------------------------
# Normalidad (Shapiro)
message("Normalidad (Shapiro) para la Tasa de Cambio (RoC)...")
safe_shapiro <- function(x) { 
  if(length(na.omit(x))>=3 && length(na.omit(x))<=5000) shapiro.test(x)$p.value else NA_real_ 
}

# Aplicar a la columna de Tasa de Cambio
shap_roc <- safe_shapiro(roc_total$RoC_2004_2021)
print(paste("p-valor de Shapiro-Wilk (RoC):", format.pval(shap_roc, digits = 4)))

### Si p(Shapiro) < 0.05 ⇒ Rechazar H0 (No es normal)

# CV e IQR de la RoC a través de los países
message("CV e IQR para la Tasa de Cambio (RoC) entre países...")
variab_roc <- roc_total %>%
  summarise(
    mean = mean(RoC_2004_2021, na.rm=TRUE),
    sd = sd(RoC_2004_2021, na.rm=TRUE),
    cv = sd/mean,
    iqr = IQR(RoC_2004_2021, na.rm=TRUE)
  )

print(variab_roc)

# --- Histograma con densidad normalizada (RoC) ----------------------------------
p_hist_roc <- ggplot(roc_total, aes(x = RoC_2004_2021)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 15, # Ajustar bins para un número menor de datos (países)
                 fill = "#2ecc71", 
                 color = "white", 
                 alpha = 0.8) +
  geom_density(color = "#3498db", 
               size = 1.2, 
               linetype = "solid", 
               alpha = 0.9) +
  labs(
    title = "Distribución de la Tasa de Cambio (RoC) 2004-2021",
    x = "Tasa de Cambio (%) de DALYs (2004 a 2021)",
    y = "Densidad"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p_hist_roc)

# --- Boxplot de la RoC ---------------------------------------------------------
# Se usa 'factor(1)' para que ggplot sepa que es una sola variable a graficar
p_box_roc <- ggplot(roc_total, aes(x = factor(1), y = RoC_2004_2021)) +
  geom_boxplot(fill="#9b59b6", color="black", outlier.shape=16, outlier.alpha=.8) +
  geom_jitter(width = 0.1, alpha = 0.6) + # Mostrar los puntos (países)
  labs(
    title = "Boxplot de la Tasa de Cambio (RoC) 2004-2021",
    x = "", 
    y = "Tasa de Cambio (%) de DALYs (2004 a 2021)"
  ) +
  scale_x_discrete(labels=NULL) + # Remover la etiqueta del eje X
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

print(p_box_roc)


# --- 1. Estadísticos básicos -------------------------------------------------
message("Descripción integral de la Tasa de Cambio (RoC) 2004-2021...")

desc_roc <- roc_total %>%
  summarise(
    n         = n(),
    media     = mean(RoC_2004_2021, na.rm=TRUE),
    mediana   = median(RoC_2004_2021, na.rm=TRUE),
    sd        = sd(RoC_2004_2021, na.rm=TRUE),
    min       = min(RoC_2004_2021, na.rm=TRUE),
    max       = max(RoC_2004_2021, na.rm=TRUE),
    q25       = quantile(RoC_2004_2021, 0.25, na.rm=TRUE),
    q75       = quantile(RoC_2004_2021, 0.75, na.rm=TRUE),
    iqr       = IQR(RoC_2004_2021, na.rm=TRUE),
    # Usamos pull y as.double para asegurar que moments lo reciba correctamente
    asimetria = moments::skewness(pull(roc_total, RoC_2004_2021), na.rm=TRUE), 
    curtosis  = moments::kurtosis(pull(roc_total, RoC_2004_2021), na.rm=TRUE)
  )

print(desc_roc)

# --- 3. QQ-plot (RoC) --------------------------------------------------------
# Crear un vector de la variable para el plot
roc_vector <- na.omit(roc_total$RoC_2004_2021) 

# Mostrar en consola/Plots
qqnorm(roc_vector, main="QQ-plot de la Tasa de Cambio (RoC)")
qqline(roc_vector, col="red")


# --- 1. Definir el vector de ISO3 de LAC (ejemplo, ajusta si es diferente) ---
# Asumo que esta lista se define en otra parte de tu script original:
LAC_ISO3 <- c("ARG", "BOL", "BRA", "CHL", "COL", "CRI", "CUB", "DOM", 
              "ECU", "SLV", "GTM", "HND", "JAM", "MEX", "NIC", "PAN", 
              "PRY", "PER", "URY", "VEN") # Añade más si aplica (e.g., HTI, TTO, BHS, etc.)

# --- 2. Definición de la Función get_lac_sf() ---
get_lac_sf <- function() {
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    st_make_valid() |>
    janitor::clean_names() |>
    mutate(iso3 = iso_a3) |>
    filter(iso3 %in% LAC_ISO3)
}

# --- 3. Definición del Tema Base para TMAP ---
base_theme <- tm_layout(
  frame = FALSE,
  fontfamily = "Times New Roman",
  title.fontface = "bold",
  title.fontfamily = "Times New Roman",
  legend.outside = TRUE,
  legend.outside.position = "right",
  legend.stack = "vertical",
  outer.margins = c(0.02, 0.25, 0.02, 0.02)
)

# --- 4. Cargar el mapa base y unir la RoC ---
# Asumo que 'roc_total' está disponible y contiene las columnas 'iso3' y 'RoC_2004_2021'.
lac_sf <- get_lac_sf()

map_df_roc <- lac_sf %>%
  left_join(roc_total, by="iso3")

tm_roc_corregido <- tm_shape(map_df_roc) +
  tm_polygons(
    fill = "RoC_2004_2021",
    # Mover el manejo de NA fuera de tm_scale_continuous()
    col.na = "lightgrey", 
    # La escala debe definirse dentro de tm_scale() en v4.0
    fill.scale = tm_scale_continuous(
      values = "RdYlBu",
      midpoint = 0
    ), 
    fill.legend = tm_legend(
      title = glue::glue("Tasa de Cambio (%) RoC\n(2004 a 2021)"),
      frame = TRUE, bg.color = "white", bg.alpha = 1
    ),
    col = "black", 
    lwd = 0.5 
  ) +
  tm_title(glue::glue("Tasa de Cambio de DALYs en LAC (2004–2021)")) +
  base_theme

print(tm_roc_corregido)


# --- CORRECCIÓN: Permitir países sin vecinos (zero.policy = TRUE) ---
# 1. Definir la matriz de pesos espaciales (vecinos, ej. Queen)
nbq <- poly2nb(lac_roc_lisa, queen = TRUE) 
# Esta vez, pasa zero.policy = TRUE
lwq <- nb2listw(nbq, style = "W", zero.policy = TRUE) 

# El error de 'Empty neighbour sets' debería desaparecer.

# 2. Calcular el LISA (Ahora debería funcionar)
lisa <- localmoran(lac_roc_lisa$RoC_2004_2021, lwq, zero.policy = TRUE)
lisa <- as.data.frame(lisa)
colnames(lisa) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr.z")

# 3. Clasificar los Cuadrantes de LISA (Clusters)
x  <- scale(lac_roc_lisa$RoC_2004_2021)[,1]
# Volver a calcular wx con zero.policy = TRUE
wx <- lag.listw(lwq, x, zero.policy = TRUE) 

lac_roc_lisa <- lac_roc_lisa %>%
  mutate(
    Zi = lisa$Z.Ii,
    Pi = lisa$Pr.z,
    quadrant = case_when(
      x >= 0 & wx >= 0 & Pi < 0.05 ~ "Alto-Alto (Aumento clusterizado)",
      x < 0 & wx < 0 & Pi < 0.05   ~ "Bajo-Bajo (Disminución clusterizada)",
      x >= 0 & wx < 0 & Pi < 0.05  ~ "Alto-Bajo (Aumento aislado)",
      x < 0 & wx >= 0 & Pi < 0.05  ~ "Bajo-Alto (Disminución aislada)",
      TRUE                         ~ "No significativo (p ≥ 0.05)"
    )
  )

# 4. Crear el Mapa de Clusters LISA (Ajustado a tmap v4)
tm_lisa_roc <- tm_shape(lac_roc_lisa) +
  tm_polygons(
    fill = "quadrant",
    # La paleta ahora va dentro de fill.scale = tm_scale_categorical() para variables factor
    fill.scale = tm_scale_categorical(
      values = c("Alto-Alto (Aumento clusterizado)" = "#CC0000",
                 "Bajo-Bajo (Disminución clusterizada)" = "#0066CC",
                 "Alto-Bajo (Aumento aislado)" = "#FF9999",
                 "Bajo-Alto (Disminución aislada)" = "#99CCFF",
                 "No significativo (p ≥ 0.05)" = "#CCCCCC")
    ),
    fill.legend = tm_legend(
      title = "Clusters LISA (p<0.05)",
      show = TRUE
    )
  ) +
  tm_title("Clusters LISA de Tasa de Cambio de DALYs (2004-2021)") +
  base_theme

print(tm_lisa_roc)

# --- 1. Cargar y preparar determinantes ---------------------------------------
message("Calculando Tasa de Cambio (RoC) 2004-2021 para cada determinante...")

# Cargar determinantes, filtrar LAC y años 2004 y 2021
# stopifnot(file.exists(PATH_DET))
det <- fread(PATH_DET) |> clean_names() |> filter(iso3 %in% LAC_ISO3)

# Identificar las variables numéricas que son determinantes (excluyendo iso3, country, year)
det_vars <- setdiff(names(det)[sapply(det, is.numeric)], c("year"))

# Función para calcular la RoC 2004-2021
calcular_roc_determinante <- function(df, var_name) {
  df |>
    # Filtrar solo 2004 y 2021
    filter(year %in% c(2004, 2021)) |>
    # Asegurarse de que hay datos para ambos años
    tidyr::drop_na(all_of(var_name)) |>
    group_by(country, iso3) |>
    # Debe haber exactamente 2 años
    filter(n() == 2) |>
    summarise(
      # Valor en 2004
      valor_2004 = .data[[var_name]][year == 2004],
      # Valor en 2021
      valor_2021 = .data[[var_name]][year == 2021],
      # Fórmula RoC: (V_final - V_inicial) / V_inicial * 100
      RoC = ((valor_2021 - valor_2004) / valor_2004) * 100,
      .groups = "drop"
    ) |>
    # Renombrar la columna RoC al formato RoC_NOMBRE
    rename_with(~glue("RoC_{var_name}"), .cols = RoC) |>
    dplyr::select(country, iso3, starts_with("RoC_"))
}

# Aplicar la función a todos los determinantes y unirlos
roc_determinantes_list <- map(det_vars, ~calcular_roc_determinante(det, .x))

# Reducir la lista a un solo dataframe unido por 'iso3'
roc_determinantes <- roc_determinantes_list |>
  reduce(full_join, by = c("country", "iso3"))

# --- 2. Crear la base de datos de modelo para RoC -----------------------------
# Asumo que 'roc_total' contiene 'iso3' y 'RoC_2004_2021'
base_modelo_roc <- roc_total |>
  dplyr::select(country, iso3, DALYs_RoC = RoC_2004_2021) |> # Renombramos para claridad
  left_join(roc_determinantes, by = c("country","iso3"))

# Seleccionar las variables candidatas (los RoC calculados)
candidatas_roc <- setdiff(names(base_modelo_roc), c("country", "iso3", "DALYs_RoC"))
message(paste("Variables RoC de Determinantes:", paste(candidatas_roc, collapse = ", ")))

# -------- BIVARIADO (Spearman con determinantes - RoC) --------------------------
message("\nCorrelaciones Spearman (LAC) con RoC de determinantes...")

cors_roc <- map_dfr(candidatas_roc, ~{
  # Usar DALYs_RoC como variable dependiente
  ct <- suppressWarnings(cor.test(base_modelo_roc$DALYs_RoC, base_modelo_roc[[.x]], method="spearman"))
  tibble(variable = .x, rho = unname(ct$estimate), p = ct$p.value)
}) |> arrange(p)

print(cors_roc)

# Necesitamos el dataframe con todas las variables RoC
vars_roc_det <- base_modelo_roc %>% dplyr::select(starts_with("RoC_")) %>% names()

# Creamos un dataframe con las correlaciones y el conteo de NAs
variable_ranking <- cors_roc %>%
  left_join(na_counts, by = c("variable")) %>%
  arrange(p) %>%
  # Añadimos una columna con el N disponible (N_total - n_na)
  mutate(N_disp = 35 - n_na)

candidatas_final_select <- c(
  "RoC_internet_usuarios",
  "RoC_vac_measles",
  "RoC_desempleo",
  "RoC_daly_hiv",
  "RoC_urbanizacion"
)

# --- 2.1. Preparar datos para OLS (solo con las 5 variables clave) -----------------
vars_ols_final <- c("DALYs_RoC", candidatas_final_select)

# Filtramos el dataframe original para tener solo las filas donde estas 6 variables están completas
ols_data_final <- base_modelo_roc |>
  dplyr::select(all_of(vars_ols_final)) |>
  na.omit() 

N_Paises_Final_Opt <- nrow(ols_data_final) # Debería ser 18
P_Predictoras_Final_Opt <- length(candidatas_final_select) # Debería ser 5

message(glue::glue("Filas usadas para OLS (N): {N_Paises_Final_Opt}"))
message(glue::glue("Variables predictoras finales (P): {P_Predictoras_Final_Opt}"))

# --- 2.2. Definir y ejecutar el modelo ------------------------------------------
form_ols_final <- as.formula(glue::glue("DALYs_RoC ~ {paste(candidatas_final_select, collapse=' + ')}"))

m_ols_final <- lm(form_ols_final, data = ols_data_final)

# --- 2.3. Selección Stepwise (para encontrar el modelo más parsimonioso) --------
message("Iniciando Selección Stepwise para OLS...")
modelo_parsimonioso_opt <- step(m_ols_final, direction = "both", trace = 0)

# --- 2.4. Reporte Final OLS -----------------------------------------------------
cat("\n======================================================\n")
cat("  REPORTE FINAL OLS (RoC DALYs vs. RoC Determinantes)  \n")
cat("======================================================\n")
print(summary(modelo_parsimonioso_opt))


# --- 3.1. Extraer coeficientes del modelo COMPLETO (m_ols_final) ---
df_ols_plot_opt_full <- broom::tidy(m_ols_final) |>
  mutate(
    # Calcular IC95% (Wald)
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    # Marcar significancia
    sig05 = case_when(
      p.value < 0.001 ~ "p < 0.001 (***)",
      p.value < 0.01  ~ "p < 0.01 (**)",
      p.value < 0.05  ~ "p < 0.05 (*)",
      TRUE            ~ "ns"
    ),
    sig05 = factor(sig05, levels = c("ns", "p < 0.05 (*)", "p < 0.01 (**)", "p < 0.001 (***)"))
  ) |>
  # Aquí filtramos el intercepto, ¡pero ahora sí tendremos variables!
  filter(term != "(Intercept)") |> 
  # Limpiar y etiquetar
  mutate(
    label = gsub("RoC_", "Tasa Cambio ", term),
    label = dplyr::case_when(
      label == "Tasa Cambio internet_usuarios" ~ "Tasa Cambio Usuarios Internet",
      label == "Tasa Cambio vac_measles" ~ "Tasa Cambio Vacunación Sarampión",
      label == "Tasa Cambio desempleo" ~ "Tasa Cambio Desempleo",
      label == "Tasa Cambio daly_hiv" ~ "Tasa Cambio DALYs VIH/SIDA",
      label == "Tasa Cambio urbanizacion" ~ "Tasa Cambio Urbanización",
      TRUE ~ label
    ),
    label = forcats::fct_rev(forcats::fct_inorder(label))
  )

# --- 3.2. Generación del gráfico (Usando el dataframe full) ---
p_forest_ols_opt_full <- ggplot(df_ols_plot_opt_full,
                                aes(x = estimate, y = label, xmin = conf.low, xmax = conf.high, color = sig05)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  geom_pointrange(linewidth = 0.6) +
  scale_color_manual(
    values = c("p < 0.001 (***)" = "#a00574", "p < 0.01 (**)" = "#238b45", "p < 0.05 (*)" = "#ff7b07", "ns" = "#7f7f7f"),
    drop = FALSE,
    name = "Significancia"
  ) +
  labs(
    title = "Forest Plot: Coeficientes del Modelo OLS (RoC DALYs vs. RoC Determinantes)",
    subtitle = glue::glue("Modelo COMPLETO con 5 variables clave. N={N_Paises_Final_Opt} Países."),
    x = "Estimación (Cambio en RoC DALYs por unidad de RoC Predictor)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top", 
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    text = element_text(family = "Times New Roman"),
    title = element_text(family = "Times New Roman", face = "bold")
  )

print(p_forest_ols_opt_full)
ggsave(here::here("output", "forest_ols_roc_roc_full_model.png"), p_forest_ols_opt_full, width = 9, height = 7, dpi = 300)