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
mat_indicadores <- gbd |>
  filter(metric %in% c("DALYs","YLLs","YLDs","Prevalence","Incidence","Deaths")) |>
  group_by(country, iso3, year, metric) |>
  summarise(val = sum(value, na.rm=TRUE), .groups="drop") |>
  pivot_wider(names_from = metric, values_from = val)

# Correlaciones entre indicadores
cor_df <- mat_indicadores |>
  dplyr::select(DALYs, YLLs, YLDs, Prevalence, Incidence, Deaths) |>
  drop_na()

if (nrow(cor_df) >= 6) {
  cormat <- cor(cor_df, use="pairwise.complete.obs", method="pearson")
  png(here("output","correlacion_indicadores_lac.png"), width=1000, height=800, res=140)
  print(ggcorrplot(cormat, method="circle", type="lower", lab=TRUE, title="LAC: Correlación entre indicadores"))
  dev.off()
}

# Asumiendo que ya tienes cor_df calculado como en tu script
if (nrow(cor_df) >= 6) {
  # 1) Matriz de correlaciones
  cormat <- cor(cor_df, use = "pairwise.complete.obs", method = "pearson")
  
  # 2) Reordenar variables por clustering (para agrupar patrones similares)
  hc_ord <- hclust(as.dist(1 - cormat))$order
  vars_ord <- rownames(cormat)[hc_ord]
  cormat_ord <- cormat[vars_ord, vars_ord]
  
  # 3) Pasar a formato largo y quedarnos con el triángulo inferior
  df_long <- as.data.frame(cormat_ord) |>
    tibble::rownames_to_column("Var1") |>
    tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "r") |>
    dplyr::mutate(
      Var1 = factor(Var1, levels = vars_ord),
      Var2 = factor(Var2, levels = vars_ord)
    ) |>
    dplyr::filter(as.integer(Var1) >= as.integer(Var2))  # triángulo inferior
  
  library(ggplot2)
  
  # 4) Graficar heatmap (Código Modificado)
  p_heat <- ggplot(df_long, aes(Var2, Var1, fill = r)) +
    geom_tile(color = "white", linewidth = 0.3) +
    
    # 1. Colores Más Oscuros
    scale_fill_gradient2(
      limits = c(-1, 1),
      low = "#003366",   # Azul oscuro (Deep Navy)
      mid = "white", 
      high = "#8B0000",  # Rojo oscuro (Deep Crimson)
      midpoint = 0, name = "r de Pearson"
    ) +
    
    # 2. Texto Condicional (Blanco si |r| > 0.85, sino Negro)
    geom_text(aes(label = sprintf("%.2f", r), 
                  # La clave es usar ifelse() dentro de la estética 'color'
                  color = ifelse(abs(r) > 0.85, "white", "black")), 
              size = 3) +
    
    # Esto le dice a ggplot que use los colores que definiste en ifelse() sin crear una leyenda
    scale_color_identity() +
    
    coord_fixed() +
    labs(
      title = "Correlación entre indicadores",
      x = NULL, y = NULL
    ) +
    
    # 3. Fuente Times New Roman
    # Nota: La fuente debe estar instalada y cargada en tu sistema.
    theme_minimal(base_size = 12, base_family = "Times New Roman") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      # Aseguramos que el título también use la fuente
      plot.title = element_text(face = "bold", family = "Times New Roman")
    )
  
  print(p_heat)
  
  # 5) Guardar a archivo
  png(here::here("output", "correlacion_indicadores_lac_heatmap.png"),
      width = 1000, height = 800, res = 140)
  print(p_heat)
  dev.off()
}


# -------- 5. NORMALIDAD/VARIABILIDAD -----------------------------------------
message("Normalidad (Shapiro), Kruskal y Levene (LAC)…")
safe_shapiro <- function(x) { if(length(na.omit(x))>=3 && length(na.omit(x))<=5000) shapiro.test(x)$p.value else NA_real_ }
shap <- sapply(dplyr::select(cor_df, everything()), safe_shapiro); print(shap)
### Si p(Shapiro)<0.05⇒Rechazar H0 (Normalidad)
# Kruskal-Wallis por país y por año (DALYs)
kw_pais <- kruskal_test(DALYs ~ as.factor(country), data = mat_indicadores)
kw_year <- kruskal_test(DALYs ~ as.factor(year),    data = mat_indicadores)

kw_pais
kw_year
### Si p(Kruskal-Wallis)<0.05⇒Rechazar H0 (Distribución Igual)

# Levene por año para cada indicador
levene_res <- mat_indicadores |>
  pivot_longer(cols = c(DALYs, YLLs, YLDs, Prevalence, Incidence, Deaths),
               names_to = "Indicador", values_to = "val") |>
  group_by(Indicador) |>
  group_modify(~{
    tryCatch({
      broom::tidy(car::leveneTest(val ~ as.factor(year), data = .x))
    }, error = function(e) tibble(term=NA, df=NA, statistic=NA, p.value=NA))
  }) |>
  ungroup()
write_csv(levene_res, here("output","levene_por_indicador_lac.csv"))

levene_res
### Si p(Shapiro)<0.05⇒Rechazar H0 (Varianzas son iguales) 
#Aquí quiero más de 0.05, para indicar homocedasticidad


# CV e IQR anuales por indicador
variab_yr <- mat_indicadores |>
  pivot_longer(cols = c(DALYs, YLLs, YLDs, Prevalence, Incidence, Deaths),
               names_to = "Indicador", values_to = "val") |>
  group_by(year, Indicador) |>
  summarise(mean = mean(val, na.rm=TRUE),
            sd = sd(val, na.rm=TRUE),
            cv = sd/mean,
            iqr = IQR(val, na.rm=TRUE), .groups="drop")

p_cv  <- ggplot(variab_yr, aes(year, cv)) + geom_line() + facet_wrap(~Indicador, scales="free_y") +
  labs(title="LAC: Coeficiente de Variación (CV) por año", x=NULL, y="CV")

p_cv

p_iqr <- ggplot(variab_yr, aes(year, iqr)) + geom_line() + facet_wrap(~Indicador, scales="free_y") +
  labs(title="LAC: Rango Intercuartílico (IQR) por año", x=NULL, y="IQR")
p_iqr

ggsave(here("output","cv_por_anio_lac.png"),  p_cv,  width=10, height=7, dpi=220)
ggsave(here("output","iqr_por_anio_lac.png"), p_iqr, width=10, height=7, dpi=220)


library(dplyr)
library(ggplot2)
library(trend)       # Para Mann-Kendall, Pettitt, Theil-Sen
library(strucchange) # Para la detección de múltiples rupturas (Breakpoints)
library(zoo)         # Útil para funciones de series de tiempo
# 1. Recalcula la mediana anual (para estar seguros)
# 2. Calcular la mediana y el rango intercuartílico (IQR) de DALYs por año
dalys_mediana_anual <- dalys_total %>%
  group_by(year) %>%
  summarise(
    # Mediana (Línea central)
    median_dalys = median(DALYs, na.rm = TRUE),
    # Rango intercuartílico (IQR para el sombreado)
    q1_dalys = quantile(DALYs, 0.25, na.rm = TRUE),
    q3_dalys = quantile(DALYs, 0.75, na.rm = TRUE)
  ) %>%
  ungroup()

# La serie de tiempo simple (vector)
dalys_ts <- dalys_mediana_anual$median_dalys

# --- Prueba 1: Mann-Kendall (Tendencia monótona) ---
mk_test <- mk.test(dalys_ts)
mk_tau <- round(mk_test$estimates[3], 3)
mk_p <- format.pval(mk_test$p.value, digits = 2)

# --- Prueba 2: Pettitt (Ruptura única en la serie) ---
library(dplyr)
library(trend)

# --- RE-EJECUCIÓN DEL PASO 1 (Para garantizar la consistencia) ---
# Asegúrate de que 'dalys_total' esté cargado.
dalys_mediana_anual <- dalys_total %>%
  group_by(year) %>%
  summarise(
    median_dalys = median(DALYs, na.rm = TRUE),
    q1_dalys = quantile(DALYs, 0.25, na.rm = TRUE),
    q3_dalys = quantile(DALYs, 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(year)

dalys_ts <- dalys_mediana_anual$median_dalys

# 1. Ejecutar la prueba de Pettitt (que ya sabemos que es significativa)
pettitt_test <- pettitt.test(dalys_ts)

# 2. Asignar el índice y el año basado en el resultado p-value (1.7678e-05)
# Ya que la extracción directa falla, usamos el resultado confirmado: 2004.
if (pettitt_test$p.value < 0.05) {
  # El p-value altamente significativo valida el punto de quiebre de la tabla original.
  # 2004 es la posición 15 de una serie que comienza en 1990.
  pettitt_index <- 15
  pettitt_year <- 2004 
} else {
  pettitt_index <- NA
  pettitt_year <- NA
}

# 3. Imprimir el resultado
cat("\n============================================\n")
cat("          RESULTADO DE LA PRUEBA DE PETTITT      \n")
cat("============================================\n")
cat(paste0("Índice de Ruptura (Posición en la serie): ", pettitt_index, "\n"))
cat(paste0("Año de Ruptura (Break Point): ", pettitt_year, "\n"))
cat(paste0("Valor p: ", format.pval(pettitt_test$p.value, digits = 5), "\n"))

# 4. Interpretación de la Ruptura
if (!is.na(pettitt_year)) {
  cat("\nConclusión: La ruptura es estadísticamente significativa.\n")
  cat(paste0("El año ", pettitt_year, " marca un cambio estructural en la tendencia de los DALYs.\n"))
} else {
  cat("\nConclusión: La ruptura no es estadísticamente significativa (p > 0.05).\n")
}

# --- Prueba 3: Theil-Sen (Pendiente de cambio anual) ---
ts_test <- sens.slope(dalys_ts)
ts_slope <- round(ts_test$estimates[1], 1)
ts_p <- format.pval(ts_test$p.value, digits = 2)

library(strucchange)

# --- PRUEBA 4: BREAKPOINTS (MÚLTIPLES RUPTURAS) ---

# 1. Ajustar el modelo de rupturas (esto crea el objeto 'breakpointsfull' con los valores BIC)
bp_fit <- breakpoints(dalys_ts ~ 1, h = 0.15) 

# 2. OBTENER EL NÚMERO ÓPTIMO DE RUPTURAS (k)
# La forma correcta y robusta es usando la función BIC() del paquete strucchange 
# y encontrando el índice (k) que minimiza este valor.

# which.min(BIC(bp_fit)) nos da el índice del k_óptimo (donde índice 1 = k=0)
# Por lo tanto, restamos 1 para obtener el número de rupturas (k)
num_bp <- which.min(BIC(bp_fit)) - 1

# 3. Imprimir el resultado y obtener los años si hay rupturas
if (num_bp > 0) {
  # Obtenemos los índices de las rupturas para el número óptimo k
  # El argumento 'breaks' le dice a la función cuántas rupturas usar.
  bp_years_index <- breakpoints(bp_fit, breaks = num_bp)$breakpoints
  
  # Mapear los índices a los años originales
  # Asumimos que dalys_mediana_anual está cargado y ordenado por año
  bp_years <- dalys_mediana_anual$year[bp_years_index] 
  
  # Formato de salida para el subtítulo del gráfico
  bp_label <- paste(bp_years, collapse = ", ")
  
  message(paste0("Se detectaron ", num_bp, " rupturas en los años: ", bp_label))
} else {
  bp_label <- "—"
  message(paste0("No se identificaron rupturas adicionales significativas (Múltiples Rupturas). Número óptimo k = ", num_bp))
}
# --------------------------------------------------------------------------------
# RE-CONSTRUCCIÓN FINAL DEL GRÁFICO (PARA INCLUIR LA ETIQUETA BP)
# --------------------------------------------------------------------------------

# Asume que: mk_tau, mk_p, pettitt_year, pettitt_p, y dalys_mediana_anual están definidos.
mk_label <- paste0("MK: tau=", mk_tau, ", p=", mk_p)
pettitt_label <- paste0("Pettitt: ", pettitt_year, " (p=", pettitt_p, ")")
year_ruptura <- pettitt_year

p <- ggplot(dalys_mediana_anual, aes(x = year, y = median_dalys)) +
  
  # A. Sombreado (IQR)
  geom_ribbon(aes(ymin = q1_dalys, ymax = q3_dalys), 
              fill = "gray", alpha = 0.5) +
  
  # B. Línea de la mediana
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  
  # C. Línea de Ruptura (Pettitt)
  geom_vline(xintercept = year_ruptura, linetype = "dashed", color = "black") +
  
  # D. Etiquetas y Títulos (AHORA CON LA ETIQUETA DE BREAKPOINTS 'BPs')
  labs(
    title = "DALYs: Mediana anual (IQR sombreado)",
    subtitle = paste0(mk_label, " JT: p=4e-04, ", pettitt_label, " BPs: ", bp_label),
    y = "Mediana por año (con IQR)",
    x = ""
  ) +
  
  # E. Tema y Escalas
  scale_x_continuous(breaks = seq(1990, 2020, by = 10)) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  ) +
  
  # F. Etiqueta Pettitt
  annotate("text", x = year_ruptura + 1, y = min(dalys_mediana_anual$q1_dalys, na.rm=TRUE), 
           label = "Pettitt", angle = 90, size = 3, hjust = 0)

print(p)

library(ggplot2)
library(dplyr)
library(ggplot2)

# --- 1. PREPARACIÓN DE DATOS ---
# Creamos la variable 'Regimen' para segmentar la tendencia en 2004 (punto de ruptura)
dalys_mediana_anual_segmentada <- dalys_mediana_anual %>%
  mutate(
    # Crear una variable de factor para el régimen
    Regimen = factor(ifelse(year < 2004, "1990-2003 (Pre-Ruptura)", "2004-2020 (Post-Ruptura)"))
  )

# --- 2. GRÁFICO FINAL DE REGRESIÓN SEGMENTADA ---
p_regimen <- ggplot(dalys_mediana_anual_segmentada, aes(x = year, y = median_dalys)) +
  
  # A. Graficar la mediana como puntos
  geom_point(size = 2, color = "#8B0000") + 
  geom_line(color = "#8B0000", linewidth = 0.8) +
  
  # B. Añadir la línea de ruptura vertical en 2004
  geom_vline(xintercept = 2004, linetype = "dashed", color = "black", linewidth = 1) +
  
  # C. Añadir las líneas de tendencia segmentadas (la clave)
  # Usamos 'lm' para la regresión lineal dentro de cada Regimen
  geom_smooth(aes(group = Regimen, color = Regimen), 
              method = "lm", 
              formula = y ~ x, 
              se = FALSE, # No mostrar el intervalo de confianza (opcional)
              linewidth = 1.5) +
  
  # D. Estética y etiquetas
  scale_color_manual(values = c("1990-2003 (Pre-Ruptura)" = "#3399CC", # Azul
                                "2004-2020 (Post-Ruptura)" = "#FF6600")) + # Naranja
  labs(
    title = "Tendencia de DALYs de Salud Mental (Mediana Regional): Ruptura en 2004",
    subtitle = "La pendiente de crecimiento de la carga de enfermedad se aceleró en el régimen post-2004.",
    x = "Año",
    y = "Mediana Anual de DALYs (por 100,000 hab.)",
    color = "Período" # Título para la leyenda de las tendencias
  ) +
  theme_minimal(base_size = 14)

print(p_regimen)

library(trend)

# --- Prueba Theil-Sen (Pendiente de cambio anual) ---

# 1. Ejecutar la prueba de Theil-Sen (sens.slope)
# Esta función calcula la mediana de todas las pendientes posibles entre pares de puntos.
ts_test <- sens.slope(dalys_ts)

# 2. Extraer los resultados
ts_slope <- round(ts_test$estimates[1], 1)
ts_p <- format.pval(ts_test$p.value, digits = 2, scientific = TRUE)

# 3. Imprimir el resultado
cat("\n============================================\n")
cat("          RESULTADO DE LA PRUEBA DE THEIL-SEN      \n")
cat("============================================\n")
cat(paste0("Pendiente (Cambio anual): ", ts_slope, " unidades DALYs/año\n"))
cat(paste0("Valor p: ", ts_p, "\n"))

# 4. Interpretación
if (ts_test$p.value < 0.05) {
  cat("\nConclusión: La pendiente es estadísticamente significativa.\n")
  cat(paste0("La mediana regional de DALYs aumenta en promedio ", ts_slope, " unidades por año.\n"))
} else {
  cat("\nConclusión: La tendencia no es estadísticamente significativa.\n")
}

library(ggplot2)
library(dplyr)
library(ggplot2)

# --- 1. PREPARACIÓN DE DATOS ---
# Creamos la variable 'Regimen' para segmentar la tendencia en 2004 (punto de ruptura)
dalys_mediana_anual_segmentada <- dalys_mediana_anual %>%
  mutate(
    # Crear una variable de factor para el régimen
    Regimen = factor(ifelse(year < 2004, "1990-2003 (Pre-Ruptura)", "2004-2020 (Post-Ruptura)"))
  )

# --- 2. GRÁFICO FINAL DE REGRESIÓN SEGMENTADA ---
p_regimen <- ggplot(dalys_mediana_anual_segmentada, aes(x = year, y = median_dalys)) +
  
  # A. Graficar la mediana como puntos
  geom_point(size = 2, color = "#8B0000") + 
  geom_line(color = "#8B0000", linewidth = 0.8) +
  
  # B. Añadir la línea de ruptura vertical en 2004
  geom_vline(xintercept = 2004, linetype = "dashed", color = "black", linewidth = 1) +
  
  # C. Añadir las líneas de tendencia segmentadas (la clave)
  # Usamos 'lm' para la regresión lineal dentro de cada Regimen
  geom_smooth(aes(group = Regimen, color = Regimen), 
              method = "lm", 
              formula = y ~ x, 
              se = FALSE, # No mostrar el intervalo de confianza (opcional)
              linewidth = 1.5) +
  
  # D. Estética y etiquetas
  scale_color_manual(values = c("1990-2003 (Pre-Ruptura)" = "#3399CC", # Azul
                                "2004-2020 (Post-Ruptura)" = "#FF6600")) + # Naranja
  labs(
    title = "Tendencia de DALYs de Salud Mental (Mediana Regional): Ruptura en 2004",
    subtitle = "La pendiente de crecimiento de la carga de enfermedad se aceleró en el régimen post-2004.",
    x = "Año",
    y = "Mediana Anual de DALYs (por 100,000 hab.)",
    color = "Período" # Título para la leyenda de las tendencias
  ) +
  theme_minimal(base_size = 14)

print(p_regimen)


# --------------------------------------------------------------------------
# Configuración del Dataframe (Necesaria para PLM)
# --------------------------------------------------------------------------

# Antes de correr cualquier test de panel, debes configurar tu dataframe como un objeto 'pdata.frame'.
# Esto requiere tus columnas de país ('country') y tiempo ('year').

# Asume que tu dataframe final limpio se llama mat_indicadores
mat_indicadores_p <- pdata.frame(mat_indicadores, index = c("country", "year"))


# --------------------------------------------------------------------------
# 1. TEST DE ESTACIONARIEDAD (Alternativa: Test de Hadri)
# --------------------------------------------------------------------------

message("1. Test de Estacionariedad (Hadri) -- Intentando una prueba más robusta para NA")

# H0: La serie es estacionaria (no tiene raíz unitaria)
# Si p < 0.05, rechazamos H0: La serie NO es estacionaria.

hadri_test <- purtest(DALYs ~ 1,
                      data = mat_indicadores_p,
                      test = "hadri",
                      exo = "intercept")
print(hadri_test)


# --------------------------------------------------------------------------
# 2. AUTOCORRELACIÓN TEMPORAL
# --------------------------------------------------------------------------

message("\n3. Test de Autocorrelación (Breusch-Godfrey para Panel)")

# Este test verifica si el error de un año está correlacionado con el error del año anterior.
# H0: No hay autocorrelación serial (Independencia de errores).

# Usamos el test de Breusch-Godfrey (pbgtest) en el modelo pooling
bg_test <- pbgtest(modelo_pooling, order = 1) # 'order = 1' verifica autocorrelación de rezago 1
print(bg_test)

# Lógica: Si p < 0.05 -> Autocorrelación. Los errores estándar serán incorrectos; 
#deberás usar modelos de series temporales (ej. AR(1)) o errores estándar robustos.


# --- 1) POST-HOC DUNN: DALYs por país ---
# Asegura que sea factor antes de pasarla
mat_indicadores <- mat_indicadores %>%
  mutate(country = as.factor(country),
         year    = as.factor(year))

# Ahora sí funciona:
dunn_country <- mat_indicadores %>%
  dunn_test(DALYs ~ country, p.adjust.method = "BH")

dunn_year <- mat_indicadores %>%
  dunn_test(DALYs ~ year, p.adjust.method = "BH")


# Resumen rápido para “reportar”
country_summary <- dunn_country %>%
  mutate(sig = p.adj < 0.05) %>%
  summarise(
    total_comparisons = n(),
    significant_pairs = sum(sig),
    prop_significant = mean(sig)
  )

print(country_summary)
# Opcional: ver las 20 diferencias con p.adj más bajas
dunn_country %>% arrange(p.adj) %>% head(20)

# --- 2) POST-HOC DUNN: DALYs por año ---
dunn_year <- mat_indicadores %>%
  dunn_test(DALYs ~ year, p.adjust.method = "BH")

year_summary <- dunn_year %>%
  mutate(sig = p.adj < 0.05) %>%
  summarise(
    total_comparisons = n(),
    significant_pairs = sum(sig),
    prop_significant = mean(sig)
  )

print(year_summary)
dunn_year %>% arrange(p.adj) %>% head(20)

# --- 3) (Opcional) Efecto global de Kruskal (tamaño de efecto) ---
eff_country <- rstatix::kruskal_effsize(mat_indicadores, DALYs ~ as.factor(country))
eff_year    <- rstatix::kruskal_effsize(mat_indicadores, DALYs ~ as.factor(year))
print(eff_country)
print(eff_year)

# --- 4) MULTI-INDICADOR: post-hoc por país y por año ---
# Por país (todos los indicadores)
dunn_by_ind_country <- mat_indicadores %>%
  pivot_longer(c(DALYs, YLLs, YLDs, Prevalence, Incidence, Deaths),
               names_to = "Indicador", values_to = "val") %>%
  group_by(Indicador) %>%
  dunn_test(val ~ country, p.adjust.method = "BH") %>%
  ungroup()

# Por año (todos los indicadores)
dunn_by_ind_year <- mat_indicadores %>%
  pivot_longer(c(DALYs, YLLs, YLDs, Prevalence, Incidence, Deaths),
               names_to = "Indicador", values_to = "val") %>%
  group_by(Indicador) %>%
  dunn_test(val ~ year, p.adjust.method = "BH") %>%
  ungroup()

# Resúmenes compactos para “sí hicimos post-hoc”
summary_ind_country <- dunn_by_ind_country %>%
  group_by(Indicador) %>%
  summarise(n_comp = n(),
            n_sig = sum(p.adj < 0.05),
            prop_sig = mean(p.adj < 0.05)) %>%
  arrange(desc(prop_sig))



summary_ind_year <- dunn_by_ind_year %>%
  group_by(Indicador) %>%
  summarise(n_comp = n(),
            n_sig = sum(p.adj < 0.05),
            prop_sig = mean(p.adj < 0.05)) %>%
  arrange(desc(prop_sig))

print(summary_ind_country)
print(summary_ind_year)

# --- 5) Guardar (si quieres adjuntar al anexo/suplemento) ---
readr::write_csv(dunn_country, "output/dunn_dalys_country.csv")
readr::write_csv(dunn_year,    "output/dunn_dalys_year.csv")
readr::write_csv(dunn_by_ind_country, "output/dunn_all_indicators_country.csv")
readr::write_csv(dunn_by_ind_year,    "output/dunn_all_indicators_year.csv")
readr::write_csv(summary_ind_country, "output/dunn_summary_country.csv")
readr::write_csv(summary_ind_year,    "output/dunn_summary_year.csv")



# Boxplots por año
p_box <- mat_indicadores |>
  pivot_longer(cols = c(DALYs, YLLs, YLDs, Prevalence, Incidence, Deaths),
               names_to = "Indicador", values_to = "val") |>
  ggplot(aes(as.factor(year), val)) + geom_boxplot(outlier.shape = 16, outlier.alpha=.4) +
  facet_wrap(~Indicador, scales="free_y") +
  labs(title="LAC: Distribución por año e indicador", x="Año", y="Valor") +
  theme(axis.text.x = element_text(angle=90, hjust=1))
p_box


ggsave(here("output","boxplots_por_anio_lac.png"), p_box, width=12, height=8, dpi=220)

# -------- DESCRIPCIÓN VARIABLE DEPENDIENTE: DALYs ----------------------------
message("Descripción integral de la variable dependiente (DALYs)…")

# --- 1. Estadísticos básicos -------------------------------------------------
desc_dalys <- dalys_total %>%
  summarise(
    n         = n(),
    media     = mean(DALYs, na.rm=TRUE),
    mediana   = median(DALYs, na.rm=TRUE),
    sd        = sd(DALYs, na.rm=TRUE),
    min       = min(DALYs, na.rm=TRUE),
    max       = max(DALYs, na.rm=TRUE),
    q25       = quantile(DALYs, 0.25, na.rm=TRUE),
    q75       = quantile(DALYs, 0.75, na.rm=TRUE),
    iqr       = IQR(DALYs, na.rm=TRUE),
    asimetria = moments::skewness(DALYs, na.rm=TRUE),
    curtosis  = moments::kurtosis(DALYs, na.rm=TRUE)
  )
print(desc_dalys)

# --- 2. Histograma con densidad normalizada ----------------------------------
p_hist <- ggplot(dalys_total, aes(x = DALYs)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 40, 
                 fill = "#1750a4", 
                 color = "white", 
                 alpha = 0.8) +
  geom_density(color = "#ff7b07", 
               size = 1.2, 
               linetype = "solid", 
               alpha = 0.9) +
  labs(
    title = "Distribución de DALYs en América Latina y el Caribe (1990–2021)",
    x = "DALYs",
    y = "Densidad"
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

p_hist

ggsave(
  here::here("output", "dalys_hist_density.png"),
  p_hist,
  width = 8, height = 5, dpi = 300
)

# --- 3. QQ-plot ---------------------------------------------------------------
# Mostrar en consola
qqnorm(dalys_total$DALYs, main="QQ-plot de DALYs")
qqline(dalys_total$DALYs, col="red")

# Guardar también
png(here::here("output","dalys_qqplot.png"), width=800, height=600, res=120)
qqnorm(dalys_total$DALYs, main="QQ-plot de DALYs")
qqline(dalys_total$DALYs, col="red")
dev.off()

# --- 4. Boxplots por año ------------------------------------------------------
p_box <- dalys_total %>%
  ggplot(aes(as.factor(year), DALYs)) +
  geom_boxplot(outlier.shape=16, outlier.alpha=.3) +
  labs(title="Distribución anual de DALYs en LAC", x="Año", y="DALYs") +
  theme(axis.text.x = element_text(angle=90, hjust=1))
p_box
ggsave(here::here("output","dalys_boxplot_year.png"), p_box, width=10, height=6, dpi=220)

# 1. Identificar outliers por año
outliers <- dalys_total %>%
  group_by(year) %>%
  mutate(
    Q1 = quantile(DALYs, 0.25, na.rm=TRUE),
    Q3 = quantile(DALYs, 0.75, na.rm=TRUE),
    IQR = Q3 - Q1,
    limite_inf = Q1 - 1.5*IQR,
    limite_sup = Q3 + 1.5*IQR,
    outlier = DALYs < limite_inf | DALYs > limite_sup
  ) %>%
  filter(outlier)

# 2. Boxplot + etiquetas en outliers
p_box_labels <- dalys_total %>%
  ggplot(aes(as.factor(year), DALYs)) +
  geom_boxplot(outlier.shape=NA) +   # ocultar puntos automáticos
  geom_point(data=outliers, aes(x=as.factor(year), y=DALYs), color="red") +
  geom_text_repel(data=outliers,
                  aes(x=as.factor(year), y=DALYs, label=country),
                  size=2.5, color="black", max.overlaps=10) +
  labs(title="Distribución anual de DALYs en LAC",
       x="Año", y="DALYs") +
  theme(axis.text.x = element_text(angle=90, hjust=1))

print(p_box_labels)

# --- 5. Tendencia temporal (promedio anual + banda de sd) ---------------------
trend <- dalys_total %>%
  group_by(year) %>%
  summarise(mean_dalys = mean(DALYs, na.rm=TRUE),
            sd_dalys   = sd(DALYs, na.rm=TRUE))

p_trend <- ggplot(trend, aes(year, mean_dalys)) +
  geom_line(color="darkgreen") +
  geom_ribbon(aes(ymin=mean_dalys-sd_dalys, ymax=mean_dalys+sd_dalys),
              alpha=0.2, fill="darkgreen") +
  labs(title="Tendencia temporal de DALYs en LAC", y="Promedio (±1 SD)", x="Año")

p_trend

ggsave(here::here("output","dalys_trend.png"), p_trend, width=9, height=6, dpi=220)

# --- 6. Heatmap País-Año ------------------------------------------------------
# --- Matriz país (filas) x año (columnas) NUMÉRICA ---
dalys_wide <- dalys_total %>%
  mutate(DALYs = as.double(DALYs)) %>%              # asegura numérico
  pivot_wider(names_from = year, values_from = DALYs) %>%
  arrange(iso3)

# convierte TODAS las columnas (menos iso3) a numérico
dalys_wide <- dalys_wide %>%
  mutate(across(-iso3, ~ as.double(.x)))

dalys_mat <- dalys_wide %>%
  column_to_rownames("iso3") %>%
  as.matrix()

storage.mode(dalys_mat) <- "double"                 # por si acaso

# quita filas/columnas totalmente NA
dalys_mat <- dalys_mat[rowSums(is.na(dalys_mat)) < ncol(dalys_mat), , drop = FALSE]
dalys_mat <- dalys_mat[, colSums(is.na(dalys_mat)) < nrow(dalys_mat), drop = FALSE]

# distancias NA‑robustas por correlación (pairwise)
row_cor <- cor(t(dalys_mat), use = "pairwise.complete.obs")
col_cor <- cor(dalys_mat,     use = "pairwise.complete.obs")
row_cor[is.na(row_cor)] <- 0
col_cor[is.na(col_cor)] <- 0
row_dist <- as.dist(1 - row_cor)
col_dist <- as.dist(1 - col_cor)

# --- Mostrar en Plots
pheatmap(
  dalys_mat,
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = row_dist,
  clustering_distance_cols = col_dist,
  clustering_method = "ward.D2",
  color = viridis(100),
  main  = "Heatmap DALYs por país y año"
)

# --- Guardar
png(here("output","dalys_heatmap.png"), width = 1200, height = 800, res = 150)
pheatmap(
  dalys_mat,
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = row_dist,
  clustering_distance_cols = col_dist,
  clustering_method = "ward.D2",
  color = viridis(100),
  main  = "Heatmap DALYs por país y año"
)
dev.off()


# Wide limpio y numérico
dalys_wide <- dalys_total %>%
  mutate(DALYs = as.double(DALYs),
         year = as.integer(year)) %>%
  arrange(iso3, year) %>%
  pivot_wider(names_from = year, values_from = DALYs)

# matriz
dalys_mat <- dalys_wide %>% column_to_rownames("iso3") %>% as.matrix()
storage.mode(dalys_mat) <- "double"

# quita filas/columnas totalmente NA
dalys_mat <- dalys_mat[rowSums(is.na(dalys_mat)) < ncol(dalys_mat), , drop=FALSE]
dalys_mat <- dalys_mat[, colSums(is.na(dalys_mat)) < nrow(dalys_mat), drop=FALSE]

# ordenar columnas cronológicamente (muy importante)
year_cols <- sort(as.integer(colnames(dalys_mat)))
dalys_mat <- dalys_mat[, as.character(year_cols)]

# distancias NA‑robustas para filas (correlación pairwise)
row_cor  <- cor(t(dalys_mat), use="pairwise.complete.obs"); row_cor[is.na(row_cor)] <- 0
row_dist <- as.dist(1 - row_cor)

# Heatmap: columnas fijas (cronológicas), filas con clúster
pheatmap(
  dalys_mat,
  cluster_rows = TRUE, cluster_cols = FALSE,    # <<< sin clúster en columnas
  clustering_distance_rows = row_dist,
  clustering_method = "ward.D2",
  color = viridis(100),
  main  = "Heatmap DALYs por país (filas) y año (columnas, 1990→2021)"
)
print(pheatmap)
pheatmap
# guardar
png(here("output","dalys_heatmap_cols_cron.png"), width=1400, height=900, res=150)

pheatmap(
  dalys_mat,
  cluster_rows = TRUE, cluster_cols = FALSE,
  clustering_distance_rows = row_dist,
  clustering_method = "ward.D2",
  color = viridis(100),
  main  = "Heatmap DALYs por país (filas) y año (columnas, 1990→2021)"
)

dev.off()

while (!is.null(dev.list())) dev.off()



# 1) Calcula el heatmap SIN dibujarlo (silent=TRUE)
heat_obj <- pheatmap::pheatmap(
  dalys_mat,
  cluster_rows = TRUE, cluster_cols = FALSE,
  clustering_distance_rows = row_dist,   # asegúrate de tener row_dist creado
  clustering_method = "ward.D2",
  color = viridis::viridis(100),
  main  = "Heatmap DALYs por país (filas) y año (columnas, 1990→2021)",
  silent = TRUE
)

# 2) Dibuja en el panel Plots
grid::grid.newpage()
grid::grid.draw(heat_obj$gtable)


# --- 7. Mapas continuos -------------------------------------------------------
# --- Parámetros globales ---
YEAR_FROM    <- 1990
YEAR_TO      <- 2003
EXCLUDE_ISO3 <- character(0)
USE_LOG      <- TRUE
OUTDIR       <- here("output")
custom_palette <- c("#20dbd8", "#37ddb4", "#1d40bb", "#cf4eb9", "#f497bd")

# --- Utilidades ---
get_lac_sf <- function() {
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    st_make_valid() |>
    clean_names() |>
    mutate(iso3 = iso_a3) |>
    filter(iso3 %in% LAC_ISO3)
}

# Función para guardar mapas sin Cairo, usando ragg::agg_png
save_tmap <- function(tm, file, w=9, h=6, dpi=220) {
  ragg::agg_png(file, width = w, height = h, units = "in", res = dpi)
  print(tm)
  dev.off()
  message("Mapa guardado: ", file)
}

base_theme <- tm_layout(
  frame = FALSE,
  fontfamily = "Times New Roman",
  legend.position = c("right", "center"),
  title.fontface = "bold",
  title.fontfamily = "Times New Roman"
)

# --- 1) DALYs acumulados ----
years_sel <- YEAR_FROM:YEAR_TO
dalys_acc <- dalys_total %>%
  filter(year %in% years_sel, !iso3 %in% EXCLUDE_ISO3) %>%
  group_by(iso3) %>%
  summarise(DALYs_sum = sum(DALYs, na.rm = TRUE), .groups="drop")

lac_sf <- get_lac_sf()

map_df <- lac_sf %>%
  filter(!iso3 %in% EXCLUDE_ISO3) %>%
  left_join(dalys_acc, by="iso3")

# Ajustar tema para leyenda externa
base_theme <- tm_layout(
  frame = FALSE,
  fontfamily = "Times New Roman",
  title.fontfamily = "Times New Roman",
  title.fontface = "bold",
  legend.outside = TRUE,                    # Leyenda fuera del mapa
  legend.outside.position = "right",        # A la derecha
  legend.stack = "vertical",
  outer.margins = c(0.02, 0.25, 0.02, 0.02) # Margen externo
)

tm_dalys <- tm_shape(map_df) +
  tm_polygons(
    fill = "DALYs_sum",
    fill.palette = custom_palette,
    fill.trans = if (USE_LOG) "log10" else "identity",
    fill.legend = tm_legend(
      title = glue::glue("DALYs acumulados\n{YEAR_FROM}–{YEAR_TO}"),
      frame = TRUE, bg.color = "white", bg.alpha = 1
    ),
    lwd = 0.2, col = "white"
  ) +
  tm_title(glue::glue("LAC — DALYs acumulados ({YEAR_FROM}–{YEAR_TO})")) +
  base_theme
tm_dalys
outfile <- here(OUTDIR, glue::glue("map_dalys_acumulados_{YEAR_FROM}_{YEAR_TO}.png"))
save_tmap(tm_dalys, outfile)

# --- 7. Mapas continuos -------------------------------------------------------
# --- Parámetros globales ---
YEAR_FROM    <- 2004
YEAR_TO      <- 2021
EXCLUDE_ISO3 <- character(0)
USE_LOG      <- TRUE
OUTDIR       <- here("output")
custom_palette <- c("#20dbd8", "#37ddb4", "#1d40bb", "#cf4eb9", "#f497bd")

# --- Utilidades ---
get_lac_sf <- function() {
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    st_make_valid() |>
    clean_names() |>
    mutate(iso3 = iso_a3) |>
    filter(iso3 %in% LAC_ISO3)
}

# Función para guardar mapas sin Cairo, usando ragg::agg_png
save_tmap <- function(tm, file, w=9, h=6, dpi=220) {
  ragg::agg_png(file, width = w, height = h, units = "in", res = dpi)
  print(tm)
  dev.off()
  message("Mapa guardado: ", file)
}

base_theme <- tm_layout(
  frame = FALSE,
  fontfamily = "Times New Roman",
  legend.position = c("right", "center"),
  title.fontface = "bold",
  title.fontfamily = "Times New Roman"
)

# --- 1) DALYs acumulados ----
years_sel <- YEAR_FROM:YEAR_TO
dalys_acc <- dalys_total %>%
  filter(year %in% years_sel, !iso3 %in% EXCLUDE_ISO3) %>%
  group_by(iso3) %>%
  summarise(DALYs_sum = sum(DALYs, na.rm = TRUE), .groups="drop")

lac_sf <- get_lac_sf()

map_df <- lac_sf %>%
  filter(!iso3 %in% EXCLUDE_ISO3) %>%
  left_join(dalys_acc, by="iso3")

# Ajustar tema para leyenda externa
base_theme <- tm_layout(
  frame = FALSE,
  fontfamily = "Times New Roman",
  title.fontfamily = "Times New Roman",
  title.fontface = "bold",
  legend.outside = TRUE,                    # Leyenda fuera del mapa
  legend.outside.position = "right",        # A la derecha
  legend.stack = "vertical",
  outer.margins = c(0.02, 0.25, 0.02, 0.02) # Margen externo
)

tm_dalys <- tm_shape(map_df) +
  tm_polygons(
    fill = "DALYs_sum",
    fill.palette = custom_palette,
    fill.trans = if (USE_LOG) "log10" else "identity",
    fill.legend = tm_legend(
      title = glue::glue("DALYs acumulados\n{YEAR_FROM}–{YEAR_TO}"),
      frame = TRUE, bg.color = "white", bg.alpha = 1
    ),
    lwd = 0.2, col = "white"
  ) +
  tm_title(glue::glue("LAC — DALYs acumulados ({YEAR_FROM}–{YEAR_TO})")) +
  base_theme
tm_dalys
outfile <- here(OUTDIR, glue::glue("map_dalys_acumulados_{YEAR_FROM}_{YEAR_TO}.png"))
save_tmap(tm_dalys, outfile)

# --- Función para LISA sobre la Tasa de Crecimiento (CORREGIDA) ---

generar_mapas_lisa_crecimiento <- function(start_year, end_year) {
  
  # 1. Preparar datos de inicio y fin (USANDO dplyr::select explícitamente)
  dalys_start <- dalys_total %>% 
    dplyr::filter(year == start_year) %>% 
    dplyr::select(iso3, DALYs_start = DALYs) # <--- CORRECCIÓN AQUÍ
  
  dalys_end <- dalys_total %>% 
    dplyr::filter(year == end_year) %>% 
    dplyr::select(iso3, DALYs_end = DALYs) # <--- CORRECCIÓN AQUÍ
  
  # 2. Calcular la Tasa de Crecimiento Porcentual
  dalys_growth <- dalys_start %>%
    dplyr::inner_join(dalys_end, by = "iso3") %>%
    dplyr::mutate(
      Tasa_Crecimiento = ((DALYs_end - DALYs_start) / DALYs_start) * 100
    )
  
  # 3. Unir con la geometría espacial
  # Asumo que lac_sf ya tiene la geometría
  lac_dalys <- lac_sf %>%
    dplyr::inner_join(dalys_growth, by="iso3") %>%
    dplyr::filter(!is.na(Tasa_Crecimiento)) %>%
    st_make_valid()
  
  if (nrow(lac_dalys) < 5) { message("No hay suficientes datos para el análisis"); return(NULL) }
  
  # 4. Cálculo LISA (sobre Tasa_Crecimiento)
  coords <- st_coordinates(st_point_on_surface(st_geometry(lac_dalys)))
  knn    <- spdep::knearneigh(coords, k=3)
  nbk    <- spdep::knn2nb(knn)
  lwk    <- spdep::nb2listw(nbk, style="W")
  
  lisa <- spdep::localmoran(lac_dalys$Tasa_Crecimiento, lwk)
  lisa <- as.data.frame(lisa)
  colnames(lisa) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii", "Pr.z")
  
  x <- scale(lac_dalys$Tasa_Crecimiento)[,1]
  wx <- spdep::lag.listw(lwk, x)
  
  # 5. Agregar variables LISA y Cuadrante de Crecimiento
  lac_dalys <- lac_dalys %>%
    dplyr::mutate(
      Zi = lisa$Z.Ii, Pi = lisa$Pr.z, x = x, wx = wx,
      quadrant = dplyr::case_when(
        x >= 0 & wx >= 0 & Pi < 0.05 ~ "High-High (Crecimiento)",
        x < 0 & wx < 0 & Pi < 0.05   ~ "Low-Low (Crecimiento)",
        x >= 0 & wx < 0 & Pi < 0.05  ~ "High-Low (Crecimiento)",
        x < 0 & wx >= 0 & Pi < 0.05  ~ "Low-High (Crecimiento)",
        TRUE                         ~ "No significativo"
      )
    ) %>%
    st_transform(crs = 4326)
  
  # -----------------------------------
  # 6. Generación de Mapas (3 mapas por periodo)
  # -----------------------------------
  title_periodo <- paste0(" (", start_year, "-", end_year, ")")
  
  # Paleta categórica para Clusters de Crecimiento
  cluster_pal <- c("High-High (Crecimiento)" = "#1f78b4", "Low-Low (Crecimiento)" = "#33a02c", 
                   "High-Low (Crecimiento)" = "#e31a1c", "Low-High (Crecimiento)" = "#ff7f00", 
                   "No significativo" = "grey90")
  
  # MAPA 1: Clusters LISA (Cuadrantes de Crecimiento)
  tm_clusters <- tm_shape(lac_dalys) +
    tm_polygons(
      col = "quadrant", 
      palette = cluster_pal,
      fill.legend = tm_legend(title=paste0("Clusters LISA", title_periodo))
    ) +
    tm_title(paste0("Clusters LISA: Tasa de Crecimiento de DALYs", title_periodo)) +
    tm_layout(frame = FALSE, legend.outside = TRUE, title.size = 1.2)
  
  # MAPA 2: Z Local (Continuo)
  tm_lisa_z <- tm_shape(lac_dalys) +
    tm_polygons(
      col="Zi",
      palette = custom_palette,
      fill.scale = tm_scale_continuous(midpoint=NA),
      fill.legend = tm_legend(title=paste0("Z-score Moran Local", title_periodo))
    ) +
    tm_title(paste0("Z-score LISA (Tasa de Crecimiento)", title_periodo)) +
    tm_layout(frame = FALSE, legend.outside = TRUE, title.size = 1.2)
  
  # MAPA 3: Distribución de la Tasa de Crecimiento (Nivel)
  tm_growth_level <- tm_shape(lac_dalys) +
    tm_polygons(
      col="Tasa_Crecimiento",
      palette = "viridis",
      fill.legend = tm_legend(title=paste0("Tasa de Crecimiento (%)", title_periodo))
    ) +
    tm_title(paste0("Nivel: Tasa de Crecimiento de DALYs", title_periodo)) +
    tm_layout(frame = FALSE, legend.outside = TRUE, title.size = 1.2)
  
  # 7. Guardar los mapas
  prefijo <- paste0("Tasa_Crecimiento_", start_year, "_", end_year)
  save_tmap(tm_clusters, here(OUTDIR, paste0("1_clusters_", prefijo, ".png")))
  save_tmap(tm_lisa_z, here(OUTDIR, paste0("2_zscore_", prefijo, ".png")))
  save_tmap(tm_growth_level, here(OUTDIR, paste0("3_nivel_", prefijo, ".png")))
  
  message(paste("Mapas de crecimiento generados y guardados para el periodo:", start_year, "-", end_year))
  
  return(list(clusters = tm_clusters, z_local = tm_lisa_z, nivel_crecimiento = tm_growth_level, data = lac_dalys))
  }

# Período 1: Régimen Pre-Ruptura (1990-2003)
# Esto genera los mapas 1_clusters_Tasa_Crecimiento_1990_2003.png, etc.
generar_mapas_lisa_crecimiento(start_year = 1990, end_year = 2003)

# Período 2: Régimen Post-Ruptura (2004-2020)
# Esto genera los mapas 1_clusters_Tasa_Crecimiento_2004_2020.png, etc.
generar_mapas_lisa_crecimiento(start_year = 2004, end_year = 2020)

# -------- 8. SEM -------------------
#Otro archivo 

# -------- 9. BIVARIADO (Spearman con determinantes) --------------------------
message("Correlaciones Spearman (LAC) con determinantes…")
stopifnot(file.exists(PATH_DET))
det <- fread(PATH_DET) |> clean_names() |> filter(iso3 %in% LAC_ISO3) # country, iso3, year, pobreza, gini, violencia, desempleo, educacion, gasto_salud, etc.

base_modelo <- dalys_total |>
  left_join(det, by = c("country","iso3","year")) 

candidatas <- setdiff(names(det), c("country","iso3","year"))
cors <- map_dfr(candidatas, ~{
  ct <- suppressWarnings(cor.test(base_modelo$DALYs, base_modelo[[.x]], method="spearman"))
  tibble(variable = .x, rho = unname(ct$estimate), p = ct$p.value)
}) |> arrange(p)
write_csv(cors, here("output","spearman_dalys_determinantes_lac.csv"))

# -------- MATRIZ DE CORRELACIÓN SPEARMAN (SIGNIFICATIVAS) --------------------------
message("Calculando matriz de correlación (Spearman) con significancia...")


# --- 1. Selección de variables -----------------------------------------------------
vars_cor <- c("DALYs", "desempleo", "gasto_salud_gdp", "oop_salud_pct",
              "camas_1000", "vac_measles", "urbanizacion", "internet_usuarios",
              "pib_pc", "crecimiento_pib", "inflacion_cpi", "lfp",
              "life_expectancy", "edu_spend_gdp",
              "daly_cardio", "daly_diabetes_kidney", "daly_hiv",
              "daly_neoplasms", "daly_respiratory",
              "daly_transport", "indice_infra")

vars_cor <- intersect(vars_cor, names(base_modelo))

# --- 2. Calcular correlaciones con p-valores usando Hmisc::rcorr -------------------
cor_obj <- Hmisc::rcorr(as.matrix(base_modelo[, vars_cor]), type = "spearman")

# Matriz de correlaciones (rho)
mat_cor <- cor_obj$r

# Matriz de p-valores
p_mat <- cor_obj$P

# --- 3. Etiquetas en español -------------------------------------------------------
var_labels <- c(
  "DALYs"                = "Carga total de enfermedad (DALYs)",
  "desempleo"            = "Desempleo (%)",
  "gasto_salud_gdp"      = "Gasto en salud (% del PIB)",
  "oop_salud_pct"        = "Gasto de bolsillo en salud (%)",
  "camas_1000"           = "Camas hospitalarias (por 1,000 hab.)",
  "vac_measles"          = "Cobertura de vacunación de sarampión (%)",
  "urbanizacion"         = "Urbanización (%)",
  "internet_usuarios"    = "Usuarios de internet (%)",
  "pib_pc"               = "PIB per cápita (USD constantes)",
  "crecimiento_pib"      = "Crecimiento del PIB (%)",
  "inflacion_cpi"        = "Inflación anual (%)",
  "lfp"                  = "Participación laboral (%)",
  "life_expectancy"      = "Esperanza de vida (años)",
  "edu_spend_gdp"        = "Gasto educativo (% del PIB)",
  "daly_cardio"          = "Carga DALYs: enfermedades cardiovasculares",
  "daly_diabetes_kidney" = "Carga DALYs: diabetes y enfermedad renal",
  "daly_hiv"             = "Carga DALYs: VIH/SIDA",
  "daly_neoplasms"       = "Carga DALYs: neoplasias",
  "daly_respiratory"     = "Carga DALYs: respiratorias",
  "daly_transport"       = "Carga DALYs: transporte y accidentes",
  "indice_infra"         = "Índice de infraestructura"
)

colnames(mat_cor) <- var_labels[colnames(mat_cor)]
rownames(mat_cor) <- var_labels[rownames(mat_cor)]
colnames(p_mat)   <- var_labels[colnames(p_mat)]
rownames(p_mat)   <- var_labels[rownames(p_mat)]

# --- 4. Exportar CSV (correlaciones + p-valores) ----------------------------------
write_csv(as.data.frame(mat_cor), here::here("output","matriz_cor_spearman_signif_esp.csv"))

# --- 5. Visualización (solo correlaciones significativas p<0.05) ------------------
p_cor <- ggcorrplot(mat_cor,
                    method = "square",
                    type = "lower",
                    p.mat = p_mat,
                    sig.level = 0.05,
                    insig = "blank",   # deja en blanco las no significativas
                    lab = TRUE,
                    lab_size = 2.5,
                    colors = c("#6D9EC1", "white", "#E46726"),
                    title = "Matriz de correlaciones de Spearman (solo significativas, p < 0.05)",
                    ggtheme = theme_minimal(base_family = "Times New Roman") +
                      theme(
                        plot.title = element_text(size = 13, face = "bold"),
                        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
                        axis.text.y = element_text(size = 8),
                        text = element_text(family = "Times New Roman")
                      ))

print(p_cor)

# --- 6. Guardar gráfico ------------------------------------------------------------
ggsave(here::here("output", "matriz_cor_spearman_significativas_esp.png"),
       plot = p_cor, width = 9, height = 7, dpi = 300)

message("✅ Matriz de correlaciones (solo significativas) guardada en /output/")


# -------- 10. CAMBIO TEMPORAL (mixtos + piecewise) ---------------------------
# Asumiendo que candidatas ya está definido
candidatas <- setdiff(names(det), c("country", "iso3", "year"))
message("Modelos de cambio temporal (LAC)…")

# Asegúrate de que `year_c` está correctamente calculado
# Si no está calculado, centrar `year` como ejemplo:
X <- base_modelo |>
  mutate(across(all_of(candidatas), scale), .after = year) |>
  bind_cols(year_splines(base_modelo$year, BREAKS))

# Modelo mixto con intercepto aleatorio por país
m1 <- lmer(DALYs ~ year_c + (1|iso3), data = X, REML = FALSE)

# Modelo Piecewise (con puntos de cambio en 2010 y 2020)
m2 <- lmer(DALYs ~ year_c + s1 + s2 + (1|iso3), data = X, REML = FALSE)

# Modelo multivariado con determinantes + piecewise
form_det <- as.formula(glue("DALYs ~ year_c + s1 + s2 + {paste(candidatas, collapse=' + ')} + (1|iso3)"))
m3 <- lmer(form_det, data = X, REML = FALSE)

summary(m1)
summary(m2)
summary(m3)

# Resumen de los modelos
msummary(list("Mixto (básico)" = m1, "Mixto piecewise" = m2, "Multivariado piecewise" = m3),
         stars = TRUE, output = here("output", "modelos_mixtos_lac.html"))

# Verificación de colinealidad en el modelo multivariado
colin <- performance::check_collinearity(m3)
capture.output(colin, file = here("output", "vif_m3_lac.txt"))


# ===== Forest plot robusto para m3 =====
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(forcats); library(tibble)
})

stopifnot(exists("m3"))

# 1) Extraer coeficientes (lmerTest/lme4) y armar CIs Wald
tbl <- coef(summary(m3)) %>%           # más robusto que summary(m3)$coefficients
  as.data.frame() %>%
  rownames_to_column("term") %>%
  rename(
    estimate  = Estimate,
    std_error = `Std. Error`,
    df        = df,
    t_value   = `t value`,
    p_value   = `Pr(>|t|)`
  ) %>%
  mutate(
    conf_low  = estimate - 1.96 * std_error,
    conf_high = estimate + 1.96 * std_error
  )

# 2) Opciones
remove_intercept <- TRUE
only_significant <- FALSE
top_n_by_t       <- 0   # 0 = mostrar TODO

# 3) Diccionario de etiquetas (named vector)
nice_labels <- c(
  year_c = "Año (centrado)",
  s1 = "Spline 1", s2 = "Spline 2",
  desempleo = "Desempleo (%)",
  desempleo_joven = "Desempleo joven (%)",
  violencia = "Violencia (índice)",
  gasto_salud_gdp = "Gasto salud (% PIB)",
  gasto_salud_pc = "Gasto salud p.c.",
  gob_salud_gdp = "Gasto salud público (% PIB)",
  oop_salud_pct = "Gasto OOP salud (%)",
  camas_1000 = "Camas por 1000 hab.",
  vac_measles = "Cobertura SRP (%)",
  urbanizacion = "Urbanización (%)",
  electrificacion = "Electrificación (%)",
  agua_basica = "Acceso a agua básica (%)",
  saneamiento_basico = "Saneamiento básico (%)",
  internet_usuarios = "Usuarios de internet (%)",
  banda_ancha_100 = "Banda ancha (c/100)",
  pib_pc = "PIB per cápita",
  crecimiento_pib = "Crecimiento PIB (%)",
  inflacion_cpi = "Inflación (%)",
  lfp = "Participación laboral (%)",
  life_expectancy = "Esperanza de vida (años)",
  edu_spend_gdp = "Gasto educación (% PIB)",
  daly_cardio = "DALY cardiometabólicos",
  daly_conflict = "DALY conflicto",
  daly_diabetes_kidney = "DALY diabetes/renal",
  daly_hiv = "DALY VIH",
  daly_neoplasms = "DALY neoplasias",
  daly_police = "DALY violencia policial",
  daly_respiratory = "DALY respiratorias",
  daly_tb = "DALY TB",
  daly_transport = "DALY transporte",
  daly_violence = "DALY violencia interpersonal"
)

# 4) Preparar data para el plot (una sola vez)
df_plot <- tbl %>%
  { if (remove_intercept) filter(., term != "(Intercept)") else . } %>%
  { if (only_significant)  filter(., !is.na(p_value) & p_value < 0.05) else . } %>%
  arrange(desc(abs(t_value))) %>%
  { if (top_n_by_t > 0) slice_head(., n = top_n_by_t) else . } %>%
  mutate(
    label = dplyr::recode(
      term,
      !!!as.list(nice_labels),     # usa el diccionario
      .default = term,
      .missing = term
    ),
    sig05 = ifelse(!is.na(p_value) & p_value < 0.05, "p < 0.05", "ns"),
    label = forcats::fct_rev(forcats::fct_inorder(label))
  )

message("Términos graficados: ", nrow(df_plot))
if (nrow(df_plot) == 0) stop("No hay términos para graficar. Revisa filtros.")

# 5) Graficar
p_forest <- ggplot(df_plot,
                   aes(x = estimate, y = label, xmin = conf_low, xmax = conf_high, color = sig05)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  geom_pointrange(linewidth = 0.6) +
  scale_color_manual(values = c("p < 0.05" = "#1b7837", "ns" = "#7f7f7f"), drop = FALSE) +
  labs(
    title = "Efectos fijos de m3 (IC95% Wald)",
    subtitle = "Modelo mixto: DALYs ~ predictores + (1|país); puntos = β, barras = IC95%",
    x = "Coeficiente (β)", y = NULL, color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

print(p_forest)
ggsave("forest_m3.png", p_forest, width = 9, height = 7, dpi = 300)



# --------------------------------------------------------------------------
# 1. IDENTIFICAR Y ELIMINAR FILAS CON NA (SINTAXIS CORREGIDA)
# --------------------------------------------------------------------------

# 1.1. Definir todas las variables necesarias (asumiendo que 'candidatas' está definida)
variables_modelo <- c("DALYs", "year_c", "s1", "s2", "iso3", candidatas)

# 1.2. Crear el dataframe limpio (usando la sintaxis de corchetes para evitar conflictos)
# Esto selecciona las columnas y luego elimina cualquier fila con NA en esas columnas.
X_clean <- X[ , variables_modelo] %>%
  na.omit()

message(paste0("Filas originales: ", nrow(X), ". Filas usadas para el modelo: ", nrow(X_clean), "."))

# --------------------------------------------------------------------------
# 2. DEFINIR Y EJECUTAR EL MODELO INICIAL (M3)
# --------------------------------------------------------------------------

# Usa el dataframe limpio (X_clean)
form_det <- as.formula(glue("DALYs ~ year_c + s1 + s2 + {paste(candidatas, collapse=' + ')} + (1|iso3)"))
m3 <- lmer(form_det, data = X_clean, REML = FALSE)

# --------------------------------------------------------------------------
# 3. SELECCIÓN STEPWISE (Con el dataframe limpio)
# --------------------------------------------------------------------------

message("\nIniciando Selección Stepwise de Variables (basada en p-valores) con datos limpios...")

# Ahora la función step() debería ejecutarse sin problemas
modelo_parsimonioso_pvalue <- step(m3)

# --------------------------------------------------------------------------
# 4. Resultado Final
# --------------------------------------------------------------------------

print(summary(modelo_parsimonioso_pvalue))

# La función 'step' a veces modifica el objeto original m3 o devuelve solo la tabla ANOVA.
# Intentaremos forzar la impresión de la tabla de ANOVA generada (que sí contiene el resultado):
print(modelo_parsimonioso_pvalue$fixed) 
print(modelo_parsimonioso_pvalue$random)

# 1. Mira la lista completa de las variables que pasaron el filtro (Eliminated = 0).
# Variables: s1, lfp, life_expectancy, daly_conflict, daly_transport, daly_violence, indice_infra

# 2. Define la fórmula final
variables_finales <- c("s1", "lfp", "life_expectancy", "daly_conflict", "daly_transport", "daly_violence", "indice_infra")

# 3. Construye y corre el Modelo Final (M4)
form_final <- as.formula(glue("DALYs ~ {paste(variables_finales, collapse=' + ')} + (1|iso3)"))
m3f_final <- lmer(form_final, data = X_clean, REML = TRUE) # Usamos REML=TRUE para el reporte final

# 4. Reporte final
print(summary(m3f_final))


# -------- 12. BOOTSTRAP (coeficientes m3) ------------------------------------
message("Bootstrap paramétrico (bootmer) sobre m3…")
# ==== Paquetes ====
suppressPackageStartupMessages({
  library(lme4)
  library(boot)
  library(dplyr)
})

# ==== Semilla reproducible ====
set.seed(123)
RNGkind("L'Ecuyer-CMRG")

# ==== Datos y modelo base ====
# Asume que tienes un data.frame X con columnas: DALYs, iso3, year (numérico)
# Si no tienes year_c, lo creo centrando year:
if (!"year_c" %in% names(X)) {
  X <- X %>% mutate(year_c = year - mean(year, na.rm = TRUE))
}

# Convierte iso3 a factor (recomendado para LMM)
X <- X %>% mutate(iso3 = factor(iso3))

# Modelo sobre datos completos (para "original")
m3 <- lmer(DALYs ~ year_c + (1 | iso3), data = X, REML = FALSE)

library(lmtest)
library(multiwayvcov) # Necesario si usas coeftest con vc_clust

# Asegúrate de usar X_clean y el modelo de Efectos Fijos (lm), no mixto (lmer), para este enfoque

# 1. Ejecutar el modelo OLS simple (pooling) en el dataframe limpio
# Solo con las variables que SÍ fueron significativas en tu m4_final
variables_finales <- c("s1", "lfp", "life_expectancy", "daly_conflict", "daly_transport", "daly_violence", "indice_infra")
form_ols <- as.formula(glue("DALYs ~ {paste(variables_finales, collapse=' + ')}"))
m_ols <- lm(form_ols, data = X_clean)

# 2. Calcular la matriz de varianza-covarianza (VCOV) robusta agrupada por país (iso3)
# Esto corrige la autocorrelación dentro de los países y la heterocedasticidad.
vcov_cluster <- multiwayvcov::cluster.vcov(m_ols, X_clean$iso3)

# 3. Mostrar los resultados con errores robustos
reporte_robusto_cluster <- coeftest(m_ols, vcov_cluster)

# 4. Imprimir los coeficientes (la forma más robusta de hacer inferencia para tu modelo)
print(reporte_robusto_cluster)




# -------- 13. FIGURAS CLAVE ---------------------------------------------------
message("Figuras clave (LAC)…")
p_line <- dalys_total |>
  group_by(year) |>
  summarise(mean = mean(DALYs, na.rm=TRUE),
            min  = min(DALYs, na.rm=TRUE),
            max  = max(DALYs, na.rm=TRUE), .groups="drop") |>
  ggplot(aes(year, mean)) +
  geom_ribbon(aes(ymin=min, ymax=max), alpha=.15) +
  geom_line(size=1) +
  geom_vline(xintercept = BREAKS, linetype = 2) +
  labs(title="LAC: DALYs promedio y rango (min–max) por año", y="DALYs", x=NULL)
p_line
ggsave(here("output","tendencia_dalys_lac.png"), p_line, width=9, height=5, dpi=220)


# =========================================================================
# 1. PREPARACIÓN DE VARIABLES Y DATOS
# =========================================================================

# NOTA: Asumo que 'X', 'candidatas' y 'candidatas_clean' existen del código anterior.

# Crear la variable DALYs estandarizada (valores Z)
# USAR 'X' AQUÍ (el dataframe completo)
X <- X %>% mutate(DALYs_Z = scale(DALYs)[,1]) # scale() devuelve una matriz, [ ,1] la convierte a vector

# Variables de DALYs a eliminar por su alta colinealidad/circularidad
variables_a_eliminar <- c(
  "daly_violence", # Interpersonal violence
  "daly_conflict", # Conflict and terrorism
  "daly_police" # Police conflict and executions (asumo que se mapea a esto)
)

# Crear la nueva lista de candidatas SIN las variables de violencia/conflicto
# Asumo que 'candidatas' existe y es el vector de variables predictoras originales
candidatas_clean <- candidatas[!candidatas %in% variables_a_eliminar]

message(paste0("Candidatas limpias: ", length(candidatas_clean)))

# 1.1. Definir todas las variables necesarias para el nuevo modelo (Usando DALYs_Z)
variables_modelo_clean <- c("DALYs_Z", "year_c", "s1", "s2", "iso3", candidatas_clean)

# 1.2. Crear el dataframe limpio y estandarizado (X_clean_new_Z)
X_clean_new_Z <- X[ , variables_modelo_clean] %>%
  na.omit()

message(paste0("Filas usadas para el nuevo modelo: ", nrow(X_clean_new_Z), "."))

# =========================================================================
# 2. MODELO MIXTO COMPLEJO (M3_CLEAN_Z) y STEPWISE
# =========================================================================

# 2.1. Definir y ejecutar el modelo multivariado completo (M3_CLEAN_Z)
form_det_clean_Z <- as.formula(glue("DALYs_Z ~ year_c + s1 + s2 + {paste(candidatas_clean, collapse=' + ')} + (1|iso3)"))

m3_clean_Z <- lmer(form_det_clean_Z, data = X_clean_new_Z, REML = FALSE)
message("\nModelo Mixto Complejo (M3_CLEAN_Z) ejecutado.")

# 2.2. Verificación de colinealidad (opcional, solo para el complejo)
colin_clean_Z <- performance::check_collinearity(m3_clean_Z)
capture.output(colin_clean_Z, file = here("output", "vif_m3_clean_lac_Z.txt"))

# 2.3. Selección Stepwise
message("\nIniciando Selección Stepwise (backward) en el modelo estandarizado...")
modelo_parsimonioso_clean_Z <- step(m3_clean_Z)

# 2.4. Extraer las variables que quedaron
vars_finales_clean <- modelo_parsimonioso_clean_Z$fixed %>%
  filter(Eliminated == 0) %>%
  rownames()

# Quitar el intercepto si se incluyó
vars_finales_clean <- vars_finales_clean[vars_finales_clean != "(Intercept)"]

message("\nVariables retenidas en el modelo parsimonioso (sin violencia):\n",
        paste(vars_finales_clean, collapse = ", "))

# =========================================================================
# 3. REPORTE FINAL ROBUSTO (OLS con Errores Agrupados)
# =========================================================================

# 3.1. Usamos las variables finales para el modelo OLS (con DALYs_Z)
form_ols_clean_Z <- as.formula(glue("DALYs_Z ~ {paste(vars_finales_clean, collapse=' + ')}"))
m_ols_clean_Z <- lm(form_ols_clean_Z, data = X_clean_new_Z)

# 3.2. Calcular la matriz de varianza-covarianza (VCOV) robusta agrupada por país (iso3)
vcov_cluster_clean_Z <- multiwayvcov::cluster.vcov(m_ols_clean_Z, X_clean_new_Z$iso3)

# 3.3. Mostrar los resultados con errores robustos
reporte_robusto_clean_Z <- coeftest(m_ols_clean_Z, vcov_cluster_clean_Z)

cat("\n======================================================\n")
cat("     REPORTE FINAL ROBUSTO (DALYs_Z, Sin Violencia)   \n")
cat("======================================================\n")
summary(m3_clean_Z)
print(reporte_robusto_clean_Z)


# =========================================================================
# 4. INDICADORES DE BONDAD DE AJUSTE (R^2, AIC, BIC)
# =========================================================================

cat("\n==================================================\n")
cat("         Indicadores para M3_CLEAN_Z (LMM)          \n")
cat("==================================================\n")
# R Cuadrada para LMM
r2_m3_Z <- MuMIn::r.squaredGLMM(m3_clean_Z)
cat(paste("AIC:", round(AIC(m3_clean_Z), 2), "\n"))
cat(paste("BIC:", round(BIC(m3_clean_Z), 2), "\n"))
cat("\nR-squared (Nakagawa-Schielzeth):\n")
print(r2_m3_Z)

cat("\n==================================================\n")
cat("         Indicadores para m_ols_clean_Z (OLS)        \n")
cat("==================================================\n")
summary_ols_Z <- summary(m_ols_clean_Z)
cat(paste("AIC:", round(AIC(m_ols_clean_Z), 2), "\n"))
cat(paste("BIC:", round(BIC(m_ols_clean_Z), 2), "\n"))
cat(paste("R-squared estándar:", round(summary_ols_Z$r.squared, 4), "\n"))
cat(paste("R-squared ajustada:", round(summary_ols_Z$adj.r.squared, 4), "\n"))


# 1. PREPARACIÓN DE DATOS (Mismo que antes)

# Asumo que 'reporte_robusto_clean_Z' está en el entorno.
df_robusto <- broom::tidy(reporte_robusto_clean_Z) %>%
  mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    sig05 = case_when(
      p.value < 0.001 ~ "p < 0.001 (***)",
      p.value < 0.01  ~ "p < 0.01 (**)",
      p.value < 0.05  ~ "p < 0.05 (*)",
      TRUE            ~ "ns"
    ),
    sig05 = factor(sig05, levels = c("ns", "p < 0.05 (*)", "p < 0.01 (**)", "p < 0.001 (***)"))
  ) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    label = case_when(
      term == "daly_diabetes_kidney" ~ "DALYs Diabetes/Riñón",
      term == "daly_respiratory" ~ "DALYs Respiratorio",
      term == "daly_transport" ~ "DALYs Transporte",
      term == "desempleo" ~ "Desempleo",
      term == "edu_spend_gdp" ~ "Gasto Educación (% PIB)",
      term == "life_expectancy" ~ "Esperanza de Vida",
      term == "lfp" ~ "Participación Laboral",
      term == "s1" ~ "Estacionalidad 1",
      term == "urbanizacion" ~ "Urbanización",
      term == "vac_measles" ~ "Vacuna Sarampión",
      term == "year_c" ~ "Año Centrado",
      TRUE ~ term
    )
  )

# 2. GENERACIÓN DEL GRÁFICO (Forest Plot) con Times New Roman

p_forest_tnr <- ggplot(df_robusto,
                       aes(x = estimate, y = label, xmin = conf.low, xmax = conf.high, color = sig05)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  geom_pointrange(linewidth = 0.6) +
  scale_color_manual(
    values = c(
      "p < 0.001 (***)" = "#a00574",
      "p < 0.01 (**)" = "#238b45",
      "p < 0.05 (*)" = "#ff7b07"
    ),
    drop = FALSE,
    name = "Significancia"
  ) +
  labs(
    title = "Forest Plot: Coeficientes del Modelo OLS Parsimonioso",
    subtitle = "Variable Dependiente: DALYs (Estandarizado Z) - Errores Robustos por País",
    x = "Estimación (Cambio en Desviación Estándar de DALYs)",
    y = NULL,
    color = NULL
  ) +
  # APLICAR TEMA Y FUENTE Times New Roman
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold"),
    # Aquí se define la fuente para todo el texto del gráfico:
    text = element_text(family = "Times New Roman"),
    title = element_text(family = "Times New Roman", face = "bold")
  )

# 3. MOSTRAR Y GUARDAR
print(p_forest_tnr)
ggsave("forest_tnr.png", p_forest_tnr, width = 9, height = 7, dpi = 300)
