# ============================================================
#  SALUD MENTAL EN AMÉRICA LATINA Y EL CARIBE — PIPELINE R
### TODO PERO SIN VIOLENCIA EN LA VARIABLE DEPENDIENTE####
### Después la quito también de la independiente###
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
BREAKS   <- c(2005)


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
unique(gbd$cause)

# -------- 4. MÉTRICAS PRINCIPALES --------------------------------------------
unique(gbd$metric)
gbd <- gbd %>% dplyr::mutate(metric = stringr::str_replace(metric, "^(DALYs|YLDs|YLLs).*", "\\1"))
unique(gbd$cause)
message("Construyendo tablas por indicador LAC…")

# 1. Definir las causas a excluir
causas_violencia <- c(
  "Interpersonal violence",
  "Police conflict and executions",
  "Conflict and terrorism"
)

gbd <- gbd %>%
  filter(!iso3 %in% c("PRI", "VIR"))


# 2. Modificar dalys_total para filtrar primero y luego sumar
dalys_total <- gbd %>%
  # Filtra solo las métricas DALYs
  filter(metric == "DALYs") %>%
  
  # 🚨 NUEVO PASO: Excluir las causas de violencia/conflicto
  filter(!cause %in% causas_violencia) %>%
  
  # Agrupar y sumar el valor (DALYs sin violencia)
  group_by(country, iso3, year) %>%
  summarise(DALYs = sum(value, na.rm=TRUE), .groups="drop")

message("DALYs totales calculados (excluyendo violencia/conflicto).")

dalys_total %>%
  summarise(
    n_obs = n(),
    n_paises = n_distinct(iso3),
    n_years = n_distinct(year)
  )

table(table(dalys_total$iso3))

dalys_total %>%
  distinct(iso3, country) %>%
  arrange(iso3) %>%
  print(n = Inf)



colnames(dalys_total)


mat_indicadores <- gbd |>
  filter(metric %in% c("DALYs","YLLs","YLDs","Prevalence","Incidence","Deaths")) |>
  group_by(country, iso3, year, metric) |>
  summarise(val = sum(value, na.rm=TRUE), .groups="drop") |>
  pivot_wider(names_from = metric, values_from = val)

mat_indicadores <- mat_indicadores %>%
  filter(year >= 1990, year <= 2021)


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
    median_dalys = median(DALYs, na.rm = TRUE),
    q1_dalys     = quantile(DALYs, 0.25, na.rm = TRUE),
    q3_dalys     = quantile(DALYs, 0.75, na.rm = TRUE)
  ) %>%
  arrange(year) %>%
  ungroup()


# La serie de tiempo simple (vector)
dalys_ts <- dalys_mediana_anual$median_dalys

# --- Prueba 1: Mann-Kendall (Tendencia monótona) ---
mk_test <- mk.test(dalys_ts)
mk_tau <- round(mk_test$estimates[3], 3)
mk_p <- mk_test$p.value  # deja el valor numérico puro


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

# ------------------------------------------------------------
# 1. Preparación de datos
# ------------------------------------------------------------

dalys_mediana_anual <- dalys_total %>%
  group_by(year) %>%
  summarise(
    median_dalys = median(DALYs, na.rm = TRUE)
  ) %>%
  arrange(year) %>%
  ungroup() %>%
  mutate(
    time = year - min(year) + 1,          # tiempo lineal
    post2005 = ifelse(year >= 2005, 1, 0), # dummy régimen
    time_post = time * post2005            # interacción (cambio de pendiente)
  )

# ------------------------------------------------------------
# 2. Modelo segmentado (cambio en nivel y pendiente)
# ------------------------------------------------------------

modelo_segmentado <- lm(
  median_dalys ~ time + post2005 + time_post,
  data = dalys_mediana_anual
)

summary(modelo_segmentado)

# ------------------------------------------------------------
# 3. Extraer pendientes
# ------------------------------------------------------------

b <- coef(modelo_segmentado)

pendiente_pre  <- b["time"]
pendiente_post <- b["time"] + b["time_post"]
cambio_pend    <- b["time_post"]

cat("\nPendiente pre-2005:", pendiente_pre)
cat("\nPendiente post-2005:", pendiente_post)
cat("\nCambio en pendiente:", cambio_pend)

library(nlme)

modelo_gls <- gls(
  median_dalys ~ time + post2005 + time_post,
  data = dalys_mediana_anual,
  correlation = corAR1()
)

summary(modelo_gls)



# --- Prueba 3: Theil-Sen (Pendiente de cambio anual) ---
ts_test <- sens.slope(dalys_ts)
ts_slope <- round(ts_test$estimates[1], 1)
ts_p <- format.pval(ts_test$p.value, digits = 2)

library(trend)

pre  <- dalys_mediana_anual %>% filter(year < 2005)
post <- dalys_mediana_anual %>% filter(year >= 2005)

ts_pre  <- sens.slope(pre$median_dalys)
ts_post <- sens.slope(post$median_dalys)

ts_pre
ts_post

library(trend)

# Theil–Sen
ts_test <- sens.slope(dalys_ts)

# Extraer estimaciones
ts_slope  <- as.numeric(ts_test$estimates)
ts_low    <- ts_test$conf.int[1]
ts_high   <- ts_test$conf.int[2]
ts_p      <- ts_test$p.value

# Redondear para reporte
ts_slope_r <- round(ts_slope, 1)
ts_low_r   <- round(ts_low, 1)
ts_high_r  <- round(ts_high, 1)
ts_p_r     <- format.pval(ts_p, digits = 3, scientific = TRUE)

cat("\n====================================\n")
cat("     RESULTADOS THEIL–SEN (GLOBAL)\n")
cat("====================================\n")
cat(paste0("Pendiente: ", ts_slope_r, " DALYs/año\n"))
cat(paste0("IC95%: ", ts_low_r, " – ", ts_high_r, "\n"))
cat(paste0("p-valor: ", ts_p_r, "\n"))



modelo_simple <- lm(median_dalys ~ time, data = dalys_mediana_anual)
anova(modelo_simple, modelo_segmentado)

###Tendencia cuadrática
modelo_quad <- lm(median_dalys ~ time + I(time^2), 
                  data = dalys_mediana_anual)

summary(modelo_quad)

#Prueba Chow
library(strucchange)

sctest(median_dalys ~ time, 
       type = "Chow", 
       point = which(dalys_mediana_anual$year == 2005),
       data = dalys_mediana_anual)

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
dalys_mediana_anual <- dalys_total %>%
  group_by(year) %>%
  summarise(
    median_dalys = median(DALYs, na.rm = TRUE),
    q1_dalys = quantile(DALYs, 0.25, na.rm = TRUE),
    q3_dalys = quantile(DALYs, 0.75, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(year)
# Extraer resultados del test de Pettitt
pettitt_year <- pettitt_test$estimate
pettitt_p <- pettitt_test$p.value

mk_label <- paste0("Mann–Kendall: tau = ", mk_tau, ", p = ", signif(mk_p, 3))
pettitt_label <- paste0("Pettitt: ", pettitt_year, " (p = ", signif(pettitt_p, 3), ")")
pettitt_year <- dalys_mediana_anual$year[pettitt_test$estimate]


year_ruptura <- pettitt_year

p <- ggplot(dalys_mediana_anual, aes(x = year, y = median_dalys)) +
  geom_ribbon(aes(ymin = q1_dalys, ymax = q3_dalys),
              fill = "gray", alpha = 0.5) +
  geom_line(linewidth = 0.8, color = "#2c3e50") +
  geom_point(size = 1.5, color = "black") +
  geom_vline(xintercept = pettitt_year,
             linetype = "dashed",
             color = "black") +
  labs(
    title = "DALYs: Mediana anual (IQR sombreado)",
    subtitle = paste0(
      "Mann–Kendall: τ = ", round(mk_tau, 3),
      ", p = ", signif(mk_p, 3)
    ),
    x = "Año",
    y = "Mediana por año (con IQR)"
  ) +
  scale_x_continuous(breaks = seq(1990, 2020, by = 5)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

print(p)


library(ggplot2)
library(dplyr)
library(ggplot2)

# --- 1. PREPARACIÓN DE DATOS ---
# Creamos la variable 'Regimen' para segmentar la tendencia en 2005 (punto de ruptura)
dalys_mediana_anual_segmentada <- dalys_mediana_anual %>%
  mutate(
    # Crear una variable de factor para el régimen
    Regimen = factor(ifelse(year < 2005, "1990-2005 (Pre-Ruptura)", "2005-2020 (Post-Ruptura)"))
  )

# --- 2. GRÁFICO FINAL DE REGRESIÓN SEGMENTADA ---
p_regimen <- ggplot(dalys_mediana_anual_segmentada, aes(x = year, y = median_dalys)) +
  
  # A. Graficar la mediana como puntos
  geom_point(size = 2, color = "#8B0000") + 
  geom_line(color = "#8B0000", linewidth = 0.8) +
  
  # B. Añadir la línea de ruptura vertical en 2005
  geom_vline(xintercept = 2005, linetype = "dashed", color = "black", linewidth = 1) +
  
  # C. Añadir las líneas de tendencia segmentadas (la clave)
  # Usamos 'lm' para la regresión lineal dentro de cada Regimen
  geom_smooth(aes(group = Regimen, color = Regimen), 
              method = "lm", 
              formula = y ~ x, 
              se = FALSE, # No mostrar el intervalo de confianza (opcional)
              linewidth = 1.5) +
  
  # D. Estética y etiquetas
  scale_color_manual(values = c("1990-2005 (Pre-Ruptura)" = "#3399CC", # Azul
                                "2005-2020 (Post-Ruptura)" = "#FF6600")) + # Naranja
  labs(
    title = "Tendencia de DALYs de Salud Mental (Mediana Regional): Ruptura en 2005",
    subtitle = "La pendiente de crecimiento de la carga de enfermedad se aceleró en el régimen post-2005.",
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
# Creamos la variable 'Regimen' para segmentar la tendencia en 2005 (punto de ruptura)
dalys_mediana_anual_segmentada <- dalys_mediana_anual %>%
  mutate(
    # Crear una variable de factor para el régimen
    Regimen = factor(ifelse(year < 2005, "1990-2004 (Pre-Ruptura)", "2005-2020 (Post-Ruptura)"))
  )

# --- 2. GRÁFICO FINAL DE REGRESIÓN SEGMENTADA ---
p_regimen <- ggplot(dalys_mediana_anual_segmentada, aes(x = year, y = median_dalys)) +
  
  # A. Graficar la mediana como puntos
  geom_point(size = 2, color = "#8B0000") + 
  geom_line(color = "#8B0000", linewidth = 0.8) +
  
  # B. Añadir la línea de ruptura vertical en 2005
  geom_vline(xintercept = 2005, linetype = "dashed", color = "black", linewidth = 1) +
  
  # C. Añadir las líneas de tendencia segmentadas (la clave)
  # Usamos 'lm' para la regresión lineal dentro de cada Regimen
  geom_smooth(aes(group = Regimen, color = Regimen), 
              method = "lm", 
              formula = y ~ x, 
              se = FALSE, # No mostrar el intervalo de confianza (opcional)
              linewidth = 1.5) +
  
  # D. Estética y etiquetas
  scale_color_manual(values = c("1990-2004 (Pre-Ruptura)" = "#3399CC", # Azul
                                "2005-2020 (Post-Ruptura)" = "#FF6600")) + # Naranja
  labs(
    title = "Tendencia de DALYs de Salud Mental (Mediana Regional): Ruptura en 2005",
    subtitle = "La pendiente de crecimiento de la carga de enfermedad se aceleró en el régimen post-2005.",
    x = "Año",
    y = "Mediana Anual de DALYs (por 100,000 hab.)",
    color = "Período" # Título para la leyenda de las tendencias
  ) +
  theme_minimal(base_size = 14)

print(p_regimen)

# --------------------------------------------------------------------------
# Configuración del Dataframe (Necesaria para PLM)
# --------------------------------------------------------------------------

library(plm)

mat_indicadores <- gbd |>
  filter(
    metric %in% c("DALYs","YLLs","YLDs","Prevalence","Incidence","Deaths"),
    year >= 1990,
    year <= 2021
  ) |>
  group_by(country, iso3, year, metric) |>
  summarise(val = sum(value, na.rm=TRUE), .groups="drop") |>
  pivot_wider(names_from = metric, values_from = val)



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
# Usa el panel que YA EXISTE
modelo_pooling <- plm(
  DALYs ~ 1,
  data  = mat_indicadores_p,
  model = "pooling"
)

bg_test <- pbgtest(modelo_pooling, order = 1)
print(bg_test)

library(lmtest)

coeftest(
  modelo_pooling,
  vcov = vcovHC(modelo_pooling,
                method = "arellano",
                type   = "HC1",
                cluster = "group")
)



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
YEAR_TO      <- 2004
EXCLUDE_ISO3 <- character(0)
USE_LOG      <- TRUE
OUTDIR       <- here("output")
custom_palette <- c(
  "#ffd03c",  # amarillo dorado
  "#f6e45e",  # amarillo suave
  "#0064a6",  # azul intenso
  "#00a9b2",  # turquesa medio
  "#00f0ca"   # verde-menta brillante
)

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
YEAR_FROM    <- 2005
YEAR_TO      <- 2021
EXCLUDE_ISO3 <- character(0)
USE_LOG      <- TRUE
OUTDIR       <- here("output")
custom_palette <- c(
  "#ffd03c",  # amarillo dorado
  "#f6e45e",  # amarillo suave
  "#0064a6",  # azul intenso
  "#00a9b2",  # turquesa medio
  "#00f0ca"   # verde-menta brillante
)
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

# Período 1: Régimen Pre-Ruptura (1990-2004)
# Esto genera los mapas 1_clusters_Tasa_Crecimiento_1990_2004.png, etc.
generar_mapas_lisa_crecimiento(start_year = 1990, end_year = 2005)

# Período 2: Régimen Post-Ruptura (2005-2020)
# Esto genera los mapas 1_clusters_Tasa_Crecimiento_2005_2020.png, etc.
generar_mapas_lisa_crecimiento(start_year = 2005, end_year = 2020)

### ¿Quién es el peor? 

ranking_paises <- dalys_total %>%
  group_by(country) %>%
  summarise(
    mean_dalys = mean(DALYs, na.rm = TRUE),
    median_dalys = median(DALYs, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_dalys))

head(ranking_paises, 5)

mediana_regional_global <- median(dalys_total$DALYs, na.rm = TRUE)

ranking_paises <- ranking_paises %>%
  mutate(
    brecha_pct = (mean_dalys - mediana_regional_global) /
      mediana_regional_global * 100
  )

quintiles_anuales <- dalys_total %>%
  group_by(year) %>%
  mutate(
    quintil = ntile(DALYs, 5)
  ) %>%
  ungroup()

persistencia_quintil <- quintiles_anuales %>%
  filter(quintil == 5) %>%
  group_by(country) %>%
  summarise(
    anos_en_quintil_superior = n(),
    prop_periodo = n() / 32
  ) %>%
  arrange(desc(prop_periodo))

top10 <- ranking_paises %>% slice_max(mean_dalys, n = 10)

ggplot(top10, aes(x = reorder(country, mean_dalys), y = mean_dalys)) +
  geom_col(fill = "#34495e") +
  coord_flip() +
  labs(
    title = "Países con mayor carga promedio (1990–2021)",
    x = NULL,
    y = "DALYs promedio"
  ) +
  theme_minimal()

ggplot(ranking_paises,
       aes(x = reorder(country, brecha_pct), y = brecha_pct)) +
  geom_col(aes(fill = brecha_pct > 0)) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#8B0000",
                               "FALSE" = "#003366"),
                    guide = "none") +
  labs(
    title = "Brecha porcentual respecto a la mediana regional",
    y = "% sobre/ bajo la mediana regional",
    x = NULL
  ) +
  theme_minimal()

ggplot(persistencia_quintil,
       aes(x = reorder(country, prop_periodo),
           y = prop_periodo)) +
  geom_col(fill = "#7f8c8d") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Persistencia en el quintil superior de carga",
    y = "Proporción del periodo 1990–2021",
    x = NULL
  ) +
  theme_minimal()

heatmap_quintiles <- quintiles_anuales %>%
  mutate(country = reorder(country, DALYs, median))

ggplot(heatmap_quintiles,
       aes(x = year, y = country, fill = as.factor(quintil))) +
  geom_tile() +
  scale_fill_manual(values = c("#003366", "#2c3e50",
                               "#7f8c8d", "#b0b0b0",
                               "#8B0000"),
                    name = "Quintil") +
  labs(
    title = "Estratificación regional de la carga de DALYs",
    x = "Año",
    y = NULL
  ) +
  theme_minimal()

# -------- 8. SEM -------------------
#Otro archivo 

# -------- 9. BIVARIADO (Spearman con determinantes) --------------------------
message("Correlaciones Spearman (LAC) con determinantes…")
stopifnot(file.exists(PATH_DET))
det <- fread(PATH_DET) |> clean_names() |> filter(iso3 %in% LAC_ISO3) # country, iso3, year, pobreza, gini, violencia, desempleo, educacion, gasto_salud, etc.

base_modelo <- dalys_total |>
  left_join(det, by = c("country","iso3","year")) 

candidatas <- setdiff(names(det), c("country","iso3","year","life_expectancy"))
cors <- map_dfr(candidatas, ~{
  ct <- suppressWarnings(cor.test(base_modelo$DALYs, base_modelo[[.x]], method="spearman"))
  tibble(variable = .x, rho = unname(ct$estimate), p = ct$p.value)
}) |> arrange(p)
write_csv(cors, here("output","spearman_dalys_determinantes_lac.csv"))

# -------- MATRIZ DE CORRELACIÓN SPEARMAN (SIGNIFICATIVAS) --------------------------
message("Calculando matriz de correlación (Spearman) con significancia...")
base_modelo

# --- 1. Selección de variables -----------------------------------------------------
vars_cor <- c("DALYs", "desempleo", "gasto_salud_gdp", "oop_salud_pct",
              "camas_1000", "vac_measles", "urbanizacion", "internet_usuarios",
              "pib_pc", "crecimiento_pib", "inflacion_cpi", "lfp",
              "edu_spend_gdp",
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

##### Bivariado sólo respecto a DALYS y considerando panel within-country

library(dplyr)

vars_biv <- vars_cor  # las que ya definiste antes

base_within <- base_modelo %>%
  group_by(country) %>%
  mutate(across(all_of(vars_biv),
                ~ . - mean(., na.rm = TRUE),
                .names = "{.col}_w")) %>%
  ungroup()

vars_within <- paste0(vars_biv, "_w")

cor_within <- Hmisc::rcorr(
  as.matrix(base_within[, vars_within]),
  type = "spearman"
)

mat_cor_w <- cor_within$r
p_mat_w   <- cor_within$P

rho_dalys_w <- mat_cor_w["DALYs_w", ]
p_dalys_w   <- p_mat_w["DALYs_w", ]

tabla_within <- tibble(
  variable = names(rho_dalys_w),
  rho = as.numeric(rho_dalys_w),
  p = as.numeric(p_dalys_w)
) %>%
  filter(variable != "DALYs_w") %>%
  mutate(
    variable = gsub("_w", "", variable),
    variable_label = var_labels[variable],
    abs_rho = abs(rho),
    signif = p < 0.05
  ) %>%
  arrange(desc(abs_rho))

tabla_within

n_obs <- nrow(base_within)  # número de observaciones efectivas

tabla_forest <- tabla_within %>%
  mutate(
    z = 0.5 * log((1 + rho) / (1 - rho)),
    se = 1 / sqrt(n_obs - 3),
    z_low = z - 1.96 * se,
    z_high = z + 1.96 * se,
    conf.low = (exp(2*z_low) - 1) / (exp(2*z_low) + 1),
    conf.high = (exp(2*z_high) - 1) / (exp(2*z_high) + 1),
    label = reorder(variable_label, rho)
  )
library("ggplot2")

p_forest_within <- ggplot(tabla_forest,
                          aes(x = rho, y = label,
                              xmin = conf.low, xmax = conf.high)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(color = "#a00574") +
  theme_minimal(base_size = 12) +
  theme(text = element_text(family = "Times New Roman"),
        axis.text.y = element_text(face = "bold")) +
  labs(
    title = "Asociaciones intra-país (Spearman within)",
    subtitle = "Correlaciones centradas por país",
    x = expression(rho),
    y = NULL
  )

print(p_forest_within)

# -------- SUPUESTOS ANTES DEL MULTIVARIADO ------------------------------------
library(dplyr)
library(plm)
library(lmtest)
library(sandwich)
library(car)

# ------------------------------------------------------------
# 1. Preparar datos panel
# ------------------------------------------------------------
panel_df <- pdata.frame(base_modelo, index = c("country", "year"))

# Seleccionar predictores numéricos
vars_pred <- base_modelo %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-DALYs, -year) %>%
  names()

formula_fe <- as.formula(
  paste("DALYs ~", paste(vars_pred, collapse = " + "))
)

# ------------------------------------------------------------
# 2. MODELO EFECTOS FIJOS
# ------------------------------------------------------------
modelo_fe <- plm(
  formula_fe,
  data  = panel_df,
  model = "within"
)

# ------------------------------------------------------------
# 3. MULTICOLINEALIDAD (VIF aproximado vía LM con dummies país)
# ------------------------------------------------------------
# Modelo sin dummies para evaluar colinealidad estructural
modelo_vif <- lm(
  as.formula(
    paste("DALYs ~", paste(vars_pred, collapse = " + "))
  ),
  data = base_modelo
)

vif_tabla <- data.frame(
  variable = names(vif(modelo_vif)),
  VIF = vif(modelo_vif)
) %>%
  arrange(desc(VIF))

cat("\n================ VIF ==================\n")
print(vif_tabla)

# ------------------------------------------------------------
# 4. VARIACIÓN WITHIN (media de varianza dentro de país)
# ------------------------------------------------------------
within_var <- base_modelo %>%
  group_by(country) %>%
  summarise(
    across(all_of(vars_pred),
           ~ var(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  summarise(
    across(everything(),
           \(x) mean(x, na.rm = TRUE))
  ) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "mean_within_variance"
  ) %>%
  arrange(mean_within_variance)

cat("\n========== VARIACIÓN WITHIN ==========\n")
print(within_var)

# ------------------------------------------------------------
# 5. AUTOCORRELACIÓN (Wooldridge / BG para panel)
# ------------------------------------------------------------
bg_test <- pbgtest(modelo_fe, order = 1)

cat("\n========== AUTOCORRELACIÓN ==========\n")
print(bg_test)

# ------------------------------------------------------------
# 6. HETEROCEDASTICIDAD (Breusch-Pagan)
# ------------------------------------------------------------
bp_test <- bptest(modelo_fe)

cat("\n========== HETEROCEDASTICIDAD ==========\n")
print(bp_test)

# ------------------------------------------------------------
# 7. DEPENDENCIA TRANSVERSAL (Pesaran CD)
# ------------------------------------------------------------
cd_test <- pcdtest(modelo_fe, test = "cd")

cat("\n========== DEPENDENCIA TRANSVERSAL ==========\n")
print(cd_test)

# ------------------------------------------------------------
# 8. COEFICIENTES ROBUSTOS (Cluster por país)
# ------------------------------------------------------------
cat("\n========== COEFICIENTES ROBUSTOS ==========\n")
coeftest(
  modelo_fe,
  vcov = vcovHC(modelo_fe,
                method = "arellano",
                type   = "HC1",
                cluster = "group")
)

coeftest(
  modelo_fe,
  vcov = vcovSCC(modelo_fe, type = "HC1", maxlag = 2)
)

# ==========================================================
# PAQUETES
# ==========================================================
library(plm)
library(lmtest)
library(sandwich)
library(dplyr)
library(performance)

# ==========================================================
# 1. DEFINICIÓN DE VARIABLES
# ==========================================================

# --- Determinantes estructurales ---
vars_det <- c(
  "desempleo", "pib_pc", "crecimiento_pib",
  "inflacion_cpi", "lfp", "edu_spend_gdp",
  "gasto_salud_gdp", "oop_salud_pct",
  "camas_1000", "vac_measles",
  "urbanizacion", "indice_infra"
)

vars_det <- intersect(vars_det, names(base_modelo))

# --- DALYs otras patologías ---
vars_epi <- c(
  "daly_cardio", "daly_diabetes_kidney",
  "daly_neoplasms", "daly_hiv",
  "daly_respiratory", "daly_transport",
  "daly_conflict", "daly_violence",
  "daly_police"
)

vars_epi <- intersect(vars_epi, names(base_modelo))

# --- Parsimonioso (coherencia teórica)
vars_pars <- c(
  "desempleo",
  "lfp",
  "indice_infra",
  "oop_salud_pct"
)

vars_pars <- intersect(vars_pars, names(base_modelo))


# ==========================================================
# 2. BASE COMÚN (MISMA PARA TODOS LOS MODELOS)
# ==========================================================

vars_all <- unique(c("DALYs", vars_det, vars_epi))

base_comun <- base_modelo %>%
  dplyr::select(country, year, all_of(vars_all)) %>%
  na.omit()

panel_df <- pdata.frame(base_comun,
                        index = c("country", "year"))


# ==========================================================
# 3. MODELO 1: Determinantes estructurales
# ==========================================================

form_m1 <- as.formula(
  paste("DALYs ~", paste(vars_det, collapse = " + "))
)

m1 <- plm(form_m1,
          data = panel_df,
          model = "within")

coefs_m1 <- coeftest(
  m1,
  vcov = vcovSCC(m1, type = "HC1", maxlag = 2)
)


# ==========================================================
# 4. MODELO 2: Ajuste epidemiológico
# ==========================================================

form_m2 <- as.formula(
  paste("DALYs ~",
        paste(c(vars_det, vars_epi), collapse = " + "))
)

m2 <- plm(form_m2,
          data = panel_df,
          model = "within")

coefs_m2 <- coeftest(
  m2,
  vcov = vcovSCC(m2, type = "HC1", maxlag = 2)
)


# ==========================================================
# 5. MODELO 3: Parsimonioso
# ==========================================================

form_m3 <- as.formula(
  paste("DALYs ~", paste(vars_pars, collapse = " + "))
)

m3 <- plm(form_m3,
          data = panel_df,
          model = "within")

coefs_m3 <- coeftest(
  m3,
  vcov = vcovSCC(m3, type = "HC1", maxlag = 2)
)


# ==========================================================
# 6. CRITERIOS DE INFORMACIÓN
# ==========================================================

extraer_metricas <- function(modelo){
  ll  <- as.numeric(logLik(modelo))
  k   <- length(coef(modelo))
  n   <- nobs(modelo)
  aic <- -2*ll + 2*k
  bic <- -2*ll + log(n)*k
  r2  <- summary(modelo)$r.squared["rsq"]
  
  data.frame(
    LogLik = ll,
    Parametros = k,
    N = n,
    AIC = aic,
    BIC = bic,
    R2_within = r2
  )
}

metricas <- dplyr::bind_rows(
  Modelo_1 = extraer_metricas(m1),
  Modelo_2 = extraer_metricas(m2),
  Modelo_3 = extraer_metricas(m3),
  .id = "Modelo"
)


# ==========================================================
# 7. OUTPUT
# ==========================================================

cat("\n================ MODELO 1 =================\n")
print(coefs_m1)

cat("\n================ MODELO 2 =================\n")
print(coefs_m2)

cat("\n================ MODELO 3 =================\n")
print(coefs_m3)

cat("\n=========== CRITERIOS DE INFORMACIÓN ===========\n")
print(metricas)


library(broom)
library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
library(kableExtra)

# ---------------------------------------------------------
# 1. Función para limpiar coeficientes
# ---------------------------------------------------------

limpiar_coefs <- function(coefs, modelo_nombre){
  
  df <- broom::tidy(coefs) %>%
    mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ ""
      ),
      coef_fmt = sprintf("%.2f%s\n(%.2f)", estimate, stars, std.error)
    ) %>%
    dplyr::select(term, coef_fmt) %>%
    rename(!!modelo_nombre := coef_fmt)
  
  return(df)
}

m1_tab <- limpiar_coefs(coefs_m1, "Modelo 1")
m2_tab <- limpiar_coefs(coefs_m2, "Modelo 2")
m3_tab <- limpiar_coefs(coefs_m3, "Modelo 3")

# ---------------------------------------------------------
# 2. Unir modelos
# ---------------------------------------------------------

tabla_modelos <- full_join(m1_tab, m2_tab, by="term") %>%
  full_join(m3_tab, by="term") %>%
  arrange(term)

# ---------------------------------------------------------
# 3. Añadir métricas al final
# ---------------------------------------------------------

metricas_tab <- metricas %>%
  dplyr::select(Modelo, N, R2_within, AIC, BIC) %>%
  pivot_longer(-Modelo) %>%
  pivot_wider(names_from = Modelo, values_from = value) %>%
  rename(term = name)

tabla_final <- bind_rows(tabla_modelos, metricas_tab)

# ---------------------------------------------------------
# 4. Mostrar en formato publicable
# ---------------------------------------------------------

kbl(tabla_final,
    align = "lccc",
    booktabs = TRUE,
    caption = "Modelos de efectos fijos (errores Driscoll–Kraay). Coeficientes β con errores estándar entre paréntesis."
) %>%
  kable_styling(full_width = FALSE)


# ==========================================================
# 11.8 ROBUSTEZ Y SENSIBILIDAD
# ==========================================================

library(plm)
library(lmtest)
library(sandwich)
library(dplyr)
library(boot)

# ----------------------------------------------------------
# BASE PANEL
# ----------------------------------------------------------

vars_pars <- c(
  "DALYs",
  "desempleo",
  "lfp",
  "indice_infra",
  "oop_salud_pct"
)

base_panel <- base_modelo %>%
  dplyr::select(country, year, all_of(vars_pars)) %>%
  na.omit()

panel_df <- pdata.frame(base_panel, index = c("country","year"))

form_m3 <- DALYs ~ desempleo + lfp + indice_infra + oop_salud_pct

# ==========================================================
# 1. MODELO BASE (Driscoll–Kraay)
# ==========================================================

m_base <- plm(
  form_m3,
  data = panel_df,
  model = "within"
)

coefs_base <- coeftest(
  m_base,
  vcov = vcovSCC(m_base, type = "HC1", maxlag = 2)
)

print(coefs_base)

# ==========================================================
# 2. BOOTSTRAP (intervalos empíricos)
# ==========================================================

boot_fun <- function(data, indices){
  
  d <- data[indices,]
  
  panel_boot <- pdata.frame(d, index = c("country","year"))
  
  m <- plm(
    form_m3,
    data = panel_boot,
    model = "within"
  )
  
  coef(m)
}

set.seed(123)

boot_res <- boot(
  data = base_panel,
  statistic = boot_fun,
  R = 1000
)

boot_ci <- apply(
  boot_res$t,
  2,
  quantile,
  probs = c(.025,.975),
  na.rm = TRUE
)

print("Intervalos bootstrap")
print(boot_ci)

# ==========================================================
# 3. COMPARACIÓN CON ERRORES ESTÁNDAR CONVENCIONALES
# ==========================================================

coefs_conv <- coeftest(
  m_base,
  vcov = vcovHC(m_base, type = "HC1")
)

print("Errores estándar convencionales")
print(coefs_conv)

# ==========================================================
# 4. UMBRALES ALTERNATIVOS DE IMPUTACIÓN
# ==========================================================

# ejemplo: eliminar observaciones con >20% NA original

base_strict <- base_modelo %>%
  dplyr::select(country, year, all_of(vars_pars)) %>%
  filter(
    rowMeans(is.na(.)) < .20
  ) %>%
  na.omit()

panel_strict <- pdata.frame(base_strict, index=c("country","year"))

m_imput <- plm(
  form_m3,
  data = panel_strict,
  model = "within"
)

coefs_imput <- coeftest(
  m_imput,
  vcov = vcovSCC(m_imput, type="HC1", maxlag=2)
)

print("Modelo con umbral alternativo de imputación")
print(coefs_imput)

# ==========================================================
# 5. EXCLUSIÓN DE PAÍSES EXTREMOS
# ==========================================================

q1 <- quantile(base_panel$DALYs, .25)
q3 <- quantile(base_panel$DALYs, .75)

iqr <- q3 - q1

lim_sup <- q3 + 1.5*iqr
lim_inf <- q1 - 1.5*iqr

paises_extremos <- base_panel %>%
  filter(DALYs > lim_sup | DALYs < lim_inf) %>%
  pull(country) %>%
  unique()

print("Países extremos detectados")
print(paises_extremos)

base_trim <- base_panel %>%
  filter(!country %in% paises_extremos)

panel_trim <- pdata.frame(base_trim, index=c("country","year"))

m_trim <- plm(
  form_m3,
  data = panel_trim,
  model = "within"
)

coefs_trim <- coeftest(
  m_trim,
  vcov = vcovSCC(m_trim, type="HC1", maxlag=2)
)

print("Modelo sin países extremos")
print(coefs_trim)



# -------- 10. CAMBIO TEMPORAL (mixtos + piecewise) ---------------------------
# Asumiendo que candidatas ya está definido
candidatas <- setdiff(names(det), c("country", "iso3", "year", "life_expectancy"))
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

library(car)
vif(m3)

# =====================================================
#   R² marginal y condicional para modelos mixtos
# =====================================================
library(performance)

# Calcular R² Nakagawa (variance explained)
r2_m1 <- performance::r2_nakagawa(m1)
r2_m2 <- performance::r2_nakagawa(m2)
r2_m3 <- performance::r2_nakagawa(m3)

# Mostrar resultados en tabla
r2_tab <- data.frame(
  Modelo = c("m1: Tendencia temporal", 
             "m2: Piecewise temporal",
             "m3: Multivariado estructural"),
  R2_Marginal = c(r2_m1$R2_marginal, r2_m2$R2_marginal, r2_m3$R2_marginal),
  R2_Condicional = c(r2_m1$R2_conditional, r2_m2$R2_conditional, r2_m3$R2_conditional)
)

print(r2_tab, digits = 4)

# Si quieres un formato más bonito
library(knitr)
kable(r2_tab, digits = 4, caption = "R² marginal y condicional (Nakagawa & Schielzeth, 2013)")

# =====================================================
#   ICC: proporción de varianza explicada por país
# =====================================================
library(performance)

# Calcular ICC para cada modelo
icc_m1 <- performance::icc(m1)
icc_m2 <- performance::icc(m2)
icc_m3 <- performance::icc(m3)

# Crear tabla resumen
icc_tab <- data.frame(
  Modelo = c("m1: Tendencia temporal",
             "m2: Piecewise temporal",
             "m3: Multivariado estructural"),
  ICC = c(icc_m1$ICC_adjusted, icc_m2$ICC_adjusted, icc_m3$ICC_adjusted)
)

print(icc_tab, digits = 4)

# Versión bonita
library(knitr)
kable(icc_tab, digits = 4, caption = "Proporción de varianza explicada por país (ICC ajustado)")


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
variables_finales <- c("s1", "lfp", "daly_conflict", "daly_transport", "daly_violence", "indice_infra")

# 3. Construye y corre el Modelo Final (M4)
form_final <- as.formula(glue("DALYs ~ {paste(variables_finales, collapse=' + ')} + (1|iso3)"))
m3f_final <- lmer(form_final, data = X_clean, REML = TRUE) # Usamos REML=TRUE para el reporte final

# 4. Reporte final
print(summary(m3f_final))

# Crear variable combinada de carga DALY metabólica
X$m_daly_metabolico <- scale((X$daly_diabetes_kidney + X$daly_cardio) / 2)

# Modelo m3b: versión depurada de m3
library(lme4)
m3b <- lmer(
  DALYs ~ 
    # Estructura temporal (puedes cambiar por poly(year_c,2) si prefieres)
    s1 + s2 +
    # Socioeconómicos
    desempleo + pib_pc + crecimiento_pib + inflacion_cpi +
    lfp + edu_spend_gdp +
    # Salud y estructura
    gasto_salud_gdp + oop_salud_pct + camas_1000 + vac_measles +
    urbanizacion + indice_infra +
    # Cargas DALY específicas
    m_daly_metabolico + daly_neoplasms + daly_conflict + daly_hiv +
    daly_respiratory + daly_transport + daly_violence + daly_police +
    # Efectos aleatorios
    (1 | iso3),
  data = X,
  REML = FALSE
)

# Revisar colinealidad
library(performance)
check_collinearity(m3b)

# Comparar con el modelo anterior
library(performance)
compare_performance(m3, m3b, rank = TRUE)

# Obtener R² marginal y condicional
r2(m3b)


# Extraer parámetros estandarizados del modelo
params_m3b <- model_parameters(m3b, standardize = "refit", ci = 0.95)


# Limpiar y preparar datos para graficar
df_plot <- params_m3b %>%
  filter(!grepl("Intercept", Parameter, ignore.case = TRUE)) %>%
  mutate(
    label = reorder(Parameter, Coefficient),
    sig05 = ifelse(p < 0.05, "p < 0.05", "ns")
  )



tabla_m3b <- df_plot %>%
  dplyr::select(label, Coefficient, CI_low, CI_high, p) %>%
  arrange(p) %>%
  mutate(across(where(is.numeric), round, 3))

# Ver en consola
print(tabla_m3b, n = Inf)

tabla_apa <- tabla_m3b %>%
  mutate(
    # Añadir estrellas según nivel de significancia
    sig = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ ""
    ),
    beta = sprintf("%.2f", Coefficient),
    CI   = sprintf("[%.2f, %.2f]", CI_low, CI_high),
    p_fmt = ifelse(p < 0.001, "< .001",
                   ifelse(p < 0.01, sprintf("%.3f", p),
                          sprintf("%.3f", p)))
  ) %>%
  dplyr::select(Variable = label, `β` = beta, `IC95%` = CI, `p` = p_fmt, `Sig.` = sig)

library(kableExtra)
# Mostrar en formato APA
kbl(tabla_apa, align = "lcccc",
    caption = "Coeficientes estandarizados (β) del modelo mixto m3b.",
    booktabs = TRUE, digits = 2) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  add_header_above(c(" " = 1, "Efectos fijos" = 4))


library(dplyr)
library(ggplot2)
library(parameters)

# Extraer coeficientes estandarizados
params_m3b <- model_parameters(m3b, standardize = "refit", ci = 0.95)

# Etiquetas legibles
etiquetas <- c(
  "s1" = "Spline 1 (tendencia temporal)",
  "desempleo" = "Desempleo (%)",
  "pib_pc" = "PIB per cápita (USD constantes)",
  "crecimiento_pib" = "Crecimiento PIB (%)",
  "inflacion_cpi" = "Inflación (%)",
  "lfp" = "Participación laboral (%)",
  "edu_spend_gdp" = "Gasto educativo (% PIB)",
  "gasto_salud_gdp" = "Gasto público en salud (% PIB)",
  "oop_salud_pct" = "Gasto de bolsillo en salud (%)",
  "camas_1000" = "Camas hospitalarias (por 1,000 hab.)",
  "vac_measles" = "Cobertura de vacunación (SRP, %)",
  "urbanizacion" = "Urbanización (%)",
  "indice_infra" = "Índice de infraestructura",
  "m_daly_metabolico" = "DALYs cardiometabólicos",
  "daly_neoplasms" = "DALYs neoplasias",
  "daly_conflict" = "DALYs conflicto",
  "daly_hiv" = "DALYs VIH/SIDA",
  "daly_respiratory" = "DALYs respiratorias",
  "daly_transport" = "DALYs transporte",
  "daly_violence" = "DALYs violencia interpersonal",
  "daly_police" = "DALYs violencia policial"
)

# Preparar dataframe
df_plot <- params_m3b %>%
  filter(!grepl("Intercept", Parameter, ignore.case = TRUE)) %>%
  transmute(
    label = etiquetas[Parameter],
    est = Coefficient,
    cil = CI_low,
    cih = CI_high,
    pval = p,
    sig = case_when(
      pval < 0.01 ~ "p < 0.01",
      pval < 0.05 ~ "p < 0.05",
      TRUE        ~ "ns"
    ),
    texto = sprintf("%.2f [%.2f, %.2f]", est, cil, cih)
  ) %>%
  # eliminar etiquetas vacías o NA
  filter(!is.na(label), !is.na(est)) %>%
  arrange(desc(abs(est)))

# Paleta de colores
colores_sig <- c("p < 0.01" = "#e31793", "p < 0.05" = "#ec855f", "ns" = "gray70")

# Plot final
p_forest <- ggplot(df_plot, aes(y = reorder(label, abs(est)))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, color = "gray60") +
  geom_pointrange(aes(x = est, xmin = cil, xmax = cih, color = sig),
                  linewidth = 0.7) +
  geom_text(aes(x = ifelse(est > 0, cih + 0.008, cil - 0.008), label = texto),
            hjust = ifelse(df_plot$est > 0, 0, 1),
            size = 2.8, color = "black") +   # tamaño de texto más pequeño
  scale_color_manual(values = colores_sig, drop = FALSE) +
  labs(
    title = "Efectos fijos del modelo m3b (IC95% Wald)",
    subtitle = "Modelo mixto: DALYs ~ predictores + (1|país); puntos = β, barras = IC95%",
    x = "Coeficiente estandarizado (β)",
    y = NULL,
    color = "Significancia estadística"   # ← aquí va
  ) +
  coord_cartesian(xlim = c(min(df_plot$cil) - 0.05, max(df_plot$cih) + 0.05)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10)
  )

print(p_forest)
ggsave("forest_m3b_sinNA.png", p_forest, width = 9, height = 7, dpi = 300)



# Versión limpia del forest plot (sin coeficientes numéricos, con Times New Roman)
p_forest_clean_tnr <- ggplot(df_plot, aes(y = reorder(label, abs(est)))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, color = "gray60") +
  geom_pointrange(aes(x = est, xmin = cil, xmax = cih, color = sig),
                  linewidth = 0.8) +
  scale_color_manual(values = colores_sig, drop = FALSE) +
  labs(
    title = "Efectos fijos del modelo m3b (IC95% Wald)",
    subtitle = "Modelo mixto: DALYs ~ predictores + (1|país); puntos = β, barras = IC95%",
    x = "Coeficiente estandarizado (β)",
    y = NULL,
    color = "Significancia estadística"
  ) +
  coord_cartesian(
    xlim = c(min(df_plot$cil) - 0.05, max(df_plot$cih) + 0.05)
  ) +
  theme_minimal(base_size = 12, base_family = "Times New Roman") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 11, family = "Times New Roman"),
    legend.title = element_text(size = 12, family = "Times New Roman", face = "bold"),
    axis.text.y = element_text(size = 10, family = "Times New Roman"),
    axis.title.x = element_text(size = 12, family = "Times New Roman"),
    plot.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    plot.subtitle = element_text(size = 11, family = "Times New Roman"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Mostrar y guardar
print(p_forest_clean_tnr)
ggsave("forest_m3b_timesnewroman.png", p_forest_clean_tnr, width = 9, height = 7, dpi = 300)


# Modelo m3b_sin_police: versión depurada sin violencia policial
library(lme4)
m3b_sin_police <- lmer(
  DALYs ~ 
    # Estructura temporal
    s1 + s2 +
    # Socioeconómicos
    desempleo + pib_pc + crecimiento_pib + inflacion_cpi +
    lfp + edu_spend_gdp +
    # Salud y estructura
    gasto_salud_gdp + oop_salud_pct + camas_1000 + vac_measles +
    urbanizacion + indice_infra +
    # Cargas DALY específicas (sin violencia policial)
    m_daly_metabolico + daly_neoplasms + daly_conflict + daly_hiv +
    daly_respiratory + daly_transport + daly_violence +
    # Efectos aleatorios
    (1 | iso3),
  data = X,
  REML = FALSE
)

# Revisar colinealidad
library(performance)
check_collinearity(m3b_sin_police)

# Comparar con el modelo anterior
compare_performance(m3b, m3b_sin_police, rank = TRUE)

# R² marginal y condicional
r2(m3b_sin_police)

# Parámetros estandarizados
library(parameters)
params_m3b_sin_police <- model_parameters(m3b_sin_police, standardize = "refit", ci = 0.95)

# Preparar tabla y gráfico como antes
df_plot <- params_m3b_sin_police %>%
  filter(!grepl("Intercept", Parameter, ignore.case = TRUE)) %>%
  mutate(
    label = reorder(Parameter, Coefficient),
    sig05 = ifelse(p < 0.05, "p < 0.05", "ns")
  )

tabla_m3b <- df_plot %>%
  dplyr::select(label, Coefficient, CI_low, CI_high, p) %>%
  arrange(p) %>%
  mutate(across(where(is.numeric), round, 3))

print(tabla_m3b, n = Inf)

# 🔹 Etiquetas legibles (ya sin violencia policial)
etiquetas <- c(
  "s1" = "Spline 1 (tendencia temporal)",
  "s2" = "Spline 2 (no linealidad temporal)",
  "desempleo" = "Desempleo (%)",
  "pib_pc" = "PIB per cápita (USD constantes)",
  "crecimiento_pib" = "Crecimiento PIB (%)",
  "inflacion_cpi" = "Inflación (%)",
  "lfp" = "Participación laboral (%)",
  "edu_spend_gdp" = "Gasto educativo (% PIB)",
  "gasto_salud_gdp" = "Gasto público en salud (% PIB)",
  "oop_salud_pct" = "Gasto de bolsillo en salud (%)",
  "camas_1000" = "Camas hospitalarias (por 1,000 hab.)",
  "vac_measles" = "Cobertura de vacunación (SRP, %)",
  "urbanizacion" = "Urbanización (%)",
  "indice_infra" = "Índice de infraestructura",
  "m_daly_metabolico" = "DALYs cardiometabólicos",
  "daly_neoplasms" = "DALYs neoplasias",
  "daly_conflict" = "DALYs conflicto armado",
  "daly_hiv" = "DALYs VIH/SIDA",
  "daly_respiratory" = "DALYs respiratorias",
  "daly_transport" = "DALYs transporte",
  "daly_violence" = "DALYs violencia interpersonal"
)

# 🔹 Preparar datos para el gráfico
df_plot <- params_m3b_sin_police %>%
  filter(!grepl("Intercept", Parameter, ignore.case = TRUE)) %>%
  transmute(
    label = etiquetas[Parameter],
    est = Coefficient,
    cil = CI_low,
    cih = CI_high,
    pval = p,
    sig = case_when(
      pval < 0.01 ~ "p < 0.01",
      pval < 0.05 ~ "p < 0.05",
      TRUE        ~ "ns"
    ),
    texto = sprintf("%.2f [%.2f, %.2f]", est, cil, cih)
  ) %>%
  filter(!is.na(label), !is.na(est)) %>%
  arrange(desc(abs(est)))

# 🔹 Paleta de colores
colores_sig <- c("p < 0.01" = "#e31793", "p < 0.05" = "#ec855f", "ns" = "gray70")

# 🔹 Forest plot
library(ggplot2)

p_forest <- ggplot(df_plot, aes(y = reorder(label, abs(est)))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4, color = "gray60") +
  geom_pointrange(aes(x = est, xmin = cil, xmax = cih, color = sig),
                  linewidth = 0.7) +
  scale_color_manual(values = colores_sig, drop = FALSE) +
  labs(
    title = "Efectos fijos del modelo m3b sin violencia policial (IC95% Wald)",
    subtitle = "Modelo mixto: DALYs ~ predictores + (1|país)\nPuntos = β estandarizado, barras = IC95%",
    x = "Coeficiente estandarizado (β)",
    y = NULL,
    color = "Significancia estadística"
  ) +
  coord_cartesian(
    xlim = c(min(df_plot$cil) - 0.05, max(df_plot$cih) + 0.05)
  ) +
  theme_minimal(base_size = 12, base_family = "Times New Roman") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 11, family = "Times New Roman"),
    legend.title = element_text(size = 12, face = "bold", family = "Times New Roman"),
    axis.text.y = element_text(size = 10, family = "Times New Roman"),
    axis.text.x = element_text(size = 10, family = "Times New Roman"),
    axis.title.x = element_text(size = 12, family = "Times New Roman"),
    plot.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    plot.subtitle = element_text(size = 11, family = "Times New Roman"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Mostrar y guardar
print(p_forest)
ggsave("forest_m3b_sin_police_simple.png", p_forest, width = 9, height = 7, dpi = 300)


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

