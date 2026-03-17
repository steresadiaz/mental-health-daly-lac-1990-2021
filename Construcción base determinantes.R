# ======================================================================
# Descarga robusta + Limpieza panel + Imputación (<=10%) + Winsorización
# ======================================================================

to_install <- c("WDI","wbstats","dplyr","janitor","countrycode","here",
                "readr","zoo","tibble","purrr","tidyr")
new_pk <- setdiff(to_install, rownames(installed.packages()))
if (length(new_pk)) install.packages(new_pk, Ncpus = parallel::detectCores())
library(dplyr)       # o bien
library(tidyverse)   # incluye dplyr, ggplot2, tidyr, etc.

library(WDI); library(wbstats)
library(dplyr); library(janitor); library(countrycode)
library(here); library(readr)
library(zoo); library(tibble); library(purrr); library(tidyr)

# --- Países de ALC (para WDI/wbstats) ---
LAC_ISO3 <- wbstats::wb_countries() %>%
  dplyr::filter(region == "Latin America & Caribbean", !is.na(iso3c)) %>%
  dplyr::pull(iso3c) %>% unique()


options(timeout = max(600, getOption("timeout", 60)))
options(download.file.method = "libcurl")

show <- function(x, n = Inf, width = Inf) {
  if (!inherits(x, "tbl_df")) x <- tibble::as_tibble(x)
  print(x, n = n, width = width)
  invisible(x)
}


# 1) Indicadores (extendidos) --------------------------------------------------
ind <- c(
  # Socioeconómico / laboral
  pobreza           = "SI.POV.NAHC",
  pobreza_215       = "SI.POV.DDAY",
  pobreza_365       = "SI.POV.LMIC",
  pobreza_685       = "SI.POV.UMIC",
  gini              = "SI.POV.GINI",
  pib_pc_const15    = "NY.GDP.PCAP.KD",
  crecimiento_pib   = "NY.GDP.MKTP.KD.ZG",
  inflacion_cpi     = "FP.CPI.TOTL.ZG",
  lfp               = "SL.TLF.CACT.ZS",
  desempleo         = "SL.UEM.TOTL.ZS",
  desempleo_joven   = "SL.UEM.1524.ZS",
  neet_joven        = "SL.UEM.NEET.ZS",     # NEET total correcto
  neet_joven_fe     = "SL.UEM.NEET.FE.ZS",  # fallback femenino
  neet_joven_ma     = "SL.UEM.NEET.MA.ZS",  # fallback masculino
  
  # Seguridad / salud
  homicidios_100k   = "VC.IHR.PSRC.P5",
  escolaridad_exp   = "SE.SCH.LIFE",
  gasto_salud_gdp   = "SH.XPD.CHEX.GD.ZS",
  gasto_salud_pc    = "SH.XPD.CHEX.PC.CD",
  gob_salud_gdp     = "SH.XPD.GHED.GD.ZS",
  oop_salud_pct     = "SH.XPD.OOPC.CH.ZS",
  camas_1000        = "SH.MED.BEDS.ZS",
  enfermeras_1000   = "SH.MED.NUMW.P3",
  uhc_indice        = "SH.UHC.SRVS.CV.XD",
  vac_measles       = "SH.IMM.MEAS",
  
  # Infraestructura / acceso
  urbanizacion      = "SP.URB.TOTL.IN.ZS",
  electrificacion   = "EG.ELC.ACCS.ZS",
  agua_basica       = "SH.H2O.BASW.ZS",
  saneamiento_basico= "SH.STA.BASS.ZS",
  internet_usuarios = "IT.NET.USER.ZS",
  banda_ancha_100   = "IT.NET.BBND.P2"
)

start_year <- 2000; end_year <- 2024

# 2) Descarga robusta ----------------------------------------------------------
# trocea países e indicadores para URLs cortas y menos fallos
chunk_vec <- function(x, k) split(x, ceiling(seq_along(x)/k))

fetch_chunk <- function(countries, indicators, start, end, tries = 4, sleep_s = 2) {
  # intenta WDI con reintentos
  for (i in seq_len(tries)) {
    res <- tryCatch(
      WDI(country = countries, indicator = indicators, start = start, end = end, extra = TRUE),
      error = function(e) NULL, warning = function(w) NULL
    )
    if (!is.null(res)) return(res)
    message(sprintf("WDI fallo intento %d/%d. Reintentando en %ds...", i, tries, sleep_s*i))
    Sys.sleep(sleep_s * i)
  }
  # fallback: wbstats
  message("Usando fallback wbstats::wb_data() para este trozo...")
  wb <- tryCatch(
    wbstats::wb_data(indicators, countries = countries,
                     start_date = start, end_date = end, return_wide = TRUE),
    error = function(e) NULL
  )
  if (is.null(wb)) return(NULL)
  # homogeniza columnas clave a estilo WDI
  wb <- wb %>%
    rename(country = country, iso2c = iso2c) %>%
    mutate(year = as.integer(date)) %>%
    select(any_of(c("country","iso2c","year", indicators)))
  wb
}

# descarga en trozos
iso_chunks <- chunk_vec(LAC_ISO3, 8)                # 8 países por llamada
ind_chunks <- chunk_vec(unname(ind), 6)             # 6 indicadores por llamada

all_list <- list()
idx <- 1L
for (ic in seq_along(iso_chunks)) {
  for (jc in seq_along(ind_chunks)) {
    message(sprintf("Descargando bloque países %d/%d, indicadores %d/%d...",
                    ic, length(iso_chunks), jc, length(ind_chunks)))
    chunk <- fetch_chunk(iso_chunks[[ic]], ind_chunks[[jc]],
                         start = start_year, end = end_year)
    if (!is.null(chunk)) {
      chunk$..block_id <- idx
      all_list[[length(all_list)+1]] <- chunk
    } else {
      message("⚠️ Bloque vacío por error persistente.")
    }
    idx <- idx + 1L
  }
}
stopifnot(length(all_list) > 0)
dat_raw <- bind_rows(all_list)

# 3) Limpieza y estandarización ------------------------------------------------
df <- dat_raw %>% clean_names()

# --- DEDUP: colapsar country-iso2c-year tomando el primer no-NA por columna ---
needed_vars <- c(
  "si_pov_nahc","si_pov_dday","si_pov_lmic","si_pov_umic","si_pov_gini",
  "ny_gdp_pcap_kd","ny_gdp_mktp_kd_zg","fp_cpi_totl_zg","sl_tlf_cact_zs",
  "sl_uem_totl_zs","sl_uem_1524_zs","sl_uem_neet_zs","sl_uem_neet_fe_zs","sl_uem_neet_ma_zs",
  "vc_ihr_psrc_p5","se_sch_life","sh_xpd_chex_gd_zs","sh_xpd_chex_pc_cd",
  "sh_xpd_ghed_gd_zs","sh_xpd_oopc_ch_zs","sh_med_beds_zs","sh_med_numw_p3",
  "sh_uhc_srvs_cv_xd","sh_imm_meas",
  "sp_urb_totl_in_zs","eg_elc_accs_zs","sh_h2o_basw_zs","sh_sta_bass_zs",
  "it_net_user_zs","it_net_bbnd_p2"
)
for (nm in setdiff(needed_vars, names(df))) df[[nm]] <- NA_real_

take_first_non_na <- function(x) {
  y <- x[!is.na(x)]
  if (length(y) == 0) NA_real_ else as.numeric(y[1])
}

df_dedup <- df %>%
  dplyr::group_by(country, iso2c, year) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(needed_vars), take_first_non_na), .groups = "drop")

dups_after <- df_dedup %>% dplyr::count(country, iso2c, year) %>% dplyr::filter(n > 1)
if (nrow(dups_after) == 0) message("✅ Duplicados resueltos (1 fila por país–año).")

# NEET total con fallback F/M
neet_mix <- with(df_dedup, {
  fm <- rowMeans(cbind(sl_uem_neet_fe_zs, sl_uem_neet_ma_zs), na.rm = TRUE); fm[is.nan(fm)] <- NA_real_
  ifelse(!is.na(sl_uem_neet_zs), sl_uem_neet_zs, fm)
})

det_raw <- df_dedup %>%
  dplyr::mutate(
    iso3 = countrycode::countrycode(iso2c, "iso2c", "iso3c"),
    pobreza_best = dplyr::coalesce(si_pov_nahc, si_pov_umic, si_pov_lmic, si_pov_dday),
    neet_joven   = neet_mix
  ) %>%
  dplyr::transmute(
    country, iso3, year,
    # socio/laboral
    pobreza = pobreza_best, gini = si_pov_gini,
    desempleo = sl_uem_totl_zs, desempleo_joven = sl_uem_1524_zs, neet_joven,
    # seguridad/salud
    violencia = vc_ihr_psrc_p5, educacion = se_sch_life,
    gasto_salud_gdp = sh_xpd_chex_gd_zs, gasto_salud_pc = sh_xpd_chex_pc_cd,
    gob_salud_gdp = sh_xpd_ghed_gd_zs, oop_salud_pct = sh_xpd_oopc_ch_zs,
    camas_1000 = sh_med_beds_zs, enfermeras_1000 = sh_med_numw_p3,
    uhc_indice = sh_uhc_srvs_cv_xd,  vac_measles = sh_imm_meas,
    # infraestructura
    urbanizacion = sp_urb_totl_in_zs, electrificacion = eg_elc_accs_zs,
    agua_basica = sh_h2o_basw_zs, saneamiento_basico = sh_sta_bass_zs,
    internet_usuarios = it_net_user_zs, banda_ancha_100 = it_net_bbnd_p2,
    # macro
    pib_pc = ny_gdp_pcap_kd, crecimiento_pib = ny_gdp_mktp_kd_zg, inflacion_cpi = fp_cpi_totl_zg,
    lfp = sl_tlf_cact_zs
  ) %>% dplyr::arrange(iso3, year)


# ==================== GINI SEDLAC/OWID → PANEL =====================
to_install <- c("readr","dplyr","janitor","stringr","countrycode","here")
new_pk <- setdiff(to_install, rownames(installed.packages()))
if (length(new_pk)) install.packages(new_pk, Ncpus = parallel::detectCores())

library(readr); library(dplyr); library(janitor)
library(stringr); library(countrycode); library(here)

# Panel base
panel <- det_raw
stopifnot(is.data.frame(panel)); panel <- panel %>% mutate(iso3 = dplyr::coalesce(.data$iso3))


# Diccionario para nombres “raros” de OWID
custom_names <- c(
  "Bahamas, The" = "BHS",
  "Venezuela, RB" = "VEN",
  "Bolivia" = "BOL",
  "Dominican Republic" = "DOM",
  "St. Kitts and Nevis" = "KNA",
  "St. Lucia" = "LCA",
  "St. Vincent and the Grenadines" = "VCT",
  "Antigua and Barbuda" = "ATG",
  "Trinidad and Tobago" = "TTO"
)

# 1) Descargar OWID
url <- "https://ourworldindata.org/grapher/income-inequality-in-latin-america.csv"
gini_raw <- readr::read_csv(url, show_col_types = FALSE) %>% clean_names()

# 2) Normalizar "(urban)" y mapear ISO3
gini_norm <- gini_raw %>%
  transmute(
    country_raw = entity,
    code = na_if(code, ""),
    year,
    gini_val = gini_coefficient_equivalized_income
  ) %>%
  mutate(
    # flag si es observación urbana
    is_urban = str_detect(country_raw, "\\(\\s*urban\\s*\\)"),
    # quitar paréntesis al FINAL (p. ej. "Argentina (urban)")
    country = str_trim(str_remove(country_raw, "\\s*\\([^)]*\\)\\s*$")),
    # ISO3: usa code si viene; si no, dedúcelo desde country
    iso3 = coalesce(
      code,
      countrycode(country, "country.name.en", "iso3c", custom_match = custom_names)
    )
  ) %>%
  # filtra ALC y ventana
  filter(!grepl("^OWID_", country_raw),
         iso3 %in% LAC_ISO3,
         dplyr::between(year, 2000, 2022)) %>%
  mutate(gini_val = ifelse(gini_val <= 1, 100*gini_val, gini_val)) %>%
  arrange(iso3, year, is_urban)  # no-urban primero

# 3) Resolver duplicados iso3×year (prefiere no-urban; si no hay, usa urban)
gini_sedlac <- gini_norm %>%
  group_by(iso3, year) %>%
  summarise(
    gini_sedlac = { x <- gini_val[!is.na(gini_val)]; if (length(x)) x[1] else NA_real_ },
    .groups = "drop"
  )

# 4) Unir al panel y construir 'gini' final (WDI → respaldo con SEDLAC)
has_wdi_gini <- "gini" %in% names(panel)
# Requisitos mínimos
stopifnot(all(c("iso3","year") %in% names(panel)),
          all(c("iso3","year","gini_sedlac") %in% names(gini_sedlac)))

# 1) Garantiza unicidad y tipos
gini_sedlac <- gini_sedlac %>%
  dplyr::group_by(iso3, year) %>%
  dplyr::summarise(gini_sedlac = dplyr::first(as.numeric(gini_sedlac)), .groups = "drop")

# 2) Une y arma columnas
# --- Une SEDLAC al panel y arma gini final (prioridad WDI -> SEDLAC) ---
panel_plus <- panel %>%
  dplyr::left_join(gini_sedlac, by = c("iso3","year")) %>%
  { df <- .;                                   # etapa intermedia para "inyectar" columnas si faltan
  if (!"gini_wdi" %in% names(df)) df$gini_wdi <- NA_real_
  if (!"gini"     %in% names(df)) df$gini     <- NA_real_  # por si tu WDI venía como otra col
  df
  } %>%
  dplyr::mutate(
    gini_wdi = dplyr::coalesce(gini_wdi, gini),            # conserva WDI en gini_wdi
    gini     = dplyr::coalesce(gini_sedlac, gini_wdi),     # PRIORIDAD: SEDLAC
    gini_source = dplyr::case_when(
      !is.na(gini_sedlac) &  is.na(gini_wdi) ~ "SEDLAC",
      is.na(gini_sedlac) & !is.na(gini_wdi) ~ "WDI",
      !is.na(gini_sedlac) & !is.na(gini_wdi) ~ "BOTH",
      TRUE                                   ~ NA_character_
    )
  ) %>%
  dplyr::relocate(gini, gini_source, dplyr::any_of(c("gini_wdi","gini_sedlac")), .after = year)

message(sprintf("Cobertura Gini final (no-NA): %.1f%%", 100*mean(!is.na(panel_plus$gini))))


# 5) (opcional) Guardar
dir.create(here("data"), showWarnings = FALSE)
write_csv(gini_sedlac, here("data","gini_sedlac_owid.csv"))
write_csv(panel_plus,  here("data","panel_con_gini.csv"))

panel_plus %>% dplyr::filter(is.na(gini)) %>% dplyr::distinct(country) %>% dplyr::arrange(country) %>% print(n=Inf)

panel_plus <- readr::read_csv(here::here("data", "panel_con_gini.csv"))

# Si tu pipeline usa det_cleanwoyears:
det_raw <- panel_plus


# =================== PANEL QA + IMPUTACIÓN + WINSORIZACIÓN ====================

message("Tratando la base como PANEL iso3 × year…")

# Utilidades
impute_series_if_lt <- function(x, time, thr = 0.10) {
  if (!is.numeric(x)) return(x)
  prop_na <- mean(is.na(x))
  if (all(is.na(x)) || sum(!is.na(x)) < 2) return(x)
  if (isTRUE(prop_na <= thr)) {
    y <- na.approx(x, x = time, na.rm = FALSE)
    y <- na.locf(y, na.rm = FALSE)
    y <- na.locf(y, fromLast = TRUE, na.rm = FALSE)
    return(y)
  }
  x
}

# Función de winsor sin tocar NA
winsorize_vec <- function(x, p = c(0.01, 0.99)) {
  if (!is.numeric(x)) return(x)
  n_obs <- sum(!is.na(x)); if (n_obs < 3) return(x)
  qs <- stats::quantile(x, probs = p, na.rm = TRUE, names = FALSE)
  if (any(is.na(qs))) return(x)
  x2 <- x
  x2[!is.na(x) & x < qs[1]] <- qs[1]
  x2[!is.na(x) & x > qs[2]] <- qs[2]
  x2
}

# Variables numéricas
measure_cols <- setdiff(names(det_raw), c("country","iso3","year"))
message("Variables numéricas detectadas (", length(measure_cols), "):")
print(measure_cols)

# Duplicados panel
panel_dups <- det_raw %>% count(iso3, year) %>% filter(n > 1)
if (nrow(panel_dups) > 0) {
  message("⚠️ Duplicados por (iso3, year) encontrados (muestra):")
  print(head(panel_dups %>% arrange(desc(n)), 20))
} else {
  message("✅ Sin duplicados por (iso3, year).")
}

# Cobertura por país (muestra top)
panel_cov <- det_raw %>%
  mutate(miss_row = rowMeans(across(all_of(measure_cols), ~is.na(.x)))) %>%
  group_by(iso3) %>%
  summarise(min_year = min(year, na.rm = TRUE),
            max_year = max(year, na.rm = TRUE),
            n_rows   = dplyr::n(),
            na_row_pct_avg = mean(miss_row)*100,
            .groups = "drop")
message("Cobertura panel por país (top 10 con más NA promedio en fila):")
print(panel_cov %>% arrange(desc(na_row_pct_avg)) %>% head(10))

# %NA antes
# % NA por AÑO
det_raw %>% dplyr::select(year, dplyr::where(is.numeric)) %>%
  tidyr::pivot_longer(-year, values_to="v") %>%
  dplyr::summarise(pct_na = 100*mean(is.na(v)), .by = year)%>%
  show()


# % NA por PAÍS (detecta iso3 o iso3c)
det_raw %>%
  dplyr::select(iso3, dplyr::where(is.numeric)) %>%
  tidyr::pivot_longer(-iso3, values_to = "v") %>%
  dplyr::summarise(pct_na = 100*mean(is.na(v)), .by = iso3) %>%
  show()

# % NA por INDICADOR
det_raw %>% dplyr::select(dplyr::where(is.numeric)) %>%
  dplyr::summarise(dplyr::across(everything(), ~100*mean(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to="indicador", values_to="pct_na") %>%
  dplyr::arrange(dplyr::desc(pct_na))%>%
  show()

#Tendencia temporal (si p<0.05, hay patrón en el tiempo ⇒ no MCAR):
det_raw %>% dplyr::select(year, dplyr::where(is.numeric), -dplyr::starts_with("block")) %>%
  tidyr::pivot_longer(-year, names_to="var", values_to="v") %>%
  dplyr::summarise(p = suppressWarnings(cor.test(as.numeric(is.na(v)), year)$p.value), .by = var) %>%
  dplyr::arrange(p)

#Dependencia por país (si p<0.05, la falta depende del país ⇒ no MCAR):
det_raw %>% dplyr::mutate(iso = factor(dplyr::coalesce(iso3))) %>%
  dplyr::select(iso, dplyr::where(is.numeric), -year, -dplyr::starts_with("block")) %>%
  tidyr::pivot_longer(-iso, names_to="var", values_to="v") %>%
  dplyr::summarise(p = suppressWarnings(chisq.test(table(iso, is.na(v)))$p.value), .by = var) %>%
  dplyr::arrange(p)


# ---------- %NA ANTES (sobre det_raw, ya deduplicado y recortado si aplica) ---
measure_cols <- setdiff(names(det_raw), c("country","iso3","year"))

missing_before_long <- det_raw %>%
  summarise(across(all_of(measure_cols), ~mean(is.na(.x))*100)) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "na_pct") %>%
  arrange(desc(na_pct))
message("Variables con >10% NA ANTES de cualquier tratamiento:")
print(missing_before_long %>% filter(na_pct > 10), n=Inf)

# ---------- (A) WINSORIZACIÓN PRIMERO (p1–p99) --------------------------------
# Calcula umbrales p1/p99 por variable usando SOLAMENTE valores observados
winsor_bounds_pre <- purrr::map(measure_cols, ~{
  x <- det_raw[[.x]]
  tryCatch(stats::quantile(x, c(0.01, 0.99), na.rm = TRUE, names = FALSE),
           error = function(e) c(NA_real_, NA_real_))
})
names(winsor_bounds_pre) <- measure_cols

# Diagnóstico: ¿cuántos valores se recortarían?
winsor_counts_pre <- purrr::map2_dfr(measure_cols, winsor_bounds_pre, function(v, b) {
  x <- det_raw[[v]]
  if (any(is.na(b))) {
    tibble(variable = v, p1 = NA_real_, p99 = NA_real_,
           n_low = NA_integer_, n_high = NA_integer_,
           n_total = sum(!is.na(x)), pct_clamped = NA_real_)
  } else {
    n_low  <- sum(!is.na(x) & x <  b[1])
    n_high <- sum(!is.na(x) & x >  b[2])
    n_tot  <- sum(!is.na(x))
    tibble(variable = v, p1 = b[1], p99 = b[2],
           n_low = n_low, n_high = n_high,
           n_total = n_tot,
           pct_clamped = 100*(n_low + n_high)/ifelse(n_tot==0, NA, n_tot))
  }
})
message("Top variables con mayor recorte por winsorización (ANTES de imputar):")
print(winsor_counts_pre %>% arrange(desc(pct_clamped)) %>% head(15))



# Aplica winsor a det_raw (NO rellena NA)
det_winsor <- det_raw %>% mutate(across(all_of(measure_cols), winsorize_vec))

# ---------- %NA TRAS WINSOR (debe ser igual a before; winsor no imputa) -------
missing_after_winsor <- det_winsor %>%
  summarise(across(all_of(measure_cols), ~mean(is.na(.x))*100)) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "na_pct") %>%
  arrange(desc(na_pct))
message("Chequeo: %NA tras winsor (debe ser igual a 'antes'):")
print(missing_after_winsor %>% filter(na_pct > 10), n=Inf)

# ---------- (B) IMPUTACIÓN DESPUÉS (por país; umbral y brecha por variable) ----------------

# 0) Detecta columnas numéricas del panel si no está `measure_cols`
if (!exists("measure_cols")) {
  measure_cols <- setdiff(names(det_winsor), c("country","iso3","year"))
  measure_cols <- measure_cols[vapply(det_winsor[measure_cols], is.numeric, logical(1))]
}

# === Clasificación data-driven de variables según volatilidad y monotonicidad ===
# Requiere: dplyr, tidyr, purrr, tibble, zoo (ya los tienes)
base_for_scores <- det_winsor   # usa tu panel antes de imputar (o el que prefieras)
stopifnot(all(c("iso3","year") %in% names(base_for_scores)))

# 1) Detecta columnas numéricas candidatas
measure_cols <- setdiff(names(base_for_scores), c("country","iso3","year"))
measure_cols <- measure_cols[vapply(base_for_scores[measure_cols], is.numeric, logical(1))]

# --- helpers "seguros" ---
safe_sd   <- function(z) if (sum(!is.na(z)) >= 2) stats::sd(z, na.rm = TRUE) else NA_real_
safe_med  <- function(z) if (all(is.na(z))) NA_real_ else stats::median(z, na.rm = TRUE)
safe_spr  <- function(y, x) {
  ok <- stats::complete.cases(y, x)
  if (sum(ok) >= 3 && length(unique(x[ok])) > 1) suppressWarnings(stats::cor(y[ok], x[ok], method="spearman"))
  else NA_real_
}

# --- métricas por país–variable (robusto a NA/constantes) ---
var_country_stats <- base_for_scores %>%
  dplyr::select(iso3, year, dplyr::all_of(measure_cols)) %>%
  tidyr::pivot_longer(-c(iso3,year), names_to="variable", values_to="x") %>%
  dplyr::group_by(iso3, variable) %>%
  dplyr::arrange(year, .by_group = TRUE) %>%
  dplyr::mutate(
    dx  = x - dplyr::lag(x),
    rel = abs(dx) / (abs(dplyr::lag(x)) + 1e-9)  # evita división por 0
  ) %>%
  dplyr::summarise(
    n_obs        = sum(!is.na(x)),
    na_prop      = mean(is.na(x)),
    sd_diff      = safe_sd(dx),
    med_rel      = safe_med(rel),
    monotone_rho = safe_spr(year, x),
    .groups = "drop"
  )

stopifnot(exists("var_country_stats"), all(c("variable","med_rel","monotone_rho") %in% names(var_country_stats)))
suppressPackageStartupMessages({library(dplyr); library(tibble)})

# 1) Resumen por variable (medianas entre países) y cortes por terciles
var_summary <- var_country_stats %>%
  group_by(variable) %>%
  summarise(
    monotone = median(monotone_rho, na.rm = TRUE),
    vol      = median(med_rel,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # fallback por si alguna queda NA
    monotone = ifelse(is.na(monotone), 0, monotone),
    vol      = ifelse(is.na(vol),      median(vol, na.rm = TRUE), vol)
  )

q_mono <- quantile(var_summary$monotone, probs = c(1/3, 2/3), na.rm = TRUE)
q_vol  <- quantile(var_summary$vol,      probs = c(1/3, 2/3), na.rm = TRUE)

class_tbl <- var_summary %>%
  mutate(class = case_when(
    vol >= q_vol[2]                              ~ "volatile",
    vol <= q_vol[1] & monotone >= q_mono[2]      ~ "slow",
    TRUE                                         ~ "mix"
  ))

# 2) Listas (por si las quieres usar después)
slow_vars <- intersect(class_tbl$variable[class_tbl$class == "slow"],     measure_cols)
mix_vars  <- intersect(class_tbl$variable[class_tbl$class == "mix"],      measure_cols)
volatile  <- intersect(class_tbl$variable[class_tbl$class == "volatile"], measure_cols)

# 3) Umbrales y límites por variable
thr_base <- 0.15
thr_by_var <- tibble(variable = measure_cols) %>%
  mutate(
    thr = case_when(
      variable %in% volatile  ~ 0.10,
      variable %in% slow_vars ~ 0.30,
      variable %in% mix_vars  ~ 0.20,
      TRUE                    ~ thr_base
    ),
    max_gap = case_when(
      variable %in% volatile  ~ 1L,
      variable %in% slow_vars ~ 3L,
      variable %in% mix_vars  ~ 2L,
      TRUE                    ~ 2L
    )
  )

# (Opcional) ver un resumen ordenado
thr_by_var %>% left_join(class_tbl, by = "variable") %>%
  arrange(match(variable, measure_cols)) %>%
  print(n = 99)

# ================= Gráfico: Volatilidad vs Monotonicidad =================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ggrepel); library(scales); library(viridis)
})

stopifnot(exists("thr_by_var"))
stopifnot(exists("var_country_stats"))  # viene del paso de métricas por país

# 1) Reconstruir resumen por variable (medianas entre países)
var_summary <- var_country_stats %>%
  group_by(variable) %>%
  summarise(
    mono_med    = median(monotone_rho, na.rm = TRUE),   # monotonicidad
    med_rel_med = median(med_rel,      na.rm = TRUE),   # volatilidad relativa
    .groups = "drop"
  ) %>%
  mutate(
    mono_med    = ifelse(is.na(mono_med),    0, mono_med),
    med_rel_med = ifelse(is.na(med_rel_med), median(med_rel_med, na.rm=TRUE), med_rel_med)
  )

# 2) Cortes (terciles) y clase (slow/mix/volatile) si no vienen ya
x_cuts <- quantile(var_summary$mono_med,    c(1/3, 2/3), na.rm = TRUE)
y_cuts <- quantile(var_summary$med_rel_med, c(1/3, 2/3), na.rm = TRUE)

class_tbl <- var_summary %>%
  mutate(class = case_when(
    med_rel_med >= y_cuts[2]                           ~ "volatile",
    med_rel_med <= y_cuts[1] & mono_med >= x_cuts[2]   ~ "slow",
    TRUE                                               ~ "mix"
  ))

# 3) Data para el gráfico: umbrales + métricas + clase
plot_df <- thr_by_var %>%
  left_join(class_tbl, by = "variable") %>%
  mutate(class = factor(class, levels = c("slow","mix","volatile")))

# 4) Plot
p <- ggplot(plot_df, aes(x = mono_med, y = med_rel_med, color = class, label = variable)) +
  geom_point(size = 3, alpha = .9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 100, seed = 123) +
  geom_vline(xintercept = x_cuts, linetype = "dashed", linewidth = .3) +
  geom_hline(yintercept = y_cuts, linetype = "dashed", linewidth = .3) +
  scale_y_continuous("Volatilidad (|Δ|/|x_{t-1}| mediana)", labels = percent_format(accuracy = 1)) +
  scale_x_continuous("Monotonicidad (ρ de Spearman por país → mediana entre países)",
                     limits = c(min(plot_df$mono_med, na.rm=TRUE)-.05, 1.05)) +
  scale_color_viridis_d(name = "Clase", end = .9) +
  labs(title = "Clasificación automática de variables (panel TSCS)",
       subtitle = "Eje Y: cambio relativo mediano |Δ|/|x_{t-1}| (volatilidad) • Eje X: monotonicidad (ρ)\nLíneas punteadas: terciles usados para clasificar → slow / mix / volatile") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "right")

print(p)


# (opcional) guardar
dir.create(here::here("output"), showWarnings = FALSE)
ggsave(here::here("output","clasificacion_volatilidad_monotonicidad.png"),
       p, width = 12, height = 8, dpi = 220)

p_counts <- plot_df %>% count(class) %>%
  ggplot(aes(x = class, y = n, fill = class)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_d(guide = "none", option = "D", end = .9) +
  labs(x = NULL, y = "Número de variables", title = "Conteo por clase") +
  theme_minimal(base_size = 12)
print(p_counts)


# 2) Umbrales y límites por variable
thr_base <- 0.15
thr_by_var <- tibble::tibble(
  variable = measure_cols,
  thr = dplyr::case_when(
    variable %in% volatile  ~ 0.10,
    variable %in% slow_vars ~ 0.30,
    variable %in% mix_vars  ~ 0.20,
    TRUE                    ~ thr_base
  ),
  max_gap = dplyr::case_when(
    variable %in% volatile  ~ 1L,
    variable %in% slow_vars ~ 3L,
    variable %in% mix_vars  ~ 2L,
    TRUE                    ~ 2L
  )
)

# 3) Función de imputación con guardas
impute_guarded <- function(x, time, thr = 0.15, max_gap = 2, min_obs = 5, allow_end = FALSE){
  if (!is.numeric(x)) return(x)
  if (sum(!is.na(x)) < min_obs) return(x)             # muy pocos datos
  if (mean(is.na(x)) > thr) return(x)                 # demasiados NA globales
  y <- zoo::na.approx(x, x = time, na.rm = FALSE, maxgap = max_gap)  # sólo huecos cortos
  if (!allow_end) {
    # evita extrapolar más de 1 año en extremos
    y <- zoo::na.locf(y, na.rm = FALSE, maxgap = 1)
    y <- zoo::na.locf(y, fromLast = TRUE, na.rm = FALSE, maxgap = 1)
  } else {
    y <- zoo::na.locf(zoo::na.locf(y, na.rm = FALSE), fromLast = TRUE, na.rm = FALSE)
  }
  y
}

# 4) Imputación por país con parámetros por variable
det_imp <- det_winsor %>%
  dplyr::group_by(iso3) %>%
  dplyr::arrange(year, .by_group = TRUE) %>%
  dplyr::mutate(dplyr::across(all_of(measure_cols), ~{
    v   <- dplyr::cur_column()
    par <- dplyr::filter(thr_by_var, variable == v)
    impute_guarded(.x, time = year, thr = par$thr, max_gap = par$max_gap, min_obs = 5, allow_end = FALSE)
  })) %>%
  dplyr::ungroup()

# 5) Resumen de imputaciones
imp_summary <- purrr::map_dfr(measure_cols, function(v) {
  before <- det_winsor[[v]]; after <- det_imp[[v]]
  n_imp  <- sum(is.na(before) & !is.na(after))
  tibble::tibble(
    variable       = v,
    n_imputed      = n_imp,
    pct_panel      = 100 * n_imp / length(before),          # % del panel imputado (tu métrica original)
    pct_de_los_NA  = 100 * n_imp / sum(is.na(before)),      # % de NA originales que se rellenaron
    pct_de_finales = 100 * n_imp / sum(!is.na(after))       # % de valores finales que provienen de imputación
  )
}) %>% dplyr::arrange(dplyr::desc(n_imputed))

message("Imputaciones por variable (winsor → imputación):")
print(imp_summary %>% dplyr::filter(n_imputed > 0) %>% dplyr::select(variable, n_imputed, pct_panel, pct_de_los_NA, pct_de_finales) %>% head(20), n=Inf)

# (Opcional) Imputaciones por país
imp_by_iso3 <- det_winsor %>%
  dplyr::select(iso3, year, dplyr::all_of(measure_cols)) %>%
  tidyr::pivot_longer(-c(iso3,year), names_to = "var", values_to = "v_before") %>%
  dplyr::bind_cols(
    det_imp %>%
      dplyr::select(iso3, year, dplyr::all_of(measure_cols)) %>%
      tidyr::pivot_longer(-c(iso3,year), names_to = "var2", values_to = "v_after") %>%
      dplyr::select(v_after)
  ) %>%
  dplyr::mutate(imp_flag = is.na(v_before) & !is.na(v_after)) %>%
  dplyr::summarise(n_imputed = sum(imp_flag), .by = iso3) %>%
  dplyr::arrange(dplyr::desc(n_imputed))

message("Países con más celdas imputadas (top 10):")
print(head(imp_by_iso3, 10))

# 6) %NA después de imputar (diagnóstico)
missing_after_imp <- det_imp %>%
  dplyr::summarise(dplyr::across(all_of(measure_cols), ~mean(is.na(.x))*100)) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "na_pct") %>%
  dplyr::arrange(dplyr::desc(na_pct))

message("Variables con >10% NA después de imputar:")
print(missing_after_imp %>% dplyr::filter(na_pct > 10), n=Inf)


det_clean <- det_imp %>% mutate(across(all_of(measure_cols), winsorize_vec))

# ---------- REPORTES RESUMEN --------------------------------------------------
# %NA final
missing_after_all <- det_clean %>%
  summarise(across(all_of(measure_cols), ~mean(is.na(.x))*100)) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "na_pct") %>%
  arrange(desc(na_pct))

message("Variables con >10% NA FINALES (winsor→impute→clip):")
print(missing_after_all %>% filter(na_pct > 10), n=Inf)

# Tabla resumen: %NA antes -> tras winsor -> tras imputar -> final
na_report <- missing_before_long %>%
  rename(na_before = na_pct) %>%
  left_join(missing_after_winsor, by = "variable") %>%
  rename(na_after_winsor = na_pct) %>%
  left_join(missing_after_imp,     by = "variable") %>%
  rename(na_after_imp = na_pct) %>%
  left_join(missing_after_all,     by = "variable") %>%
  rename(na_final = na_pct) %>%
  arrange(desc(na_before))
message("Resumen %NA por variable (antes -> winsor -> imputar -> final):")
print(head(na_report, 30), n= Inf)

# ---------- GUARDADOS ---------------------------------------------------------
dir.create(here::here("data"), showWarnings = FALSE)
readr::write_csv(det_raw,   here::here("data","determinantes_wb_raw.csv"))

# 5) Chequeos finales ----------------------------------------------------------
message("Cobertura temporal final por país:")
det_clean %>%
  group_by(iso3) %>%
  summarise(min_year = min(year, na.rm=TRUE),
            max_year = max(year, na.rm=TRUE),
            n_rows = dplyr::n(),
            .groups = "drop") %>%
  arrange(iso3) %>%
  print(n = 99)

na_report <- missing_before_long %>%
  rename(na_before = na_pct) %>%
  left_join(missing_after_imp, by = "variable") %>%
  rename(na_after_imp = na_pct) %>%
  left_join(missing_after_all, by = "variable") %>%
  rename(na_final = na_pct) %>%
  arrange(desc(na_before))

message("Resumen %NA por variable (antes -> después imputación -> final):")
print(head(na_report, 30))

# Mantén solo indicadores con cobertura global ≥ 40%
VAR_COV_THR <- 0.40
keep_vars <- det_clean %>% dplyr::select(dplyr::where(is.numeric), -year) %>%
  dplyr::summarise(dplyr::across(everything(), ~mean(!is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to="var", values_to="cov") %>%
  dplyr::filter(cov >= VAR_COV_THR) %>% dplyr::pull(var)

det_cleanwoind <- det_clean %>% dplyr::select(country, dplyr::any_of(c("iso3","iso3c")), year, dplyr::any_of(keep_vars))

# Elimina años con >35% de faltantes (en columnas numéricas)
THR_YEAR_NA <- 0.35

num_cols <- det_cleanwoind %>%
  dplyr::select(dplyr::where(is.numeric)) %>%
  dplyr::select(-dplyr::any_of(c("year")), -dplyr::starts_with("block")) %>%
  names()

miss_year <- det_cleanwoind %>%
  dplyr::select(year, dplyr::all_of(num_cols)) %>%
  tidyr::pivot_longer(-year, values_to="v") %>%
  dplyr::summarise(pct_na = mean(is.na(v)), .by = year)

years_drop <- miss_year %>% dplyr::filter(pct_na > THR_YEAR_NA) %>% dplyr::pull(year)

message("Años a eliminar (>", THR_YEAR_NA*100, "% NA): ",
        ifelse(length(years_drop)==0, "ninguno", paste(sort(years_drop), collapse=", ")))

det_cleanwoyears <- det_cleanwoind %>% dplyr::filter(!year %in% years_drop)


# %NA después
# % NA por AÑO
det_cleanwoyears %>% dplyr::select(year, dplyr::where(is.numeric)) %>%
  tidyr::pivot_longer(-year, values_to="v") %>%
  dplyr::summarise(pct_na = 100*mean(is.na(v)), .by = year)%>%
  show()


# % NA por PAÍS (detecta iso3 o iso3c)
det_cleanwoyears %>% 
  dplyr::select(iso3, dplyr::where(is.numeric)) %>%
  tidyr::pivot_longer(-iso3, values_to="v") %>%
  dplyr::summarise(pct_na = 100*mean(is.na(v)), .by = iso3)%>%
  show()

# % NA por INDICADOR
det_cleanwoyears %>% dplyr::select(dplyr::where(is.numeric)) %>%
  dplyr::summarise(dplyr::across(everything(), ~100*mean(is.na(.)))) %>%
  tidyr::pivot_longer(everything(), names_to="indicador", values_to="pct_na") %>%
  dplyr::arrange(dplyr::desc(pct_na))%>%
  show()




# ===== HEATMAP DE COBERTURA (denominador = años del país) =====
library(dplyr); library(tidyr); library(ggplot2); library(scales); library(viridis)

base <- det_cleanwoyears
iso_cols <- intersect(c("iso3","iso3c"), names(base)); stopifnot(length(iso_cols)>=1)
base <- base %>% mutate(iso = dplyr::coalesce(!!!rlang::syms(iso_cols))) %>%
  relocate(country, iso, year, .before = 1)

# Solo indicadores numéricos, sin year/aux, y descarta 100% NA
measures <- base %>% dplyr::select(where(is.numeric)) %>%
  dplyr::select(-any_of(c("year","block_id"))) %>% names()
measures <- measures[sapply(measures, \(v) any(!is.na(base[[v]])))]

# Denominador estable = años por país; cobertura = años con dato / años del país
cov_tbl <- base %>%
  pivot_longer(all_of(measures), names_to = "variable", values_to = "val") %>%
  group_by(iso, country, variable) %>%
  summarise(
    n_years_non_na = n_distinct(year[!is.na(val)]),
    n_years_country = n_distinct(year),
    cobertura = ifelse(n_years_country > 0, n_years_non_na / n_years_country, NA_real_),
    .groups = "drop"
  )

# Orden para el plot
country_order <- cov_tbl %>% group_by(country) %>%
  summarise(cov_mean = mean(cobertura, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(cov_mean)) %>% pull(country)
var_order <- cov_tbl %>% group_by(variable) %>%
  summarise(cov_mean = mean(cobertura, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(cov_mean)) %>% pull(variable)

cov_tbl <- cov_tbl %>%
  mutate(country = factor(country, levels = country_order),
         variable = factor(variable, levels = var_order))

# Heatmap
p_cov <- ggplot(cov_tbl, aes(variable, country, fill = cobertura)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis(name = "Cobertura", limits = c(0,1), oob = squish,
                     labels = percent_format(1)) +
  labs(
    title = "Cobertura temporal por indicador y país",
    subtitle = paste0("Proporción de años con dato no-NA | Ventana: ",
                      min(base$year, na.rm=TRUE), "–", max(base$year, na.rm=TRUE)),
    x = "Indicador", y = "País"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cov)

# 1) Conjuntos completos de países e indicadores
all_countries <- det_cleanwoyears %>% dplyr::distinct(iso3, country)
all_vars      <- cov_tbl %>% dplyr::distinct(variable)

# 2) Completar la malla y rellenar faltantes con 0
cov_tbl2_full <- tidyr::expand_grid(all_countries, all_vars) %>%
  dplyr::left_join(cov_tbl, by = c("country","variable")) %>%
  dplyr::mutate(cobertura = dplyr::coalesce(cobertura, 0))

# (opcional) ordenar ejes
cov_tbl2_full <- cov_tbl2_full %>%
  dplyr::arrange(country, variable) %>%
  dplyr::mutate(
    country  = factor(country, levels = unique(country)),
    variable = factor(variable, levels = unique(variable))
  )

# 3) Graficar con la tabla completa
ggplot(cov_tbl2_full, aes(variable, country, fill = cobertura)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis(limits = c(0,1), oob = scales::squish,
                     labels = scales::percent_format(accuracy = 1),
                     name = "Cobertura") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cobertura temporal por indicador y país",
       subtitle = paste0("Proporción de años con dato no-NA | Ventana: ",
                         min(det_cleanwoyears$year, na.rm = TRUE), "–",
                         max(det_cleanwoyears$year, na.rm = TRUE)),
       x = "Indicador", y = "País")

# ---------- GUARDADOS FINALES ----------
dir.create(here::here("data"),  showWarnings = FALSE)
dir.create(here::here("output"), showWarnings = FALSE)

# Panel final para modelar (tras filtrar variables y años)
readr::write_csv(det_cleanwoyears, here::here("data","determinantes_wb.csv"))

# Metadatos útiles
saveRDS(list(
  keep_vars   = keep_vars,
  thr_by_var  = thr_by_var,
  class_tbl   = class_tbl,
  na_report   = na_report,
  cov_tbl     = cov_tbl,
  cov_tbl2_full = cov_tbl2_full
), file = here::here("data","determinantes_meta.rds"))

message("✅ Pipeline cerrado. Archivos en 'data/':",
        "\n  - determinantes_wb_raw.csv",
        "\n  - panel_con_gini.csv",
        "\n  - determinantes_wb.csv",
        "\n  - determinantes_panel_final.csv",
        "\n  - determinantes_meta.rds")

colnames(det_cleanwoyears)

######## Análisis exploratorio de las variables independientes #######
library(dplyr)
library(tidyr)
library(psych) # Para la estadística descriptiva (describe)
library(ggplot2) # Para histogramas

# Identificar las variables ID (no numéricas para el análisis)
id_vars <- c("country", "iso3", "year")

# Crear un dataframe solo con las variables numéricas para el análisis
det_numeric <- det_cleanwoyears %>%
  select(-all_of(id_vars))

# Obtener los nombres de las variables numéricas
numeric_cols <- colnames(det_numeric)

# Convertir el dataframe numérico a formato largo para usar ggplot eficientemente
data_long_hist <- det_numeric %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Generar y mostrar los histogramas en un solo gráfico (usando facetas)
histogramas <- ggplot(data_long_hist, aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1F77B4", color = "white") +
  facet_wrap(~ variable, scales = "free") + # 'scales = "free"' permite que cada gráfico tenga su propio eje X e Y
  labs(
    title = "Histogramas de Distribución para Variables Numéricas",
    x = "Valor",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 7)) # Ajustar el tamaño del texto de las facetas para que quepa

print(histogramas)

# Generar la estadística descriptiva
# se = error estándar de la media
estadistica_descriptiva <- psych::describe(det_numeric)

# Mostrar el resultado (se presenta por filas)
print(estadistica_descriptiva)

# Sugerencia de interpretación:
# - Un 'n' bajo indica alta tasa de NA para esa variable.
# - 'skew' > 1 o < -1 indica un sesgo significativo (la distribución no es normal).

######## Análisis exploratorio de las variables independientes #######
library(dplyr)
library(tidyr)
library(psych) # Para la estadística descriptiva (describe)
library(ggplot2) # Para histogramas

# Identificar las variables ID (no numéricas para el análisis)
id_vars <- c("country", "iso3", "year")

# Crear un dataframe solo con las variables numéricas para el análisis
det_numeric <- det_z %>%
  select(-all_of(id_vars))

# Obtener los nombres de las variables numéricas
numeric_cols <- colnames(det_numeric)

# Convertir el dataframe numérico a formato largo para usar ggplot eficientemente
data_long_hist <- det_numeric %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

# Generar y mostrar los histogramas en un solo gráfico (usando facetas)
histogramas <- ggplot(data_long_hist, aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1F77B4", color = "white") +
  facet_wrap(~ variable, scales = "free") + # 'scales = "free"' permite que cada gráfico tenga su propio eje X e Y
  labs(
    title = "Histogramas de Distribución para Variables Numéricas",
    x = "Valor",
    y = "Frecuencia"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 7)) # Ajustar el tamaño del texto de las facetas para que quepa

print(histogramas)

# Generar la estadística descriptiva
# se = error estándar de la media
estadistica_descriptiva <- psych::describe(det_numeric)

# Mostrar el resultado (se presenta por filas)
print(estadistica_descriptiva)

# Sugerencia de interpretación:
# - Un 'n' bajo indica alta tasa de NA para esa variable.
# - 'skew' > 1 o < -1 indica un sesgo significativo (la distribución no es normal).


# Calcular la matriz de correlación de Pearson
# 'use = "pairwise.complete.obs"' calcula la correlación para cada par de variables
# utilizando sólo las observaciones que no son NA para *ese par* específico.
cor_matrix <- cor(det_numeric, use = "pairwise.complete.obs", method = "pearson")

# Mostrar la matriz completa (será grande)
print(cor_matrix)

# -------------------------------------------------------------
# Alternativa: Visualización de la Correlación (Muy Recomendable)
# El paquete 'corrplot' es excelente para esto.

# Primero, instala y carga el paquete si no lo tienes:
# install.packages("corrplot")
library(corrplot)

# Visualizar la matriz de correlación
# Método 'circle' y orden 'hclust' para agrupar variables similares
corrplot(cor_matrix,
         method = "circle",
         type = "upper", # Muestra solo la parte superior (la matriz es simétrica)
         order = "hclust", # Agrupa visualmente las variables correlacionadas
         tl.cex = 0.5, # Tamaño del texto de las etiquetas
         p.mat = cor_matrix, # Usar la misma matriz para el valor p si se quisiera
         sig.level = 0.05,
         insig = "blank", # Dejar en blanco las correlaciones no significativas (si se usa p.mat)
         title = "Matriz de Correlación de Pearson entre Variables Numéricas",
         mar = c(0, 0, 1, 0)) # Ajuste de márgenes

# Sugerencia de interpretación:
# - Valores cercanos a +1: Correlación positiva fuerte (se mueven en la misma dirección).
# - Valores cercanos a -1: Correlación negativa fuerte (se mueven en direcciones opuestas).
# - Valores cercanos a 0: Poca o ninguna relación lineal.

library(dplyr)
library(tidyr)

# Definir el umbral de correlación
R_UMBRAL <- 0.75

# --- PASO 1: Preparar los datos de Correlación (ya calculada como cor_matrix) ---

# 1. Convertir la matriz de correlación a un formato largo (tidy)
cor_df <- cor_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable1") %>%
  pivot_longer(
    cols = -variable1,
    names_to = "variable2",
    values_to = "r"
  )

# 2. Filtrar y limpiar las correlaciones
cor_filtrada <- cor_df %>%
  # Eliminar correlaciones de una variable consigo misma (r=1)
  filter(variable1 != variable2) %>%
  # Eliminar duplicados (ej. A-B vs. B-A). Mantiene solo donde variable1 < variable2
  rowwise() %>%
  filter(variable1 < variable2) %>%
  ungroup() %>%
  # Filtrar por el umbral absoluto
  filter(abs(r) >= R_UMBRAL)

# --- PASO 2: Extraer los conteos 'n' de la estadística descriptiva ---

# --- PASO 2: Crear el dataframe de conteos 'n' (Esto es lo que faltaba) ---

# Crear un dataframe simple de variable y n a partir del resumen descriptivo
n_counts <- estadistica_descriptiva %>%
  # Seleccionar solo las columnas de la estadística descriptiva que contienen el nombre de la variable y 'n'
  select(n) %>%
  # Mover los nombres de las filas (nombres de las variables) a una columna
  tibble::rownames_to_column(var = "variable")

# --- PASO 3: Unir Correlaciones y Conteos 'n' ---

# 1. Unir el conteo 'n' para la variable1
tabla_colinealidad <- cor_filtrada %>%
  left_join(n_counts, by = c("variable1" = "variable")) %>%
  rename(n_var1 = n)

# 2. Unir el conteo 'n' para la variable2
tabla_colinealidad <- tabla_colinealidad %>%
  left_join(n_counts, by = c("variable2" = "variable")) %>%
  rename(n_var2 = n)

# 3. Reordenar y formatear la tabla final
tabla_colinealidad_final <- tabla_colinealidad %>%
  select(
    variable1, n_var1,
    variable2, n_var2,
    r
  ) %>%
  # Ordenar por el valor absoluto de la correlación de mayor a menor
  arrange(desc(abs(r)))

# Mostrar el resultado final
print(tabla_colinealidad_final)

# Lista de variables a ELIMINAR directamente
variables_a_eliminar <- c(
  "desempleo_joven",
  "gasto_salud_pc",
  "gob_salud_gdp",
  "violencia",
  "banda_ancha_100",
  "DALY_tb" # Se elimina porque las variables de infraestructura son más causales/predictoras
)

# Crear un nuevo dataframe sin estas columnas
det_prelim_clean <- det_cleanwoyears %>%
  select(-all_of(variables_a_eliminar))

# Definir las variables de infraestructura
vars_infra <- c("electrificacion", "agua_basica", "saneamiento_basico")

# Crear el índice de infraestructura y eliminar las tres originales
det_final_clean <- det_prelim_clean %>%
  # Calcular la media del índice (promedio de las 3 variables por fila/país-año)
  # rowMeans(na.rm = TRUE) asegura que se calcule el índice incluso si hay NA en alguna de las 3
  rowwise() %>%
  mutate(indice_infra = mean(c_across(all_of(vars_infra)), na.rm = TRUE)) %>%
  ungroup() %>%
  # Eliminar las variables originales
  select(-all_of(vars_infra))

# Resumen de la nueva variable para verificar su creación
print(det_final_clean %>% select(indice_infra) %>% summary())

# Revisa las columnas del dataframe final
print(colnames(det_final_clean))


# Crear un dataframe solo con las variables numéricas para el análisis
det_numeric <- det_final_clean %>%
  select(-all_of(id_vars))

# Obtener los nombres de las variables numéricas
numeric_cols <- colnames(det_numeric)
# Calcular la matriz de correlación de Pearson
# 'use = "pairwise.complete.obs"' calcula la correlación para cada par de variables
# utilizando sólo las observaciones que no son NA para *ese par* específico.
cor_matrix <- cor(det_numeric, use = "pairwise.complete.obs", method = "pearson")

# Mostrar la matriz completa (será grande)
print(cor_matrix)

# 1. Convertir la matriz de correlación a un formato largo (tidy)
cor_df <- cor_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "variable1") %>%
  pivot_longer(
    cols = -variable1,
    names_to = "variable2",
    values_to = "r"
  )

# 2. Filtrar y limpiar las correlaciones
cor_filtrada <- cor_df %>%
  # Eliminar correlaciones de una variable consigo misma (r=1)
  filter(variable1 != variable2) %>%
  # Eliminar duplicados (ej. A-B vs. B-A). Mantiene solo donde variable1 < variable2
  rowwise() %>%
  filter(variable1 < variable2) %>%
  ungroup() %>%
  # Filtrar por el umbral absoluto
  filter(abs(r) >= R_UMBRAL)
# --- PASO 2: Crear el dataframe de conteos 'n' (Esto es lo que faltaba) ---

# Crear un dataframe simple de variable y n a partir del resumen descriptivo
n_counts <- estadistica_descriptiva %>%
  # Seleccionar solo las columnas de la estadística descriptiva que contienen el nombre de la variable y 'n'
  select(n) %>%
  # Mover los nombres de las filas (nombres de las variables) a una columna
  tibble::rownames_to_column(var = "variable")

# --- PASO 3: Unir Correlaciones y Conteos 'n' ---

# 1. Unir el conteo 'n' para la variable1
tabla_colinealidad <- cor_filtrada %>%
  left_join(n_counts, by = c("variable1" = "variable")) %>%
  rename(n_var1 = n)

# 2. Unir el conteo 'n' para la variable2
tabla_colinealidad <- tabla_colinealidad %>%
  left_join(n_counts, by = c("variable2" = "variable")) %>%
  rename(n_var2 = n)

# 3. Reordenar y formatear la tabla final
tabla_colinealidad_final <- tabla_colinealidad %>%
  select(
    variable1, n_var1,
    variable2, n_var2,
    r
  ) %>%
  # Ordenar por el valor absoluto de la correlación de mayor a menor
  arrange(desc(abs(r)))

# Mostrar el resultado final
print(tabla_colinealidad_final)

# Generar la estadística descriptiva
# se = error estándar de la media
estadistica_descriptiva <- psych::describe(det_numeric)

# Mostrar el resultado (se presenta por filas)
print(estadistica_descriptiva)

# Calcular la matriz de correlación de Pearson
# 'use = "pairwise.complete.obs"' calcula la correlación para cada par de variables
# utilizando sólo las observaciones que no son NA para *ese par* específico.
cor_matrix <- cor(det_numeric, use = "pairwise.complete.obs", method = "pearson")

# Mostrar la matriz completa (será grande)
print(cor_matrix)

# -------------------------------------------------------------
# Alternativa: Visualización de la Correlación (Muy Recomendable)
# El paquete 'corrplot' es excelente para esto.

# Primero, instala y carga el paquete si no lo tienes:
# install.packages("corrplot")
library(corrplot)

# Visualizar la matriz de correlación
# Método 'circle' y orden 'hclust' para agrupar variables similares
corrplot(cor_matrix,
         method = "circle",
         type = "upper", # Muestra solo la parte superior (la matriz es simétrica)
         order = "hclust", # Agrupa visualmente las variables correlacionadas
         tl.cex = 0.5, # Tamaño del texto de las etiquetas
         p.mat = cor_matrix, # Usar la misma matriz para el valor p si se quisiera
         sig.level = 0.05,
         insig = "blank", # Dejar en blanco las correlaciones no significativas (si se usa p.mat)
         title = "Matriz de Correlación de Pearson entre Variables Numéricas",
         mar = c(0, 0, 1, 0)) # Ajuste de márgenes

####### Estandarizar valores Z#####
det_z <- det_final_clean %>%
  mutate(across(
    where(is.numeric) & !any_of(c("year")),
    ~ scale(.x)[,1]
  ))

# Guardar versión estandarizada
write_csv(det_z, here::here("data","determinantes_wb_z.csv"))


grViz("
digraph variables {

  graph [
    layout = dot,
    rankdir = TB,
    fontsize = 10,
    fontname = Helvetica,
    labelloc = t
  ]

  node [
    shape = box,
    style = filled,
    fillcolor = \"#F7F7F7\",
    color = \"#444444\",
    fontname = Helvetica,
    fontsize = 10
  ]

  edge [
    color = \"#444444\",
    arrowsize = 0.7,
    fontname = Helvetica,
    fontsize = 9
  ]

  # =========================
  # UNIVERSO INICIAL
  # =========================
  A [label = \"Universo inicial\\n40 indicadores candidatos\\nDeterminantes sociales + indicadores epidemiológicos\\n(WDI + GBD + OWID)\"]

  # =========================
  # ARMONIZACIÓN
  # =========================
  B [label = \"Armonización y limpieza\\n• Unificación país–año\\n• Resolución de duplicados\\n• Homologación de escalas\\n• Prioridad de fuentes (SEDLAC > WDI)\"]

  # =========================
  # COBERTURA
  # =========================
  C [label = \"Evaluación de cobertura\\n• % NA por variable\\n• % NA por país\\n• % NA por año\"]

  D1 [label = \"Exclusión temporal\\nAños con >35% de NA global\\n(eliminación de años inestables)\",
      fillcolor = \"#EFEFEF\"]

  D2 [label = \"Criterio de cobertura mínima por indicador\\n≥ 40% de observaciones válidas\\n(panel país–año)\"]

  # =========================
  # CLASIFICACIÓN EMPÍRICA
  # =========================
  F [label = \"Clasificación empírica de variables\\n• Lentas\\n• Mixtas\\n• Volátiles\\n(según monotonicidad y cambio relativo)\"]

  # =========================
  # DATOS FALTANTES
  # =========================
  E [label = \"Tratamiento de datos faltantes\\n• Winsorización p1–p99\\n• Imputación restringida por país\\n• Umbrales diferenciados según volatilidad\\n• Sin extrapolación estructural\"]

  # =========================
  # COLINEALIDAD
  # =========================
  G [label = \"Diagnóstico de colinealidad\\n|r| ≥ 0.75\\n+ comparación de cobertura efectiva (n)\"]

  # =========================
  # DECISIONES TEÓRICO-EMPÍRICAS
  # =========================
  H [label = \"Decisiones conceptuales\\n• Eliminación de redundantes\\n• Priorización causal\\n• Preservación interpretativa\\n• Construcción de índices compuestos\"]

  I [label = \"Variables eliminadas\\n• Redundancia empírica\\n• Baja cobertura longitudinal\\n• Colinealidad estructural sin valor teórico añadido\",
      fillcolor = \"#EFEFEF\"]

  J [label = \"Construcción de índice\\nÍndice de infraestructura\\n(electricidad + agua + saneamiento)\",
      fillcolor = \"#EFEFEF\"]

  # =========================
  # RESULTADO FINAL
  # =========================
  K [label = \"Conjunto analítico final\\n33 variables independientes\\nPanel país–año balanceado\",
      fillcolor = \"#E8F0FE\"]

  # =========================
  # FLUJOS
  # =========================
  A -> B
  B -> C
  C -> D1
  D1 -> D2
  D2 -> F
  F -> E
  E -> G
  G -> H
  H -> I
  H -> J
  I -> K
  J -> K
}
")

