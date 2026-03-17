# File: R/pca_causas_latam_fix.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(psych)
  library(ggrepel)
  library(patchwork)
  library(glue)
  library(janitor)
})

# ----------------------- 1) Datos -------------------------------------------
raw <- read_csv(
  "/Users/sofiadiaz/Library/Mobile Documents/com~apple~CloudDocs/tesis-anahuac/data/dalys_causas_latam.csv",
  show_col_types = FALSE
) |>
  clean_names()

trad <- c(
  "Depressive disorders" = "Trastornos depresivos",
  "Anxiety disorders" = "Trastornos de ansiedad",
  "Bipolar disorder" = "Trastorno bipolar",
  "Autism spectrum disorders" = "Trastornos del espectro autista",
  "Attention-deficit/hyperactivity disorder" = "TDAH",
  "Idiopathic developmental intellectual disability" = "Discapacidad intelectual",
  "Conduct disorder" = "Trastorno de la conducta",
  "Schizophrenia" = "Esquizofrenia",
  "Alcohol use disorders" = "Trastornos por consumo de alcohol",
  "Drug use disorders" = "Trastornos por consumo de drogas",
  "Eating disorders" = "Trastornos de la conducta alimentaria",
  "Self-harm" = "Autolesiones",
  "Other mental disorders" = "Otros trastornos mentales",
  "Conflict and terrorism" = "Conflicto y terrorismo",
  "Police conflict and executions" = "Violencia policial",
  "Interpersonal violence" = "Violencia interpersonal"
)

stopifnot(all(c("cause","year","dalys") %in% names(raw)))

raw <- raw |> mutate(cause = recode(cause, !!!trad))

# ----------------------- 2) Matriz cause x year -----------------------------
years_all <- sort(unique(raw$year))  # fija orden
mat_tbl <- raw |>
  group_by(cause, year) |>
  summarise(dalys = sum(dalys, na.rm = TRUE), .groups = "drop") |>
  complete(cause, year = years_all)  # por qué: asegura mismo set de años en todas las causas

# ⚠️ Decide cómo tratar ausencias: NA (recomendado para correlaciones pairwise) o 0.
USE_ZERO <- FALSE
mat_tbl <- if (USE_ZERO) mutate(mat_tbl, dalys = replace_na(dalys, 0)) else mat_tbl

mat <- mat_tbl |>
  arrange(cause, year) |>
  pivot_wider(names_from = year, values_from = dalys) |>
  column_to_rownames("cause") |>
  as.matrix()

# Elimina causas con varianza ~0 (no aportan al PCA; evitan NaN en escalado)
var_ok <- apply(mat, 1, function(v) sd(v, na.rm = TRUE)) > 0
mat <- mat[var_ok, , drop = FALSE]

# ----------------------- 3) Modo de PCA -------------------------------------
MODE <- "causas_as_variables"  # "anhos_as_variables" para tu enfoque original

if (MODE == "causas_as_variables") {
  # Variables = causas; Observaciones = años
  X <- scale(t(mat), center = TRUE, scale = TRUE)    # columnas=causas
  R <- cor(X, use = "pairwise.complete.obs")
  R_smooth <- psych::cor.smooth(R)
  KMO_res <- psych::KMO(R_smooth)
  bartlett_try <- try(psych::cortest.bartlett(R_smooth, n = nrow(X)), silent = TRUE)
  bartlett_res <- if (inherits(bartlett_try, "try-error")) list(statistic = NA, p.value = NA) else bartlett_try
  
  set.seed(123)
  pca <- psych::principal(X, nfactors = 2, rotate = "varimax", scores = FALSE)
  
  df_coord <- as.data.frame(unclass(pca$loadings[, 1:2])) |>
    rownames_to_column("causa") |>
    rename(PC1 = RC1, PC2 = RC2)
} else {
  # Variables = años; Observaciones = causas
  X <- scale(mat, center = TRUE, scale = TRUE)
  R <- cor(X, use = "pairwise.complete.obs")
  R_smooth <- psych::cor.smooth(R)
  KMO_res <- psych::KMO(R_smooth)
  bartlett_try <- try(psych::cortest.bartlett(R_smooth, n = nrow(X)), silent = TRUE)
  bartlett_res <- if (inherits(bartlett_try, "try-error")) list(statistic = NA, p.value = NA) else bartlett_try
  
  set.seed(123)
  pca <- psych::principal(X, nfactors = 2, rotate = "varimax", scores = TRUE)
  
  df_coord <- as.data.frame(pca$scores[, 1:2]) |>
    setNames(c("PC1","PC2")) |>
    rownames_to_column("causa")
}

# ----------------------- 4) Indicadores -------------------------------------
var_exp_total <- sum(pca$Vaccounted["Proportion Var", ]) * 100
bartlett_val <- ifelse(is.numeric(bartlett_res$statistic), round(bartlett_res$statistic, 1), NA)
bartlett_p   <- ifelse(is.numeric(bartlett_res$p.value), formatC(bartlett_res$p.value, format="e", digits=2), NA)

indicadores <- glue(
  "KMO: {round(KMO_res$MSA, 2)} (adecuado)\n",
  "Bartlett: χ² = {bartlett_val}, p = {bartlett_p}\n",
  "Varianza total explicada: {round(var_exp_total,1)} %\n",
  "PC1: {round(pca$Vaccounted['Proportion Var',1]*100,1)} % | ",
  "PC2: {round(pca$Vaccounted['Proportion Var',2]*100,1)} %\n",
  "Rotación: Varimax | Criterio: Autovalores > 1"
)

# ----------------------- 5) Clustering + pesos ------------------------------
set.seed(123)
k <- 3
km <- kmeans(df_coord[, c("PC1","PC2")], centers = k, nstart = 50)  # por qué: estabilidad
df_coord$cluster <- factor(km$cluster)

pesos <- df_coord |>
  group_by(cluster) |>
  mutate(
    peso_raw = abs(PC1) + abs(PC2),
    peso_pct = 100 * peso_raw / sum(peso_raw)
  ) |>
  ungroup()

# ----------------------- 6) Plots -------------------------------------------
theme_set(theme_minimal(base_family = "Times New Roman"))

# 6.1 PCA plot
p_pca <- ggplot(df_coord, aes(PC1, PC2, fill = cluster)) +
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.9) +
  geom_text_repel(aes(label = causa), size = 3, family = "Times New Roman",
                  color = "black", segment.color = "gray70", max.overlaps = 50) +
  stat_ellipse(data = df_coord |> group_by(cluster) |> filter(n() > 2),
               geom = "polygon", alpha = 0.12, color = "black") +
  scale_fill_manual(
    values = c("1" = "#e31793", "2" = "#ec855f", "3" = "#698eff"),
    labels = c("Afectivo–Conductual", "Disruptivo–Estructural", "Relacional–Moralizado")
  ) +
  labs(
    title    = "Estructura y composición interna de las constelaciones del sufrimiento mental",
    subtitle = indicadores,
    x = "PC1", y = "PC2"
  ) +
  theme(
    plot.subtitle = element_text(size = 8, hjust = 0, color = "gray20", lineheight = 1.08),
    axis.title   = element_text(size = 9),
    axis.text    = element_text(size = 8),
    plot.title   = element_text(face = "bold", size = 12),
    panel.grid   = element_line(color = "gray90", linewidth = 0.25)
  )

# 6.2 Pesos plot
p_pesos <- pesos |>
  mutate(
    causa = str_wrap(causa, 25),
    cluster = factor(
      cluster,
      levels = c("1", "2", "3"),
      labels = c("Afectivo–Conductual", "Disruptivo–Estructural", "Relacional–Moralizado")
    )
  ) |>
  ggplot(aes(x = peso_pct, y = reorder(causa, peso_pct), fill = cluster)) +
  geom_col(alpha = 0.88, width = 0.6, show.legend = FALSE) +
  geom_text(
    aes(label = paste0(round(peso_pct, 1), "%")),
    hjust = -0.1, size = 2.7, family = "Times New Roman"
  ) +
  facet_wrap(~cluster, scales = "free_y", nrow = 1) +
  labs(x = "Peso (%)", y = NULL) +
  scale_fill_manual(
    values = c("Afectivo–Conductual" = "#e31793",
               "Disruptivo–Estructural" = "#ec855f",
               "Relacional–Moralizado" = "#698eff")
  ) +
  theme_minimal(base_size = 10, base_family = "Times New Roman") +
  theme(
    strip.text       = element_text(face = "bold", size = 9),
    axis.text.y      = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.25)
  ) +
  expand_limits(x = 105)

# 6.3 Composición final
final_plot <- (p_pca + p_pesos + plot_layout(widths = c(1.05, 1.35))) +
  plot_annotation(
    title = "Estructura y composición interna de las constelaciones del sufrimiento mental",
    theme = theme(plot.title = element_text(family = "Times New Roman",
                                            face = "bold", size = 11, hjust = 0.5))
  )

ggsave("Figura_PCA_Causas_FINAL.png", final_plot, width = 11, height = 6, dpi = 400)

final_plot
