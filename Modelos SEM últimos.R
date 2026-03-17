# file: scripts/sem_modelos_asc_todos.R
# Objetivo: ajustar TODOS los modelos (CFA/SEM) y ordenarlos por complejidad ascendente.
# Nota: requiere un data.frame `sem_data` con las variables indicadas abajo.

# --- Paquetes -----------------------------------------------------------------
ensure_pkgs <- function(pkgs) {
  miss <- setdiff(pkgs, rownames(installed.packages()))
  if (length(miss)) install.packages(miss, Ncpus = max(1, parallel::detectCores() - 1))
  invisible(lapply(pkgs, require, character.only = TRUE))
}
ensure_pkgs(c(
  "lavaan","semPlot","ggplot2","dplyr","tibble","tidyr","here","glue","psych","purrr"
))

# --- 0) Paquetes ---------------------------------------------------------------
pkgs <- c("dplyr","tidyr","stringr","janitor","readr","purrr")
invisible(lapply(setdiff(pkgs, rownames(installed.packages())), install.packages))
invisible(lapply(pkgs, library, character.only = TRUE))

# --- 1) Lee tu base ------------------------------------------------------------
# Ajusta la ruta si usas otra
raw <- readr::read_csv("/Users/sofiadiaz/Library/Mobile Documents/com~apple~CloudDocs/tesis-anahuac/data/dalys_causas_latam.csv", show_col_types = FALSE) |>
  janitor::clean_names()
colnames(raw)

raw %>%
  group_by(cause) %>%
  summarise(total_dalys = sum(dalys, na.rm = TRUE)) %>%
  ggplot(aes(x = reorder(cause, total_dalys), y = total_dalys / 1e6)) +
  geom_col(fill = "#3B82F6") +
  coord_flip() +
  labs(
    title = "Carga total de DALYs por trastorno (1990–2021)",
    x = "Causa",
    y = "DALYs totales (millones)"
  ) +
  theme_minimal(base_size = 13)

raw %>%
  group_by(year) %>%
  summarise(total_dalys = sum(dalys, na.rm = TRUE)) %>%
  ggplot(aes(x = year, y = total_dalys / 1e6)) +
  geom_line(color = "#1D4ED8", size = 1) +
  geom_point(color = "#1D4ED8") +
  labs(
    title = "Evolución regional de la carga total de DALYs (1990–2021)",
    x = "Año",
    y = "DALYs totales (millones)"
  ) +
  theme_minimal(base_size = 13)

raw %>%
  group_by(country) %>%
  summarise(mean_dalys = mean(dalys, na.rm = TRUE)) %>%
  ggplot(aes(x = reorder(country, mean_dalys), y = mean_dalys / 1e3)) +
  geom_col(fill = "#2563EB") +
  coord_flip() +
  labs(
    title = "Carga promedio de DALYs por país (1990–2021)",
    x = "País",
    y = "DALYs promedio (miles)"
  ) +
  theme_minimal(base_size = 13)

library(tidyverse)

# Matriz de causas (años como filas, causas como columnas)
mat_causas <- raw %>%
  group_by(year, cause) %>%
  summarise(dalys = sum(dalys, na.rm = TRUE)) %>%
  pivot_wider(names_from = cause, values_from = dalys)

# Escalar y aplicar PCA
pca <- prcomp(mat_causas[,-1], scale. = TRUE)
summary(pca)

# Visualizar componentes principales
biplot(pca, main = "Estructura interna de la carga por causa (PCA)")


mat_causas <- raw %>%
  group_by(cause, year) %>%
  summarise(dalys = sum(dalys, na.rm = TRUE)) %>%
  pivot_wider(names_from = year, values_from = dalys)

# 🔹 Asegurar que los nombres de fila sean las causas
mat_causas <- mat_causas %>%
  column_to_rownames("cause") %>%
  as.matrix() %>%
  scale()


# --- 2. PCA
pca <- prcomp(mat_causas, scale. = TRUE)
summary(pca)  # para revisar la varianza explicada

# --- 3. Filtrar variables dentro del 80% de varianza acumulada
# (primeros componentes hasta PC2 en tu caso)
var_exp <- summary(pca)$importance[3, ]  # proporción acumulada
selected_components <- which(var_exp <= 0.80)
selected_components <- max(selected_components)  # última dentro del 80%

# --- 4. Agrupar causas por categoría DSM-5
grupos <- tibble(
  cause = colnames(mat_causas),
  grupo = case_when(
    cause %in% c("Depressive disorders", "Anxiety disorders", "Bipolar disorder") ~ "Ánimo y ansiedad",
    cause %in% c("Autism spectrum disorders", "Attention-deficit/hyperactivity disorder",
                 "Idiopathic developmental intellectual disability") ~ "Neurodesarrollo",
    cause %in% c("Conduct disorder") ~ "Conducta e impulsos",
    cause %in% c("Schizophrenia") ~ "Psicóticos",
    cause %in% c("Alcohol use disorders", "Drug use disorders") ~ "Sustancias",
    cause %in% c("Eating disorders") ~ "Conducta alimentaria",
    cause %in% c("Self-harm") ~ "Autolesivas",
    cause %in% c("Conflict and terrorism", "Interpersonal violence",
                 "Police conflict and executions") ~ "Contextuales y violencia",
    TRUE ~ "Otros"
  )
)

# --- 5. Visualización publicación-ready
fviz_pca_biplot(pca,
                axes = c(1, 2),
                repel = TRUE,
                col.var = grupos$grupo,   # color por grupo DSM-5
                col.ind = "grey80",       # puntos grises para años
                addEllipses = TRUE,       # agrupa visualmente por categoría
                label = "var",
                pointshape = 21,
                fill.ind = "white",
                geom.ind = "point") +
  labs(
    title = "Estructura interna de la carga por causa (PCA, varianza >80%)",
    subtitle = "Agrupado según categorías DSM-5",
    x = "Componente principal 1 (≈70%)",
    y = "Componente principal 2 (≈14%)",
    color = "Grupo DSM-5"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")


# Clustering jerárquico sobre los scores del PCA
clusters <- kmeans(pca$x[, 1:2], centers = 3)  # prueba con 3–4 clusters
fviz_cluster(list(data = pca$x[, 1:2], cluster = clusters$cluster),
             geom = "text",
             repel = TRUE,
             labelsize = 4) +
  labs(title = "Agrupación de causas según estructura PCA",
       x = "PC1 (Clínico ↔ Social)",
       y = "PC2 (Etiológico ↔ Moralizado)")

library(tidyverse)

pesos %>%
  mutate(
    causa = str_wrap(causa, 20),
    cluster = factor(cluster,
                     levels = c("1", "2", "3"),
                     labels = c("Afectivo–Conductual",
                                "Disruptivo–Estructural",
                                "Relacional–Moralizado"))
  ) %>%
  ggplot(aes(x = peso_pct, 
             y = reorder(causa, peso_pct), 
             fill = cluster)) +
  geom_col(alpha = 0.85, width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(round(peso_pct, 1), "%")),
            hjust = -0.1, size = 3, color = "gray15") +
  facet_wrap(~cluster, scales = "free_y", nrow = 1) +
  labs(
    title = "Peso relativo de cada causa dentro de su cluster",
    subtitle = "Distancia normalizada al centroide del grupo (PCA)",
    x = "Peso (%)",
    y = NULL
  ) +
  scale_fill_manual(values = c("#EF4444", "#2563EB", "#10B981")) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    axis.text.y = element_text(size = 9),
    panel.spacing = unit(0.8, "lines"),
    plot.margin = margin(5, 10, 5, 5)
  ) +
  expand_limits(x = 105)

library(ggplot2)
library(ggrepel)
library(patchwork)

# --- Tema base
theme_set(theme_minimal(base_family = "Times New Roman"))

# --- PCA limpio, con tres clusters en escala de grises
p_pca <- ggplot(df_pca, aes(PC1, PC2, fill = cluster)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, color = "black", show.legend = FALSE) +
  geom_point(size = 2, shape = 21, color = "black", alpha = 0.7) +
  geom_text_repel(aes(label = causa),
                  size = 2.6, family = "Times New Roman",
                  max.overlaps = 25, color = "black", segment.color = "gray70") +
  annotate("text", x = -9, y = 0.25, label = "Afectivo–Conductual",
           family = "Times New Roman", fontface = "bold", size = 3.2) +
  annotate("text", x = 3.8, y = 0.3, label = "Disruptivo–Estructural",
           family = "Times New Roman", fontface = "bold", size = 3.2) +
  annotate("text", x = -14, y = 0.6, label = "Relacional–Moralizado",
           family = "Times New Roman", fontface = "bold", size = 3.2) +
  scale_fill_manual(values = c("gray60", "gray80", "gray40")) +
  labs(
    x = "PC1 (Clínico ↔ Social)",
    y = "PC2 (Etiológico ↔ Moralizado)"
  ) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid = element_line(color = "gray90", size = 0.25),
    plot.margin = margin(5, 5, 5, 5)
  )

# --- Pesos (sin cluster 3)
p_pesos <- pesos %>%
  filter(cluster != 3) %>%
  mutate(
    causa = str_wrap(causa, 25),
    cluster = factor(cluster,
                     levels = c("1", "2"),
                     labels = c("Afectivo–Conductual",
                                "Disruptivo–Estructural"))
  ) %>%
  ggplot(aes(x = peso_pct, y = reorder(causa, peso_pct))) +
  geom_col(fill = "gray40", alpha = 0.8, width = 0.55) +
  geom_text(aes(label = ifelse(peso_pct == 0 | is.na(peso_pct), "",
                               paste0(round(peso_pct, 1), "%"))),
            hjust = -0.15, size = 2.5, family = "Times New Roman", color = "black") +
  facet_wrap(~cluster, scales = "free_y", nrow = 2) +
  labs(
    x = "Peso (%)",
    y = NULL
  ) +
  theme_minimal(base_size = 9, base_family = "Times New Roman") +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    axis.text.y = element_text(size = 7.5),
    axis.text.x = element_text(size = 7.5),
    axis.title = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90", size = 0.25),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  expand_limits(x = 105)

# --- Combina directamente con patchwork y un título global centrado
final_plot <- (p_pca + p_pesos + plot_layout(widths = c(1.1, 1.2))) +
  plot_annotation(
    title = "Estructura y composición interna de las constelaciones del sufrimiento mental",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", face = "bold",
                                size = 11, hjust = 0.5)
    )
  )

# --- Guarda figura final
ggsave("Figura_PCA_Clusters_Pesos_FINAL.png", final_plot, width = 10.5, height = 5.8, dpi = 400)
final_plot


library(ggplot2)
library(ggrepel)
library(patchwork)

theme_set(theme_minimal(base_family = "Times New Roman"))

# --- PCA colorido
p_pca <- ggplot(df_pca, aes(PC1, PC2, fill = cluster)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, color = "black", show.legend = FALSE,
               data = df_pca %>% group_by(cluster) %>% filter(n() > 2)) +
  geom_point(size = 2.5, shape = 21, color = "black", alpha = 0.8) +
  geom_text_repel(
    data = df_pca %>%
      filter(!(cluster == 2 & !(causa %in% c("Conduct disorder",
                                             "Other mental disorders",
                                             "Conflict and terrorism")))),
    aes(label = causa),
    size = 2.6, family = "Times New Roman",
    max.overlaps = 25, color = "black", segment.color = "gray70"
  ) +
  annotate("text", x = -9, y = 0.25, label = "Afectivo–Conductual",
           family = "Times New Roman", fontface = "bold", size = 3.2, color = "#e31793") +
  annotate("text", x = 3.8, y = 0.3, label = "Disruptivo–Estructural",
           family = "Times New Roman", fontface = "bold", size = 3.2, color = "#ec855f") +
  annotate("text", x = -14, y = 0.6, label = "Relacional–Moralizado",
           family = "Times New Roman", fontface = "bold", size = 3.2, color = "#698eff") +
  scale_fill_manual(values = c("#e31793", "#ec855f", "#698eff")) +
  labs(
    x = "PC1 (Clínico ↔ Social)",
    y = "PC2 (Etiológico ↔ Moralizado)"
  ) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid = element_line(color = "gray90", size = 0.25),
    plot.margin = margin(5, 5, 5, 5)
  )

# --- Pesos con colores por cluster (idénticos al PCA) y sin etiquetas
p_pesos <- pesos %>%
  filter(cluster != 3) %>%
  mutate(
    causa = str_wrap(causa, 25),
    cluster = factor(cluster,
                     levels = c("1", "2"),
                     labels = c("Afectivo–Conductual",
                                "Disruptivo–Estructural"))
  ) %>%
  ggplot(aes(x = peso_pct, y = reorder(causa, peso_pct), fill = cluster)) +
  geom_col(alpha = 0.9, width = 0.55, show.legend = FALSE) +
  facet_wrap(~cluster, scales = "free_y", nrow = 2) +
  labs(
    x = "Peso (%)",
    y = NULL
  ) +
  # Colores iguales a los del PCA 💗🧡
  scale_fill_manual(values = c(
    "Afectivo–Conductual" = "#e31793",
    "Disruptivo–Estructural" = "#ec855f"
  )) +
  theme_minimal(base_size = 9, base_family = "Times New Roman") +
  theme(
    strip.text = element_text(face = "bold", size = 8),
    axis.text.y = element_text(size = 7.5, family = "Times New Roman", color = "black"),
    axis.text.x = element_text(size = 7.5),
    axis.title = element_text(size = 8),
    panel.grid.major = element_line(color = "gray90", size = 0.25),
    panel.spacing = unit(0.6, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  expand_limits(x = 105)



# --- Combinación final
final_plot <- (p_pca + p_pesos + plot_layout(widths = c(1.1, 1.2))) +
  plot_annotation(
    title = "Estructura y composición interna de las constelaciones del sufrimiento mental",
    theme = theme(
      plot.title = element_text(family = "Times New Roman", face = "bold",
                                size = 11, hjust = 0.5)
    )
  )

ggsave("Figura_PCA_Clusters_Pesos_COLOR_FILTRADO.png", final_plot,
       width = 10.5, height = 5.8, dpi = 400)

final_plot


raw %>%
  group_by(country, year) %>%
  summarise(total_dalys = sum(dalys, na.rm = TRUE)) %>%
  ggplot(aes(x = year, y = reorder(country, -total_dalys), fill = total_dalys)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "Evolución espacio–temporal de la carga de enfermedad mental",
    x = "Año", y = "País", fill = "DALYs"
  ) +
  theme_minimal(base_size = 12)


vars_needed <- c(
  "depresivos","ansiedad","bipolar","esquizofrenia","tea","tdah","conducta","otros_trs",
  "intelectual","autolesion","alcohol","drogas","tca","conflicto_terrorismo",
  "violencia_interpersonal","policia_ejecuciones"
)

# 1) Tabla de mapeo cause -> target (exactamente tus 16 etiquetas)
map_tbl <- tibble::tibble(
  cause = c(
    "Alcohol use disorders",
    "Anxiety disorders",
    "Attention-deficit/hyperactivity disorder",
    "Autism spectrum disorders",
    "Bipolar disorder",
    "Conduct disorder",
    "Conflict and terrorism",
    "Depressive disorders",
    "Drug use disorders",
    "Eating disorders",
    "Idiopathic developmental intellectual disability",
    "Interpersonal violence",
    "Other mental disorders",
    "Police conflict and executions",
    "Schizophrenia",
    "Self-harm"
  ),
  target = c(
    "alcohol",
    "ansiedad",
    "tdah",
    "tea",
    "bipolar",
    "conducta",
    "conflicto_terrorismo",
    "depresivos",
    "drogas",
    "tca",
    "intelectual",
    "violencia_interpersonal",
    "otros_trs",
    "policia_ejecuciones",
    "esquizofrenia",
    "autolesion"
  )
)

# 2) Construir sem_data via join (sin usar select/recode que se pisan)
sem_data <- raw %>%
  dplyr::inner_join(map_tbl, by = "cause") %>%        # agrega columna 'target'
  dplyr::group_by(country, iso3, year, target) %>%
  dplyr::summarise(val = sum(dalys, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = target, values_from = val) %>%
  dplyr::arrange(country, year)

# 3) Asegurar columnas esperadas por tu script
vars_needed <- c(
  "depresivos","ansiedad","bipolar","esquizofrenia","tea","tdah","conducta","otros_trs",
  "intelectual","autolesion","alcohol","drogas","tca","conflicto_terrorismo",
  "violencia_interpersonal","policia_ejecuciones"
)
faltan <- setdiff(vars_needed, names(sem_data))
if (length(faltan)) sem_data[faltan] <- NA_real_

sem_data <- sem_data %>%
  dplyr::select(country, iso3, year, dplyr::all_of(vars_needed)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(vars_needed), as.numeric))

stopifnot(all(vars_needed %in% names(sem_data)))

# --- Heatmaps 2019 por causa (color relativo dentro de cada causa) -----------
library(dplyr); library(tidyr); library(ggplot2); library(viridisLite)

causas_all <- c(
  "depresivos","ansiedad","bipolar","esquizofrenia",
  "tea","tdah","conducta","otros_trs","intelectual",
  "autolesion","alcohol","drogas","tca",
  "conflicto_terrorismo","violencia_interpersonal","policia_ejecuciones"
)

# Etiqueta compacta para valores crudos
compact_si <- function(x) {
  out <- ifelse(is.na(x), "",
                ifelse(abs(x) >= 1e9, sprintf("%.1fB", x/1e9),
                       ifelse(abs(x) >= 1e6, sprintf("%.1fM", x/1e6),
                              ifelse(abs(x) >= 1e3, sprintf("%.1fk", x/1e3),
                                     sprintf("%.0f", x)))))
  sub("\\.0([kMB])$", "\\1", out)
}

# Base 2019
df19 <- sem_data %>%
  dplyr::filter(year == 2019) %>%
  dplyr::select(country, iso3, year, dplyr::all_of(causas_all))

# Orden de países: fijo según patrón z para que ambas figuras tengan la misma fila
mat_z <- scale(as.matrix(df19 %>% dplyr::select(dplyr::all_of(causas_all))))
rownames(mat_z) <- df19$country
countries_ord <- rownames(mat_z)[stats::hclust(stats::dist(mat_z))$order]

# Contraste de texto según luminancia de viridis (usa el fill que dibujamos)
text_contrast <- function(fill_scaled) {
  s <- pmax(0, pmin(1, fill_scaled))
  hex <- viridisLite::viridis(256)[pmax(1, pmin(256, round(s*255)+1))]
  R <- strtoi(substr(hex, 2, 3), 16) / 255
  G <- strtoi(substr(hex, 4, 5), 16) / 255
  B <- strtoi(substr(hex, 6, 7), 16) / 255
  L <- 0.2126*R + 0.7152*G + 0.0722*B
  ifelse(L < 0.50, "white", "black")
}

# Botón: "raw" (color min–max por causa, etiqueta = crudo) o "z" (z por causa)
make_heatmap_2019 <- function(modo = c("raw","z")) {
  modo <- match.arg(modo)
  
  if (modo == "raw") {
    # Largo con valor crudo y reescala 0–1 *por causa* para el color
    df_plot <- df19 %>%
      tidyr::pivot_longer(cols = dplyr::all_of(causas_all),
                          names_to = "causa", values_to = "val") %>%
      group_by(causa) %>%
      mutate(fill_scaled = (val - min(val, na.rm = TRUE)) /
               (max(val, na.rm = TRUE) - min(val, na.rm = TRUE))) %>%
      ungroup() %>%
      mutate(
        country = factor(country, levels = countries_ord),
        causa   = factor(causa, levels = causas_all),
        lbl     = compact_si(val),
        txt_col = text_contrast(fill_scaled)
      )
    
    gg <- ggplot(df_plot, aes(causa, country, fill = fill_scaled)) +
      geom_tile(color = "white", linewidth = 0.15) +
      geom_text(aes(label = lbl), size = 3.1, color = "black", alpha = 0.18) +
      geom_text(aes(label = lbl, color = txt_col), size = 3) +
      scale_color_identity() +
      scale_fill_viridis_c(name = "Relativo dentro de causa",
                           limits = c(0,1), breaks = c(0,0.5,1),
                           labels = c("Mín","Medio","Máx"), na.value = "grey90") +
      labs(title = "Causas (16) en 2019 — color relativo por causa (etiqueta = DALYs)",
           x = "Causa", y = "País")
    
  } else { # modo == "z"
    # z por causa (columna) para el color; etiqueta = z
    df_z <- df19 %>%
      mutate(across(all_of(causas_all), ~ as.numeric(scale(.x))))
    
    df_plot <- df_z %>%
      tidyr::pivot_longer(cols = dplyr::all_of(causas_all),
                          names_to = "causa", values_to = "z") %>%
      mutate(
        country = factor(country, levels = countries_ord),
        causa   = factor(causa, levels = causas_all),
        lbl     = ifelse(is.na(z), "", sprintf("%.1f", z)),
        # decidimos contraste con el *mismo* fill (z reescalado a 0–1)
        z01     = (z - min(z, na.rm = TRUE)) / (max(z, na.rm = TRUE) - min(z, na.rm = TRUE)),
        txt_col = text_contrast(z01)
      )
    
    gg <- ggplot(df_plot, aes(causa, country, fill = z)) +
      geom_tile(color = "white", linewidth = 0.15) +
      geom_text(aes(label = lbl), size = 3.1, color = "black", alpha = 0.18) +
      geom_text(aes(label = lbl, color = txt_col), size = 3) +
      scale_color_identity() +
      scale_fill_gradient2(name = "z-score (por causa)", limits = c(-2.5, 2.5),
                           low = "#2166ac", mid = "white", high = "#b2182b",
                           na.value = "grey90") +
      labs(title = "Causas (16) en 2019 — z por causa (etiqueta = z)",
           x = "Causa", y = "País")
  }
  
  gg + theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
}

# Render y guardar (mismo orden de países)
p_raw <- make_heatmap_2019("raw")
p_z   <- make_heatmap_2019("z")

if (!exists("outdir")) outdir <- here::here("output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(outdir, "heatmap2019_raw_by_cause.png"), p_raw, width = 12, height = 9, dpi = 300)
ggsave(file.path(outdir, "heatmap2019_z_by_cause.png"),   p_z,   width = 12, height = 9, dpi = 300)

p_raw; p_z

make_heatmap_2019 <- function(modo = c("raw","z")) {
  modo <- match.arg(modo)
  
  # ---- base larga (solo países x 16 causas) ----
  base_long <- df19 |>
    tidyr::pivot_longer(cols = dplyr::all_of(causas_all),
                        names_to = "causa", values_to = "val") |>
    mutate(
      country = as.character(country),
      causa   = as.character(causa)
    )
  
  # Totales por país (suma por fila) -> columna extra
  tot_fila <- base_long |>
    dplyr::group_by(country) |>
    dplyr::summarise(val = sum(val, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(causa = "Total_fila")
  
  # Totales por causa (suma por columna) -> fila extra
  tot_col <- base_long |>
    dplyr::group_by(causa) |>
    dplyr::summarise(val = sum(val, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(country = "Total_columna")
  
  # Tabla de min–max por causa (para NO alterar escalas de las 16 causas)
  rng_by_cause <- base_long |>
    dplyr::group_by(causa) |>
    dplyr::summarise(cmin = min(val, na.rm = TRUE),
                     cmax = max(val, na.rm = TRUE),
                     .groups = "drop")
  
  if (modo == "raw") {
    # Escala 0–1 SOLO con la base; luego aplicamos a los totales
    base_scaled <- base_long |>
      dplyr::left_join(rng_by_cause, by = "causa") |>
      dplyr::mutate(fill_scaled = (val - cmin) / (cmax - cmin))
    
    # Totales por columna (fila “Total (causa)”) se escalan con el rng de cada causa
    tot_col_scaled <- tot_col |>
      dplyr::left_join(rng_by_cause, by = "causa") |>
      dplyr::mutate(fill_scaled = (val - cmin) / (cmax - cmin))
    
    # Totales por fila (columna “Total (país)”) requieren su propio rango
    rmin <- min(tot_fila$val, na.rm = TRUE); rmax <- max(tot_fila$val, na.rm = TRUE)
    tot_fila_scaled <- tot_fila |>
      dplyr::mutate(fill_scaled = (val - rmin) / (rmax - rmin))
    
    df_plot <- dplyr::bind_rows(base_scaled, tot_col_scaled, tot_fila_scaled) |>
      dplyr::mutate(
        country = factor(country, levels = c(countries_ord, "Total_columna")),
        causa   = factor(causa,   levels = c(causas_all, "Total_fila")),
        lbl     = compact_si(val),
        txt_col = text_contrast(pmax(0, pmin(1, fill_scaled))),
        is_total = (country == "Total_columna" | causa == "Total_fila")
      )
    
    gg <- ggplot(df_plot, aes(causa, country, fill = fill_scaled)) +
      geom_tile(color = "white", linewidth = 0.15) +
      geom_text(aes(label = lbl, color = txt_col,
                    fontface = ifelse(is_total, "bold", "plain")),
                size = 3) +
      scale_color_identity() +
      scale_fill_viridis_c(name = "Relativo dentro de causa",
                           limits = c(0,1), breaks = c(0,0.5,1),
                           labels = c("Mín","Medio","Máx"), na.value = "grey90") +
      labs(title = "Causas (16) en 2019 — con totales (etiqueta = DALYs)",
           x = "Causa", y = "País")
    
  } else { # modo == "z"
    # z por causa con la base; luego calculamos z de los totales usando
    # media y sd de la base (no alteramos z de las celdas originales)
    z_stats <- base_long |>
      dplyr::group_by(causa) |>
      dplyr::summarise(m = mean(val, na.rm = TRUE),
                       s = sd(val, na.rm = TRUE), .groups = "drop")
    
    base_z <- base_long |>
      dplyr::left_join(z_stats, by = "causa") |>
      dplyr::mutate(z = (val - m) / s)
    
    tot_col_z <- tot_col |>
      dplyr::left_join(z_stats, by = "causa") |>
      dplyr::mutate(z = (val - m) / s)
    
    # Para “Total_fila” (col extra) z interno a esa columna de totales
    r_m <- mean(tot_fila$val, na.rm = TRUE); r_s <- sd(tot_fila$val, na.rm = TRUE)
    tot_fila_z <- tot_fila |>
      dplyr::mutate(z = (val - r_m) / r_s)
    
    df_plot <- dplyr::bind_rows(base_z, tot_col_z, tot_fila_z) |>
      dplyr::mutate(
        country = factor(country, levels = c(countries_ord, "Total_columna")),
        causa   = factor(causa,   levels = c(causas_all, "Total_fila")),
        lbl     = ifelse(is.na(z), "", sprintf("%.1f", z)),
        z01     = (z - min(z, na.rm = TRUE)) / (max(z, na.rm = TRUE) - min(z, na.rm = TRUE)),
        txt_col = text_contrast(z01),
        is_total = (country == "Total_columna" | causa == "Total_fila")
      )
    
    gg <- ggplot(df_plot, aes(causa, country, fill = z)) +
      geom_tile(color = "white", linewidth = 0.15) +
      geom_text(aes(label = lbl, color = txt_col,
                    fontface = ifelse(is_total, "bold", "plain")),
                size = 3) +
      scale_color_identity() +
      scale_fill_gradient2(name = "z-score (por causa)", limits = c(-2.5, 2.5),
                           low = "#2166ac", mid = "white", high = "#b2182b",
                           na.value = "grey90") +
      labs(title = "Causas (16) en 2019 — con totales (etiqueta = z)",
           x = "Causa", y = "País")
  }
  
  gg + theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank()
    )
}


p_raw_tot <- make_heatmap_2019("raw")
p_raw_tot   # se abre en la pestaña Plots

p_z_tot <- make_heatmap_2019("z")
p_z_tot


# Asegura columnas que falten (por si algún país/año no tuvo cierta causa)
faltan <- setdiff(vars_needed, names(sem_data))
if (length(faltan)) sem_data[faltan] <- NA_real_

# Orden y tipos
# Reordenar y tipar, forzando dplyr/tidyr
sem_data <- sem_data %>%
  dplyr::select(country, iso3, year, dplyr::all_of(vars_needed)) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(vars_needed), as.numeric)) %>%
  dplyr::arrange(country, year)


# Check que tu script pedía
stopifnot(all(vars_needed %in% names(sem_data)))


# --- Datos requeridos ---------------------------------------------------------
vars_needed <- c(
  "depresivos","ansiedad","bipolar","esquizofrenia","tea","tdah","conducta","otros_trs",
  "intelectual","autolesion","alcohol","drogas","tca","conflicto_terrorismo",
  "violencia_interpersonal","policia_ejecuciones"
)
stopifnot(exists("sem_data"))
stopifnot(all(vars_needed %in% names(sem_data)))

# --- Salidas ------------------------------------------------------------------
outdir <- here::here("output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# --- Helpers (solo lo esencial) ----------------------------------------------
count_pat <- function(txt, pat) {
  m <- gregexpr(pat, txt, perl = TRUE)[[1]]
  sum(m > 0)
}
model_complexity <- function(s) {
  s_clean <- gsub("#.*$", "", s)
  load_defs <- count_pat(s_clean, "=~")
  n_plus    <- count_pat(s_clean, "\\+")
  paths     <- count_pat(s_clean, "(?<!~)~(?!~)")
  covs      <- count_pat(s_clean, "~~")
  score <- 3*paths + 2*covs + 1*(load_defs + n_plus) # prioriza estructura
  c(load_defs = load_defs, indicators_extra = n_plus, paths = paths, covs = covs, score = score)
}

safe_fit <- function(name, syntax, stdlv = TRUE, estimator = "MLR", data = sem_data) {
  tryCatch(sem(syntax, data = data, estimator = estimator, missing = "fiml", std.lv = stdlv),
           error = function(e) NULL)
}

has_converged <- function(f) !is.null(f) && isTRUE(tryCatch(lavInspect(f, "converged"), error = function(e) FALSE))

get_measures_safe <- function(f) {
  if (!has_converged(f)) return(tibble(CFI=NA_real_,TLI=NA_real_,RMSEA=NA_real_,SRMR=NA_real_,AIC=NA_real_,BIC=NA_real_))
  vals <- suppressWarnings(tryCatch(fitmeasures(f, c("cfi","tli","rmsea","srmr","aic","bic")),
                                    error=function(e) setNames(rep(NA_real_,6), c("cfi","tli","rmsea","srmr","aic","bic"))))
  tibble(CFI=as.numeric(vals["cfi"]),TLI=as.numeric(vals["tli"]),RMSEA=as.numeric(vals["rmsea"]),
         SRMR=as.numeric(vals["srmr"]),AIC=as.numeric(vals["aic"]),BIC=as.numeric(vals["bic"]))
}

# --- Modelos: medición CFA (1F→6F) ------------------------------------------
m1 <- '
  F =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual + autolesion + alcohol + drogas + tca + conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
'

m2 <- '
  TSM  =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual
  CONS =~ autolesion + alcohol + drogas + tca + conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  TSM ~~ CONS
'

m3 <- '
  TSM      =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual
  BEHAVIOR =~ autolesion + alcohol + drogas + tca
  VIOLENCE =~ conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  TSM ~~ BEHAVIOR
  TSM ~~ VIOLENCE
  BEHAVIOR ~~ VIOLENCE
'

m4 <- '
  TSM_CONGENITO  =~ tea + tdah + conducta + intelectual
  TSM_ADQUIRIDO  =~ depresivos + ansiedad + bipolar + esquizofrenia + otros_trs
  BEHAVIOR       =~ autolesion + alcohol + drogas + tca
  VIOLENCE       =~ conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  TSM_CONGENITO ~~ TSM_ADQUIRIDO + BEHAVIOR + VIOLENCE
  TSM_ADQUIRIDO ~~ BEHAVIOR + VIOLENCE
  BEHAVIOR ~~ VIOLENCE
'

m5 <- '
  F_Neurodesarrollo  =~ tea + tdah + intelectual + conducta
  F_Emocional        =~ depresivos + ansiedad + bipolar + esquizofrenia + otros_trs
  F_Comportamiento   =~ autolesion + tca
  F_Adicciones       =~ alcohol + drogas
  F_Violencia        =~ conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  F_Neurodesarrollo ~~ F_Emocional + F_Comportamiento + F_Adicciones + F_Violencia
  F_Emocional       ~~ F_Comportamiento + F_Adicciones + F_Violencia
  F_Comportamiento  ~~ F_Adicciones + F_Violencia
  F_Adicciones      ~~ F_Violencia
'

m6 <- '
  F_Neurodesarrollo  =~ tea + tdah + intelectual
  F_Emocional        =~ depresivos + ansiedad + bipolar
  F_Psicosis         =~ esquizofrenia + otros_trs
  F_Comportamiento   =~ autolesion + tca + conducta
  F_Adicciones       =~ alcohol + drogas
  F_Violencia        =~ conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  F_Neurodesarrollo ~~ F_Emocional + F_Psicosis + F_Comportamiento + F_Adicciones + F_Violencia
  F_Emocional       ~~ F_Psicosis + F_Comportamiento + F_Adicciones + F_Violencia
  F_Psicosis        ~~ F_Comportamiento + F_Adicciones + F_Violencia
  F_Comportamiento  ~~ F_Adicciones + F_Violencia
  F_Adicciones      ~~ F_Violencia
'

# --- Nuevos modelos sugeridos -----------------------------------------------
# 2F alternativo: Internalizante vs Externalizante
m2_int_ext <- '
  INT =~ depresivos + ansiedad + bipolar + tca + autolesion
  EXT =~ conducta + tdah + alcohol + drogas
  INT ~~ EXT
'

# 3F alternativo: INT/EXT/Thought
m3_int_ext_thought <- '
  INT     =~ depresivos + ansiedad + tca + autolesion + bipolar
  EXT     =~ conducta + tdah + alcohol + drogas
  THOUGHT =~ esquizofrenia + otros_trs + intelectual
  INT ~~ EXT
  INT ~~ THOUGHT
  EXT ~~ THOUGHT
'

# Orden superior 3F → GLOBAL
m3_ho <- '
  TSM      =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual
  BEHAVIOR =~ autolesion + alcohol + drogas + tca
  VIOLENCE =~ conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  F_GLOBAL =~ TSM + BEHAVIOR + VIOLENCE
'

# Bifactor en TSM: general + grupos ortogonales
m_bifact_tsm <- '
  G       =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual
  S_INT   =~ depresivos + ansiedad + bipolar
  S_NEURO =~ tea + tdah + intelectual + conducta
  S_PSYCH =~ esquizofrenia + otros_trs
  G ~~ 0*S_INT
  G ~~ 0*S_NEURO
  G ~~ 0*S_PSYCH
  S_INT ~~ 0*S_NEURO
  S_INT ~~ 0*S_PSYCH
  S_NEURO ~~ 0*S_PSYCH
'

# --- Orden superior, variantes CFA/SEM, fisiopatológicos --------------------
m6_ho <- '
  F_Neurodesarrollo  =~ tea + tdah + intelectual
  F_Emocional        =~ depresivos + ansiedad + bipolar
  F_Psicosis         =~ esquizofrenia + otros_trs
  F_Comportamiento   =~ autolesion + tca + conducta
  F_Adicciones       =~ alcohol + drogas
  F_Violencia        =~ conflicto_terrorismo + violencia_interpersonal + policia_ejecuciones
  F_GLOBAL_DALY =~ F_Neurodesarrollo + F_Emocional + F_Psicosis + F_Comportamiento + F_Adicciones + F_Violencia
'

m_final <- '
  F_Neurodesarrollo  =~ tdah + intelectual
  F_Emocional        =~ depresivos + ansiedad + bipolar + esquizofrenia
  F_Comportamiento   =~ autolesion + tca + conducta
  F_Adicciones       =~ alcohol + drogas
  F_Violencia        =~ conflicto_terrorismo + violencia_interpersonal
  F_Neurodesarrollo ~~ F_Emocional + F_Comportamiento + F_Adicciones + F_Violencia
  F_Emocional       ~~ F_Comportamiento + F_Adicciones + F_Violencia
  F_Comportamiento  ~~ F_Adicciones + F_Violencia
  F_Adicciones      ~~ F_Violencia
'

m_pca <- '
  F1_Mixed =~ conducta + esquizofrenia + bipolar
  F2_Violence_Addiction =~ alcohol + violencia_interpersonal + conflicto_terrorismo
  F3_Mood_Selfharm =~ depresivos + autolesion
  F4_Neurodevelopmental =~ tdah + intelectual
  F5_Anxiety_Other =~ ansiedad + tca
  F1_Mixed ~~ F2_Violence_Addiction + F3_Mood_Selfharm + F4_Neurodevelopmental + F5_Anxiety_Other
  F2_Violence_Addiction ~~ F3_Mood_Selfharm + F4_Neurodevelopmental + F5_Anxiety_Other
  F3_Mood_Selfharm ~~ F4_Neurodevelopmental + F5_Anxiety_Other
  F4_Neurodevelopmental ~~ F5_Anxiety_Other
'

model_sem_2F <- '
  TSM  =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual
  CONS =~ autolesion + alcohol + drogas + tca
  CONS ~ TSM
'

model_sem_3F <- '
  TSM    =~ depresivos + ansiedad + bipolar + esquizofrenia + tea + tdah + conducta + otros_trs + intelectual
  AUTTCA =~ autolesion + tca
  SUD    =~ drogas
  AUTTCA ~ TSM
  SUD    ~ TSM
  AUTTCA ~~ SUD
'

model_physio <- '
  TSM_N =~ tea + tdah + intelectual + esquizofrenia + conducta
  TSM_A =~ depresivos + ansiedad + bipolar + otros_trs
  AUTTCA =~ autolesion + tca
  AUTTCA ~ g1*TSM_A + g2*TSM_N
  drogas  ~ b1*TSM_A + b2*TSM_N
  TSM_A ~~ TSM_N
  AUTTCA ~~ drogas
'

model_physio_clean <- '
  TSM_N =~ tea + intelectual + esquizofrenia
  TSM_A =~ depresivos + ansiedad + bipolar
  AUTTCA =~ autolesion + tca
  AUTTCA ~ TSM_A + TSM_N
  drogas  ~ TSM_A + TSM_N
  TSM_A ~~ TSM_N
  AUTTCA ~~ drogas
'

model_final <- '
  TSM =~ depresivos + ansiedad + bipolar + esquizofrenia + intelectual
  autolesion ~ TSM
  tca        ~ TSM
  drogas     ~ TSM
  autolesion ~~ tca
  autolesion ~~ drogas
  tca ~~ drogas
'

model_twof <- '
  TSM_A =~ depresivos + ansiedad + bipolar
  TSM_N =~ esquizofrenia + intelectual
  AUTTCA =~ autolesion + tca
  AUTTCA ~ TSM_A + TSM_N
  drogas ~ TSM_A + TSM_N
  TSM_A ~~ TSM_N
  AUTTCA ~~ drogas
'

model_twof_obs <- '
  TSM_A =~ depresivos + ansiedad + bipolar
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N
  tca        ~ TSM_A + TSM_N
  drogas     ~ TSM_A + TSM_N
  TSM_A ~~ TSM_N
  autolesion ~~ tca + drogas
  tca ~~ drogas
'

model_A <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N
  tca        ~ TSM_A + TSM_N
  drogas     ~ TSM_A + TSM_N
  TSM_A ~~ TSM_N
  autolesion ~~ tca + drogas
  tca ~~ drogas
'

model_B <- '
  TSM_A =~ depresivos + ansiedad + bipolar
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N
  tca        ~ TSM_A + TSM_N
  drogas     ~ TSM_A + TSM_N
  TSM_A ~~ TSM_N
  autolesion ~~ tca + drogas
  tca ~~ drogas
  ansiedad ~~ depresivos
  bipolar ~~ esquizofrenia
'

model_C <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N + bipolar
  tca        ~ TSM_A + TSM_N + bipolar
  drogas     ~ TSM_A + TSM_N + bipolar
  TSM_A ~~ TSM_N
  bipolar ~~ TSM_A + TSM_N
  depresivos ~~ ansiedad
  bipolar ~~ esquizofrenia
  autolesion ~~ tca + drogas
  tca ~~ drogas
'

model_D <- '
  TSM_A =~ depresivos + ansiedad + bipolar
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N
  tca        ~ TSM_A + TSM_N
  drogas     ~ TSM_A + TSM_N
  depresivos ~~ ansiedad
  bipolar ~~ esquizofrenia
  TSM_A ~~ TSM_N
  autolesion ~~ tca + drogas
  tca ~~ drogas
  bipolar ~~ 0.05*bipolar
'

model_C0 <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N + bipolar
  tca        ~ TSM_A + TSM_N + bipolar
  drogas     ~ TSM_A + TSM_N + bipolar
  TSM_A ~~ TSM_N
  depresivos ~~ ansiedad
  bipolar ~~ esquizofrenia
  autolesion ~~ tca + drogas
  tca ~~ drogas
'

r <- 0.80; sqrt_r <- sqrt(r)
model_E <- glue::glue("
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  BIP =~ {round(sqrt_r, 6)}*bipolar
  bipolar ~~ {round(1 - r, 6)}*bipolar
  autolesion ~ TSM_A + TSM_N + BIP
  tca        ~ TSM_A + TSM_N + BIP
  drogas     ~ TSM_A + TSM_N + BIP
  TSM_A ~~ TSM_N + BIP
  TSM_N ~~ BIP
  depresivos ~~ ansiedad
  bipolar ~~ esquizofrenia
  autolesion ~~ tca + drogas
  tca ~~ drogas
")

model_C_corr <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N + bipolar
  tca        ~ TSM_A + TSM_N + bipolar
  drogas     ~ TSM_A + TSM_N + bipolar
  TSM_A ~~ TSM_N
  bipolar ~~ TSM_A + TSM_N
  bipolar ~~ esquizofrenia
  autolesion ~~ tca + drogas
  tca ~~ drogas
  depresivos ~~ autolesion
  intelectual ~~ tca
'

model_C_AUT_safe <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  AUTTCA =~ autolesion + tca
  AUTTCA ~ TSM_A + TSM_N + bipolar
  drogas ~ TSM_A + TSM_N + bipolar
  TSM_A ~~ TSM_N
  bipolar ~~ TSM_A + TSM_N
  bipolar ~~ esquizofrenia
'

model_C_corr_fixH <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N + bipolar
  tca        ~ TSM_A + TSM_N + bipolar
  drogas     ~ TSM_A + TSM_N + bipolar
  TSM_A ~~ TSM_N
  bipolar ~~ TSM_A + TSM_N
  bipolar ~~ esquizofrenia
  autolesion ~~ tca + drogas
  tca ~~ drogas
  depresivos ~~ autolesion
  intelectual ~~ tca
  esquizofrenia ~~ 0.05*esquizofrenia
'

model_C_corr_marker <- '
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ 1*esquizofrenia + intelectual
  autolesion ~ TSM_A + TSM_N + bipolar
  tca        ~ TSM_A + TSM_N + bipolar
  drogas     ~ TSM_A + TSM_N + bipolar
  TSM_A ~~ TSM_N
  bipolar ~~ TSM_A + TSM_N
  bipolar ~~ esquizofrenia
  autolesion ~~ tca + drogas
  tca ~~ drogas
  depresivos ~~ autolesion
  intelectual ~~ tca
'

# --- Catálogo y ajuste --------------------------------------------------------
models <- tibble::tibble(
  name = c(
    "CFA_1F",
    "CFA_2F_TSM_vs_CONS",
    "CFA_2F_INT_EXT",
    "CFA_3F_TSM_BEHAV_VIOL",
    "CFA_3F_INT_EXT_THOUGHT",
    "CFA_4F_CONG_ADQ_BEH_VIOL",
    "CFA_5F_PCA",
    "CFA_6F_Granular",
    "HO_6F_to_GLOBAL",
    "HO_3F_to_GLOBAL",
    "CFA_5F_Clean",
    "CFA_PCA_Driven",
    "SEM_2F_CONS_on_TSM",
    "SEM_3F_AUTTCA_SUD_on_TSM",
    "PHYSIO_Full",
    "PHYSIO_Clean",
    "SEM_1F_TSM_to_Outcomes",
    "SEM_2F_A_N_to_AUT_SUD",
    "SEM_2F_A_N_ObsOut",
    "SEM_A_Minimal",
    "SEM_B_WithClinCovs",
    "SEM_C_Bip_Exogenous",
    "SEM_D_Bip_In_TSM_A_HeywoodFix",
    "SEM_C0_NoCovsWithTSM",
    "SEM_E_BIP_Latent_r80",
    "SEM_C_Corrected",
    "SEM_C_AUT_Safe",
    "SEM_C_Corrected_Heywood",
    "SEM_C_Corrected_Marker",
    "BIFACTOR_TSM_G_Specifics"
  ),
  label = c(
    "CFA 1F: Factor global",
    "CFA 2F: TSM vs Consecuencias",
    "CFA 2F: Internalizante vs Externalizante",
    "CFA 3F: TSM / Conducta-Adicciones / Violencia",
    "CFA 3F: INT / EXT / Trastorno del pensamiento",
    "CFA 4F: Congénito / Adquirido / Conducta-Add / Violencia",
    "CFA 5F: Basado en PCA (5 factores)",
    "CFA 6F: Granular (Neurodev, Emocional, Psicosis, Conducta, Adicciones, Violencia)",
    "Orden superior: 6F → F_GLOBAL_DALY",
    "Orden superior: 3F → F_GLOBAL",
    "CFA 5F limpio (sin tea/otros_trs/policía)",
    "CFA data-driven (PCA)",
    "SEM simple: CONS ~ TSM",
    "SEM 3F: AUTTCA, SUD ~ TSM",
    "SEM fisiopatológico completo",
    "SEM fisiopatológico limpio",
    "SEM final: 1F TSM → outcomes",
    "SEM 2F: TSM_A & TSM_N → AUTTCA, drogas",
    "SEM 2F: outcomes observados",
    "SEM A: TSM_A/TSM_N mínimos",
    "SEM B: +covs clínicas (ans~dep; bip~esq)",
    "SEM C: bipolar exógeno observado",
    "SEM D: bipolar en TSM_A (var resid 0.05)",
    "SEM C′: sin covs con TSM",
    "SEM E: BIP latente 1 indicador (r=0.80)",
    "SEM C corregido: covs residuales plausibles",
    "SEM C con AUTTCA latente (seguro)",
    "SEM C corregido + Heywood fix (esq)",
    "SEM C corregido con marcador en TSM_N",
    "Bifactor TSM: G + INT/NEURO/PSYCH (ortogonales)"
  ),
  note = c(
    "Unidimensionalidad total; línea base parsimoniosa.",
    "Separa causas TSM de consecuencias/violencia; útil para SEM posterior.",
    "Clásico eje INT/EXT en psicopatología.",
    "Añade violencia como factor aparte además de conducta/adicciones.",
    "Incluye un tercer factor de pensamiento/psicosis.",
    "Divide TSM en congénito/adquirido + consecuencias separadas.",
    "Estructura de 5 dominios sugeridos por PCA.",
    "Seis dominios más finos; mayor interpretabilidad clínica.",
    "Relación de segundo orden que resume seis factores en uno global.",
    "Segundo orden con 3 dominios (TSM/BEHAV/VIOL) hacia un global.",
    "Versión depurada para estabilidad/identificación.",
    "Cargas según PCA, enfoque puramente empírico.",
    "Ruta simple de CONS sobre TSM.",
    "Dos outcomes latentes/observados regresados en TSM.",
    "Causas upstream diferenciadas (neuro vs afectivo).",
    "Versión depurada del fisiopatológico.",
    "Un solo factor TSM prediciendo outcomes.",
    "Dos factores TSM_A/TSM_N hacia AUTTCA y SUD.",
    "Igual que anterior pero outcomes como observados.",
    "Modelo mínimo con 2 indicadores por factor.",
    "Añade comorbilidades residuales plausibles.",
    "Tratamiento de bipolar como predictor exógeno.",
    "Bipolar en el factor con varianza residual fijada (Heywood fix).",
    "Elimina covarianzas de bipolar con TSM para sensibilidad.",
    "Latente de 1 indicador para bipolar con fiabilidad fijada.",
    "Incluye residuales clínicamente plausibles (dep↔aut, int↔tca).",
    "AUTTCA como latente; evita problemas de starts.",
    "Fija var(e) en esquizofrenia para evitar Heywood.",
    "Escalado por marcador en TSM_N.",
    "Bifactor: G captura covarianza común, específicos ortogonales."
  ),
  syntax = c(
    m1, m2, m2_int_ext, m3, m3_int_ext_thought, m4, m5, m6,
    m6_ho, m3_ho,
    m_final, m_pca,
    model_sem_2F, model_sem_3F,
    model_physio, model_physio_clean, model_final, model_twof, model_twof_obs,
    model_A, model_B, model_C, model_D, model_C0, model_E,
    model_C_corr, model_C_AUT_safe, model_C_corr_fixH, model_C_corr_marker,
    m_bifact_tsm
  )
) |> 
  dplyr::mutate(stdlv = ifelse(name == "SEM_C_Corrected_Marker", FALSE, TRUE))

# Complejidad (robusta)
comp <- t(vapply(models$syntax, model_complexity, numeric(5)))
models <- dplyr::bind_cols(models, as.data.frame(comp))

fits <- purrr::pmap(list(models$name, models$syntax, models$stdlv), ~ safe_fit(..1, ..2, stdlv = ..3))

models$converged <- vapply(fits, has_converged, logical(1))
meas_df <- purrr::map_dfr(fits, get_measures_safe)
models <- dplyr::bind_cols(models, meas_df)

models_ord <- models |>
  dplyr::arrange(score, paths, covs, load_defs, dplyr::desc(CFI)) |>
  dplyr::mutate(idx = dplyr::row_number())

# --- Exportar tabla y avisos --------------------------------------------------
readr::write_csv(models_ord |>
                   dplyr::select(idx, name, label, note, load_defs, indicators_extra, paths, covs, score, converged, CFI, TLI, RMSEA, SRMR, AIC, BIC),
                 file = file.path(outdir, "models_table.csv"))

if (any(!models_ord$converged)) {
  message("Modelos NO convergentes: ", paste(models_ord$name[!models_ord$converged], collapse = ", "))
}

fits_named <- purrr::set_names(fits, models$name)
saveRDS(list(models=models_ord, fits=fits_named), file = file.path(outdir, "models_fits.rds"))
print(outdir)

# --- Diagramas del más simple y del más complejo -----------------------------
if (nrow(models_ord) > 0) {
  first <- models_ord$name[[1]]
  last  <- models_ord$name[[nrow(models_ord)]]
  f1 <- fits_named[[first]]; f2 <- fits_named[[last]]
  if (!is.null(f1)) {
    png(file.path(outdir, paste0("diagram_", first, ".png")), width=1400, height=900, res=150)
    semPlot::semPaths(f1, what="std", whatLabels="std", style="lisrel", rotation=2, residuals=FALSE, edge.label.cex=.8, layout="tree", nCharNodes=0)
    dev.off()
  }
  if (!is.null(f2)) {
    png(file.path(outdir, paste0("diagram_", last, ".png")), width=1400, height=900, res=150)
    semPlot::semPaths(f2, what="std", whatLabels="std", style="lisrel", rotation=2, residuals=FALSE, edge.label.cex=.8, layout="tree", nCharNodes=0)
    dev.off()
  }
}

#### secuencial secuestration 
# Requisitos
library(dplyr)
library(tidyr)

# 1) Verifica que estén las columnas necesarias
need <- c("country","iso3","year",
          "autolesion","violencia_interpersonal","policia_ejecuciones")
stopifnot(all(need %in% names(sem_data)))

# 2) Construye outcome_compuesto (suma) y su z-score
X <- sem_data %>%
  mutate(
    # Suma "efectos" (maneja NA como 0 o ignóralos: elige una de las dos líneas)
    efectos_suma = autolesion + violencia_interpersonal + policia_ejecuciones,
    # Si prefieres ignorar NA en la suma, usa esta en su lugar:
    # efectos_suma = rowSums(across(c(autolesion, violencia_interpersonal, policia_ejecuciones)), na.rm = TRUE),
    
    # z-score global (media 0, sd 1)
    outcome_compuesto = as.numeric(scale(efectos_suma))
  ) %>%
  # Deja sólo lo necesario para el flujo de “sequential sequestration”
  dplyr::select(country, iso3, year,
         autolesion, violencia_interpersonal, policia_ejecuciones,
         depresivos, ansiedad, bipolar, esquizofrenia, tea, tdah, conducta,
         otros_trs, intelectual, tca, conflicto_terrorismo,  # si existen
         outcome_compuesto) 

X <- sem_data %>%
  mutate(efectos_suma = rowSums(across(c(autolesion, violencia_interpersonal, policia_ejecuciones)), na.rm = TRUE)) %>%
  group_by(year) %>%
  mutate(outcome_compuesto = as.numeric(scale(efectos_suma))) %>%
  ungroup()

X <- sem_data %>%
  mutate(
    autolesion_z              = as.numeric(scale(autolesion)),
    violencia_interpersonal_z = as.numeric(scale(violencia_interpersonal)),
    policia_ejecuciones_z     = as.numeric(scale(policia_ejecuciones)),
    efectos_suma              = autolesion_z + violencia_interpersonal_z + policia_ejecuciones_z,
    outcome_compuesto         = as.numeric(scale(efectos_suma))
  )

library(dplyr)
install.packages("pROC")
library(pROC)


vars_cand <- c("autolesion","violencia_interpersonal","policia_ejecuciones",
               "depresivos","ansiedad","bipolar","esquizofrenia","tdah","tea")

# Función que arma outcome excluyendo v_si_existe (si está en el set original)
mk_outcome_loo <- function(df, v){
  base <- c("autolesion","violencia_interpersonal","policia_ejecuciones")
  keep <- setdiff(base, v)                   # excluye la candidata si está en base
  # si quitas dos y queda uno, va ese; si quitas las tres (p.ej. candidata no está en base), usa las tres
  if (length(keep) == 0) keep <- base
  rowSums(df[keep], na.rm = TRUE)
}

# 1) rank1 sin circularidad
rank1 <- lapply(vars_cand, function(v){
  y_raw <- mk_outcome_loo(X, v)
  y <- as.numeric(scale(y_raw)) > median(as.numeric(scale(y_raw)), na.rm = TRUE)
  x <- X[[v]]
  auc <- tryCatch(pROC::roc(y, x, quiet = TRUE)$auc, error=function(e) NA_real_)
  tibble(var = v, auc = as.numeric(auc))
}) %>% bind_rows() %>% arrange(desc(auc), var)

best1 <- rank1$var[1]
q1 <- quantile(X[[best1]], 0.75, na.rm = TRUE)
G3 <- X[[best1]] >= q1

# 2) rank2 en el remanente, también sin circularidad
resto <- X[!G3, ]
rank2 <- lapply(setdiff(vars_cand, best1), function(v){
  y_raw <- mk_outcome_loo(resto, v)
  y <- as.numeric(scale(y_raw)) > median(as.numeric(scale(y_raw)), na.rm = TRUE)
  x <- resto[[v]]
  auc <- tryCatch(pROC::roc(y, x, quiet = TRUE)$auc, error=function(e) NA_real_)
  tibble(var = v, auc = as.numeric(auc))
}) %>% bind_rows() %>% arrange(desc(auc), var)

best2 <- rank2$var[1]
q2 <- quantile(resto[[best2]], 0.75, na.rm = TRUE)
G2 <- (!G3) & (X[[best2]] >= q2)

X$grado_ss <- factor(ifelse(G3, 3, ifelse(G2, 2, 1)),
                     levels = c(1,2,3), labels = c("Mejor","Intermedio","Peor"))

# 5) Checa separación
X %>%
  group_by(grado_ss) %>%
  summarise(n = n(),
            mean_out = mean(outcome_compuesto, na.rm = TRUE),
            sd_out   = sd(outcome_compuesto, na.rm = TRUE)) %>%
  arrange(grado_ss) -> resumen_grados

print(rank1)
print(rank2)
print(resumen_grados)

library(dplyr)
library(tidyr)

vars_all <- c(
  "depresivos","ansiedad","bipolar","esquizofrenia",
  "tea","tdah","conducta","otros_trs","intelectual",
  "autolesion","alcohol","drogas","tca",
  "conflicto_terrorismo","violencia_interpersonal","policia_ejecuciones"
)

df <- sem_data %>% dplyr::select(country, iso3, year, dplyr::all_of(vars_all))

# Definición de bloques
blocks <- list(
  INT  = c("depresivos","ansiedad","bipolar","esquizofrenia"),
  EXT  = c("conducta","alcohol","drogas"),
  EFEC = c("autolesion","violencia_interpersonal","policia_ejecuciones")
)

# --- correlación var vs. promedio de su bloque (excluyéndose a sí misma) -----
corr_block <- dplyr::bind_rows(lapply(names(blocks), function(b){
  vset <- blocks[[b]]
  dplyr::bind_rows(lapply(vset, function(v){
    others <- setdiff(vset, v)
    if (length(others) == 0) {
      tibble(var = v, block = b, r_block = NA_real_)
    } else {
      # promedio por fila de las "otras" variables del bloque
      block_mean <- rowMeans(as.matrix(df %>% dplyr::select(dplyr::all_of(others))),
                             na.rm = TRUE)
      r <- suppressWarnings(cor(df[[v]], block_mean, use = "pairwise.complete.obs"))
      tibble(var = v, block = b, r_block = as.numeric(r))
    }
  }))
}))

print(corr_block)

screen <- screen %>% left_join(corr_block, by = "var")

# 3) sugerencia (umbrales razonables; ajusta si quieres)
screen <- screen %>%
  mutate(suggest_drop = (pct_na > 25) | (pct_zero > 25) | (is.na(r_block) | r_block < 0.25))

print(screen %>% arrange(desc(suggest_drop), desc(pct_zero), desc(pct_na)))

vars_core10 <- c(
  "depresivos","ansiedad","bipolar","esquizofrenia",   # INT
  "conducta","alcohol","drogas",                       # EXT
  "autolesion","violencia_interpersonal","policia_ejecuciones"  # EFECTOS
)

sem_data_core10 <- sem_data %>%
  dplyr::select(country, iso3, year, all_of(vars_core10)) %>%
  mutate(across(all_of(vars_core10), as.numeric)) %>%
  arrange(country, year)

# Guarda por si quieres replicar figuras/modelos
if (!exists("outdir")) outdir <- here::here("output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(sem_data_core10, file.path(outdir, "sem_data_core10.csv"))

library(lavaan)

# CFA 3 factores (INT / EXT / EFECTOS)
mod_cfa3 <- '
  INT     =~ depresivos + ansiedad + bipolar + esquizofrenia
  EXT     =~ conducta + alcohol + drogas
  EFECTOS =~ autolesion + violencia_interpersonal + policia_ejecuciones
  INT ~~ EXT
'

fit_cfa3 <- cfa(mod_cfa3, data = sem_data_core10, std.lv = TRUE, estimator = "MLR", missing = "fiml")
summary(fit_cfa3, fit.measures = TRUE, standardized = TRUE)

# SEM con factor de MÉTODO para posibles sesgos de reporte (subregistro)
mod_sem_method <- '
  INT     =~ depresivos + ansiedad + bipolar + esquizofrenia
  EXT     =~ conducta + alcohol + drogas
  EFECTOS =~ autolesion + violencia_interpersonal + policia_ejecuciones

  # Factor de método: común a ítems susceptibles a subregistro
  METH    =~ autolesion + violencia_interpersonal + policia_ejecuciones
  INT ~~ EXT
  # Hacer METH ortogonal a INT/EXT por identificación (opcional, recomendado)
  METH ~~ 0*INT
  METH ~~ 0*EXT

  # Estructural
  EFECTOS ~ INT + EXT
'

fit_sem_m <- sem(mod_sem_method, data = sem_data_core10, std.lv = TRUE, estimator = "MLR", missing = "fiml")
summary(fit_sem_m, fit.measures = TRUE, standardized = TRUE)

core <- c("depresivos","ansiedad","bipolar","esquizofrenia",
          "conducta","alcohol","drogas",
          "autolesion","violencia_interpersonal","policia_ejecuciones")

sem_data_std <- sem_data %>%
  mutate(across(all_of(core), ~ scale(log1p(.x))[,1]))   # log1p + z


mod_clean <- '
  INT   =~ depresivos + ansiedad + bipolar + esquizofrenia
  SUBS  =~ alcohol + drogas
  COND  =~ conducta

  # Outcomes observados (no latentes)
  autolesion ~ INT + SUBS + COND
  violencia_interpersonal ~ INT + SUBS + COND
  policia_ejecuciones ~ INT + SUBS + COND

  # correlaciones plausibles
  INT ~~ SUBS + COND
  SUBS ~~ COND
  autolesion ~~ violencia_interpersonal + policia_ejecuciones
  violencia_interpersonal ~~ policia_ejecuciones
'

fit_clean <- sem(mod_clean, data = sem_data_std, std.lv = TRUE,
                 estimator = "MLR", missing = "fiml")
summary(fit_clean, fit.measures = TRUE, standardized = TRUE)

mod_alt <- '
  INT   =~ depresivos + ansiedad + bipolar + esquizofrenia
  SUBS  =~ alcohol + drogas
  COND  =~ conducta
  EFEC  =~ 1*violencia_interpersonal + policia_ejecuciones   # fija la escala

  EFEC ~ INT + SUBS + COND
  autolesion ~ INT + SUBS + COND

  INT ~~ SUBS + COND
  SUBS ~~ COND
  autolesion ~~ EFEC
'
fit_alt <- sem(mod_alt, data = sem_data_std, std.lv = TRUE,
               estimator = "MLR", missing = "fiml")
summary(fit_alt, fit.measures = TRUE, standardized = TRUE)

# p.ej., por grado_ss ya calculado
fit_mg <- sem(mod_clean, data = sem_data_std, group = "grado_ss",
              std.lv = TRUE, estimator = "MLR", missing = "fiml")

# ------------------- RESCATE: pipeline mínimo y robusto ----------------------
library(dplyr)
library(lavaan)
library(glue)
library(tibble)
library(readr)

# ---- A) Dataset: asegúrate de tener SUBS_idx y (si existe) grado_ss ----
core <- c("depresivos","ansiedad","bipolar","esquizofrenia",
          "conducta","alcohol","drogas",
          "autolesion","violencia_interpersonal","policia_ejecuciones")

sem_data_std <- sem_data %>%
  mutate(across(all_of(core), ~ scale(log1p(.x))[,1]))

# pega grado_ss si tienes X
if (exists("X") && all(c("country","iso3","year","grado_ss") %in% names(X))) {
  sem_data_std <- sem_data_std %>%
    left_join(X %>% dplyr::select(country, iso3, year, grado_ss),
              by = c("country","iso3","year"))
}

# índice observado de sustancias
sem_data_std <- sem_data_std %>%
  mutate(SUBS_idx = scale(rowMeans(cbind(alcohol, drogas), na.rm = TRUE))[,1])

stopifnot("SUBS_idx" %in% names(sem_data_std))

# ---- B) Modelos ----
M0_cfa3 <- '
  INT     =~ depresivos + ansiedad + bipolar + esquizofrenia
  EXT     =~ conducta + alcohol + drogas
  EFECTOS =~ autolesion + violencia_interpersonal + policia_ejecuciones
  INT ~~ EXT
'

r_COND <- 0.70; sqrt_r <- sqrt(r_COND)

M1_clean_comp <- glue('
  INT   =~ 1*depresivos + ansiedad + bipolar + esquizofrenia
  COND  =~ {sqrt_r}*conducta
  conducta ~~ {1 - r_COND}*conducta

  autolesion ~ INT + SUBS_idx + COND
  violencia_interpersonal ~ INT + SUBS_idx + COND
  policia_ejecuciones ~ INT + SUBS_idx + COND

  INT ~~ COND
  autolesion ~~ violencia_interpersonal + policia_ejecuciones
  violencia_interpersonal ~~ policia_ejecuciones
')

M2_efec2 <- glue('
  INT   =~ 1*depresivos + ansiedad + bipolar + esquizofrenia
  COND  =~ {sqrt_r}*conducta
  conducta ~~ {1 - r_COND}*conducta

  EFEC  =~ 1*violencia_interpersonal + policia_ejecuciones
  EFEC  ~ INT + SUBS_idx + COND
  autolesion ~ INT + SUBS_idx + COND

  INT ~~ COND
  violencia_interpersonal ~~ policia_ejecuciones
  autolesion ~~ EFEC
')

# ---- C) Ajuste "safe" (sin muffleWarning) ----
fit_safe <- function(expr) {
  out <- tryCatch(suppressWarnings(eval(expr)),
                  error = function(e) NULL)
  out
}

fit_M0 <- fit_safe(quote(cfa(M0_cfa3, data=sem_data_std, std.lv=FALSE,
                             estimator="MLR", missing="fiml")))
# fallback si M0 truena
if (is.null(fit_M0) || !isTRUE(tryCatch(lavInspect(fit_M0, "converged"), error=function(e) FALSE))) {
  M0_alt <- '
    INT =~ depresivos + ansiedad + bipolar + esquizofrenia
    EXT =~ conducta + alcohol + drogas
    INT ~~ EXT
  '
  fit_M0 <- fit_safe(quote(cfa(M0_alt, data=sem_data_std, std.lv=FALSE,
                               estimator="MLR", missing="fiml")))
}

fit_M1 <- fit_safe(quote(sem(M1_clean_comp, data=sem_data_std, std.lv=FALSE,
                             estimator="MLR", missing="fiml")))
fit_M2 <- fit_safe(quote(sem(M2_efec2,     data=sem_data_std, std.lv=FALSE,
                             estimator="MLR", missing="fiml")))

# ---- D) Tabla de métricas (solo lo que esté disponible) ----
pick <- c("cfi","tli","rmsea","srmr","aic","bic")

get_measures <- function(fit, pick) {
  ok <- !is.null(fit) && isTRUE(tryCatch(lavInspect(fit, "converged"), error=function(e) FALSE))
  if (!ok) return(setNames(rep(NA_real_, length(pick)), pick))
  as.numeric(fitmeasures(fit, pick))
}

fits <- list(M0_ref = fit_M0, M1_clean_comp = fit_M1, M2_efec2 = fit_M2)
meas <- lapply(fits, get_measures, pick = pick)
tab <- tibble(modelo = names(fits)) %>%
  bind_cols(as.data.frame(do.call(rbind, meas)))
names(tab)[-1] <- pick

print(tab)

# ---- E) Guarda (si quieres) ----
if (!exists("outdir")) outdir <- here::here("output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(tab, file.path(outdir, "models_final_shortlist.csv"))
saveRDS(fits, file.path(outdir, "models_final_fits.rds"))
message("Listo. Shortlist guardado en: ", outdir)

# --- BIVARIADOS: correlaciones y regresiones simples (fix de select) ----------
library(dplyr); library(tidyr); library(purrr)
library(broom); library(lmtest); library(sandwich); library(MASS)

# 1) define predictores y outcomes
preds <- c("depresivos","ansiedad","bipolar","esquizofrenia",
           "conducta","alcohol","drogas")
outs  <- c("autolesion","violencia_interpersonal","policia_ejecuciones")

# 2) std dataset (log1p + z) + índice de efectos opcional
df_std <- sem_data %>%
  dplyr::mutate(across(all_of(c(preds, outs)), ~ scale(log1p(.x))[,1])) %>%
  dplyr::mutate(EFEC_idx = scale(rowMeans(cbind(autolesion,
                                                violencia_interpersonal,
                                                policia_ejecuciones),
                                          na.rm = TRUE))[,1])

outs2 <- c(outs, "EFEC_idx")

# 3) helpers (con namespaces explícitos)
fit_ols_hc3 <- function(y, x, data){
  f <- as.formula(paste(y, "~", x))
  m <- lm(f, data = data)
  ct <- lmtest::coeftest(m, vcov = sandwich::vcovHC(m, type = "HC3"))
  tibble::tibble(
    term   = x,
    beta   = unname(ct[2,1]),
    se     = unname(ct[2,2]),
    t      = unname(ct[2,3]),
    p      = unname(ct[2,4]),
    r2     = summary(m)$r.squared,
    n      = sum(complete.cases(model.frame(f, data)))
  )
}

fit_rlm <- function(y, x, data){
  f <- as.formula(paste(y, "~", x))
  m <- tryCatch(MASS::rlm(f, data = data), error = function(e) NULL)
  if (is.null(m)) return(tibble::tibble(term=x, beta_rlm=NA_real_, se_rlm=NA_real_, t_rlm=NA_real_))
  tt <- tryCatch(broom::tidy(m), error=function(e) NULL)
  if (is.null(tt)) return(tibble::tibble(term=x, beta_rlm=NA_real_, se_rlm=NA_real_, t_rlm=NA_real_))
  tt <- tt %>% dplyr::filter(term == x) %>%
    dplyr::transmute(term, beta_rlm = estimate, se_rlm = std.error, t_rlm = statistic)
  tt
}

corr_xy <- function(y, x, data){
  v <- suppressWarnings(cor(data[[y]], data[[x]], use = "pairwise.complete.obs"))
  tibble::tibble(term = x, r = as.numeric(v))
}

# 4) barrido bivariado: para cada outcome y cada predictor
biv_results <- purrr::map_dfr(outs2, function(y){
  purrr::map_dfr(preds, function(x){
    ols <- fit_ols_hc3(y, x, df_std)
    rl  <- fit_rlm(y, x, df_std)
    rr  <- corr_xy(y, x, df_std)
    dplyr::bind_cols(tibble::tibble(outcome = y),
                     rr %>% dplyr::select(r),
                     ols,
                     rl %>% dplyr::select(beta_rlm, se_rlm, t_rlm))
  })
}) %>%
  dplyr::arrange(outcome, dplyr::desc(abs(beta)))  # orden por |beta| OLS

print(biv_results)

# 5) tabla bonita y guardado
biv_pretty <- biv_results %>%
  dplyr::mutate(sig = dplyr::case_when(
    p < .001 ~ "***",
    p < .01  ~ "**",
    p < .05  ~ "*",
    TRUE     ~ ""
  )) %>%
  dplyr::select(outcome, predictor = term, r,
                beta_ols = beta, se_ols = se, t_ols = t, p_ols = p, sig,
                beta_rlm, se_rlm, t_rlm, r2, n)

if (!exists("outdir")) outdir <- here::here("output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(biv_pretty, file.path(outdir, "bivariados_resumen.csv"))
message("Bivariados guardados en: ", file.path(outdir, "bivariados_resumen.csv"))

library(lavaan)
r_COND <- 0.75  # fiabilidad asumida de 'conducta' como 1-indicador

mod_sem_final <- glue::glue('
  # Medición
  INT  =~ 1*depresivos + ansiedad + bipolar + esquizofrenia
  EXT1 =~ {sqrt(r_COND)}*conducta
  conducta ~~ {1 - r_COND}*conducta
  INT ~~ EXT1

  # EFEC latente de 2 indicadores, cargas iguales para estabilidad
  EFEC =~ l*violencia_interpersonal + l*policia_ejecuciones

  # Estructural (alcohol observado, drogas fuera del principal)
  EFEC      ~ INT + EXT1 + alcohol
  autolesion~ INT + EXT1 + alcohol

  # Correlaciones residuales entre outcomes si usas observados
  violencia_interpersonal ~~ policia_ejecuciones
')

fit_sem_final <- sem(mod_sem_final, data = sem_data_std,
                     std.lv = FALSE, estimator = "MLR", missing = "fiml")
summary(fit_sem_final, fit.measures = TRUE, standardized = TRUE)

# SIN EXT latente; conducta y alcohol observadas
mod_sem_min <- '
  # Medición INT
  INT =~ 1*depresivos + ansiedad + bipolar + esquizofrenia

  # EFEC latente a 2 ítems (cargas iguales)
  EFEC =~ l*violencia_interpersonal + l*policia_ejecuciones

  # Estructural
  autolesion ~ INT                    # lo que dictan los bivariados
  EFEC       ~ INT + alcohol          # alcohol es fuerte y estable; INT opcional

  # (Opcional) mete conducta como covariable observada si quieres:
  # autolesion ~ INT + conducta
  # EFEC       ~ INT + alcohol + conducta

  # (Opcional) cov residuales, si EFEC se reemplaza por observados:
  violencia_interpersonal ~~ policia_ejecuciones
'

fit_min <- sem(mod_sem_min, data = sem_data_std, std.lv = FALSE,
               estimator = "MLR", missing = "fiml")
summary(fit_min, fit.measures = TRUE, standardized = TRUE)

mod_dep_aut <- '
  INT  =~ 1*depresivos + ansiedad + bipolar + esquizofrenia
  EFEC =~ l*violencia_interpersonal + l*policia_ejecuciones

  autolesion ~ depresivos        # directo, como mostraron bivariados
  EFEC       ~ alcohol + INT     # alcohol fijo; INT opcional

  violencia_interpersonal ~~ policia_ejecuciones
'
fit_dep_aut <- sem(mod_dep_aut, data=sem_data_std, std.lv=FALSE,
                   estimator="MLR", missing="fiml")
summary(fit_dep_aut, fit.measures=TRUE, standardized=TRUE)

mod_int_cov <- '
  INT  =~ 1*depresivos + ansiedad + bipolar + esquizofrenia
  EFEC =~ l*violencia_interpersonal + l*policia_ejecuciones

  autolesion ~ INT
  EFEC       ~ alcohol + INT

  depresivos ~~ autolesion       # comorbilidad/medición
  violencia_interpersonal ~~ policia_ejecuciones
'
fit_dep_aut <- sem(mod_int_cov, data=sem_data_std, std.lv=FALSE,
                   estimator="MLR", missing="fiml")
summary(fit_dep_aut, fit.measures=TRUE, standardized=TRUE)

mod_final_stable <- '
  # Medición
  INT  =~ 1*depresivos + ansiedad + bipolar + esquizofrenia
  EFEC =~ 1*violencia_interpersonal + policia_ejecuciones   # una fija, otra libre

  # Estructural (hallazgos clave)
  autolesion ~ depresivos
  EFEC       ~ alcohol + INT

  # Residuos plausibles
  violencia_interpersonal ~~ policia_ejecuciones
  # Si no da problemas numéricos, puedes conservar:
  # EFEC ~~ autolesion
'

fit_final_stable <- sem(mod_final_stable, data = sem_data_std,
                        std.lv = TRUE, estimator = "MLR", missing = "fiml")
summary(fit_final_stable, fit.measures = TRUE, standardized = TRUE)

sem_data_std <- sem_data_std %>%
  mutate(EFEC_idx = scale(rowMeans(cbind(violencia_interpersonal,
                                         policia_ejecuciones), na.rm=TRUE))[,1])
# 1) Checa que existan las columnas
stopifnot(all(c("violencia_interpersonal","policia_ejecuciones") %in% names(sem_data_std)))

# 2) Construye el índice (promedio de los dos; maneja NA)
#    Si tus columnas ya están z (como antes), promediar está bien;
#    luego lo re–estandarizamos para que quede con media 0, sd 1.
sem_data_std <- sem_data_std %>%
  dplyr::mutate(
    EFEC_idx = rowMeans(dplyr::across(c(violencia_interpersonal, policia_ejecuciones)),
                        na.rm = TRUE),
    EFEC_idx = as.numeric(scale(EFEC_idx))
  )

# (Opcional) confirma que quedó
summary(sem_data_std$EFEC_idx)

mod_A <- '
  INT =~ depresivos + ansiedad + bipolar + esquizofrenia
  autolesion ~ depresivos
  EFEC_idx  ~ alcohol + INT
'
fit_A <- sem(mod_A, data=sem_data_std, std.lv=TRUE, estimator="MLR", missing="fiml")

summary(fit_A, fit.measures = TRUE, standardized = TRUE)

mod_B <- '
  INT  =~ depresivos + ansiedad + bipolar + esquizofrenia
  EFEC =~ 1*violencia_interpersonal + policia_ejecuciones

  autolesion ~ depresivos
  EFEC       ~ alcohol + INT

  # Fix mínimo para evitar Heywood en violencia
  violencia_interpersonal ~~ 0.05*violencia_interpersonal
  # (opcional) remover esta si estabiliza mejor:
  # violencia_interpersonal ~~ 0*policia_ejecuciones
'
fit_B <- sem(mod_B, data=sem_data_std, std.lv=TRUE, estimator="MLR", missing="fiml")

summary(fit_B, fit.measures = TRUE, standardized = TRUE)

mod_C <- '
  INT =~ depresivos + ansiedad + bipolar + esquizofrenia

  autolesion             ~ depresivos
  violencia_interpersonal ~ alcohol + INT
  policia_ejecuciones     ~ alcohol + INT

  violencia_interpersonal ~~ policia_ejecuciones
'
fit_C <- sem(mod_C, data=sem_data_std, std.lv=TRUE, estimator="MLR", missing="fiml")
summary(fit_C, fit.measures = TRUE, standardized = TRUE)

# --- Gráfico de SEM (publicación) --------------------------------------------
library(lavaan)
library(semPlot)
library(qgraph)

plot_sem_clean <- function(fit, file_base = "sem_clean",
                           hide_below = 0.10,
                           p_sig_only = FALSE,
                           preview = TRUE,          # <- mostrar en RStudio
                           save_png = TRUE,
                           save_svg = FALSE,        # <- requiere svglite
                           width_px = 2000, height_px = 1300, dpi = 220) {
  stopifnot(inherits(fit, "lavaan"))
  
  # --- Qué edges mostrar (por p y |est.std|) ---
  ss <- standardizedSolution(fit)
  ss <- ss[ss$op %in% c("~","=~","~~"), ]
  ss$plot_keep <- if (p_sig_only) ss$pvalue < .05 else TRUE
  ss$plot_keep <- ss$plot_keep & (abs(ss$est.std) >= hide_below | ss$op == "=~")
  
  # --- Modelo semPlot + etiquetas/anchos ---
  sp  <- semPlot::semPlotModel(fit, what = "std")
  lab <- sprintf("%.2f", sp@Pars$est)
  idx <- match(paste(sp@Pars$lhs, sp@Pars$op, sp@Pars$rhs),
               paste(ss$lhs,       ss$op,     ss$rhs))
  show_edge <- !is.na(idx) & ss$plot_keep[idx]
  lab[!show_edge] <- ""
  w <- pmax(1, 8*abs(sp@Pars$est)); w[!show_edge] <- 0.5
  
  # --- Etiquetas de nodos ---
  nice <- c(
    depresivos="Depresivos", ansiedad="Ansiedad", bipolar="Bipolar",
    esquizofrenia="Esquizofrenia",
    autolesion="Autolesión",
    violencia_interpersonal="Violencia\ninterpersonal",
    policia_ejecuciones="Policía/\nejecuciones",
    INT="Internalizante"
  )
  repl <- setNames(sp@Vars$name, sp@Vars$name)
  common <- intersect(names(nice), names(repl))
  repl[common] <- nice[common]
  sp@Vars$label <- unname(repl[sp@Vars$name])
  
  col_lat  <- "#2E86AB"; col_man <- "#7D8C8D"; col_edge <- "#333333"
  
  draw_it <- function() semPlot::semPaths(
    sp,
    what="std", whatLabels="no",
    style="lisrel", layout="tree", rotation=2,
    residuals=FALSE, intercepts=FALSE,
    edgeLabels = lab,
    edge.width = w,
    edge.color = col_edge,
    curvePivot=TRUE, curveAdjacent=TRUE,
    color = list(lat=col_lat, man=col_man),
    nCharNodes=0, sizeLat=10, sizeMan=8, mar=c(6,6,6,6)
  )
  
  # 1) Vista previa en el panel Plots
  if (isTRUE(preview)) draw_it()
  
  # 2) Guardar PNG
  if (isTRUE(save_png) && !is.null(file_base)) {
    png(paste0(file_base, ".png"), width = width_px, height = height_px, res = dpi)
    draw_it(); dev.off()
    message("PNG: ", normalizePath(paste0(file_base, ".png")))
  }
  
  # 3) Guardar SVG (con svglite)
  if (isTRUE(save_svg) && !is.null(file_base)) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      message("Instala 'svglite' para SVG: install.packages('svglite')")
    } else {
      svglite::svglite(paste0(file_base, ".svg"),
                       width = width_px/96, height = height_px/96)
      draw_it(); dev.off()
      message("SVG: ", normalizePath(paste0(file_base, ".svg")))
    }
  }
}


# Define la sintaxis del modelo SOLO con los comandos de lavaan
model_C_corr_syntax_FINAL <- '
  # Factores Latentes de TSM (CFA)
  TSM_A =~ depresivos + ansiedad
  TSM_N =~ esquizofrenia + intelectual

  # Rutas de Regresión (SEM)
  autolesion ~ TSM_A + TSM_N + bipolar
  tca ~ TSM_A + TSM_N + bipolar
  drogas ~ TSM_A + TSM_N + bipolar

  # Covarianzas
  TSM_A ~~ TSM_N
  bipolar ~~ TSM_A + TSM_N
  bipolar ~~ esquizofrenia

  # Covarianzas Residuales
  autolesion ~~ tca + drogas
  tca ~~ drogas
  depresivos ~~ autolesion
  intelectual ~~ tca
'

# Ahora, ejecuta el ajuste con esta sintaxis limpia
fit_final <- sem(model_C_corr_syntax_FINAL,
                 data = sem_data,
                 estimator = "MLR",
                 missing = "fiml",
                 std.lv = TRUE)

# Verifica la convergencia
if (!lavaan::lavInspect(fit_final, "converged")) {
  stop("El modelo no ha convergido. Revisa los datos o la sintaxis.")
}
message("¡El modelo ha ajustado con éxito! 🎉")

# Extraer métricas de ajuste
fit_measures <- lavaan::fitmeasures(fit_final, c("cfi", "tli", "rmsea", "srmr"))

print("--- Métricas de Ajuste (SEM_C_Corrected) ---")
print(fit_measures)

# 1. Extraer la solución estandarizada
results_std <- lavaan::standardizedSolution(fit_final) %>%
  # 2. CREAR LA COLUMNA DE RELACIÓN PRIMERO
  dplyr::mutate(
    Relacion = paste(lhs, op, rhs),
    # 3. Calcular la significancia
    p_sig = ifelse(pvalue < 0.001, "***", 
                   ifelse(pvalue < 0.01, "**", 
                          ifelse(pvalue < 0.05, "*", ""))),
    # 4. Redondear valores
    est.std = round(est.std, 3),
    pvalue = round(pvalue, 3)
  ) %>%
  # 5. Seleccionar las columnas finales (filtrando solo cargas y regresiones)
  dplyr::filter(op %in% c("=~", "~")) %>% 
  dplyr::select(Relacion, Beta_std = est.std, P_val = pvalue, Signif = p_sig)


print("--- Rutas de Regresión Clave (Beta Estandarizado) ---")
regression_paths <- results_std %>% dplyr::filter(grepl("~", Relacion))
print(knitr::kable(regression_paths, caption = "Rutas de Regresión Clave", format = "pipe"))

print("--- Cargas Factoriales Clave (Beta Estandarizado) ---")
loading_paths <- results_std %>% dplyr::filter(grepl("=~", Relacion))
print(knitr::kable(loading_paths, caption = "Cargas Factoriales", format = "pipe"))

# Asegura que 'here' esté cargado para la ruta del archivo
library(here) 

library(semPlot)

print("--- Diagrama del Modelo SEM_C_Corrected (Viewer) ---")

library(lavaan)
library(dplyr)

# Asumiendo que 'fit_final' es el objeto de tu modelo SEM ajustado

# 1. Extraer la solución estandarizada
results_std <- lavaan::standardizedSolution(fit_final) %>%
  # 2. CREAR LA COLUMNA DE RELACIÓN
  dplyr::mutate(
    Relacion = paste(lhs, op, rhs),
    # 3. Calcular la significancia (p_sig)
    p_sig = ifelse(pvalue < 0.001, "***", 
                   ifelse(pvalue < 0.01, "**", 
                          ifelse(pvalue < 0.05, "*", ""))),
    # 4. Redondear valores
    est.std = round(est.std, 3),
    pvalue = round(pvalue, 3)
  ) %>%
  # 5. Seleccionar las columnas finales (filtrando solo cargas y regresiones)
  dplyr::filter(op %in% c("=~", "~")) %>% 
  dplyr::select(Relacion, Beta_std = est.std, P_val = pvalue, Signif = p_sig)

# 6. Mostrar Rutas de Regresión Clave
# Esta tabla contiene tu hallazgo (autolesion ~ TSM_A) y las rutas anómalas
regression_paths <- results_std %>% dplyr::filter(grepl("~", Relacion))
print(knitr::kable(regression_paths, caption = "Rutas de Regresión Clave (Beta Estandarizado)", format = "pipe"))


library(semPlot)

# Código para generar el diagrama en el Viewer
semPlot::semPaths(fit_final, 
                  what = "est",           # Mostrar estimaciones (no estandarizadas, como en la imagen que subiste)
                  whatLabels = "est",     
                  layout = "tree",        
                  edge.label.cex = 0.5,   
                  fade = FALSE,           
                  residuals = FALSE,      
                  style = "lisrel",       
                  label.cex = 1,        
                  esize = 0.03)


# --- Paquetes base ---
library(dplyr)
library(tidyr)
library(lavaan)
library(semPlot)
library(here)

# --- 1️⃣ Cargar datos DALYs base ---
# Ajusta la ruta a donde esté tu CSV (usa la misma que en tus scripts)
raw <- readr::read_csv(
  "/Users/sofiadiaz/Library/Mobile Documents/com~apple~CloudDocs/tesis-anahuac/data/dalys_causas_latam.csv",
  show_col_types = FALSE
) %>%
  janitor::clean_names()

# --- 2️⃣ Mapeo de variables al formato de análisis ---
map_tbl <- tibble::tibble(
  cause = c(
    "Alcohol use disorders",
    "Anxiety disorders",
    "Attention-deficit/hyperactivity disorder",
    "Autism spectrum disorders",
    "Bipolar disorder",
    "Conduct disorder",
    "Depressive disorders",
    "Drug use disorders",
    "Eating disorders",
    "Idiopathic developmental intellectual disability",
    "Other mental disorders",
    "Schizophrenia",
    "Self-harm"
  ),
  target = c(
    "alcohol", "ansiedad", "tdah", "tea", "bipolar", "conducta",
    "depresivos", "drogas", "tca", "intelectual", "otros_trs",
    "esquizofrenia", "autolesion"
  )
)

sem_data <- raw %>%
  inner_join(map_tbl, by = "cause") %>%
  group_by(country, year, target) %>%
  summarise(val = sum(dalys, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = target, values_from = val) %>%
  arrange(country, year)

# --- 3️⃣ Escalamiento (log + z) ---
sem_data_std <- sem_data %>%
  mutate(across(where(is.numeric), ~ scale(log1p(.x))[, 1]))

# --- 4️⃣ Crear carpeta de salida ---
outdir <- here::here("output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- 5️⃣ Definir y ajustar modelo sin violencia ---
modelo_sin_violencia <- '
  # Factores latentes según PCA (sin violencia)
  Afectivo_Conductual =~ autolesion + ansiedad + depresivos + alcohol
  Disruptivo_Estructural =~ conducta + otros_trs + drogas + tca +
                            tdah + tea + esquizofrenia + intelectual + bipolar

  # Relación estructural principal
  Afectivo_Conductual ~ Disruptivo_Estructural

  # Correlación residual
  Afectivo_Conductual ~~ Disruptivo_Estructural
'

fit_sin_violencia <- sem(
  modelo_sin_violencia,
  data = sem_data_std,
  std.lv = TRUE,
  estimator = "MLR",
  missing = "fiml"
)

# --- 6️⃣ Resultados y guardado ---
summary(fit_sin_violencia, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
fitMeasures(fit_sin_violencia, c("cfi", "tli", "rmsea", "srmr"))

# Guardar coeficientes estandarizados del modelo
params_sin_violencia <- parameterEstimates(fit_sin_violencia, standardized = TRUE) %>%
  dplyr::filter(op %in% c("~", "=~")) %>%
  dplyr::select(
    lhs, op, rhs,
    Beta_std = std.all,   # <-- columna correcta universal
    se, pvalue
  )

# Revisa el resultado
print(head(params_sin_violencia))

# Guardar CSV
readr::write_csv(params_sin_violencia,
                 file.path(outdir, "params_sin_violencia.csv"))


# --- 7️⃣ Diagrama (colores del PCA) ---
png(file.path(outdir, "diagram_sin_violencia.png"),
    width = 1600, height = 1000, res = 180)
semPlot::semPaths(
  fit_sin_violencia,
  what = "std", whatLabels = "std",
  style = "lisrel", layout = "tree",
  residuals = FALSE, edge.label.cex = 0.9,
  fade = FALSE,
  color = list(lat = "#e31793", man = "#ec855f"),
  sizeLat = 10, sizeMan = 7, nCharNodes = 0
)
dev.off()

message("✅ Modelo sin violencia ajustado y guardado en: ", outdir)
