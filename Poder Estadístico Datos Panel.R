library(fixest)
library(purrr)
library(dplyr)
library(tibble)


# ---------- 1) Simulación de poder con FE país + año (versión mejorada) --------------------------
power_panel_sim <- function(G = 33, T = 23,
                            beta1 = 0.10,
                            sigma_eps = 1, sigma_i = .6, sigma_t = .3,
                            rho_x = 0.5,
                            R = 1000, alpha = 0.05, two_sided = TRUE,
                            seed = 123){
  set.seed(seed)
  id <- rep(seq_len(G), each = T)
  tt <- rep(seq_len(T), times = G)
  mu_i <- rnorm(G, 0, sigma_i)
  tau_t <- rnorm(T, 0, sigma_t)
  mu_i_rep <- mu_i[id]
  tau_t_rep <- tau_t[tt]
  x_i <- rnorm(G, 0, 1)
  x_it <- rho_x * x_i[id] + sqrt(1 - rho_x^2) * rnorm(G*T, 0, 1)
  x_it <- scale(x_it)[,1]
  
  reject <- logical(R)
  for (r in seq_len(R)){
    eps <- rnorm(G*T, 0, sigma_eps)
    y <- mu_i_rep + tau_t_rep + beta1 * x_it + eps
    
    data_sim <- data.frame(y = y, x_it = x_it, id = id, tt = tt)
    
    est <- feols(y ~ x_it | id + tt, data = data_sim, cluster = ~id, warn = FALSE)
    
    # === MODIFICACIÓN CLAVE: Extraer el p-valor directamente del resumen ===
    # Esto es más robusto y preciso que el cálculo manual.
    p <- summary(est, se = "cluster")$coeftable["x_it", "Pr(>|t|)"]
    
    reject[r] <- (!is.na(p)) && p < alpha
  }
  power <- mean(reject)
  se <- sqrt(power*(1-power)/R)
  list(power = power, se = se,
       design = list(G=G, T=T, beta1=beta1, sigma_eps=sigma_eps, sigma_i=sigma_i, sigma_t=sigma_t, rho_x=rho_x, R=R, alpha=alpha))
}

# ---------- 2) Búsqueda del T mínimo (o G mínimo) para una potencia objetivo ---
find_T_for_power <- function(target = .80, G = 33, T_grid = 8:30, ...){
  res <- map_df(T_grid, ~{
    out <- power_panel_sim(G = G, T = .x, ...)
    tibble(T = .x, power = out$power)
  })
  list(table = res, T_min = res %>% filter(power >= target) %>% slice_head(n=1) %>% pull(T))
}

find_G_for_power <- function(target = .80, T = 20, G_grid = 10:50, ...){
  res <- map_df(G_grid, ~{
    out <- power_panel_sim(G = .x, T = T, ...)
    tibble(G = .x, power = out$power)
  })
  list(table = res, G_min = res %>% filter(power >= target) %>% slice_head(n=1) %>% pull(G))
}

# ================== EJEMPLOS ==========================================
# 1) Poder con 33 países x 23 años, efecto estandarizado beta1 = 0.10
pwr1 <- power_panel_sim(G = 33, T = 23, beta1 = 0.10, R = 1000)
print(pwr1)

# 2) ¿Cuántos años necesito (con 33 países) para alcanzar 80% de poder?
gridT <- find_T_for_power(target = .80, G = 33, beta1 = 0.10, R = 600)
print(gridT)

# 3) ¿Cuántos países necesito si T = 15?
gridG <- find_G_for_power(target = .80, T = 15, beta1 = 0.10, R = 600)
print(gridG)
