# ==============================================================================
# SCRIPT: PP.R (Procedimento Padrão)
# ==============================================================================

source("./funcoes.R")

# 1. PARÂMETROS ORIGINAIS -----------------------------------------------------
arquivo_estacao_unica <- file.choose()
D_min <- c(seq(10,50,10), seq(60,1440,60))
D_h <- D_min / 60 
D_nomes <- sprintf("%.6f", D_h)
Tr <- c(2, 10, 30, 50, 100)
distribuicao <- "GEV" #GEV ou GUMBEL

# Parâmetros DEoptim e nloptr
de_control <- list(strategy=4, itermax=10000, trace=200) 
local_methods <- c("NLOPT_LN_COBYLA", "NLOPT_LN_PRAXIS", "NLOPT_LN_BOBYQA") 
# 2. DADOS --------------------------------------------------------------------
dados_base <- carregar_dados(arquivo_estacao_unica)
cod_estacao_atual <- as.character(dados_base$dados_estacoes$codEstacao)
cidade <- as.character(dados_base$dados_estacoes$municipio)[1]

dados_prep <- preparar_dados(dados_base$dados_limpos)
df_completos <- completar_dados(dados_prep, mes_inicio = 9)
pcp_10min <- data.frame(
  datetime = df_completos$datahora,
  P.mm. = df_completos$valorMedida,
  waterYear = df_completos$waterYear
)
annualMax <- calcular_maximos_anuais(pcp_10min, D_min)
#Ou
annualMax <- readRDS(file.choose())
colnames(annualMax) <- D_nomes 

# 3. GEV E CURVA ALVO ---------------------------------------------------------
res_gev <- ajustar_modelo_gev(annualMax, modelo = distribuicao)
IDFe <- calcular_curva_idf(res_gev$parameters, Tr, D_nomes)

# Vetores para Otimização
target_vec <- unlist(IDFe)
grid <- expand.grid(Tr = Tr, D = D_h)
valid <- is.finite(target_vec)
I_opt <- target_vec[valid]
R_opt <- grid$Tr[valid]
D_opt <- grid$D[valid]

# Função Objetivo
obj_fun <- function(p) {
  # p: a, b, c, d
  calc <- IDF_type_III(R_opt, D_opt, p[1], p[2], p[3], p[4])
  if(any(!is.finite(calc))) return(1e10)
  sum((calc - I_opt)^2) + penalidade_suavidade(calc)
}

lower <- c(0, 0, 0, 0); upper <- c(10000, 500, 20, 1)

# 4. OTIMIZAÇÃO  ----------------------------------------------
cat("\n>>> OTIMIZAÇÃO GLOBAL (DEoptim)...\n")
res_de <- DEoptim(obj_fun, lower, upper, control = do.call(DEoptim.control, de_control))
par_de <- res_de$optim$bestmem
val_de <- res_de$optim$bestval

cat("\n>>> REFINAMENTO LOCAL (Loop 3 Algoritmos)...\n")
melhor_res <- list(solution = par_de, objective = val_de, method = "DEoptim")

for(metodo in local_methods) {
  cat(" Refinando com", metodo, "...\n")
  tryCatch({
    fit_loc <- nloptr(
      x0 = par_de,
      eval_f = obj_fun,
      lb = lower, ub = upper,
      opts = list(algorithm = metodo, maxeval = 10000, xtol_rel = 1e-100, print_level = 0)
    )
    
    if(fit_loc$objective < melhor_res$objective) {
      melhor_res <- list(solution = fit_loc$solution, objective = fit_loc$objective, method = metodo)
      cat("  -> Novo melhor encontrado!\n")
    }
  }, error = function(e) cat("  Erro no método:", metodo, "\n"))
}

par_final <- melhor_res$solution
names(par_final) <- c("a", "b", "c", "d")
cat("Melhor método:", melhor_res$method, "\n")
print(par_final)

# 5. RESULTADOS ---------------------------------------------------------------
IDF_final <- matrix(NA, length(Tr), length(D_min))
for(i in seq_along(Tr)) {
  IDF_final[i, ] <- IDF_type_III(Tr[i], D_h, par_final[1], par_final[2], par_final[3], par_final[4])
}
# Métricas do Ajuste (Equação vs GEV Teórica)
metricas_ajuste <- calcular_metricas_erro(IDFe, IDF_final)
cat("Métricas do Ajuste (vs GEV):\n"); print(metricas_ajuste)

# --- [INSERÇÃO] CÁLCULO DAS MÉTRICAS GLOBAIS (EQUAÇÃO vs DADOS OBSERVADOS) ---
# 1. Prepara dados observados (Ordenação + Tr Empírico)
annualMax_sorted <- annualMax
for (j in 1:ncol(annualMax_sorted)) {
  val <- na.omit(annualMax[, j])
  annualMax_sorted[, j] <- c(sort(val, decreasing = TRUE), rep(NA, nrow(annualMax) - length(val)))
}

n_anos <- nrow(annualMax_sorted)
rank_mat <- matrix(1:n_anos, nrow=n_anos, ncol=ncol(annualMax_sorted))
Tjl_emp <- (n_anos + 0.12) / (rank_mat - 0.44)

# 2. Consolida em dataframe único para validação
df_obs <- as.data.frame(annualMax_sorted) %>% mutate(id = 1:n()) %>% pivot_longer(-id, names_to="D_h", values_to="I_obs")
df_tr  <- as.data.frame(Tjl_emp) %>% mutate(id = 1:n()) %>% pivot_longer(-id, names_to="D_h", values_to="Tr_emp")

dados_globais <- left_join(df_obs, df_tr, by=c("id", "D_h")) %>%
  mutate(D_h = as.numeric(D_h)) %>% # Converte nome da coluna (horas) para numérico
  na.omit()

# 3. Calcula o valor predito pela equação final para cada ponto real observado
dados_globais$I_pred <- IDF_type_III(dados_globais$Tr_emp, dados_globais$D_h, 
                                     par_final["a"], par_final["b"], par_final["c"], par_final["d"])

# 4. Métricas Globais (Realidade vs Modelo)
metricas_globais <- calcular_metricas_erro(dados_globais$I_obs, dados_globais$I_pred)
cat("\nMétricas Globais (vs Observado):\n"); print(metricas_globais)
# -----------------------------------------------------------------------------

dir_out <- file.path("./Resultados/PP", cidade, distribuicao)
if (!dir.exists(dir_out)) {dir.create(dir_out, recursive = TRUE)}

# [AJUSTE] Salva lista completa com as novas métricas e dados de validação
saveRDS(list(
  params = par_final, 
  IDF = IDF_final,
  metrics_fit = metricas_ajuste,      # Quão bem a equação colou na GEV
  metrics_global = metricas_globais,  # Quão bem a equação representa a chuva real
  validation_data = dados_globais     # Dados ponto a ponto para diagnósticos
), file.path(dir_out, "resultados_pp.rds"))

plotar_curvas_idf(IDF_final, Tr, D_min, cidade, cod_estacao_atual, dir_out, paste("PP -", distribuicao))
