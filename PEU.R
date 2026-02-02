# ==============================================================================
# SCRIPT: PEU.R (Procedimento de Etapa Única)
# ==============================================================================
source("./funcoes.R")

# 1. PARÂMETROS ORIGINAIS -----------------------------------------------------
arquivo_estacao_unica <- file.choose()
distribuicao <- "gumbel"      # "gumbel" ou "gev"
metodo_empirico <- "gringorten"
D_IDF <- c(seq(10, 50, 10), seq(60, 1440, 60)) # Minutos
D_IDF <- D_IDF/60
Tr <- c(2, 10, 30, 50, 100)

# Parâmetros de Otimização 
itermax_global <- 2000
itermax_local  <- 10000
populacao_deoptim_mult <- 30
estrategia_deoptim <- 4  
metodos_busca_local <- c("NLOPT_LN_COBYLA", "NLOPT_LN_PRAXIS", "NLOPT_LN_BOBYQA")

# 2. CARGA E PREPARAÇÃO -------------------------------------------------------
dados_base <- carregar_dados(arquivo_estacao_unica)
cod_estacao_atual <- as.character(dados_base$dados_estacoes$codEstacao[1])
cidade <- as.character(dados_base$dados_estacoes$municipio)[1]

dados_prep <- preparar_dados(dados_base$dados_limpos)
df_completo <- completar_dados(dados_prep)
pcp_10min <- df_completo %>% select(datahora, valorMedida, waterYear)

annualMax <- calcular_maximos_anuais(pcp_10min, D_IDF)
annualMax <- readRDS(file.choose()) # Caso queira carregar o arquivo já gerado

# Limpeza e Ordenação
for (j in seq_len(ncol(annualMax))) {
  col <- na.omit(annualMax[, j])
  annualMax[, j] <- c(sort(col, decreasing = TRUE), rep(NA, nrow(annualMax) - length(col)))
}
annualMax <- annualMax[rowSums(is.na(annualMax)) < ncol(annualMax), , drop = FALSE]

# Tempo de Retorno Empírico
n <- nrow(annualMax)
rank_mat <- matrix(1:n, nrow=n, ncol=ncol(annualMax))
Tjl <- (n + 0.12) / (rank_mat - 0.44)
Tjl[is.na(annualMax)] <- NA

cols_validas <- colSums(!is.na(annualMax)) >= 2
annualMax_otim <- annualMax[, cols_validas, drop=FALSE]
Tjl_otim <- Tjl[, cols_validas, drop=FALSE]
D_IDF_ativos <- D_IDF[cols_validas]

# 3. OTIMIZAÇÃO GLOBAL (DEoptim) ----------------------------------------------
cat("\n>>> INICIANDO OTIMIZAÇÃO GLOBAL (DEoptim)...\n")
if (distribuicao == "gumbel") {
  lower <- c(a=0, b=0, c=0, d=0); upper <- c(a=10000, b=500, c=100, d=1)
  n_vars <- 4
} else {
  lower <- c(a=0, b=0, k=0, c=0, d=0); upper <- c(a=10000, b=500, k=100, c=100, d=1)
  n_vars <- 5
}

otim_global <- DEoptim(
  fn = function(p) calcular_erro_absoluto_total(p, as.matrix(annualMax_otim), Tjl_otim, D_IDF_ativos, distribuicao),
  lower = lower, upper = upper,
  control = DEoptim.control(
    strategy = estrategia_deoptim,
    NP = max(50, populacao_deoptim_mult * n_vars),
    itermax = itermax_global,
    trace = 50,
    parallelType = 0
  )
)
par_global <- otim_global$optim$bestmem
obj_global <- otim_global$optim$bestval
cat("DEoptim Obj:", obj_global, "\n")

# 4. REFINAMENTO LOCAL  ----------------------------------------
cat("\n>>> INICIANDO REFINAMENTO LOCAL (nloptr - 3 Algoritmos)...\n")

melhor_res <- list(solution = par_global, objective = obj_global, method = "DEoptim")

for (algo in metodos_busca_local) {
  cat(" Testando algoritmo:", algo, "...\n")
  tryCatch({
    res <- nloptr(
      x0 = par_global,
      eval_f = function(x) calcular_erro_absoluto_total(x, as.matrix(annualMax_otim), Tjl_otim, D_IDF_ativos, distribuicao),
      lb = lower, ub = upper,
      opts = list(algorithm = algo, maxeval = itermax_local, xtol_rel = 1e-100, print_level = 0)
    )
    
    if (res$objective < melhor_res$objective) {
      melhor_res <- list(solution = res$solution, objective = res$objective, method = algo)
      cat("  -> Novo melhor encontrado!\n")
    }
  }, error = function(e) cat("  Erro no algoritmo:", algo, "\n"))
}

par_final <- melhor_res$solution
names(par_final) <- names(lower)
cat("Melhor método final:", melhor_res$method, "\n")
print(par_final)

# 5. RESULTADOS ---------------------------------------------------------------
IDF_final <- matrix(NA, length(Tr), length(D_IDF))
for(i in seq_along(Tr)) {
  for(j in seq_along(D_IDF)) {
    IDF_final[i, j] <- calcular_intensidade_idf_modificada(Tr[i], D_IDF[j], par_final, distribuicao)
  }
}

# Métricas
IDF_calc_otim <- matrix(NA, nrow(annualMax_otim), ncol(annualMax_otim))
for(j in 1:ncol(annualMax_otim)) {
  for(i in 1:nrow(annualMax_otim)) {
    if(!is.na(Tjl_otim[i,j])) {
      IDF_calc_otim[i,j] <- calcular_intensidade_idf_modificada(Tjl_otim[i,j], D_IDF_ativos[j], par_final, distribuicao)
    }
  }
}
metrics <- calcular_metricas_erro(annualMax_otim, IDF_calc_otim)
print(metrics)

dir_out <- file.path("./Resultados/PEU", cidade, distribuicao)
if (!dir.exists(dir_out)) {dir.create(dir_out, recursive = TRUE)}
saveRDS(list(params=par_final, metrics=metrics, IDF=IDF_final), file.path(dir_out, "resultados_peu.rds"))
plotar_curvas_idf(IDF_final, Tr, D_IDF, cidade, cod_estacao_atual, dir_out, "PEU")
