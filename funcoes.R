# ==============================================================================
# ARQUIVO: funcoes.R
# DESCRIÇÃO: Bibliotecas e funções matemáticas/estatísticas.
# ==============================================================================

# Carregamento de pacotes 
install.packages("pacman")
library(pacman)
pacman::p_load("tidyverse", "lubridate", "zoo", "evd", "DEoptim", "nloptr", "viridis", "gridExtra")

# ------------------------------------------------------------------------------
# 1. DADOS
# ------------------------------------------------------------------------------

carregar_dados <- function(caminho_arquivo) {
  dados <- read.csv(caminho_arquivo, header = TRUE, sep = ",")
  dados_estacoes <- dados %>%
    group_by(municipio, codEstacao, latitude, longitude) %>%
    summarize(.groups = 'drop')
  dados_limpos <- dados %>%
    select(any_of(c('datahora','valorMedida','data','hora')))
  return(list(dados = dados, dados_limpos = dados_limpos, dados_estacoes = dados_estacoes))
}

preparar_dados <- function(dados) {
  dados$datahora <- ifelse(nchar(dados$datahora) == 10, 
                           paste(dados$data, dados$hora, sep = " "),
                           dados$datahora)  
  dados$datahora <- lubridate::ymd_hms(dados$datahora)
  return(dados)
}

completar_dados <- function(dados, mes_inicio = "01") {
  dados <- dados %>% mutate(datahora_truncada = floor_date(datahora, unit = "hour"))
  datahora_completa <- seq(min(dados$datahora_truncada), max(dados$datahora_truncada), by = "hour")
  
  df <- dados %>%
    tidyr::complete(datahora_truncada = datahora_completa) %>%
    mutate(datahora = dplyr::if_else(is.na(datahora), datahora_truncada, datahora)) %>%
    complete(datahora = seq(min(datahora), max(datahora), by = "10 mins")) %>%
    mutate(valorMedida = replace_na(valorMedida, 0))
  
  mes_inicio <- as.integer(mes_inicio)
  df$waterYear <- lubridate::year(df$datahora)
  mask <- lubridate::month(df$datahora) < mes_inicio
  df$waterYear[mask] <- df$waterYear[mask] - 1
  
  ano_inicial <- min(year(df$datahora))
  if(any(year(df$datahora) < ano_inicial)) df$waterYear[year(df$datahora) < ano_inicial] <- ano_inicial - 1
  
  return(df)
}

calcular_maximos_anuais <- function(pcp_10min, D_IDF) {
  # Obter o primeiro e último ano dos dados
  primeiro_ano <- min(pcp_10min$waterYear)
  ultimo_ano <- max(pcp_10min$waterYear)
  
  # Crie um data frame vazio para armazenar os máximos anuais
  annualMax <- data.frame(row.names = as.character(seq(primeiro_ano, ultimo_ano)))
  
  # Criando máximos anuais
  for (d in D_IDF) {
    # Calculando a intensidade média para a duração da tempestade atual
    int_d <- zoo::rollapply(pcp_10min$P.mm., width = d/10, FUN = sum, align = "center", fill = NA) / (d / 10) 
    int_d <- int_d*6 #conversão de mm/10min para mm/h.
    
    # Agrupando por ano e calculando o máximo anual
    annualMax[, as.character(d)] <- tapply(int_d, pcp_10min$waterYear, max, na.rm = TRUE)
  }
  
  # Renomeando as colunas
  colnames(annualMax) <- as.character(D_IDF)
  
  return(annualMax)
}

# ------------------------------------------------------------------------------
# 2. MODELOS ESTATÍSTICOS
# ------------------------------------------------------------------------------

ajustar_modelo_gev <- function(annualMax, modelo = c("GEV", "GUMBEL")) {
  modelo <- match.arg(modelo)
  fit_list <- apply(annualMax, 2, function(x) {
    x <- na.omit(x)
    if(length(x) < 2) return(NULL)
    if (modelo == "GUMBEL") evd::fgev(x, shape = 0, std.err = FALSE) 
    else tryCatch(evd::fgev(x, std.err = TRUE), error = function(e) evd::fgev(x, std.err = FALSE))
  })
  parameters <- t(sapply(fit_list, function(fit) if(is.null(fit)) c(NA,NA,NA) else fit$estimate))
  return(list(fit_list = fit_list, parameters = parameters))
}

calcular_curva_idf <- function(parameters, Tr, D_IDF) {
  IDFe <- matrix(NA, nrow = length(Tr), ncol = length(D_IDF))
  colnames(IDFe) <- D_IDF; rownames(IDFe) <- Tr
  
  for (i in seq_along(D_IDF)) {
    p <- parameters[i,]
    if(any(is.na(p))) next
    prob <- 1 - 1/Tr
    if (length(na.omit(p)) == 3) IDFe[, i] <- evd::qgev(prob, loc=p[1], scale=p[2], shape=p[3])
    else IDFe[, i] <- evd::qgumbel(prob, loc=p[1], scale=p[2])
  }
  return(IDFe)
}

# ------------------------------------------------------------------------------
# 3. EQUAÇÕES MATEMÁTICAS E OBJETIVOS
# ------------------------------------------------------------------------------

IDF_type_III <- function(R, D, a, b, c, d) { (a * R ^ b) / (D + c)^d }

# Funções auxiliares PEU
calcular_b_d <- function(D, c, d) { (D + c)^d }
calcular_a_T_gumb <- function(Tr, a, b) { a * (b - log(-log(1 - 1/Tr))) }
calcular_a_T_gev <- function(Tr, a, b, k) { a * (b + (-log(1 - 1/Tr))^-k) }

calcular_intensidade_idf_modificada <- function(Tr_val, D_val, parametros, distribuicao="gumbel") {
  if (!is.null(names(parametros))) {
    a <- parametros["a"]; b <- parametros["b"]
    if(distribuicao == "gumbel") { c <- parametros["c"]; d <- parametros["d"]; k <- NULL }
    else { k <- parametros["k"]; c <- parametros["c"]; d <- parametros["d"] }
  } else {
    a <- parametros[1]; b <- parametros[2]
    if(distribuicao == "gumbel") { c <- parametros[3]; d <- parametros[4]; k <- NULL }
    else { k <- parametros[3]; c <- parametros[4]; d <- parametros[5] }
  }
  
  if (distribuicao == "gumbel") a_T <- calcular_a_T_gumb(Tr_val, a, b)
  else a_T <- calcular_a_T_gev(Tr_val, a, b, k)
  
  b_d <- calcular_b_d(D_val, c, d)
  val <- a_T / b_d
  return(if(!is.finite(val) || val < 0) NA_real_ else val)
}

# Função Objetivo do PEU
calcular_erro_absoluto_total <- function(parametros, annualMax, Tjl, D_IDF_ativos, distribuicao="gumbel") {
  parametros <- as.numeric(parametros)
  if (any(!is.finite(parametros))) return(1e12)
  
  erro_total_abs <- 0.0
  n_pontos_validos <- 0
  
  # Penalidade se 'd' for muito pequeno (instabilidade) ou parametros negativos
  if(parametros[length(parametros)] < 0.001) return(1e12) 
  
  for (j in seq_len(ncol(annualMax))) {
    dur_atual <- D_IDF_ativos[j]
    obs_col <- annualMax[, j]
    tr_col <- Tjl[, j]
    
    valid <- !is.na(obs_col) & !is.na(tr_col) & tr_col > 0
    if(!any(valid)) next
    
    obs <- obs_col[valid]
    tr <- tr_col[valid]
    
    mod <- numeric(length(obs))
    for(k in seq_along(obs)) {
      mod[k] <- calcular_intensidade_idf_modificada(tr[k], dur_atual, parametros, distribuicao)
    }
    
    if (any(!is.finite(mod) | mod < 0)) {
      erro_total_abs <- erro_total_abs + 1e8 
    } else {
      erro_total_abs <- erro_total_abs + sum(abs(obs - mod))
      n_pontos_validos <- n_pontos_validos + length(obs)
    }
  }
  return(if(n_pontos_validos == 0) 1e12 else erro_total_abs)
}

penalidade_suavidade <- function(I_model) {
  if (any(is.nan(I_model)) || any(is.infinite(I_model))) return(1e10)
  difs <- abs(diff(I_model))
  return(sum(difs > 1000) * 1e6)
}

# ------------------------------------------------------------------------------
# 4. EXPORTAÇÃO E PLOTS
# ------------------------------------------------------------------------------

calcular_metricas_erro <- function(observado, predito) {
  obs <- as.vector(as.matrix(observado)); pred <- as.vector(as.matrix(predito))
  valid <- is.finite(obs) & is.finite(pred)
  obs <- obs[valid]; pred <- pred[valid]
  
  res <- obs - pred
  list(
    RMSE = sqrt(mean(res^2)),
    MAE = mean(abs(res)),
    NSE = 1 - (sum(res^2) / sum((obs - mean(obs))^2)),
    R2 = cor(obs, pred)^2
  )
}

plotar_curvas_idf <- function(IDF_mat, Tr, D_plot, cidade, cod_estacao, dir_saida, prefixo) {
  df <- data.frame(
    D = rep(D_plot, times = length(Tr)),
    Intensidade = as.vector(t(IDF_mat)),
    Tr = factor(rep(Tr, each = length(D_plot)))
  )
  
  g <- ggplot(df, aes(x = D, y = Intensidade, color = Tr)) +
    geom_line(linewidth = 1) +
    scale_color_viridis_d(name = "TR (anos)") +
    labs(title = paste(cidade, "-", cod_estacao), subtitle = prefixo, x = "Duração (min)", y = "Intensidade (mm/h)") +
    theme_minimal() + theme(legend.position = "bottom")
  
  dir.create(dir_saida, recursive = T, showWarnings = F)
  ggsave(file.path(dir_saida, paste0(prefixo, "_linear.png")), g, width=8, height=6, bg="white")
  ggsave(file.path(dir_saida, paste0(prefixo, "_loglog.png")), g + scale_x_log10() + scale_y_log10(), width=8, height=6, bg="white")
}