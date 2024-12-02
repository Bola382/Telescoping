# ==============================================================================
# Aplica "Relabel ECR File.R" em paralelo
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("Relabel ECR File.R")
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
library(foreach)
library(doParallel)

R = 10 # numero de replicacoes
file = "reescrita n100type1"

load("Geracao de dados/n100type1.RData")

dados = data[,-1]
y = dados[,1]
X = cbind(1,dados[,-1])

cl <- makeForkCluster(15)
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)

result = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 relabelECR(file,re, y, X)
};stopCluster(cl)

