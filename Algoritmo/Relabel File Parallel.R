setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("Re-labeling File.R")
library(foreach)
library(doParallel)

R = 5 # numero de replicacoes
file = "test1"

cl <- makeForkCluster(12) # linux
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)

# arrumar funcao
result = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 relabel(file,re)
}
stopCluster(cl)
