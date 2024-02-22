# ==============================================================================
# Aplica "Re-labeling File.R" em paralelo
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source("Re-labeling File.R")
library(foreach)
library(doParallel)

R = 10 # numero de replicacoes
file = "test2"

cl <- makeCluster(4) 
registerDoParallel(cl)
clusterSetRNGStream(cl, 1)

result = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 relabel(file,re)
}
stopCluster(cl)

