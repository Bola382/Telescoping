setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
load("Geracao de dados/dados1cov.RData")
source("Funcoes auxiliares/logbnb.R")
source("Funcao Telescoping File.R")
library(foreach)
library(doParallel)

dados = as.matrix(data[,-1]) # primeira coluna contem grupos reais

R = 10 # numero de replicacoes
Q =  260000 # iteracoes por replica
burn = 10000 # burnin
thin = 50 # salto
Gmax = 100 # maximo de componentes possivel
lprio_G = logbnb(1:Gmax) # valor da log priori de G para 1:G.max
nburn = ceiling((Q-burn)/thin) # numero de amostras pos burn

file = "test2"
dir.create(paste0("Outputs/Paralelo/",file))
sapply(1:R, function(a) dir.create(paste0("Outputs/Paralelo/",file,"/",a),showWarnings = F))
cl <- makeForkCluster(15) # linux
registerDoParallel(cl)
clusterSetRNGStream(cl, 7481)

tempo = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 telescope(file,re,dados,ncomps=10,nclusters=10,G.max=Gmax,lprio_G,
           phi=.3,phigamma=2.5,Q,burn,thin)
}
stopCluster(cl)
