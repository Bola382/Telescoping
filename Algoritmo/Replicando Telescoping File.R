# ==============================================================================
# Replicacoes em paralelo do Telescoping utlizando o "Funcao Telescoping File.R"
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
load("Geracao de dados/n100type1.RData")
source("Funcoes auxiliares/logbnb.R")
source("Funcao Telescoping File reescrita.R")
library(foreach)
library(doParallel)

dados = as.matrix(data[,-1]) # primeira coluna contem grupos reais

R = 10 # numero de replicacoes
Q =  155000 # iteracoes por replica
burn = 5000 # burnin
thin = 50 # salto
Gmax = 100 # maximo de componentes possivel
lprio_G = logbnb(1:Gmax) # valor da log priori de G para 1:G.max
nburn = ceiling((Q-burn)/thin) # numero de amostras pos burn

# hiperparametros
c = 10 # de beta 
eta = 0; omega = 10 # da Delta
r = 2.1
s = 1.1 # da tau2
loc = 1.1

file = "reescrita n100type1"
dir.create(paste0("Outputs/Paralelo/",file))
sapply(1:R, function(a) dir.create(paste0("Outputs/Paralelo/",file,"/",a),showWarnings = F))
cl <- makeForkCluster(15) # linux
registerDoParallel(cl)
clusterSetRNGStream(cl, 7481)

tempo = foreach(re = 1:R, .packages = "compiler",.verbose=T) %dopar% {
 telescope(file,re,dados,ncomps=10,G.max=Gmax,lprio_G,
           c, eta, omega, r, s, loc,
           phi=.3,phigamma=2.5,Q,burn,thin)
};stopCluster(cl)
