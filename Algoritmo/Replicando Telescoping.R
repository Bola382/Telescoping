setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
load("Geracao de dados/dados.Rdata")
source("Funcoes auxiliares/logbnb.R")
source("Telescoping + Re-labeling.R")
library(compiler)
library(parallel)

data = data[,-1] # primeira coluna contem grupos reais

R = 5 # numero de replicacoes
Q =  260000 # iteracoes por replica
burn = 10000 # burnin
thin = 50 # salto
Gmax = 100 # maximo de componentes possivel
lprio_G = logbnb(1:Gmax) # valor da log priori de G para 1:G.max
nburn = ceiling((Q-burn)/thin) # numero de amostras pos burn

cl <- makeForkCluster(R) # linux
clusterExport(cl, varlist = names(Filter(is.function, mget(ls(all=T)))))
clusterSetRNGStream(cl, 1)

resul = clusterCall(cl,telescope,data=data,ncomps=10,nclusters=10,G.max=Gmax,
            lprio_G=lprio_G,phi=.3,phigamma=2.5,
            Q=Q,burn=burn,thin=thin)
stopCluster(cl)
