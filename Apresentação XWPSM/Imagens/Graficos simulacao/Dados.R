setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
wd = getwd()
local = paste0(stringr::str_sub(wd,end=-46),"Algoritmo/Geracao de dados")
setwd(local)

load("dados1cov.RData")

beta_data1 = beta
sigma_data1 = sigma
lambda_data1 = lambda
nu_data1 = nu
probs_data1 = probs
data1 = data

rm(list = setdiff(ls(),c(grep("_data1$",ls(),value = T),"data1")))

load("dados1cov2.RData")

beta_data2 = beta
sigma_data2 = sigma
lambda_data2 = lambda
nu_data2 = nu
probs_data2 = probs
data2 = data

rm(list = setdiff(ls(),c(grep("_data1$",ls(),value = T),c(grep("_data2$",ls(),value = T),"data1","data2"))))

par(mfrow=c(1,2),mar = c(5.1, 4.1, 4.1, 2.1))
plot(data1$x1,data1$resp, col=data1$comp, pch = 16, cex = .6, 
     ylab = "Resposta", xlab = "Covariável")
for(i in 1:length(probs_data1)){
 abline(beta_data1[,i], col=i)
}
par(mar = c(5.1,1.5, 4.1, 2.1))
plot(data2$x1,data2$resp, col=data2$comp, pch=16, cex = .6, 
     ylab = "", xlab = "Covariável",ylim = c(-22,30))
for(i in 1:length(probs_data2)){
 abline(beta_data2[,i], col=i)
}

################################################################################
# estimativas
################################################################################

# dados 1
wd = getwd()
local = paste0(stringr::str_sub(wd,end=-17))
setwd(local)

moda = function(vec){
 unique(sort(vec))[which.max(table(vec))]
}

file = "test2"

beta_list1 = list()
G_list1 = list()

R = 10
for(repl in 1:R){
 caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
 beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
 beta_list1[[repl]] = matrix(colMeans(beta.samp),nrow = 2, ncol = ncol(beta.samp)/2)
 param = read.table(paste0(caminho,"param.txt"), h =T)
 G.samp = param$G; rm("param")
 G_list1[[repl]] = moda(G.samp)
}

Gplus_list1 = lapply(beta_list1, function(a) ncol(a))

# dados 2

file = "test3"

beta_list2 = list()
G_list2 = list()

R = 10
for(repl in 1:R){
 caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
 beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
 beta_list2[[repl]] = matrix(colMeans(beta.samp),nrow = 2, ncol = ncol(beta.samp)/2)
 param = read.table(paste0(caminho,"param.txt"), h =T)
 G.samp = param$G; rm("param")
 G_list2[[repl]] = moda(G.samp)
}

Gplus_list2 = lapply(beta_list2, function(a) ncol(a))

# plot
par(mfrow=c(1,2),mar = c(5.1, 4.1, 4.1, 2.1))
plot(data1$x1,data1$resp, col="darkgray", pch = 16, cex = .6, 
     ylab = "Resposta", xlab = "Covariável")
set.seed(1)
for(repl in 1:R){
 for(i in 1:Gplus_list1[[repl]]){
  abline(beta_list1[[repl]][,i] + rnorm(1,sd=.1), col=i)
 }
}
par(mar = c(5.1,1.5, 4.1, 2.1))
plot(data2$x1,data2$resp, col="darkgray", pch=16, cex = .6, 
     ylab = "", xlab = "Covariável",ylim = c(-22,30))
set.seed(1)
for(repl in 1:R){
 for(i in 1:Gplus_list2[[repl]]){
  abline(beta_list2[[repl]][,i] + rnorm(1,sd=.1), col=i)
 }
}
