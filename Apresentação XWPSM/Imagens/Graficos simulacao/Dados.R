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

col1 = t(sapply(data1$comp, function(a) if(a==3){c(1,0,0)}else{if(a==2){c(0,1,0)}else{c(0,0,1)}}))
rgb1 = apply(col1,1, function(a){rgb(a[1],a[2],a[3])})
reorder1 = c(3,2,1)
recolor1 = c("#FF0000","#00FF00","#0000FF")

col2 = t(sapply(data2$comp, function(a) if(a==3){c(1,0,0)}else{if(a==2){c(0,0,1)}else{c(0,1,0)}}))
rgb2 = apply(col2,1, function(a){rgb(a[1],a[2],a[3])})
reorder2 = c(2,3,1)
recolor2 = c("#FF0000","#00FF00","#0000FF")

par(mfrow=c(1,2),mar = c(5.1, 4.1, 4.1, 2.1))
plot(data1$x1,data1$resp, col=rgb1, pch = 16, cex = .6, 
     ylab = "Resposta", xlab = "Covariável")
for(i in 1:length(probs_data1)){
 abline(beta_data1[,i], col=recolor1[reorder1[i]])
}
par(mar = c(5.1,1.5, 4.1, 2.1))
plot(data2$x1,data2$resp, col=rgb2, pch=16, cex = .6, 
     ylab = "", xlab = "Covariável",ylim = c(-22,30))
for(i in 1:length(probs_data2)){
 abline(beta_data2[,i], col=recolor2[reorder2[i]])
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

contagem = function(vec,categ){
 sapply(1:categ, function(a) sum(vec==a))
}

file = "test2"
R = 10

beta_list1 = list()
z_list1 = matrix(NA, nrow = R, ncol = nrow(data1))
G_list1 = list()


for(repl in 1:R){
 # recuperando parametros por rep
 caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
 beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
 z.samp = read.table(paste0(caminho,"z.txt"), h =T)
 z_tmp = apply(z.samp,2,moda) # classificacao pelo maximo a posteriori
 beta_tmp = matrix(colMeans(beta.samp),nrow = 2, ncol = ncol(beta.samp)/2)
 
 # rerotulando baseado nas retas de reg
 ordem = crossprod(c(1,max(data1$x1)),beta_tmp)
 novaordem = order(ordem, decreasing = T)
 
 z_list1[repl,] = sapply(z_tmp, function(a) which(novaordem==a)) # reordena
 beta_list1[[repl]] = sapply(1:(ncol(beta.samp)/2), function(a) beta_tmp[,novaordem[a]])
 param = read.table(paste0(caminho,"param.txt"), h =T)
 G.samp = param$G; rm("param")
 G_list1[[repl]] = moda(G.samp)
}

Gplus_list1 = lapply(beta_list1, function(a) ncol(a))

# assume que o numero de Gplus nao ultrapassa 3
col1 = apply(z_list1, 2,function(a) contagem(a,Gplus_list1[[1]])/R)
col1 = apply(col1,2,function(a) rgb(a[1],a[2],a[3]))

# dados 2

file = "test4"
quais = read.table(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/test4quais.txt"))
quais = which(quais[,1]>=.5)

R = 100

beta_list2 = list()
z_list2 = matrix(NA, nrow = length(quais), ncol = nrow(data1))
G_list2 = list()
i = 1

for(repl in quais){
 # recuperando parametros por rep
 caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
 beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
 z.samp = read.table(paste0(caminho,"z.txt"), h =T)
 z_tmp = apply(z.samp,2,moda) # classificacao pelo maximo a posteriori
 beta_tmp = matrix(colMeans(beta.samp),nrow = 2, ncol = ncol(beta.samp)/2)
 
 # rerotulando baseado nas retas de reg
 ordem = crossprod(c(1,max(data1$x1)),beta_tmp)
 novaordem = order(ordem, decreasing = T)
 
 z_list2[i,] = sapply(z_tmp, function(a) which(novaordem==a)) # reordena
 beta_list2[[i]] = sapply(1:(ncol(beta.samp)/2), function(a) beta_tmp[,novaordem[a]])
 param = read.table(paste0(caminho,"param.txt"), h =T)
 G.samp = param$G; rm("param")
 G_list2[[i]] = moda(G.samp)
 i = i+1
}

Gplus_list2 = lapply(beta_list2, function(a) ncol(a))

# assume que o numero de Gplus nao ultrapassa 3
col2 = apply(z_list2, 2,function(a) contagem(a,3)/R)
col2 = apply(col2,2,function(a) rgb(a[1],a[2],a[3]))

# plot
par(mfrow=c(1,2),mar = c(5.1, 4.1, 4.1, 2.1))
plot(data1$x1,data1$resp, col=col1, pch = 16, cex = .6, 
     ylab = "Resposta", xlab = "Covariável")
par(mar = c(5.1,1.5, 4.1, 2.1))
plot(data2$x1,data2$resp, col=col2, pch=16, cex = .6, 
     ylab = "", xlab = "Covariável",ylim = c(-22,30))

