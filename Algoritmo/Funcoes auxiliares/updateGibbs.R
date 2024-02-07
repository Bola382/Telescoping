# atualiza somente parametros que dependem da componente
# param: lista com valores na ultima atualizacao dos parametros
#        deve estar nomeada na ordem c("beta","tau2","Delta")
# arg: tipo do parametro
# comp: qual componente atualizar
# y: vetor de respostas
# X: matriz de covariaveis, primeira coluna de uns
# z: variavel indicadora de componentes
# u,t: variaveis latentes da representacao da ST
# M: numero de observacoes por componente
# c: desvio padrao da priori de beta
# r,s: hiperparametros da priori de tau2
# omega: desvio padrao da priori de Delta

# lembrando que beta matriz pXGplus, tau2 e Delta vetores de tamanho Gplus

updateGibbs = function(param,arg,comp,y,X,z,u,t,M,b,c,r,s,omega){
 paramnames = names(param)
 index = which(z==comp)
 n = nrow(X)
 p = ncol(X)
 Xj = X[index,]
 uj = u[index]
 Uj = diag(uj, nrow = M[comp], ncol = M[comp]) 
 yj = y[index] 
 tj = t[index]
 
 novo = switch(arg,
               beta = full_beta.TS(comp,b,c,n,p,Xj,Uj,yj,tj,M,param[["tau2"]],param[["Delta"]]),
               tau2 = full_tau2.TS(comp,b,r,s,Xj,uj,yj,tj,M,param[["beta"]],param[["Delta"]]),
               Delta = full_Delta.TS(comp,b,omega,Xj,uj,yj,tj,M,param[["beta"]],param[["tau2"]]))
 if(arg == "beta"){
  param[[arg]][,comp] = novo
 }else{
  param[[arg]][comp] = novo
 }
 return(param) 
}
