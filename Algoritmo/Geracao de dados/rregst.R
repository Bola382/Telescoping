# Gera observacoes do modelo de regressão t assimetrico
# covs: matiz com covariaveis do modelo, para ter intercepto a primeira coluna deve ser de uns
# beta: vetor de coeficientes de regressão, um pra cada coluna de covs
# sigma: parametro de escala, comum a todas as observacoes
# lambda: parametro de forma
# nu: graus de liberdade

source("rst.R")
rregst = function(covs, beta, sigma, lambda, nu){
 covs = as.matrix(covs)
 n = nrow(covs) # numero de amostras
 
 # centralizando a media em x*beta
 b = -sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
 Delta = sigma*(lambda/sqrt(1+lambda^2))
 
 mu = covs%*%beta + b*Delta
 
 rst(n,mu,sigma,lambda,nu)
}


# covs = cbind(1,rnorm(1))
# beta = c(1,2)
# sigma = 2
# lambda = -8
# nu = 5.5
# 
# b = -sqrt(nu/pi)*exp(lgamma(nu-1/2)-lgamma(nu/2)) #exp(log(.))
# Delta = sigma*(lambda/sqrt(1+lambda^2))
# 
# y = rregst(covs,beta,sigma,lambda,nu)
# plot(covs[,2],y)
# hist(y-covs%*%beta-b*Delta,breaks=30, freq=F)
# curve(dst(x,0,sigma,lambda,nu),add=T)
# lines(density(y-covs%*%beta-b*Delta), col = 2)
