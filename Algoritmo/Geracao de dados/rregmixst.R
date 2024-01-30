# Gera observacoes do modelo de misturas de regressão t assimetrico
# covs: matiz com covariaveis do modelo, para ter intercepto a primeira coluna deve ser de uns
# probs: pesos das componentes, determina o numero de componentes
# beta: matriz com coeficientes de regressão
        # colunas correspondem as componentes enquanto as linhas correspondem 
        # as colunas de covs
# sigma: vetor com parametros de escala, um para cada componente
# lambda: vetor com parametros de forma, um para cada componente
# nu: graus de liberdade, e comum entre as componentes, apenas um escalar

source("rst.R")
rregmixst = function(covs, probs, beta, sigma, lambda, nu){
 covs = as.matrix(covs)
 n = nrow(covs)
 G = length(probs)
 
 Z = sample(1:G, n, replace = T, prob = probs)
 
 # centralizando a media em x*beta
 b = -sqrt(nu/pi)*gamma((nu-1)/2)/gamma(nu/2)
 Delta = sigma*(lambda/sqrt(1+lambda^2))
 
 mu = covs%*%beta + matrix(rep(b*Delta,n), nrow = n, ncol = G, byrow = T)
 
 resp = sapply(1:n, function(a) rst(1,mu[a,Z[a]], sigma[Z[a]], lambda[Z[a]], nu))
 
 return(list(resposta = resp, grupo = Z))
}


# covs = cbind(1,rnorm(10000))
# probs = c(.4,.6)
# beta = cbind(c(5,2), c(-5,-1))
# sigma = c(2,4)
# lambda = c(3,-8)
# nu = 5.5
# 
# mu = covs%*%beta + b*Delta
# 
# b = -sqrt(nu/pi)*exp(lgamma(nu-1/2)-lgamma(nu/2)) #exp(log(.))
# Delta = sigma*(lambda/sqrt(1+lambda^2))
# 
# resul = rregmixst(covs,probs,beta,sigma,lambda,nu)
# y = resul$resposta
# 
# plot(covs[,2],y, pch = 16, col = resul$grupo)
# abline(beta[,1])
# abline(beta[,2], col=2)
# 
# 
# hist(y-covs%*%beta-b*Delta,breaks=30, freq=F)
# curve(dst(x,0,sigma,lambda,nu),add=T)
# lines(density(y-covs%*%beta-b*Delta), col = 2)

