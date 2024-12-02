# logposteriori do modelo de misturas de regressoes ST
# precisa carregar dst.R da pasta Geracao de dados,
# contagem.R e logbnb.R

# y: vetor de amostras
# z: vetor de alocacoes
# X: matriz de planejamento, primeira coluna de uns, demais colunas covariaveis
# p: vetor de pesos das componentes
# beta: matriz de coefs de regressao, componentes por coluna
# tau2: vetor de parametros de escala por componente
# Delta: vetor de parametros de forma por componente
# nu: graus de liberdade
# alpha: hiperparametro de nu
# gama: hiperparametro de p
# c: desvio padrao a priori de beta
# r,s: hiperparametros de tau2
# eta, omega: media e desvio padrao a priori de Delta
# loc: parametro de localizacao da priori de nu

logpost = function(y, z, X, p, beta, tau2, Delta, nu, alpha, gama, c, r, s, eta, omega, loc){
  tau2 = as.numeric(tau2)
  Delta = as.numeric(Delta)
  
  n = length(y) # numero de amostras
  G = length(p) # numero de comps
  m = contagem(z,G) # vetor com o numero de amostras por componente
  Gplus = sum(m!=0)
  
  # recuperando parametros originais
  sigma = sqrt(tau2 + Delta^2)
  lambda = Delta/sqrt(tau2)
  
  mu = X%*%beta # matriz n x G contendo x'beta_j em cada coluna
  b = ifelse(nu > 1, -exp(log(nu)/2 - log(pi)/2 + lgamma((nu-1)/2) - lgamma(nu/2)), 
             stop("nu <= 1"))
  
  # termo no interior da segunda soma da log-verossimilhanca
  loglik_aux = sapply(1:G, function(j) 
    sapply(1:n, function(i) dst(y[i], mu[i,j] + b*Delta[j], sigma[j], lambda[j], nu, log=T)))
  
  # log-verossimilhanca (menos constantes)
  loglik = sum(m[1:Gplus]*log(p[1:Gplus])) + sum(sapply(1:G, function(j) sum(loglik_aux[which(z==j),j])))
  
  # log-priori (menos constantes), priori de alpha omitida por ser constante
  logprio = lgamma(gama) - G*lgamma(gama/G) + 
    ifelse(any(p == 0),0,sum((m + gama/G - 1)*log(p))) + # z e p
    2*log(gama) - 4.5*log(1+2*gama) + # gama
    log(alpha) - alpha*(nu-loc) - # nu
    1/(2*c^2)*sum(sapply(1:G, function(j) sum(beta[,j]^2))) + # beta
    (r-1)*sum(log(tau2)) - s*sum(1/tau2) - # tau2
    1/(2*omega^2)*sum((Delta-eta)^2) + # Delta
    logbnb(G) # G
  
  
  # log-posteriori (menos constantes)
  loglik + logprio
}
