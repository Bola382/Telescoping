# logaritmo da condicional completa para os graus de liberdade com priori 
# ExpLoc(loc, alpha)
# mu eh uma matriz que ja calcula xbeta_j+bDelta por conta disso ja estar
# calculado no meio do amostrador(linhas amostras, colunas componentes)
# prob, sigma, lambda vetores com um elemento para cada componente
# resp: a var resposta
# alpha hiperparametro de nu
# loc: parametro de localizacao
# A INDICADORA DE NU > loc FOI OMITIDA, ESTOU ASSUMINDO QUE A PROPOSTA TEM SUPORTE 
# EM (loc, Inf)
fullnuExpLoc = function(nu,alpha, prob, resp, mu, sigma, lambda,loc){
  resp_mu = resp-mu
  n = length(resp)
  G = length(prob)
  -alpha*(nu-loc) + sum(log(rowSums(exp(matrix(rep(log(prob),n),byrow = T,nrow = n, ncol = G)+sapply(1:G, function(a) dst(resp_mu[,a], 0, sigma[a], lambda[a], nu,log=T))))))
}
