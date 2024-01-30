# logaritmo da de uma parte da condicional completa para os graus de liberdade
# mu eh uma matriz que ja calcula xbeta_j+bDelta por conta disso ja estar
# calculado no meio do gibbs sampler(linhas amostras, colunas componentes)
# prob, sigma, lambda vetores com um elemento para cada componente
# resp e a var resposta
Phinu = function(nu, prob, resp, mu, sigma, lambda){
 resp_mu = resp-mu
 n = length(resp)
 G = length(prob)
 sum(log(rowSums(exp(matrix(rep(log(prob),n),byrow = T,nrow = n, ncol = G)+sapply(1:G, function(a) dst(resp_mu[,a], 0, sigma[a], lambda[a], nu,log=T))))))
}