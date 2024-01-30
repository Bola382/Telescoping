# gera amostras de uma distribuicao t assimetrica
# n: numero de amostras
# mu: parametro de localizacao
# sigma: parametro de escala (raiz de sigma2)
# lambda: parametro de forma
# nu: graus de liberdade

rst = function(n,mu,sigma,lambda,nu){
 # reparametrizacao
 delta = lambda/sqrt(1+lambda^2)
 Delta = sigma*delta
 tau = sigma*sqrt(1-delta^2)
 
 u = rgamma(n, shape = nu/2, rate = nu/2)
 t1 = rnorm(n)
 t2 = rnorm(n)
 
 mu + (Delta*abs(t1)+tau*t2)/sqrt(u)
}