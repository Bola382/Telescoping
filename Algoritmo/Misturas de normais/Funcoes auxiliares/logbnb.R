# calcula o logaritmo da pmf de uma beta-binomial negativa em K-1
# K: valor da variavel aleatoria
# demais argumentos sao parametros
# por padrao valores recomendados pela Sylvia Fruhwirth-Schnatter para o telescoping

logbnb = function(K,alpha.lam=1,a.pi=4,b.pi=3){
 lgamma(alpha.lam + K - 1) + lbeta(alpha.lam + a.pi, K - 1 + b.pi) - 
  lgamma(alpha.lam) - lgamma(K) - lbeta(a.pi,b.pi)
}
