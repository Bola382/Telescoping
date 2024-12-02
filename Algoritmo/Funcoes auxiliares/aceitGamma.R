# Calcula a probabilidade de aceitacao do passo de MH para gammaProb
# considerando priori F(6,3) em gammaProb
# gammaP: valor no passo anterior
# gammaP.new: valor proposto]
# n: tamanho amostral
# M: tamanho de cada componente
# G: numero de componentes
# Gplus: numero de clusters

aceitGamma = function(gammaP,gammaP.new,n,M,G,Gplus){
 # termo dentro do somatorio
 aux = sapply(1:Gplus, function(a) lgamma(M[a] + gammaP.new/G) - lgamma(M[a] + gammaP/G) -
               lgamma(1+gammaP.new/G) + lgamma(1+gammaP/G))
 
 # probabilidade de aceitacao
 lprob = (Gplus+1)*log(gammaP.new/gammaP) - 4.5*(log(1+2*gammaP.new)-log(1+2*gammaP)) +
  lgamma(gammaP.new)-lgamma(gammaP) + lgamma(n+gammaP) - lgamma(n+gammaP.new) + sum(aux)
 
 exp(lprob)
}
