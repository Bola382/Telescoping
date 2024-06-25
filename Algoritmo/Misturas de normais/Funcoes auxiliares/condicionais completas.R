full_prob.TS = function(z,G,gammaP){
 gtools::rdirichlet(1, alpha = gammaP/G + contagem(z,G))
}

full_K.TS = function(Gplus,Gmax,M,gammaProb, lpriori){
 options(digits=10)
 lpriori = if(any(lpriori == "unif")){rep(0,Gmax)}else{lpriori} # vetor com lprioris de K
 lprob_K = NULL
 Mmax = sapply(M+gammaProb/Gplus, lgamma) # evitar Inf
 # log probabilidade
 
 for(G in Gplus:Gmax){
  gammaRatio = gammaProb/G
  
  # para calcular o somatorio
  aux = sapply(1:Gplus, function(a) lgamma(M[a] + gammaRatio) - lgamma(1 + gammaRatio))
  
  lprob_K[G-Gplus+1] = Gplus*log(gammaProb) - Gplus*log(G) + lfactorial(G) - 
   lfactorial(G-Gplus) + sum(aux)
 }
 lprob_K = lprob_K - sum(Mmax)
 
 maior = max(abs(lprob_K))
 
 # evitar Inf
 prob_K = exp(lpriori[Gplus:Gmax]) * exp((lprob_K-max(lprob_K))/maior)^maior 
 if(is.nan(sum(prob_K)) | is.infinite(sum(prob_K)) | is.na(sum(prob_K)) | sum(prob_K) == 0){return("erro")}
 
 sample(Gplus:Gmax,1,prob = prob_K)
}


full_gammaProb.TS = function(gammaP,phigamma,n,M,G,Gplus){
 prop = rlnorm(1, log(gammaP), phigamma)
 
 aceit = min(1,aceitGamma(gammaP,prop,n,M,G,Gplus))
 if(aceit > 1 | aceit <0){stop("Problema no passo MH")}
 if(runif(1)<=aceit){
  out = prop
 }else{
  out = gammaP
 }
 return(out)
}

full_Z.TS = function(G,prob,mu,sigma2,n,y){
 out = NULL
 for(k in 1:n){
  aux_prob = exp(log(prob)+sapply(1:G,function(a) dnorm(y[k], mu[a], sqrt(sigma2[a]),log=T))) 
  if(any(is.na(aux_prob))){stop("erro na condicional de Z")}
  out[k] = sample(1:G,1, prob = aux_prob) # prob eh normalizado internamente
 }
 return(out)
}

full_theta.TS = function(j,M,S1,S2, eta, omega, a, b){
 param1 = (omega*eta+S1)/(omega+M[j])
 param2 = omega+M[j]
 param3 = M[j]/2 + a
 param4 = b + S2/2 + omega*eta^2/2 - (omega*eta+S1)^2/(2*(omega+M[j]))
 
 sig2 = 1/rgamma(1, shape = param3, rate = param4)
 med = rnorm(1, mean = param1, sd = sqrt(sig2/param2))
 
 return(c(med,sig2))
}