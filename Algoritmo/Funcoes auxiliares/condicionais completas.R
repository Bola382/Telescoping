# condicionais completas para cada parametro.

full_alpha = function(nu){
 Runuran::urgamma(1, shape = 2, scale = 1/(nu), lb = .02, ub = .5)
}

full_alphaExpLoc = function(nu,loc){
  Runuran::urgamma(1, shape = 2, scale = 1/(nu-loc), lb = .02, ub = .5)
}

full_nu = function(prob,b,n,y,X,phi,beta,tau2,Delta,alpha,nu){ 
 prob = c(prob)
 G = length(prob)
 aux_mu = X%*%beta+matrix(rep(b*Delta,n), nrow = n, ncol = G, byrow = T)
 aux_sig = sqrt(tau2+Delta^2)
 aux_lam = Delta/sqrt(tau2)
 
 prop = rlnorm(1, log(nu), phi)
 
 aceit = min(1,exp(fullnu(prop,alpha,prob,y,aux_mu,aux_sig,aux_lam) + log(prop) - fullnu(nu,alpha,prob,y,aux_mu,aux_sig,aux_lam) - log(nu)))
 if(aceit > 1 | aceit <0){stop("Problema no passo MH")}
 if(runif(1)<=aceit){
  out = prop
 }else{
  out = nu
 }
 return(out)
}

# A INDICADORA DE NU > loc FOI OMITIDA, ESTOU ASSUMINDO QUE A PROPOSTA TEM SUPORTE 
# EM (loc, Inf)
full_nuExpLoc = function(prob,b,n,y,X,phi,beta,tau2,Delta,alpha,nu,loc){ 
  prob = c(prob)
  G = length(prob)
  aux_mu = X%*%beta+matrix(rep(b*Delta,n), nrow = n, ncol = G, byrow = T)
  aux_sig = sqrt(tau2+Delta^2)
  aux_lam = Delta/sqrt(tau2)
  
  prop = loc+rexp(1,rate = phi)
  
  aceit = min(1,exp(fullnuExpLoc(prop,alpha,prob,y,aux_mu,aux_sig,aux_lam,loc) - fullnuExpLoc(nu,alpha,prob,y,aux_mu,aux_sig,aux_lam,loc) + phi*(prop - nu)))
  if(aceit > 1 | aceit <0){stop("Problema no passo MH")}
  if(runif(1)<=aceit){
    out = prop
  }else{
    out = nu
  }
  return(out)
}

# ==============================================================================
#                                  Telescoping
# ==============================================================================

# G: numero de componentes
# prob: vetor de pesos das componentes
# b: centralizacao da media
# y: vetor de amostras
# X: matriz de planejamento
# beta: matriz de coefs de regressao por componente, dim = k x G
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
# nu: escalar contendo nu
full_Z.TS = function(G,prob,b,y,X,beta,tau2,Delta,nu){
 if(!abs(sum(prob) - 1) < 1e-5){stop("pesos das componentes no passo anterior nao somam 1")}
 n = length(y)
 aux_mu = X%*%beta+matrix(rep(b*Delta,n), nrow = n, ncol = G, byrow = T)
 aux_sig = sqrt(tau2+Delta^2)
 aux_lam = Delta/sqrt(tau2)
 out = NULL
 for(ind in 1:n){
  aux_prob = exp(log(prob)+sapply(1:G,function(a) dst(y[ind], aux_mu[ind,a], aux_sig[a], aux_lam[a], nu,log=T))) # (eta1,eta2,eta3) x (st1,st2,st3)
  if(any(is.na(aux_prob))){return("erro")}
  out[ind] = sample(1:G,1, prob = aux_prob) # prob eh normalizado internamente
 }
 return(out)
}

# j: indica componente
# b: centralizacao da media
# c: hiperparametro de beta
# Xj, Uj, yj, tj: X, U, y e t restritos aos z[i] == j
# M: vetor de tamanhos das componentes
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
full_beta.TS = function(j,b,c,Xj,Uj,yj,tj,M,tau2,Delta){
 if(M[j]==1){ # o R trata uma unica linha de uma matriz como um vetor coluna
  XU = as.matrix(Xj)%*%Uj # t(Xj)%*%Uj, dim = k x mj
 }else{
  XU = crossprod(Xj,Uj) # t(Xj)%*%Uj, dim = k x mj
 } 
 k = nrow(XU)
 aux_y = yj-(b+tj)*Delta[j] 
 sigma_inv = XU%*%Xj + diag(tau2[j]/c^2,nrow = k)
 
 # ----------------------
 # Atualizando beta
 # ----------------------
 s_sigma = chol2inv(chol(sigma_inv)) 
 mmm = solve(sigma_inv,XU%*%aux_y)  
 ppp = tau2[j]*s_sigma      
 # novo valor pro vetor beta_j
 MASS::mvrnorm(1, mu = mmm,Sigma = ppp)
}

# j: indica componente
# b: centralizacao da media
# r,s: hiperparametros de tau2
# Xj, Uj, yj, tj: X, U, y e t restritos aos z[i] == j
# M: vetor de tamanhos das componentes
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
full_tau2.TS = function(j,b,r,s,Xj,uj,yj,tj,M,beta,Delta){
 beta = as.matrix(beta) # caso tenha so 1 grupo
 S3 = sum(uj*(yj-Xj%*%beta[,j]-(b+tj)*Delta[j])^2)
 
 1/rgamma(1, shape = M[j]/2+r, rate = (S3/2+s))
}

# hiperparametro eta oculto aqui, lembrar de alterar caso eta != 0
# j: indica componente
# b: centralizacao da media
# omega: hiperparametro de Delta
# Xj, Uj, yj, tj: X, U, y e t restritos aos z[i] == j
# M: vetor de tamanhos das componentes
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
full_Delta.TS = function(j,b,omega,Xj,uj,yj,tj,M,beta,tau2){
 beta=as.matrix(beta) # caso tenha so 1 grupo
 S1 = sum(uj*(yj-Xj%*%beta[,j])*(b+tj))
 S2 = sum(uj*(b+tj)^2)
 denom = omega^2*S2 + tau2[j]
 
 # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
 # alterar a media caso eta != 0.
 rnorm(1, omega^2*S1/denom, omega*sqrt(tau2[j]/denom))
}

# aux_mu: media por individuo por componente, ou seja, X%*%beta + bDelta, dim = n x G
# aux_yxbeta: resposta - media por individuo por componente, dim  = n x G
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
# nu: escalar contendo nu
# z: vetor de alocacoes
# u: vetor da variavel latente U
# t: vetor da variavel latente T
full_TU.TS = function(aux_mu,aux_yxbeta,tau2,Delta, nu,z,u,t){
  n = nrow(aux_mu)
  out = matrix(NA, nrow = 2, ncol = n) # primeira linha t segunda linha u
  for(ind in 1:n){
    aux_denom = Delta[z[ind]]^2 + tau2[z[ind]]
    aux_mean = (aux_yxbeta[ind,z[ind]]*Delta[z[ind]])/aux_denom
    aux_sd = sqrt(tau2[z[ind]]/(u[ind]*aux_denom))
    
    aux_truncnorm = rtruncnorm(1,-aux_mean/aux_sd)$amostra
    
    out[1,ind] = aux_mean+aux_sd*aux_truncnorm
    D1 =(aux_yxbeta[ind,z[ind]]-out[1,ind]*Delta[z[ind]])^2/(2*tau2[z[ind]]) 
    
    out[2,ind] = rgamma(1, shape = nu/2+1, rate = (out[1,ind]^2+nu)/2 + D1)
  }
  return(out)
}

# aux_mu: media por individuo por componente, ou seja, X%*%beta + bDelta, dim = n x G
# aux_yxbeta: resposta - media por individuo por componente, dim  = n x G
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
# nu: escalar contendo nu
# t: vetor da variavel latente T
# z: vetor de alocacoes
full_U.TS = function(aux_mu,aux_yxbeta,tau2,Delta,nu,t,z){
 n = nrow(aux_mu)
 out = NULL
 for(ind in 1:n){
  # atualizando U
  D1 =(aux_yxbeta[ind,z[ind]]-t[ind]*Delta[z[ind]])^2/(2*tau2[z[ind]]) 
  
  out[ind] = rgamma(1, shape = nu/2+1, rate = (t[ind]^2+nu)/2 + D1)
 }
 return(out)
}

# aux_mu: media por individuo por componente, ou seja, X%*%beta + bDelta, dim = n x G
# aux_yxbeta: resposta - media por individuo por componente, dim  = n x G
# tau2: vetor com tau2 por componente
# Delta:--------- Delta --------------
# nu: escalar contendo nu
# u: vetor da variavel latente U
# z: vetor de alocacoes
full_T.TS = function(aux_mu,aux_yxbeta,tau2,Delta,u,z){
 n = nrow(aux_mu)
 out = NULL
 for(ind in 1:n){
  aux_denom = Delta[z[ind]]^2 + tau2[z[ind]]
  aux_mean = (aux_yxbeta[ind,z[ind]]*Delta[z[ind]])/aux_denom
  aux_sd = sqrt(tau2[z[ind]]/(u[ind]*aux_denom))
  
  aux_truncnorm = rtruncnorm(1,-aux_mean/aux_sd)$amostra
  
  out[ind] = aux_mean+aux_sd*aux_truncnorm
 }
 return(out)
}


# Gplus: numero de grupos
# Gmax: numero maximo de componentes
# M: vetor de tamanhos das componentes
# gammaProb: hiperparametro da priori dirichlet
# lpriori: logpriori de G
full_K.TS = function(Gplus,Gmax,M,gammaProb, lpriori){
 options(digits=10)
 lpriori = if(any(lpriori == "unif")){rep(0,Gmax)}else{lpriori} # vetor com lprioris de K
 lprob_K = NULL
 
 for(G in Gplus:Gmax){
  gammaRatio = gammaProb/G
  
  # para calcular o somatorio
  aux = sapply(1:Gplus, function(a) lgamma(M[a] + gammaRatio) - lgamma(1 + gammaRatio))
  
  lprob_K[G-Gplus+1] = Gplus*log(gammaProb) - Gplus*log(G) + lfactorial(G) - 
   lfactorial(G-Gplus) + sum(aux)
 }
 
 maior = max(abs(lprob_K))
 
 # evitar Inf
 prob_K = exp(lpriori[Gplus:Gmax]) * exp((lprob_K-maior)/maior)^maior 
 if(is.nan(sum(prob_K)) | is.infinite(sum(prob_K)) | is.na(sum(prob_K)) | sum(prob_K) == 0){return("erro")}
 
 sample(Gplus:Gmax,1,prob = prob_K)
}

# gammaP: valor anterior de gama
# phigamma: logvariancia da distribuicao proposta
# n: tamanho amostral
# M: vetor de tamanhos das componentes
# G: numero de componentes
# Gplus: numero de grupos
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

# z: vetor de alocacoes
# G: numero de componentes
# gammaP: hiperparametro da priori dirichlet
full_prob.TS = function(z,G,gammaP){
 gtools::rdirichlet(1, alpha = gammaP/G + contagem(z,G))
}
