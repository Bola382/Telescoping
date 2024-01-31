# condicionais completas para cada parametro.

full_alpha = function(nu){
 Runuran::urgamma(1, shape = 2, scale = 1/(nu), lb = .02, ub = .5)
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

# ==============================================================================
#                                  Telescoping
# ==============================================================================

full_Z.TS = function(G,prob,b,n,y,X,beta,tau2,Delta,nu){
 aux_mu = X%*%beta+matrix(rep(b*Delta,n), nrow = n, ncol = G, byrow = T)
 aux_sig = sqrt(tau2+Delta^2)
 aux_lam = Delta/sqrt(tau2)
 out = NULL
 for(k in 1:n){
  aux_prob = exp(log(prob)+sapply(1:G,function(a) dst(y[k], aux_mu[k,a], aux_sig[a], aux_lam[a], nu,log=T))) # (eta1,eta2,eta3) x (st1,st2,st3)
  out[k] = sample(1:G,1, prob = aux_prob) # prob eh normalizado internamente
 }
 return(out)
}

full_beta.TS = function(j,b,c,n,p,Xj,Uj,yj,tj,M,tau2,Delta){
 if(M[j]==1){ # o R trata uma unica linha de uma matriz como um vetor coluna
  XU = as.matrix(Xj)%*%Uj # t(Xj)%*%Uj
 }else{
  XU = crossprod(Xj,Uj) # t(Xj)%*%Uj
 } 
 aux_y = yj-(b+tj)*Delta[j] 
 sigma_inv = XU%*%Xj + diag(tau2[j]/c^2,nrow = p)
 
 # ----------------------
 # Atualizando beta
 # ----------------------
 s_sigma = chol2inv(chol(sigma_inv)) 
 mmm = solve(sigma_inv,XU%*%aux_y)  
 ppp = tau2[j]*s_sigma      
 # novo valor pro vetor beta_j
 MASS::mvrnorm(1, mu = mmm,Sigma = ppp)
}

full_tau2.TS = function(j,b,r,s,Xj,uj,yj,tj,M,beta,Delta){
 beta = as.matrix(beta) # caso tenha so 1 grupo
 S3 = sum(uj*(yj-Xj%*%beta[,j]-(b+tj)*Delta[j])^2)
 
 1/rgamma(1, shape = M[j]/2+r, rate = (S3+s)/2)
}

full_Delta.TS = function(j,b,omega,Xj,uj,yj,tj,M,beta,tau2){
 beta=as.matrix(beta) # caso tenha so 1 grupo
 S1 = sum(uj*(yj-Xj%*%beta[,j])*(b+tj))
 S2 = sum(uj*(b+tj)^2)
 denom = omega^2*S2 + tau2[j]
 
 # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
 # alterar a media caso eta != 0.
 rnorm(1, omega^2*S1/denom, omega*sqrt(tau2[j]/denom))
}

full_U.TS = function(aux_mu,aux_yxbeta,n,tau2,Delta,nu,t,z){
 out = NULL
 for(k in 1:n){
  # atualizando U
  D1 =(aux_yxbeta[k,z[k]]-t[k]*Delta[z[k]])^2/(2*tau2[z[k]]) 
  
  out[k] = rgamma(1, shape = nu/2+1, rate = D1 + (t[k]^2+nu)/2)
 }
 return(out)
}

full_T.TS = function(aux_mu,aux_yxbeta,n,tau2,Delta,u,z){
 out = NULL
 for(k in 1:n){
  aux_denom = Delta[z[k]]^2 + tau2[z[k]]
  aux_mean = (aux_yxbeta[k,z[k]]*Delta[z[k]])/aux_denom
  aux_sd = sqrt(tau2[z[k]]/(u[k]*aux_denom))
  
  aux_truncnorm = rtruncnorm(1,-aux_mean/aux_sd)$amostra
  
  out[k] = aux_mean+aux_sd*aux_truncnorm
 }
 return(out)
}

full_K.TS = function(Gplus,Gmax,M,gammaProb, lpriori){
 options(digits=10)
 # funcao ainda tem um bug que nao consegui encontrar
 lpriori = if(any(lpriori == "unif")){rep(0,Gmax)}else{lpriori} # vetor com lprioris de K
 lprob_K = NULL
 Mmax = sapply(M, lgamma) # evitar Inf
 
 # log probabilidade
 
 for(G in Gplus:Gmax){
  gammaRatio = gammaProb/G
  
  # para calcular o somatorio
  aux = sapply(1:Gplus, function(a) lgamma(M[a] + gammaRatio) - lgamma(1 + gammaRatio))
  
  lprob_K[G-Gplus+1] = Gplus*log(gammaProb) - Gplus*log(G) + lfactorial(G) - 
   lfactorial(G-Gplus) + sum(aux) - sum(Mmax) # evitar Inf
 }
 
 prob_K = exp(lpriori[Gplus:Gmax]+lprob_K)
 if(is.infinite(sum(prob_K)) | is.na(sum(prob_K))){print(lprob_K);stop("Problema na condicional de G")}
 
 sample(Gplus:Gmax,1,prob = prob_K)
}


full_gammaProb.TS = function(gammaP,phigamma,n,M,G,Gplus){
 #phigamma ja deve estar especificado
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

full_prob.TS = function(z,G,gammaP){
 gtools::rdirichlet(1, alpha = gammaP/G + contagem(z,G))
}
