# ==================================================================
# Trata o problema de label switching (Sylvia Fruhwirth-Schnatter)
# as cadeias geradas nao estao disponiveis por conta de seu tamanho
# ===================================================================
# file: nome da pasta em Outputs/Paralelo onde estao salvos os resultados
# repl: qual replicacao escolher
# try: numero de tentativas maximas do kmeans

relabel = function(file,repl,try=10){
  caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/")
  caminho2 = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
  dir.create(caminho2,showWarnings = F)
  
  beta.samp = read.table(paste0(caminho,"beta.txt"), h =T)
  tau2.samp = read.table(paste0(caminho,"tau2.txt"), h =T)
  Delta.samp = read.table(paste0(caminho,"Delta.txt"), h =T)
  prob.samp = read.table(paste0(caminho,"prob.txt"), h =T)
  z.samp = read.table(paste0(caminho,"z.txt"), h =T)
  ut = read.table(paste0(caminho,"ut.txt"), h =T)
  umean.samp = ut$u.bar; tmean.samp = ut$t.bar; rm("ut")
  param = read.table(paste0(caminho,"param.txt"), h =T)
  nu.samp = param$nu; gammaProb = param$gammaP; alpha.samp = param$alpha
  Gplus.samp = param$Gplus; G.samp = param$G; rm("param")
  
  # organizando beta.samp em um array
  Q = nrow(beta.samp)
  Gmax = ncol(tau2.samp)
  n = ncol(z.samp)
  p = ncol(beta.samp)/Gmax # numero de betas (covs + 1)
  beta.aux = array(NA, dim = c(Q,p,Gmax), dimnames = list(1:Q,1:p,1:Gmax))
  beta.samp = as.matrix(beta.samp)
  enableJIT(3)
  for(i in 1:Q){
    qual = c(1:(p-1),0)
    for(j in 1:p){
      quais = which(1:ncol(beta.samp) %% p == qual[j])
      beta.aux[i,j,] = beta.samp[i,quais]
    }
  }
  beta.samp = beta.aux; rm("beta.aux") 
  
  # ==================================================================
  #                       Tratando label-switching
  # ==================================================================
  
  # estimando G+
  GplusHat = unique(sort(Gplus.samp))[which.max(table(Gplus.samp))]
  
  # filtrando amostras com G+ igual ao estimado
  index2 = which(Gplus.samp==GplusHat)
  
  beta.samp = beta.samp[index2,,]
  tau2.samp = tau2.samp[index2,]
  Delta.samp = Delta.samp[index2,]
  nu.samp = nu.samp[index2]
  prob.samp = prob.samp[index2,]
  gammaProb = gammaProb[index2]
  alpha.samp = alpha.samp[index2]
  umean.samp = umean.samp[index2]
  tmean.samp = tmean.samp[index2]
  z.samp = z.samp[index2,]
  G.samp = G.samp[index2]
  Gplus.samp = Gplus.samp[index2]
  
  if(GplusHat==1){ # se temos apenas 1 cluster o processo nao e realizado
    prob_ok = prob.samp[,1]
    beta_ok = beta.samp[,,1]
    tau2_ok = tau2.samp[,1]
    Delta_ok = Delta.samp[,1]
    t_ok = tmean.samp
    u_ok = umean.samp
    nu_ok = nu.samp
    
    alpha_ok = alpha.samp
    gammaProb_ok = gammaProb
    G_ok = G.samp
    
    params = cbind(beta_ok,tau2_ok,Delta_ok,nu_ok,prob_ok,
                   u_ok,t_ok,alpha_ok,gammaProb_ok,G_ok)
    colnames(params)[1:p] = paste0("beta",1:p)
    
    write.table(params,file=paste0(caminho2,"1cluster.txt"), row.names = F)
    
    return("1 cluster")
  }
  
  # montando uma matriz de parametros
  nsamp = length(index2)
  
  betamat = matrix(NA,nrow = nsamp*GplusHat, ncol = p) # cada coluna vai ser um dos coefs
  # misturando os coefs de comps diferentes, ou seja, coluna 1 beta1, qual comp? sim.
  for(j in 1:p){
    betamat[,j] = c(beta.samp[,j,1:GplusHat])
  }
  
  tau2vec = c(as.matrix(tau2.samp[,1:GplusHat]))
  Deltavec = c(as.matrix(Delta.samp[,1:GplusHat]))
  # ficam organizados por iteracao X cluster
  
  # matriz de parametros
  theta = cbind(betamat,tau2vec,Deltavec)
  
  # novos rotulos
  aux = kmeans(theta,centers=GplusHat,iter.max = 5000)
  labels = aux$cluster
  new_label = matrix(labels,nrow=nsamp,ncol=GplusHat)
  
  # verificando quais rotulos sao uma permutacao valida de {1,...,G+}
  idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:GplusHat)){a}))
  npermu = length(idpermu)
  if(npermu == 0){
    cc = 0
    repeat{
      aux = kmeans(theta,centers=GplusHat,iter.max = 5000)
      labels = aux$cluster
      new_label = matrix(labels,nrow=nsamp,ncol=GplusHat)
      
      # verificando quais rotulos sao uma permutacao valida de {1,...,G+}
      idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:GplusHat)){a}))
      npermu = length(idpermu)
      cc = cc+1
      if(npermu != 0){break}
      if(npermu == 0 & cc == try){return("exception")}
    }
  }
  
  taxa_permu = npermu/nsamp # das amostras com GplusHat comps, quantas tem permutacoes validas 
  # ajustando os rotulos dos parametros em cada permutacao valida
  
  prob_ok = matrix(NA, nrow = npermu, ncol = GplusHat, dimnames = list(1:npermu,paste0("p_",1:GplusHat)))
  beta_ok = array(NA, dim = c(npermu,p,GplusHat), 
                  dimnames = list(1:npermu,1:p,1:GplusHat))
  tau2_ok = matrix(NA, nrow = npermu, ncol = GplusHat, dimnames = list(1:npermu,paste0("comp",1:GplusHat)))
  Delta_ok = matrix(NA, nrow = npermu, ncol = GplusHat, dimnames = list(1:npermu,paste0("comp",1:GplusHat)))
  z_ok = matrix(NA, nrow = npermu, ncol = n, dimnames = list(1:npermu,1:n))
  t_ok = tmean.samp[idpermu]
  u_ok = umean.samp[idpermu]
  nu_ok = nu.samp[idpermu]
  
  alpha_ok = alpha.samp[idpermu]
  gammaProb_ok = gammaProb[idpermu]
  G_ok = G.samp[idpermu]
  param = cbind(nu_ok, gammaProb_ok, alpha_ok, G_ok); colnames(param) = c("nu", "gammaP", "alpha", "G")
  ut_ok = cbind(u_ok,t_ok); colnames(ut_ok) = c("u.bar","t.bar")
  
  for(i in 1:npermu){
    for(j in 1:GplusHat){
      prob_ok[i,new_label[idpermu[i],j]] = prob.samp[idpermu[i],j]
      tau2_ok[i,new_label[idpermu[i],j]] = tau2.samp[idpermu[i],j]
      Delta_ok[i,new_label[idpermu[i],j]] = Delta.samp[idpermu[i],j]
      for(k in 1:p){
        beta_ok[i,k,new_label[idpermu[i],j]] = beta.samp[idpermu[i],k,j]
      }
      for(l in 1:n){
        z_ok[i,l] = new_label[idpermu[i],z.samp[idpermu[i],l]]
      }
    }
  }
  
  # convertendo beta_ok em uma matriz npermu x (p*GplusHat)
  beta.aux = matrix(NA, nrow = npermu, ncol = p*GplusHat, dimnames = list(1:npermu,paste0(rep(paste0("beta",1:p),GplusHat),"_",rep(1:GplusHat,each=p))))
  
  for(ii in 1:npermu){
    beta.aux[ii,] = c(beta_ok[ii,,])
  }
  
  # salvando
  write.table(beta.aux, file = paste0(caminho2,"beta.txt"),row.names = F)
  write.table(Delta_ok, file = paste0(caminho2,"Delta.txt"),row.names = F)
  write.table(param, file = paste0(caminho2,"param.txt"),row.names = F)
  write.table(prob_ok, file = paste0(caminho2,"prob.txt"),row.names = F)
  write.table(tau2_ok, file = paste0(caminho2,"tau2.txt"),row.names = F)
  write.table(ut_ok, file = paste0(caminho2,"ut.txt"),row.names = F)
  write.table(z_ok, file = paste0(caminho2,"z.txt"),row.names = F)

  return(taxa_permu)
}
