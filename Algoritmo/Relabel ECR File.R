# =====================================================================
# Trata o problema de label switching com o ECR
# =====================================================================
# estima o numero de grupos pelo MAP e dado isso aplica o algoritmo ECR 
# Panagiotis Papastamoulis & George Iliopoulos (2010)
# precisa carregar logpost.R e do pacote label.switching

# file: nome da pasta em Outputs/Paralelo onde estao salvos os resultados
# repl: qual replicacao escolher
# y: vetor de amostras
# X: matriz de planejamento, primeira coluna de uns, demais colunas covariaveis

relabelECR = function(file,repl,y,X){
  caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/")
  caminho2 = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
  dir.create(caminho2,showWarnings = F)
  
  beta.samp = as.matrix(read.table(paste0(caminho,"beta.txt"), h =T))
  tau2.samp = read.table(paste0(caminho,"tau2.txt"), h =T)
  Delta.samp = read.table(paste0(caminho,"Delta.txt"), h =T)
  prob.samp = read.table(paste0(caminho,"prob.txt"), h =T)
  z.samp = read.table(paste0(caminho,"z.txt"), h =T)
  param = read.table(paste0(caminho,"param.txt"), h =T)
  nu.samp = param$nu; gammaProb = param$gammaP; alpha.samp = param$alpha
  Gplus.samp = param$Gplus; G.samp = param$G; rm("param")
  ut = read.table(paste0(caminho,"ut.txt"), h =T)
  umean.samp = ut$u.bar; tmean.samp = ut$t.bar; rm("ut")
  
  # estimando G+
  GplusHat = unique(sort(Gplus.samp))[which.max(table(Gplus.samp))]
  
  Q = sum(Gplus.samp==GplusHat) # amostras iguais ao MAP
  index = which(Gplus.samp==GplusHat) # indices iguais ao MAP
  n = length(y) # tamanho amostral
  k = ncol(X) # numero de betas
  
  # filtrando amostras com G+ igual ao estimado
  beta.samp = beta.samp[index,]
  tau2.samp = tau2.samp[index,]
  Delta.samp = Delta.samp[index,]
  nu.samp = nu.samp[index]
  prob.samp = prob.samp[index,]
  gammaProb = gammaProb[index]
  alpha.samp = alpha.samp[index]
  umean.samp = umean.samp[index]
  tmean.samp = tmean.samp[index]
  z.samp = z.samp[index,]
  G.samp = G.samp[index]
  Gplus.samp = Gplus.samp[index]
  
  if(GplusHat==1){ # se temos apenas 1 cluster o processo nao e realizado
    prob_ok = prob.samp[,1]
    beta_ok = beta.samp[,1:k]
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
    colnames(params) = c(paste0("beta",1:k), 
                         "tau2", "Delta", "nu", "prob", "umean", "tmean",
                         "alpha", "gama", "G")
    
    write.table(params,file=paste0(caminho2,"1cluster.txt"), row.names = F)
    
    return("1 cluster")
  }
  
  loglikMCMC = NULL
  enableJIT(3)
  for(i in 1:Q){
    beta_aux = matrix(beta.samp[i,1:(k*GplusHat)], k, GplusHat)
    loglikMCMC[i] = loglik(y, z.samp[i,], X, prob.samp[i, 1:GplusHat], 
                     beta_aux, tau2.samp[i, 1:GplusHat], Delta.samp[i, 1:GplusHat],
                     nu.samp[i]) # note que loglik = 0 em grupos vazios
  }
  
  # determinando EMV dado G+
  idEMV = which.max(loglikMCMC)
  zEMV = z.samp[idEMV,]
  
  # resolvendo label switching
  L = label.switching::label.switching(method = "ECR", zpivot = as.numeric(zEMV), z = as.matrix(z.samp), K = GplusHat)
  new_label = L$permutations[[1]]
  
  # para guardar valores pos correcao
  prob_ok = matrix(NA, nrow = Q, ncol = GplusHat, dimnames = list(1:Q,paste0("p_",1:GplusHat)))
  beta_ok = matrix(NA, nrow = Q, ncol = k*GplusHat, dimnames = list(1:Q,paste0(rep(paste0("beta",1:k),GplusHat),"_",rep(1:GplusHat,each=k))))
  tau2_ok = matrix(NA, nrow = Q, ncol = GplusHat, dimnames = list(1:Q,paste0("tau2_",1:GplusHat)))
  Delta_ok = matrix(NA, nrow = Q, ncol = GplusHat, dimnames = list(1:Q,paste0("Delta_",1:GplusHat)))
  z_ok = matrix(NA, nrow = Q, ncol = n, dimnames = list(1:Q,1:n))
  nu_ok = nu.samp
  
  alpha_ok = alpha.samp
  gammaProb_ok = gammaProb
  G_ok = G.samp
  param = cbind(nu_ok, gammaProb_ok, alpha_ok, G_ok); colnames(param) = c("nu", "gammaP", "alpha", "G")
  ut_ok = cbind(umean.samp,tmean.samp); colnames(ut_ok) = c("u.bar","t.bar")
  
  idbeta = rep(1:GplusHat, each = k)
  
  for(q in 1:Q){
    for(j in 1:GplusHat){
      prob_ok[q,new_label[q,j]] = prob.samp[q,j]
      tau2_ok[q,new_label[q,j]] = tau2.samp[q,j]
      Delta_ok[q,new_label[q,j]] = Delta.samp[q,j]
      
      newidbeta = which(idbeta == new_label[q,j])
      oldidbeta = which(idbeta == j)
      beta_ok[q,newidbeta] = beta.samp[q,oldidbeta]
    }
    for(i in 1:n){
      z_ok[q,i] = new_label[q,z.samp[q,i]]
    }
  }
  
  # salvando
  write.table(beta_ok, file = paste0(caminho2,"beta.txt"),row.names = F)
  write.table(Delta_ok, file = paste0(caminho2,"Delta.txt"),row.names = F)
  write.table(param, file = paste0(caminho2,"param.txt"),row.names = F)
  write.table(prob_ok, file = paste0(caminho2,"prob.txt"),row.names = F)
  write.table(tau2_ok, file = paste0(caminho2,"tau2.txt"),row.names = F)
  write.table(ut_ok, file = paste0(caminho2,"ut.txt"),row.names = F)
  write.table(z_ok, file = paste0(caminho2,"z.txt"),row.names = F)
}


