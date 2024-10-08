# ==================================================================
# Trata o problema de label switching (Sylvia Fruhwirth-Schnatter)
# ===================================================================
# file: nome da pasta em Outputs/Paralelo onde estao salvos os resultados
# try: numero de tentativas maximas do kmeans

relabel = function(mu,sigma2,prob,z,gammaP,Gplus,G,try=10, test,saida = ""){
 caminho2 = paste0(paste0("Outputs",saida,"/test",test,"/corrigido/"))
 # estimando G+
 GplusHat = unique(sort(Gplus))[which.max(table(Gplus))]
 
 # filtrando amostras com G+ igual ao estimado
 index2 = which(Gplus==GplusHat)
 
 mu.samp = mu[index2,]
 sigma2.samp = sigma2[index2,]
 prob.samp = prob[index2,]
 gammaProb = gammaP[index2]
 z.samp = z[index2,]
 G.samp = G[index2]
 Gplus.samp = Gplus[index2]
 
 if(GplusHat==1){ # se temos apenas 1 cluster o processo nao e realizado
  prob_ok = prob.samp[,1]
  mu_ok = mu.samp[,1]
  sigma2_ok = sigma2.samp[,1]
  
  alpha_ok = alpha.samp
  gammaProb_ok = gammaProb
  G_ok = G.samp
  
  params = cbind(mu_ok,sigma2_ok,prob_ok,
                 gammaProb_ok,G_ok)
  
  write.table(params,file=paste0(caminho2,"1cluster.txt"), row.names = F)
  
  return("1 cluster")
 }
 
 # montando uma matriz de parametros
 nsamp = length(index2)
 
 muvec = c(as.matrix(mu.samp[,1:GplusHat]))
 sigma2vec = c(as.matrix(sigma2.samp[,1:GplusHat]))
 # ficam organizados por iteracao X cluster
 
 # matriz de parametros
 theta = cbind(muvec,sigma2vec)
 theta = scale(theta) # media 0 var 1
 
 # novos rotulos
 aux = kmeans(theta,centers=GplusHat,iter.max = 5000)
 labels = aux$cluster
 new_label = matrix(labels,nrow=nsamp,ncol=GplusHat)
 
 # verificando quais rotulos sao uma permutacao valida de {1,...,G+}
 idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:GplusHat)){a}))
 npermu = length(idpermu)
 taxa_permu = npermu/nsamp # das amostras com GplusHat comps, quantas tem permutacoes validas
 if(taxa_permu <= .5){
  cc = 0
  repeat{
   aux = kmeans(theta,centers=GplusHat,iter.max = 5000)
   labels = aux$cluster
   new_label = matrix(labels,nrow=nsamp,ncol=GplusHat)
   
   # verificando quais rotulos sao uma permutacao valida de {1,...,G+}
   idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:GplusHat)){a}))
   npermu = length(idpermu)
   taxa_permu = npermu/nsamp
   cc = cc+1
   if(taxa_permu > .5){break}
   if(taxa_permu <= .5 & cc == try){return(taxa_permu)}
  }
 }
 
 # ajustando os rotulos dos parametros em cada permutacao valida
 
 prob_ok = matrix(NA, nrow = npermu, ncol = GplusHat, dimnames = list(1:npermu,paste0("p_",1:GplusHat)))
 mu_ok = matrix(NA, nrow = npermu, ncol = GplusHat, dimnames = list(1:npermu,paste0("comp",1:GplusHat)))
 sigma2_ok = matrix(NA, nrow = npermu, ncol = GplusHat, dimnames = list(1:npermu,paste0("comp",1:GplusHat)))
 z_ok = matrix(NA, nrow = npermu, ncol = n, dimnames = list(1:npermu,1:n))
 gammaProb_ok = gammaProb[idpermu]
 G_ok = G.samp[idpermu]
 param = cbind(gammaProb_ok, G_ok); colnames(param) = c("gammaP", "G")
 
 for(i in 1:npermu){
  for(j in 1:GplusHat){
   prob_ok[i,new_label[idpermu[i],j]] = prob.samp[idpermu[i],j]
   mu_ok[i,new_label[idpermu[i],j]] = mu.samp[idpermu[i],j]
   sigma2_ok[i,new_label[idpermu[i],j]] = sigma2.samp[idpermu[i],j]
   for(l in 1:n){
    z_ok[i,l] = new_label[idpermu[i],z.samp[idpermu[i],l]]
   }
  }
 }
 
 # salvando
 write.table(mu_ok, file = paste0(caminho2,"mu.txt"),row.names = F)
 write.table(param, file = paste0(caminho2,"param.txt"),row.names = F)
 write.table(prob_ok, file = paste0(caminho2,"prob.txt"),row.names = F)
 write.table(sigma2_ok, file = paste0(caminho2,"sigma2.txt"),row.names = F)
 write.table(z_ok, file = paste0(caminho2,"z.txt"),row.names = F)
 
 return(taxa_permu)
}
