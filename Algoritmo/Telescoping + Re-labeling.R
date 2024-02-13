# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping e solucionando label-switching
# ==============================================================================
# necessita carregar a pasta "Funcoes auxiliares" e a funcao "dst.R"
# e do pacote "compiler"

# valores iniciais exceto G e G_+ gerados das respectivas prioris

# data: matriz de dados, primeira coluna respostas, demais valores das covs
# ncomps: valor inicial de G
# nclusters: valor inicial de G_+
# G.max: numero maximo de componentes
# lprio_G: vetor com a log probabilidade para cada G entre 1 e G.max
# phi: "desvio padrao" da proposta de MH para nu
# phigamma: "desvio padrao" da proposta de MH para gammaProb
# Q: numero de iteracoes do algoritmo
# burn: numero de amostras descartadas
# thin: tamanho dos saltos apos o burn



telescope = function(data,ncomps,nclusters,G.max,lprio_G,phi=.3,phigamma=2.5,Q,burn,thin){
 library(compiler)
 n = nrow(data)
 p = ncol(data) # numero de betas
 
 y = data[,1] # resp
 X = as.matrix(cbind(1,data[,-1])) # covariaveis com intercepto
 
 # ~~~~~~~~~~~~~~~~
 # valores iniciais
 # ~~~~~~~~~~~~~~~~
 
 G.samp = ncomps
 Gplus.samp = nclusters # numero de componentes e clusters
 
 # para guardar amostras
 prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = Q, ncol = G.max)
 u.samp = t.samp = z.samp = matrix(NA, nrow = Q, ncol = n)
 beta.samp = array(NA, dim = c(Q,p,G.max), dimnames = list(1:Q, 1:p, 1:G.max))
 
 # hiperparametros
 gammaProb = rf(1,6,3)
 xi = rep(gammaProb/G.samp,G.samp) # prob
 c = 10 # da beta (c <- sqrt(c) das minhas contas)
 eta = 0; omega = 10 # da Delta
 r = s = .1 # da tau2
 
 prob.samp[1,1:G.samp] = gtools::rdirichlet(1, alpha = xi) # pesos
 beta.samp[1,,1:G.samp] = matrix(rnorm(p*G.samp,sd=c),ncol=G.samp) # coef reg
 tau2.samp[1,1:G.samp] = 1/rgamma(G.samp,r,s)# escala
 Delta.samp[1,1:G.samp] = rnorm(G.samp,eta,omega)# forma
 alpha.samp = runif(1,.02,.5)# hiperparametro da priori de nu
 nu.samp = rexp(1,rate=alpha.samp)# gl
 
 u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
 t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
 z.samp[1,] = sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T) # latente que indica os grupos
 
 # -------------------------------------------------------------------------------
 #                                    Amostrador
 # -------------------------------------------------------------------------------
 
 cont = 0 # contador de aceites de MH para nu
 contgamma = 0 # contador de aceites de MH para gammaProb
 
 
 enableJIT(3)
 t_tmp = Sys.time()
 for(i in 2:Q){
  # centralizacao da media
  b = -sqrt(nu.samp[i-1]/pi)*(gamma((nu.samp[i-1]-1)/2)/gamma(nu.samp[i-1]/2))
  
  # ===========================================================================
  #                                  Passo 1
  # ===========================================================================
  
  # ----------------------
  # a) Atualizando Z
  # ----------------------
  z.samp[i,] = full_Z.TS(G.samp[i-1],prob.samp[i-1,1:G.samp[i-1]],
                         b,n,y,X,
                         beta.samp[i-1,,1:G.samp[i-1]],
                         tau2.samp[i-1,1:G.samp[i-1]],
                         Delta.samp[i-1,1:G.samp[i-1]],
                         nu.samp[i-1])
  
  # ------------------------
  # b) determinando Gplus e 
  # M: tamanho de cada comp
  # ------------------------
  
  # re-rotulando
  z.samp[i,] = rotulador(z.samp[i,])
  
  # contando
  M = contagem(z.samp[i,],G.samp[i-1])
  
  # atualizando Gplus
  Gplus.samp[i] = sum(M!=0)
  
  # removendo custers vazios da contagem
  M = M[1:Gplus.samp[i]]
  
  # ===========================================================================
  #                                  Passo 2
  # ===========================================================================
  
  # ----------------------------
  # a) atualizando parametros 
  # das componentes preenchidas
  # ----------------------------
  for(j in 1:Gplus.samp[i]){
   index = which(z.samp[i,]==j)
   Xj = X[index,]
   uj = u.samp[i-1,index]
   Uj = diag(uj, nrow = M[j], ncol = M[j]) 
   yj = y[index] 
   tj = t.samp[i-1,index]
   
   # Atualizando beta
   beta.samp[i,,j] = full_beta.TS(j,b,c,n,p,Xj,Uj,yj,tj,M,
                                  tau2.samp[i-1,1:G.samp[i-1]],
                                  Delta.samp[i-1,1:G.samp[i-1]])
   
   # Atualizando tau2
   tau2.samp[i,j] = full_tau2.TS(j,b,r,s,Xj,uj,yj,tj,M,
                                 beta.samp[i,,1:G.samp[i-1]],
                                 Delta.samp[i-1,1:G.samp[i-1]])
   
   # Atualizando Delta
   Delta.samp[i,j] = full_Delta.TS(j,b,omega,Xj,uj,yj,tj,M,
                                   beta.samp[i,,1:G.samp[i-1]],
                                   tau2.samp[i,1:G.samp[i-1]])
   # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
   # alterar a media na funcao full_delta.TS caso eta != 0.
  }
  
  # Atualizando U e T
  aux_mu = X%*%beta.samp[i,,1:Gplus.samp[i]]+
   matrix(rep(b*Delta.samp[i,1:Gplus.samp[i]],n),
          nrow = n, ncol = Gplus.samp[i], byrow = T)
  aux_yxbeta = y-aux_mu
  
  # atualizando U
  u.samp[i,] = full_U.TS(aux_mu,aux_yxbeta,n,
                         tau2.samp[i,1:Gplus.samp[i]],
                         Delta.samp[i,1:Gplus.samp[i]],
                         nu.samp[i-1], t.samp[i-1,], z.samp[i,])
  
  # atualizando T
  t.samp[i,] = full_T.TS(aux_mu,aux_yxbeta,n,
                         tau2.samp[i,1:Gplus.samp[i]],
                         Delta.samp[i,1:Gplus.samp[i]],
                         u.samp[i,], z.samp[i,])
  
  
  # -------------------------------
  # b) atualizando hiperparametros
  # -------------------------------
  alpha.samp[i] = full_alpha(nu.samp[i-1])
  
  # ===========================================================================
  #                                  Passo 3
  # ===========================================================================
  
  # ----------------------
  # a) atualizando G
  # ----------------------
  G.samp[i] = full_K.TS(Gplus.samp[i],G.max,M,gammaProb[i-1],lprio_G)
  
  # ----------------------
  # b) atualizando gammaP
  # ----------------------
  gammaProb[i] = full_gammaProb.TS(gammaProb[i-1],phigamma,n,M,G.samp[i],Gplus.samp[i])
  contgamma = ifelse(gammaProb[i]==gammaProb[i-1],contgamma,contgamma + 1)
  
  
  # ===========================================================================
  #                                  Passo 4
  # ===========================================================================
  
  # ---------------------------
  # a) adicionando componentes 
  # vazios
  # ---------------------------
  dif = G.samp[i]-Gplus.samp[i]
  if(dif > 0){
   M = c(M,rep(0,dif))
   
   mtyindex = (Gplus.samp[i]+1):G.samp[i]
   
   # Atualizando beta
   beta.samp[i,,mtyindex] = matrix(rnorm(p*dif,sd=c),ncol=dif)
   
   # Atualizando tau2
   tau2.samp[i,mtyindex] = 1/rgamma(dif,r,s)
   
   # Atualizando Delta
   Delta.samp[i,mtyindex] = rnorm(dif,eta,omega)
  }
  
  # ----------------------
  # Atualizando prob
  # ----------------------
  compindex = 1:G.samp[i]
  prob.samp[i,compindex] = full_prob.TS(z.samp[i,],G.samp[i],gammaProb[i])
  
  # ----------------------
  # Atualizando nu
  # ----------------------
  nu.samp[i] = full_nu(prob.samp[i,compindex],b,n,y,X,phi,beta.samp[i,,compindex],
                       tau2.samp[i,compindex],Delta.samp[i,compindex],
                       alpha.samp[i],nu.samp[i-1])
  
  cont = ifelse(nu.samp[i]==nu.samp[i-1],cont,cont + 1)
  
 };time_ok = Sys.time()-t_tmp
 
 # ==================================================================
 # Trata o problema de label switching (Sylvia Fruhwirth-Schnatter)
 # as cadeias geradas nao estao disponiveis por conta de seu tamanho
 # ===================================================================
 
 # aplicando burn-in e thin
 index = seq(floor(burn)+1,Q,by=thin)
 
 beta.samp = beta.samp[index,,]
 tau2.samp = tau2.samp[index,]
 Delta.samp = Delta.samp[index,]
 nu.samp = nu.samp[index]
 prob.samp = prob.samp[index,]
 gammaProb = gammaProb[index]
 alpha.samp = alpha.samp[index]
 u.samp = u.samp[index,]
 t.samp = t.samp[index,]
 z.samp = z.samp[index,]
 G.samp = G.samp[index]
 Gplus.samp = Gplus.samp[index]
 
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
 u.samp = u.samp[index2,]
 t.samp = t.samp[index2,]
 z.samp = z.samp[index2,]
 G.samp = G.samp[index2]
 Gplus.samp = Gplus.samp[index2]
 
 if(GplusHat==1){ # se temos apenas 1 cluster o processo nao e realizado
  prob_ok = prob.samp[,1]
  beta_ok = beta.samp[,,1]
  tau2_ok = tau2.samp[,1]
  Delta_ok = Delta.samp[,1]
  z_ok = "todos no cluster 1"
  t_ok = t.samp
  u_ok = u.samp
  nu_ok = nu.samp
  
  alpha_ok = alpha.samp
  gammaProb_ok = gammaProb
  G_ok = G.samp
  
  out = list(time = time_ok, taxa_nu = cont/Q, GplusHat = GplusHat,
             taxa_gama = contgamma/Q, taxa_permu = taxa_permu,
             beta = beta_ok, tau2 = tau2_ok, Delta = Delta_ok,
             nu = nu_ok, prob = prob_ok, gammaProb = gammaProb_ok, alpha = alpha_ok,
             z = z_ok, t = t_ok, u = u_ok, G = G_ok)
  return(out)
 }
 
 # montando uma matriz de parametros
 nsamp = length(index2)
 
 betamat = matrix(NA,nrow = nsamp*GplusHat, ncol = ncol(X)) # cada coluna vai ser um dos coefs
 # misturando os coefs de comps diferentes, ou seja, coluna 1 beta1, qual comp? sim.
 for(j in 1:ncol(X)){
  betamat[,j] = c(beta.samp[,j,1:GplusHat])
 }
 
 tau2vec = c(tau2.samp[,1:GplusHat])
 Deltavec = c(Delta.samp[,1:GplusHat])
 # ficam organizados por iteracao X cluster
 
 # matriz de parametros
 theta = cbind(betamat,tau2vec,Deltavec)
 
 # novos rotulos
 aux = kmeans(theta,centers=GplusHat,iter.max = 300)
 labels = aux$cluster
 new_label = matrix(labels,nrow=nsamp,ncol=GplusHat)
 
 # verificando quais rotulos sao uma permutacao valida de {1,...,G+}
 idpermu = unlist(sapply(1:nsamp, function(a) if(identical(sort(new_label[a,]),1:GplusHat)){a}))
 npermu = length(idpermu)
 if(npermu == 0){return("exception")}
 
 taxa_permu = npermu/nsamp # das amostras com GplusHat comps, quantas tem permutacoes validas 
 # ajustando os rotulos dos parametros em cada permutacao valida
 
 prob_ok = matrix(NA, nrow = npermu, ncol = GplusHat)
 beta_ok = array(NA, dim = c(npermu,ncol(X),GplusHat), 
                 dimnames = list(1:npermu,1:ncol(X),1:GplusHat))
 tau2_ok = matrix(NA, nrow = npermu, ncol = GplusHat)
 Delta_ok = matrix(NA, nrow = npermu, ncol = GplusHat)
 z_ok = matrix(NA, nrow = npermu, ncol = n)
 t_ok = t.samp[idpermu,]
 u_ok = u.samp[idpermu,]
 nu_ok = nu.samp[idpermu]
 
 alpha_ok = alpha.samp[idpermu]
 gammaProb_ok = gammaProb[idpermu]
 G_ok = G.samp[idpermu]
 
 for(i in 1:npermu){
  for(j in 1:GplusHat){
   prob_ok[i,new_label[idpermu[i],j]] = prob.samp[idpermu[i],j]
   tau2_ok[i,new_label[idpermu[i],j]] = tau2.samp[idpermu[i],j]
   Delta_ok[i,new_label[idpermu[i],j]] = Delta.samp[idpermu[i],j]
   for(k in 1:ncol(X)){
    beta_ok[i,k,new_label[idpermu[i],j]] = beta.samp[idpermu[i],k,j]
   }
   for(l in 1:n){
    z_ok[i,l] = new_label[idpermu[i],z.samp[idpermu[i],l]]
   }
  }
 }
 # iteracoes MCMC, (betas,tau2,Delta,nu,prob), componente
 paramComp = array(NA, dim = c(npermu, ncol(X) + 3, GplusHat), 
               dimnames = list(1:npermu,
                               c(paste0("beta",1:ncol(X)),"tau2","Delta","prob"),
                               1:GplusHat))
 paramComp[1:npermu,1:ncol(X),] = beta_ok
 paramComp[1:npermu,1:(ncol(X)+1),] = tau2_ok
 paramComp[1:npermu,1:(ncol(X)+2),] = Delta_ok
 paramComp[1:npermu,1:(ncol(X)+3),] = prob_ok
 
 # iteracoes MCMC, (nu,gammaProb,alpha,G)
 param = array(NA, dim = c(npermu,4), 
               dimnames = list(1:npermu,c("nu","gammaProb","alpha","G")))
 param[1:npermu,1] = nu_ok
 param[1:npermu,2] = gammaProb_ok
 param[1:npermu,3] = alpha_ok
 param[1:npermu,4] = G_ok
 
 # iteracoes MCMC, qual latente, individuos
 latent = array(NA, dim = c(npermu,3,n), 
                dimnames = list(1:npermu,
                                c("z","u","t"),
                                1:n))
 latent[1:npermu,1,] = z_ok
 latent[1:npermu,2,] = u_ok
 latent[1:npermu,3,] = t_ok
 
 out = list(time = time_ok, taxa_nu = cont/Q,taxa_gama = contgamma/Q, 
            taxa_permu = taxa_permu, paramComp = paramComp, param = param,
            latent = latent)
 return(out)
}
