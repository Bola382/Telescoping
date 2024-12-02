# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# salva resultados em um .txt
# ==============================================================================
# necessita carregar a pasta "Funcoes auxiliares" e a funcao "dst.R"
# e do pacote "compiler"

# valores iniciais exceto G e G_+ gerados das respectivas prioris

# folder: nome da pasta em qual serao salvos os resultados de cada replicacao
# repl: numero da replicacao atual
# data: matriz de dados, primeira coluna respostas, demais valores das covs
# ncomps: valor inicial de G
# nclusters: valor inicial de G_+
# G.max: numero maximo de componentes
# lprio_G: vetor com a log probabilidade para cada G entre 1 e G.max
# phi: "desvio padrao" da proposta de MH para nu
# phigamma: "desvio padrao" da proposta de MH para gammaProb
# Q: numero de iteracoes do algoritmo
# burn: amostras a serem descartadas
# thin: tamanho dos saltos pos burnin

telescope = function(folder,repl,dados,ncomps,nclusters,G.max,lprio_G,phi=.3,phigamma=2.5,Q,burn,thin){
 n = nrow(dados)
 p = ncol(dados) # numero de betas
 
 y = dados[,1] # resp
 X = as.matrix(cbind(1,dados[,-1])) # covariaveis com intercepto
 
 # ~~~~~~~~~~~~~~~~
 # valores iniciais
 # ~~~~~~~~~~~~~~~~
 
 G.samp = ncomps
 Gplus.samp = nclusters # numero de componentes e clusters
 
 # para guardar amostras
 prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = 2, ncol = G.max)
 u.samp = t.samp = z.samp = matrix(NA, nrow = 2, ncol = n)
 beta.samp = array(NA, dim = c(2,p,G.max), dimnames = list(c("old","new"), 1:p, 1:G.max))
 
 # hiperparametros
 gammaProb = rf(1,6,3)
 xi = rep(gammaProb/G.samp,G.samp) # prob
 c = 10 # da beta (c <- sqrt(c) das minhas contas)
 eta = 0; omega = 10 # da Delta
 r = 2.1
 s = 1.1 # da tau2
 
 prob.samp[1,1:G.samp] = gtools::rdirichlet(1, alpha = xi) # pesos
 beta.samp[1,,1:G.samp] = matrix(rnorm(p*G.samp,sd=c),ncol=G.samp) # coef reg
 tau2.samp[1,1:G.samp] = 1/rgamma(G.samp,r,s)# escala
 Delta.samp[1,1:G.samp] = rnorm(G.samp,eta,omega)# forma
 alpha.samp = runif(1,.02,.5)# hiperparametro da priori de nu
 nu.samp = rexp(1,rate=alpha.samp)# gl
 
 # centralizacao da media
 b = -sqrt(nu.samp/pi)*(gamma((nu.samp-1)/2)/gamma(nu.samp/2))
 
 u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
 t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
 z.samp[1,] = sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T) # latente que indica os grupos
 
 # ~~~~~~~~~~~~~~~~~~~
 # Config output
 # ~~~~~~~~~~~~~~~~~~~
 cat(paste0(rep(paste0("beta",1:p),G.max),"_",rep(1:G.max,each=p)),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/beta.txt"))
 cat(paste0("tau2_",1:G.max),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2.txt"))
 cat(paste0("Delta_",1:G.max),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Delta.txt"))
 cat(paste0("p_",1:G.max),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/prob.txt"))
 cat("nu","gammaP","alpha","Gplus","G","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/param.txt"))
 cat(paste0("z_",1:n),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/z.txt"))
 cat("u.bar","t.bar","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/ut.txt"))
 
 # -------------------------------------------------------------------------------
 #                                    Amostrador
 # -------------------------------------------------------------------------------
 
 cont = 0 # contador de aceites de MH para nu
 contgamma = 0 # contador de aceites de MH para gammaProb
 
 enableJIT(3)
 t_tmp = Sys.time()
 for(i in 2:Q){
  # ===========================================================================
  #                                  Passo 1
  # ===========================================================================
  
  # ----------------------
  # a) Atualizando Z
  # ----------------------
  z.samp[2,] = full_Z.TS(G.samp[1],prob.samp[1,1:G.samp[1]],
                         b,y,X,
                         beta.samp[1,,1:G.samp[1]],
                         tau2.samp[1,1:G.samp[1]],
                         Delta.samp[1,1:G.samp[1]],
                         nu.samp[1])
  if(z.samp[2,1]=="erro"){
   cat(c(beta.samp[1,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/betaerro.txt"))
   cat(tau2.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2erro.txt"))
   cat(Delta.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Deltaerro.txt"))
   cat(prob.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/proberro.txt"))
   cat(nu.samp[1],gammaProb[1],alpha.samp[1],Gplus.samp[1],G.samp[1],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/paramerro.txt"))
   cat(z.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/zerro.txt"))
   cat(u.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/uerro.txt"))
   cat(t.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/terro.txt"))
   
   cat(c(beta.samp[2,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/betaerro.txt"), append=T)
   cat(tau2.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2erro.txt"), append=T)
   cat(Delta.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Deltaerro.txt"), append=T)
   cat(prob.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/proberro.txt"), append=T)
   cat(nu.samp[2],gammaProb[2],alpha.samp[2],Gplus.samp[2],G.samp[2],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/paramerro.txt"), append=T)
   cat(z.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/zerro.txt"), append=T)
   cat(u.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/uerro.txt"), append=T)
   cat(t.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/terro.txt"), append=T)
   stop("erro em z")
  }
  
  # ------------------------
  # b) determinando Gplus e 
  # M: tamanho de cada comp
  # ------------------------
  
  # re-rotulando
  z.samp[2,] = rotulador(z.samp[2,])
  
  # contando
  M = contagem(z.samp[2,],G.samp[1])
  
  # atualizando Gplus
  Gplus.samp[2] = sum(M!=0)
  
  # removendo custers vazios da contagem
  M = M[1:Gplus.samp[2]]
  
  # ===========================================================================
  #                                  Passo 2
  # ===========================================================================
  
  # ----------------------------
  # a) atualizando parametros 
  # das componentes preenchidas
  # ----------------------------
  for(j in 1:Gplus.samp[2]){
   index = which(z.samp[2,]==j)
   Xj = X[index,]
   uj = u.samp[1,index]
   Uj = diag(uj, nrow = M[j], ncol = M[j]) 
   yj = y[index] 
   tj = t.samp[1,index]
   
   # Atualizando beta
   beta.samp[2,,j] = full_beta.TS(j,b,c,Xj,Uj,yj,tj,M,
                                  tau2.samp[1,1:G.samp[1]],
                                  Delta.samp[1,1:G.samp[1]])
   
   # Atualizando tau2
   tau2.samp[2,j] = full_tau2.TS(j,b,r,s,Xj,uj,yj,tj,M,
                                 beta.samp[2,,1:G.samp[1]],
                                 Delta.samp[1,1:G.samp[1]])
   
   # Atualizando Delta
   Delta.samp[2,j] = full_Delta.TS(j,b,omega,Xj,uj,yj,tj,M,
                                   beta.samp[2,,1:G.samp[1]],
                                   tau2.samp[2,1:G.samp[1]])
   # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
   # alterar a media na funcao full_delta.TS caso eta != 0.
  }
  
  # Atualizando U e T
  aux_mu = X%*%beta.samp[2,,1:Gplus.samp[2]]+
   matrix(rep(b*Delta.samp[2,1:Gplus.samp[2]],n),
          nrow = n, ncol = Gplus.samp[2], byrow = T)
  aux_yxbeta = y-aux_mu
  
  # atualizando U
  u.samp[2,] = full_U.TS(aux_mu,aux_yxbeta,
                         tau2.samp[2,1:Gplus.samp[2]],
                         Delta.samp[2,1:Gplus.samp[2]],
                         nu.samp[1], t.samp[1,], z.samp[2,])
  
  # atualizando T
  t.samp[2,] = full_T.TS(aux_mu,aux_yxbeta,
                         tau2.samp[2,1:Gplus.samp[2]],
                         Delta.samp[2,1:Gplus.samp[2]],
                         u.samp[2,], z.samp[2,])
  
  
  # -------------------------------
  # b) atualizando hiperparametros
  # -------------------------------
  alpha.samp[2] = full_alpha(nu.samp[1])
  
  # ===========================================================================
  #                                  Passo 3
  # ===========================================================================
  
  # ----------------------
  # a) atualizando G
  # ----------------------
  G.samp[2] = full_K.TS(Gplus.samp[2],G.max,M,gammaProb[1],lprio_G)
  if(G.samp[2]=="erro"){
   cat(c(beta.samp[1,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/betaerro.txt"))
   cat(tau2.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2erro.txt"))
   cat(Delta.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Deltaerro.txt"))
   cat(prob.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/proberro.txt"))
   cat(nu.samp[1],gammaProb[1],alpha.samp[1],Gplus.samp[1],G.samp[1],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/paramerro.txt"))
   cat(z.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/zerro.txt"))
   cat(u.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/uerro.txt"))
   cat(t.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/terro.txt"))
   
   cat(c(beta.samp[2,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/betaerro.txt"), append=T)
   cat(tau2.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2erro.txt"), append=T)
   cat(Delta.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Deltaerro.txt"), append=T)
   cat(prob.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/proberro.txt"), append=T)
   cat(nu.samp[2],gammaProb[2],alpha.samp[2],Gplus.samp[2],G.samp[2],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/paramerro.txt"), append=T)
   cat(z.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/zerro.txt"), append=T)
   cat(u.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/uerro.txt"), append=T)
   cat(t.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/terro.txt"), append=T)
   stop("erro em G")
  }
  # ----------------------
  # b) atualizando gammaP
  # ----------------------
  gammaProb[2] = full_gammaProb.TS(gammaProb[1],phigamma,n,M,G.samp[2],Gplus.samp[2])
  contgamma = ifelse(gammaProb[2]==gammaProb[1],contgamma,contgamma + 1)
  
  
  # ===========================================================================
  #                                  Passo 4
  # ===========================================================================
  
  # ---------------------------
  # a) adicionando componentes 
  # vazios
  # ---------------------------
  dif = G.samp[2]-Gplus.samp[2]
  if(dif > 0){
   M = c(M,rep(0,dif))
   
   mtyindex = (Gplus.samp[2]+1):G.samp[2]
   
   # Atualizando beta
   beta.samp[2,,mtyindex] = matrix(rnorm(p*dif,sd=c),ncol=dif)
   
   # Atualizando tau2
   tau2.samp[2,mtyindex] = 1/rgamma(dif,r,s)
   
   # Atualizando Delta
   Delta.samp[2,mtyindex] = rnorm(dif,eta,omega)
  }
  
  # ----------------------
  # Atualizando prob
  # ----------------------
  compindex = 1:G.samp[2]
  prob.samp[2,compindex] = full_prob.TS(z.samp[2,],G.samp[2],gammaProb[2])
  
  # ----------------------
  # Atualizando nu
  # ----------------------
  nu.samp[2] = full_nu(prob.samp[2,compindex],b,n,y,X,phi,beta.samp[2,,compindex],
                       tau2.samp[2,compindex],Delta.samp[2,compindex],
                       alpha.samp[2],nu.samp[1])
  
  cont = ifelse(nu.samp[2]==nu.samp[1],cont,cont + 1)
  
  # centralizacao da media
  b = ifelse(nu.samp[2] > 1,-exp(log(nu.samp[2])/2 - log(pi)/2 + lgamma((nu.samp[2]-1)/2) - lgamma(nu.samp[2]/2)),-sqrt(nu.samp[2]/pi)*(gamma((nu.samp[2]-1)/2)/gamma(nu.samp[2]/2)))
  
  
  # output
  if(i>burn & i%%thin==0){
   cat(c(beta.samp[2,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/beta.txt"), append=T)
   cat(tau2.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2.txt"), append=T)
   cat(Delta.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Delta.txt"), append=T)
   cat(prob.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/prob.txt"), append=T)
   cat(nu.samp[2],gammaProb[2],alpha.samp[2],Gplus.samp[2],G.samp[2],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/param.txt"), append=T)
   cat(z.samp[2,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/z.txt"), append=T)
   cat(colMeans(cbind(u.samp[2,],t.samp[2,])),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/ut.txt"), append=T)
  }
  # o novo se torna antigo na proxima iteracao
  z.samp[1,] = z.samp[2,]
  beta.samp[1,,] = beta.samp[2,,]
  tau2.samp[1,] = tau2.samp[2,]
  Delta.samp[1,] = Delta.samp[2,]
  u.samp[1,] = u.samp[2,]
  t.samp[1,] = t.samp[2,]
  alpha.samp[1] = alpha.samp[2]
  Gplus.samp[1] = Gplus.samp[2]
  G.samp[1] = G.samp[2]
  gammaProb[1] = gammaProb[2]
  prob.samp[1,] = prob.samp[2,]
  nu.samp[1] = nu.samp[2]
  
  # limpando o novo
  z.samp[2,] = u.samp[2,] = t.samp[2,] = rep(NA,n)
  prob.samp[2,] = tau2.samp[2,] = Delta.samp[2,] = rep(NA, G.max)
  alpha.samp[2] = Gplus.samp[2] = G.samp[2] = gammaProb[2] = nu.samp[2] = NA
  beta.samp[2,,] = matrix(NA,p,G.max)
  
 };time_ok = Sys.time()-t_tmp
}
