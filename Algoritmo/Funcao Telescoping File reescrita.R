# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# salva resultados em um .txt
# ==============================================================================
# necessita carregar a pasta "Funcoes auxiliares" e a funcao "dst.R"
# e do pacote "compiler"

# valores iniciais exceto G, p, vetores de alocacao e G_+ 
# gerados das respectivas prioris

# sao dados pesos iguais as G componentes especificadas para a geracao de Z
# com base no vetor de alocacoes obtido calcula-se G_+

# folder: nome da pasta em qual serao salvos os resultados de cada replicacao
# repl: numero da replicacao atual
# data: matriz de dados, primeira coluna respostas, demais valores das covs
# ncomps: valor inicial de G
# G.max: numero maximo de componentes
# lprio_G: vetor com a log probabilidade para cada G entre 1 e G.max
# c: hiperparametro da priori de beta
# eta, omega: hiperparametros da priori de Delta
# r, s: hiperparametros da priori de tau2 
# loc: parametro de localizacao da priori de nu
# phi: parametro de taxa da proposta de MH para nu
# phigamma: "desvio padrao" da proposta de MH para gammaProb
# Q: numero de iteracoes do algoritmo
# burn: amostras a serem descartadas
# thin: tamanho dos saltos pos burnin

telescope = function(folder,repl,dados,ncomps,G.max,lprio_G, c, eta, omega, r, s, loc,phi=.3,phigamma=2.5,Q,burn,thin){
  n = nrow(dados)
  k = ncol(dados) # numero de betas
  
  y = dados[,1] # resp
  X = as.matrix(cbind(1,dados[,-1])) # covariaveis com intercepto
  
  # ~~~~~~~~~~~~~~~~
  # valores iniciais
  # ~~~~~~~~~~~~~~~~
  
  G.samp = ncomps # numero de componentes
  
  # para guardar amostras
  prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = 2, ncol = G.max)
  u.samp = t.samp = z.samp = matrix(NA, nrow = 2, ncol = n)
  beta.samp = array(NA, dim = c(2,k,G.max), dimnames = list(c("old","new"), 1:k, 1:G.max))
  
  # hiperparametros
  gammaProb = rf(1,6,3)
  xi = rep(gammaProb/G.samp,G.samp) # prob
  
  prob.samp[1,1:G.samp] = rep(1/G.samp,G.samp) # pesos
  beta.samp[1,,1:G.samp] = matrix(rnorm(k*G.samp,sd=c),ncol=G.samp) # coef reg
  tau2.samp[1,1:G.samp] = 1/rgamma(G.samp,r,s)# escala
  Delta.samp[1,1:G.samp] = rnorm(G.samp,eta,omega)# forma
  alpha.samp = runif(1,.02,.5)# hiperparametro da priori de nu
  nu.samp = loc + rexp(1,rate=alpha.samp)# gl
  
  u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
  t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
  z.samp[1,] = sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T)
  z.samp[1,] = rotulador(z.samp[1,])
  
  Gplus.samp = sum(contagem(z.samp[1,],G.samp) > 0) # numero de clusters
  
  # ~~~~~~~~~~~~~~~~~~~
  # Config output
  # ~~~~~~~~~~~~~~~~~~~
  cat(paste0(rep(paste0("beta",1:k),G.max),"_",rep(1:G.max,each=k)),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/beta.txt"))
  cat(paste0("tau2_",1:G.max),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2.txt"))
  cat(paste0("Delta_",1:G.max),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Delta.txt"))
  cat(paste0("p_",1:G.max),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/prob.txt"))
  cat("nu","gammaP","alpha","Gplus","G","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/param.txt"))
  cat(paste0("z_",1:n),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/z.txt"))
  cat("u.bar","t.bar","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/ut.txt"))
  
  if(burn == 0){
    cat(c(beta.samp[1,,]),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/beta.txt"), append=T)
    cat(tau2.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/tau2.txt"), append=T)
    cat(Delta.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/Delta.txt"), append=T)
    cat(prob.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/prob.txt"), append=T)
    cat(nu.samp[1],gammaProb[1],alpha.samp[1],Gplus.samp[1],G.samp[1],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/param.txt"), append=T)
    cat(z.samp[1,],"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/z.txt"), append=T)
    cat(colMeans(cbind(u.samp[1,],t.samp[1,])),"\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/ut.txt"), append=T)
  }
  
  # -------------------------------------------------------------------------------
  #                                    Amostrador
  # -------------------------------------------------------------------------------
  
  cont = 0 # contador de aceites de MH para nu
  contgamma = 0 # contador de aceites de MH para gammaProb
  
  enableJIT(3)
  t_tmp = Sys.time()
  
  for(i in 2:Q){
    # centralizacao da media
    b = -exp(log(nu.samp[1])/2 - log(pi)/2 + lgamma((nu.samp[1]-1)/2) - lgamma(nu.samp[1]/2))
    
    # ===========================================================================
    #                                  Passo 2
    # ===========================================================================
    
    U = diag(u.samp[1,])
    M = contagem(z.samp[1,],G.samp[1])
    
    # ----------------------------
    # a) atualizando parametros 
    # das componentes preenchidas
    # ----------------------------
    
    for(j in 1:G.samp[1]){
      if(M[j] != 0){
        index = which(z.samp[1,]==j)
        Xj = X[index,]
        uj = u.samp[1,index]
        Uj = U[index,index] 
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
      }else{ # atualizando da priori em caso de componente vazia
        dif = G.samp[1] - Gplus.samp[1]
        id.mty = (Gplus.samp[1]+1):G.samp[1]
        
        # Atualizando beta
        beta.samp[2,,id.mty] = matrix(rnorm(k*dif,sd=c),ncol=dif)
        
        # Atualizando tau2
        tau2.samp[2,id.mty] = 1/rgamma(dif,r,s)
        
        # Atualizando Delta
        Delta.samp[2,id.mty] = rnorm(dif,eta,omega)
      }
    }
    
    
    # ===========================================================================
    #                                  Passo 1
    # ===========================================================================
    
    # ----------------------
    # a) Atualizando Z
    # ----------------------
    z.samp[2,] = full_Z.TS(G.samp[1],prob.samp[1,1:G.samp[1]],
                           b,y,X,
                           beta.samp[2,,1:G.samp[1]],
                           tau2.samp[2,1:G.samp[1]],
                           Delta.samp[2,1:G.samp[1]],
                           nu.samp[1])
    if(any(z.samp[2,]=="erro")){
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
    
    
    # ------------------------------
    # atualizando variaveis latentes
    # ------------------------------
    aux_mu = X%*%beta.samp[2,,1:G.samp[1]]+ # media por individuo por componente
      matrix(rep(b*Delta.samp[2,1:G.samp[1]],n),
             nrow = n, ncol = G.samp[1], byrow = T) 
    aux_yxbeta = y-aux_mu # resposta - media por individuo por componente
   
    
    TU.samp = full_TU.TS(aux_mu, aux_yxbeta, tau2.samp[2,1:G.samp[1]], Delta.samp[2,1:G.samp[1]], nu.samp[1], z.samp[2,], u.samp[1,], t.samp[1,])
    
    # Atualizando T
    
    t.samp[2,] = TU.samp[1,]
    
    # Atualizando U
    
    u.samp[2,] = TU.samp[2,]
    
    # ------------------------
    # b) determinando Gplus e 
    # M: tamanho de cada comp
    # ------------------------
    
    # ------------
    # re-rotulando
    # ------------
    zcopy = z.samp[2,]
    z.samp[2,] = rotulador(z.samp[2,])
    
    # contando
    M = contagem(z.samp[2,],G.samp[1])
    
    # atualizando Gplus
    Gplus.samp[2] = sum(M!=0)
    
    tab = table(zcopy, z.samp[2,])
    new_label = sapply(1:Gplus.samp[2], function(a) which(tab[a,] > 0))
    id.old = as.numeric(rownames(tab))
    id.new = as.numeric(colnames(tab))[new_label]
    id.mty = which(M == 0)
    
    
    # removendo custers vazios da contagem
    M = M[1:Gplus.samp[2]]
    
    # permutando beta
    beta.samp[2,,id.new] = beta.samp[2,,id.old]
    beta.samp[2,,id.mty] = NA
    
    # permutando tau2
    tau2.samp[2, id.new] = tau2.samp[2, id.old]
    tau2.samp[2, id.mty] =  NA 
    
    # permutando Delta
    Delta.samp[2, id.new] = Delta.samp[2, id.old]
    Delta.samp[2, id.mty] = NA 
    
    # -------------------------------
    # b) atualizando hiperparametros
    # -------------------------------
    alpha.samp[2] = full_alphaExpLoc(nu.samp[1],loc)
    
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
      id.mty = (Gplus.samp[2]+1):G.samp[2]
      
      # Atualizando beta
      beta.samp[2,,id.mty] = matrix(rnorm(k*dif,sd=c),ncol=dif)
      
      # Atualizando tau2
      tau2.samp[2,id.mty] = 1/rgamma(dif,r,s)
      
      # Atualizando Delta
      Delta.samp[2,id.mty] = rnorm(dif,eta,omega)
    }
    
    # ----------------------
    # Atualizando prob
    # ----------------------
    compindex = 1:G.samp[2]
    prob.samp[2,compindex] = full_prob.TS(z.samp[2,],G.samp[2],gammaProb[2])
    
    # ----------------------
    # Atualizando nu
    # ----------------------
    nu.samp[2] = full_nuExpLoc(prob.samp[2,compindex],b,n,y,X,phi,beta.samp[2,,compindex],
                         tau2.samp[2,compindex],Delta.samp[2,compindex],
                         alpha.samp[2],nu.samp[1],loc)
    
    cont = ifelse(nu.samp[2]==nu.samp[1],cont,cont + 1)
    
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
    beta.samp[2,,] = matrix(NA,k,G.max)
    
  };time_ok = Sys.time()-t_tmp
  cat("nu.rate","gammaP.rate","\n", file = paste0("Outputs/Paralelo/",file,"/",repl,"/taxas.txt"))
  cat(cont/Q, contgamma/Q, file = paste0("Outputs/Paralelo/",file,"/",repl,"/taxas.txt"), append = T)
  return(time_ok)
}
