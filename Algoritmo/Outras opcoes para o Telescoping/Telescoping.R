# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                  ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
load("Geracao de dados/dados.Rdata")

head(data)

n = nrow(data)
p = ncol(data)-1 # -1 pois a primeira coluna contem os grupos
#obs: p = k+1

y = data[,2] # resp
X = as.matrix(cbind(1,data[,-(1:2)])) # covariaveis com intercepto

Q = 100000 # numero de iteracoes

# ~~~~~~~~~~~~~~~~
# valores iniciais
# ~~~~~~~~~~~~~~~~

G.max = 100 # maximo de componentes a ser verificado
G.samp = Gplus.samp = 10 # numero de componentes e clusters

lprio_G = logbnb(1:G.max) # valor da log priori de G para 1:G.max

# propostas geradas a partir das prioris
set.seed(2)
# set.seed(19) test 2

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

# ~~~~~~~~~~~~~~~~~~~
# Barra de progresso
# ~~~~~~~~~~~~~~~~~~~
library(progress)
format = "(:spin) [:bar] :percent [Decorrido: :elapsedfull || Estimado: :eta] Taxa gl :taxa || Taxa gammaP :prop"
pb = progress_bar$new(format, clear = FALSE, total = Q, complete = "=", incomplete = "-", width = 100)

# -------------------------------------------------------------------------------
#                                    Amostrador
# -------------------------------------------------------------------------------

phi = 0.3 # "desvio padrao" da proposta de MH para nu
cont = 0 # contador de aceites de MH para nu
phigamma = 2.5 # "desvio padrao" da proposta de MH para gammaProb
contgamma = 0 # contador de aceites de MH para gammaProb

library(compiler)
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
 
 
 pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%"),
                       prop = paste0(formatC(contgamma/i * 100,2,format="f"),"%")))
};beepr::beep();time_ok = Sys.time()-t_tmp
