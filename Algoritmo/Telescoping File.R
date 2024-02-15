# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# ==============================================================================
# salva resultados em um .txt

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))
source("Geracao de dados/dst.R")
load("Geracao de dados/dados1cov.Rdata")

head(data)

n = nrow(data)
p = ncol(data)-1 # -1 pois a primeira coluna contem os grupos
#obs: p = k+1

y = data[,2] # resp
X = as.matrix(cbind(1,data$x1)) # covariaveis com intercepto

Q = 100 # numero de iteracoes

# ~~~~~~~~~~~~~~~~
# valores iniciais
# ~~~~~~~~~~~~~~~~

G.max = 100 # maximo de componentes a ser verificado
G.samp = Gplus.samp = 10 # numero de componentes e clusters

lprio_G = logbnb(1:G.max) # valor da log priori de G para 1:G.max

# propostas geradas a partir das prioris
set.seed(2)

# para guardar amostras
prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = 2, ncol = G.max)
u.samp = t.samp = z.samp = matrix(NA, nrow = 2, ncol = n)
beta.samp = array(NA, dim = c(2,p,G.max), dimnames = list(c("old","new"), 1:p, 1:G.max))

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

# centralizacao da media
b = -sqrt(nu.samp/pi)*(gamma((nu.samp-1)/2)/gamma(nu.samp/2))

u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
z.samp[1,] = sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T) # latente que indica os grupos

# ~~~~~~~~~~~~~~~~~~~
# Barra de progresso
# ~~~~~~~~~~~~~~~~~~~
library(progress)
format = "(:spin) [:bar] :percent [Decorrido: :elapsedfull || Estimado: :eta] Taxa gl :taxa || Taxa gammaP :prop"
pb = progress_bar$new(format, clear = FALSE, total = Q, complete = "=", incomplete = "-", width = 100)

# ~~~~~~~~~~~~~~~~~~~
# Config output
# ~~~~~~~~~~~~~~~~~~~
cat(paste0(rep(paste0("beta",1:p),G.max),"_",rep(1:G.max,each=p)),"\n", file = "Outputs/test1/beta.txt")
cat(paste0("tau2_",1:G.max),"\n", file = "Outputs/test1/tau2.txt")
cat(paste0("Delta_",1:G.max),"\n", file = "Outputs/test1/Delta.txt")
cat(paste0("p_",1:G.max),"\n", file = "Outputs/test1/prob.txt")
cat("nu","gammaP","alpha","Gplus","G","\n", file = "Outputs/test1/param.txt")
cat(paste0("z_",1:n),"\n", file = "Outputs/test1/z.txt")
cat(paste0("u_",1:n),"\n", file = "Outputs/test1/u.txt")
cat(paste0("t_",1:n),"\n", file = "Outputs/test1/t.txt")

cat(c(beta.samp[1,,]),"\n", file = "Outputs/test1/beta.txt", append=T)
cat(tau2.samp[1,],"\n", file = "Outputs/test1/tau2.txt", append=T)
cat(Delta.samp[1,],"\n", file = "Outputs/test1/Delta.txt", append=T)
cat(prob.samp[1,],"\n", file = "Outputs/test1/prob.txt", append=T)
cat(nu.samp[1],gammaProb[1],alpha.samp[1],Gplus.samp[1],G.samp[1],"\n", file = "Outputs/test1/param.txt", append=T)
cat(z.samp[1,],"\n", file = "Outputs/test1/z.txt", append=T)
cat(u.samp[1,],"\n", file = "Outputs/test1/u.txt", append=T)
cat(t.samp[1,],"\n", file = "Outputs/test1/t.txt", append=T)

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
 # ===========================================================================
 #                                  Passo 1
 # ===========================================================================
 
 # ----------------------
 # a) Atualizando Z
 # ----------------------
 z.samp[2,] = full_Z.TS(G.samp[1],prob.samp[1,1:G.samp[1]],
                        b,n,y,X,
                        beta.samp[1,,1:G.samp[1]],
                        tau2.samp[1,1:G.samp[1]],
                        Delta.samp[1,1:G.samp[1]],
                        nu.samp[1])
 
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
  beta.samp[2,,j] = full_beta.TS(j,b,c,n,p,Xj,Uj,yj,tj,M,
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
 u.samp[2,] = full_U.TS(aux_mu,aux_yxbeta,n,
                        tau2.samp[2,1:Gplus.samp[2]],
                        Delta.samp[2,1:Gplus.samp[2]],
                        nu.samp[1], t.samp[1,], z.samp[2,])
 
 # atualizando T
 t.samp[2,] = full_T.TS(aux_mu,aux_yxbeta,n,
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
 b = -sqrt(nu.samp[2]/pi)*(gamma((nu.samp[2]-1)/2)/gamma(nu.samp[2]/2))
 
 # output
 cat(c(beta.samp[2,,]),"\n", file = "Outputs/test1/beta.txt", append=T)
 cat(tau2.samp[2,],"\n", file = "Outputs/test1/tau2.txt", append=T)
 cat(Delta.samp[2,],"\n", file = "Outputs/test1/Delta.txt", append=T)
 cat(prob.samp[2,],"\n", file = "Outputs/test1/prob.txt", append=T)
 cat(nu.samp[2],gammaProb[2],alpha.samp[2],Gplus.samp[2],G.samp[2],"\n", file = "Outputs/test1/param.txt", append=T)
 cat(z.samp[2,],"\n", file = "Outputs/test1/z.txt", append=T)
 cat(u.samp[2,],"\n", file = "Outputs/test1/u.txt", append=T)
 cat(t.samp[2,],"\n", file = "Outputs/test1/t.txt", append=T)
 
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
 
 pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%"),
                       prop = paste0(formatC(contgamma/i * 100,2,format="f"),"%")))
};time_ok = Sys.time()-t_tmp
