# ==============================================================================
#                     TELESCOPING REESCRITO
# ==============================================================================
rm(list=ls())
source("Geracao de dados/dst.R")
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))

load("Geracao de dados/n100type1.RData")

dados = data[,-1]

n = nrow(dados)
k = ncol(dados) # numero de betas

y = dados[,1] # resp
X = as.matrix(cbind(1,dados[,-1])) # covariaveis com intercepto

rm(beta); rm(sigma); rm(lambda); rm(nu); rm(probs)
# ~~~~~~~~~~~~~~~~
# valores iniciais
# ~~~~~~~~~~~~~~~~

Gmax = 100
G.samp = 10
burn = 0 # burnin
thin = 1 # salto
Q = 500000
loc = 1.1 # parametro de localizacao da priori de nu

lprio_G = logbnb(1:Gmax) # valor da log priori de G para 1:Gmax

# para guardar amostras
prob.samp = tau2.samp = Delta.samp = matrix(NA, nrow = Q, ncol = Gmax)
u.samp = t.samp = z.samp = matrix(NA, nrow = Q, ncol = n)
beta.samp = array(NA, dim = c(Q,k,Gmax), dimnames = list(1:Q, 1:k, 1:Gmax))

# hiperparametros
c = 10 # de beta 
eta = 0; omega = 10 # da Delta
r = 2.1
s = 1.1 # da tau2

#deltinha = lambda/sqrt(1+lambda^2)

set.seed(1)
gammaProb = rf(1,6,3) # hiperparametro dos pesos das componentes
prob.samp[1,1:G.samp] = rep(1/G.samp,G.samp)#gtools::rdirichlet(1, alpha = rep(gammaProb,G.samp)/G.samp) # pesos
beta.samp[1,,1:G.samp] = matrix(rnorm(k*G.samp,sd=c),ncol=G.samp) # coef reg
tau2.samp[1,1:G.samp] = 1/rgamma(G.samp,r,s)# escala 
Delta.samp[1,1:G.samp] = rnorm(G.samp,eta,omega) # forma
alpha.samp = runif(1,.02,.5) # hiperparametro da priori de nu
nu.samp = loc + rexp(1,rate=alpha.samp) # gl 

u.samp[1,] = rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
t.samp[1,] = abs(rnorm(n))/sqrt(u.samp[1,])
zcopy = sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T)
z.samp[1,] = rotulador(zcopy)

Gplus.samp = sum(contagem(z.samp[1,],G.samp) > 0)

# ~~~~~~~~~~~~~~~~~~~
# Barra de progresso
# ~~~~~~~~~~~~~~~~~~~
library(progress)
format = "(:spin) [:bar] :percent [Decorrido: :elapsedfull || Estimado: :eta] Taxa gl :taxa || Taxa gammaP :prop"
pb = progress_bar$new(format, clear = FALSE, total = Q, complete = "=", incomplete = "-", width = 100)


# -------------------------------------------------------------------------------
#                                    Amostrador
# -------------------------------------------------------------------------------

enableJIT(3)
t_tmp = Sys.time()
cont = contgamma = 0
phi=.3
phigamma = 2.5
zcopy.samp = matrix(NA, nrow = Q, ncol = n)
set.seed(1)
for(i in 2:Q){
  # centralizacao da media
  b = -exp(log(nu.samp[i-1])/2 - log(pi)/2 + lgamma((nu.samp[i-1]-1)/2) - lgamma(nu.samp[i-1]/2))
  
  # ------------------------------------
  # atualizando pametros por componente
  # ------------------------------------
  
  U = diag(u.samp[i-1,])
  m = contagem(z.samp[i-1,],G.samp[i-1]) # n de amostras por comp
  aux_cluster = sum(m>0)
  for(j in 1:G.samp[i-1]){
    if(m[j] != 0){
      index = which(z.samp[i-1,]==j)
      Xj = X[index,] # dim = mjXp
      uj = u.samp[i-1, index]
      Uj = U[index,index] # dim = mjXmj
      yj = y[index] # dim = mjX1
      tj = t.samp[i-1, index] # dim = mjX1
      
      # ----------------------
      # Atualizando beta
      # ----------------------
      
      beta.samp[i,,j] = full_beta.TS(j, b, c, Xj, Uj, yj, tj, m, tau2.samp[i-1,1:G.samp[i-1]], Delta.samp[i-1,1:G.samp[i-1]])
      
      # ----------------------
      # Atualizando tau2
      # ----------------------
      tau2.samp[i,j] = full_tau2.TS(j, b, r, s, Xj, uj, yj, tj, m, beta.samp[i,,1:G.samp[i-1]], Delta.samp[i-1,1:G.samp[i-1]])
      
      # ----------------------
      # Atualizando Delta
      # ----------------------
      
      # aqui tomamos eta = 0 para evitar mais conta, lembrar de 
      # alterar a media caso eta != 0.
      Delta.samp[i,j] = full_Delta.TS(j, b, omega, Xj, uj, yj, tj, m, beta.samp[i,,1:G.samp[i-1]], tau2.samp[i,1:G.samp[i-1]])
    }else{
      dif = G.samp[i-1] - aux_cluster
      id.mty = (aux_cluster+1):G.samp[i-1]
      
      # Atualizando beta
      beta.samp[i,,id.mty] = matrix(rnorm(k*dif,sd=c),ncol=dif)
      
      # Atualizando tau2
      tau2.samp[i,id.mty] = 1/rgamma(dif,r,s)
      
      # Atualizando Delta
      Delta.samp[i,id.mty] = rnorm(dif,eta,omega)
    }
  }
  
  # ----------------------
  # Atualizando Z
  # ----------------------
  
  z.samp[i,] = full_Z.TS(G.samp[i-1], prob.samp[i-1,1:G.samp[i-1]], b, y, X, beta.samp[i,,1:G.samp[i-1]], tau2.samp[i,1:G.samp[i-1]], Delta.samp[i,1:G.samp[i-1]], nu.samp[i-1])
  # z.samp[i,z.samp[i,]==2] = 1
  # plot(X[,-1], y, col = z.samp[i,])
  
  # ------------------------------
  # atualizando variaveis latentes
  # ------------------------------
  
  # media por individuo por componente
  aux_mu = X%*%beta.samp[i,,1:G.samp[i-1]]+matrix(rep(b*Delta.samp[i,1:G.samp[i-1]],n),nrow = n, ncol = G.samp[i-1], byrow = T)
  aux_yxbeta = y-aux_mu # resposta - media por individuo por componente
  
  TU.samp = full_TU.TS(aux_mu, aux_yxbeta, tau2.samp[i,1:G.samp[i-1]], Delta.samp[i,1:G.samp[i-1]], nu.samp[i-1], z.samp[i,], u.samp[i-1,], t.samp[i-1,])
  
  # Atualizando T
  
  t.samp[i,] = TU.samp[1,]
  
  # Atualizando U
  
  u.samp[i,] = TU.samp[2,]
  
  # ----------------------
  # re-rotulando
  # ----------------------
  
  zcopy = z.samp[i,]
  zcopy.samp[i,] = zcopy
  z.samp[i,] = rotulador(z.samp[i,])
  
  # contando
  m = contagem(z.samp[i,],G.samp[i-1])
  Gplus.samp[i] = sum(m != 0) # atualizando Gplus
  
  tab = table(zcopy, z.samp[i,])
  new_label = sapply(1:Gplus.samp[i], function(a) which(tab[a,] > 0))
  id.old = as.numeric(rownames(tab))
  id.new = as.numeric(colnames(tab))[new_label]
  id.mty = which(m == 0)
  
  m = m[1:Gplus.samp[i]] # removendo componentes vazias
  
  # permutando beta
  beta.samp[i,,id.new] = beta.samp[i,,id.old]
  beta.samp[i,,id.mty] = NA
  
  # permutando tau2
  tau2.samp[i, id.new] = tau2.samp[i, id.old]
  tau2.samp[i, id.mty] =  NA 
  
  # permutando Delta
  Delta.samp[i, id.new] = Delta.samp[i, id.old]
  Delta.samp[i, id.mty] = NA 
  
  # ----------------------
  # Atualizando alpha
  # ----------------------
  
  # gama truncada em (0.02,0.5)
  alpha.samp[i] = full_alphaExpLoc(nu.samp[i-1],loc)
  
  # ----------------------
  # atualizando G
  # ----------------------
  G.samp[i] = full_K.TS(Gplus.samp[i],Gmax,m,gammaProb[i-1],lprio_G)
  
  # ----------------------
  # atualizando gammaP
  # ----------------------
  gammaProb[i] = full_gammaProb.TS(gammaProb[i-1],phigamma,n,m,G.samp[i],Gplus.samp[i])
  contgamma = ifelse(gammaProb[i]==gammaProb[i-1],contgamma,contgamma + 1)
  
  # -------------------------------
  # atualizando componentes vazios
  # -------------------------------
  dif = G.samp[i]-Gplus.samp[i]
  if(dif > 0){
    m = c(m, rep(0, dif)) # adiciona componentes vazias
    id.mty = (Gplus.samp[i]+1):G.samp[i]
    
    # Atualizando beta
    beta.samp[i,,id.mty] = matrix(rnorm(k*dif,sd=c),ncol=dif)
    
    # Atualizando tau2
    tau2.samp[i,id.mty] = 1/rgamma(dif,r,s)
    
    # Atualizando Delta
    Delta.samp[i,id.mty] = rnorm(dif,eta,omega)
  }
  
  # ----------------------
  # Atualizando prob
  # ----------------------
  
  prob.samp[i,1:G.samp[i]] = full_prob.TS(z.samp[i,], G.samp[i], gammaProb[i])
  
  # ----------------------
  # Atualizando nu
  # ----------------------
  
  # Passo de MH
  nu.samp[i] = full_nuExpLoc(prob.samp[i,1:G.samp[i]], b, n, y, X, phi, beta.samp[i,,1:G.samp[i]], tau2.samp[i,1:G.samp[i]], Delta.samp[i,1:G.samp[i]], alpha.samp[i], nu.samp[i-1], loc)
  cont = ifelse(nu.samp[i]==nu.samp[i-1],cont,cont + 1)
  
  # A INDICADORA DE NU > loc FOI OMITIDA, ESTOU ASSUMINDO QUE A PROPOSTA TEM SUPORTE 
  # EM (loc, Inf)
  
  if(nu.samp[i] < 1){print(i)}
  pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%"),
                        prop = paste0(formatC(contgamma/i * 100,2,format="f"),"%")))
}

rm(list=ls())
load("~/Marcus/Algoritmo/Outputs/testes reescrito outra ordem 500k loc 1.1.RData")

plot.ts(Gplus.samp)

# ========================
# G
# ========================

barplot(prop.table(table(G.samp[5000:6000])))
barplot(prop.table(table(G.samp[(Q-1000):Q])), col = rgb(1, 0, 0, alpha=0.5), add = T)

# ========================
# G+
# ========================

barplot(prop.table(table(Gplus.samp[5000:6000])))
barplot(prop.table(table(Gplus.samp[(Q-1000):Q])), col = rgb(1, 0, 0, alpha=0.5), add = T)

# ========================
# gammaP
# ========================

plot.ts(cumsum(gammaProb)/1:Q)
plot.ts(gammaProb)
acf(gammaProb)

# ========================
# beta
# ========================

plot.ts(beta.samp[,1,1]) 
plot.ts(beta.samp[,2,1])

plot.ts(beta.samp[,1,2]) 
plot.ts(beta.samp[,2,2]) 

plot.ts(beta.samp[,1,3]) 
plot.ts(beta.samp[,2,3]) 

# ========================
# tau2
# ========================

plot.ts(cumsum(tau2.samp[,1])/1:Q)
plot.ts(tau2.samp[,1])
acf(tau2.samp[,1])

plot.ts(log(tau2.samp[,1]))
plot.ts(cumsum(log(tau2.samp[,1]))/1:Q)

plot.ts(cumsum(tau2.samp[,2])/1:Q)
plot.ts(tau2.samp[,2])

plot.ts(log(tau2.samp[,2]))
plot.ts(cumsum(log(tau2.samp[,2]))/1:Q)

plot.ts(cumsum(tau2.samp[,3])/1:Q)
plot.ts(tau2.samp[,3])

plot.ts(log(tau2.samp[,3]))

# ========================
# Delta
# ========================

plot.ts(cumsum(Delta.samp[,1])/1:Q)
plot.ts(Delta.samp[,1])
acf(Delta.samp[,1])

plot.ts(Delta.samp[seq(7e+4,8e+4, by = 20),1])

plot.ts(cumsum(Delta.samp[,2])/1:Q)
plot.ts(Delta.samp[,2])
acf(Delta.samp[,2])

plot.ts(cumsum(Delta.samp[,3])/1:Q)
plot.ts(Delta.samp[,3])
acf(Delta.samp[is.na(Delta.samp[,3]) == F,3])

# ========================
# nu
# ========================

plot.ts(cumsum(nu.samp)/1:Q)
plot.ts(nu.samp)
plot.ts(log(nu.samp))
acf(nu.samp[seq(5000,Q, by = 15)])

# ========================
# zmean
# ========================

plot.ts(rowMeans(z.samp))
plot.ts(cumsum(rowMeans(z.samp))/1:Q)
acf(rowMeans(z.samp))

# ========================
# tmean
# ========================

plot.ts(rowMeans(t.samp))
plot.ts(cumsum(rowMeans(t.samp))/1:Q)
acf(rowMeans(t.samp))

# ========================
# umean
# ========================

plot.ts(rowMeans(u.samp))
plot.ts(cumsum(rowMeans(u.samp))/1:Q)
acf(rowMeans(u.samp))


for(ii in 1:Q){
  plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
  for(jj in 1:G.samp[ii]){
    abline(beta.samp[ii,,jj], col = jj)
  }
  Sys.sleep(.5)
}

rm(list=ls())
load("~/Marcus/Algoritmo/Outputs/testes reescrito outra ordem 500k loc 1.1.RData")

burn = 5000
thin = 15

index0 = seq(burn,Q, by = thin)

# estimando G+
GplusHat = unique(sort(Gplus.samp[index0]))[which.max(table(Gplus.samp[index0]))]

Q = sum(Gplus.samp[index0]==GplusHat) # amostras iguais ao MAP
index = which(Gplus.samp[index0]==GplusHat) # indices iguais ao MAP
n = length(y) # tamanho amostral
k = ncol(X) # numero de betas

# filtrando amostras com G+ igual ao estimado
beta.samp = beta.samp[index0,,][index,,]
tau2.samp = tau2.samp[index0,][index,]
Delta.samp = Delta.samp[index0,][index,]
nu.samp = nu.samp[index0][index]
prob.samp = prob.samp[index0,][index,]
gammaProb = gammaProb[index0][index]
alpha.samp = alpha.samp[index0][index]
z.samp = z.samp[index0,][index,]
G.samp = G.samp[index0][index]
Gplus.samp = Gplus.samp[index0][index]

L = label.switching::label.switching(method = "ECR", zpivot = data[,1], z = as.matrix(z.samp), K = GplusHat)
#L = label.switching::label.switching(method = "STEPHENS", p = paux, z = as.matrix(z.samp))
new_label = L$permutations[[1]]

# para guardar valores pos correcao
prob_ok = matrix(NA, nrow = Q, ncol = GplusHat, dimnames = list(1:Q,paste0("p_",1:GplusHat)))
beta_ok = array(NA, dim = c(Q,k,GplusHat), dimnames = list(1:Q,1:k,1:GplusHat))
tau2_ok = matrix(NA, nrow = Q, ncol = GplusHat, dimnames = list(1:Q,paste0("tau2_",1:GplusHat)))
Delta_ok = matrix(NA, nrow = Q, ncol = GplusHat, dimnames = list(1:Q,paste0("Delta_",1:GplusHat)))
z_ok = matrix(NA, nrow = Q, ncol = n, dimnames = list(1:Q,1:n))
nu_ok = nu.samp

for(q in 1:Q){
  for(j in 1:GplusHat){
    prob_ok[q,new_label[q,j]] = prob.samp[q,j]
    tau2_ok[q,new_label[q,j]] = tau2.samp[q,j]
    Delta_ok[q,new_label[q,j]] = Delta.samp[q,j]
    beta_ok[q,,new_label[q,j]] = beta.samp[q,,j]
  }
  for(i in 1:n){
    z_ok[q,i] = new_label[q,z.samp[q,i]]
  }
}

# plot.ts(beta_ok[,1,1])
# plot.ts(beta_ok[,2,1])
# 
# plot.ts(beta_ok[,1,2])
# plot.ts(beta_ok[,2,2])
# 
# plot.ts(beta_ok[,1,3])
# plot.ts(beta_ok[,2,3])

plot.ts(rowMeans(t.samp[index0,][index,])[1:3000])

plot.ts(Delta_ok[1:3000,3])
plot.ts(cumsum(Delta_ok[,3])/1:Q)

acf(Delta_ok[,3])

plot.ts(Delta_ok[(Q-2000):Q,3])

ii = which.max(beta_ok[,1,1])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

