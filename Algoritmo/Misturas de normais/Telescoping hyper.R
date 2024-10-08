# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# salva resultados em um .txt
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))

load("dados/dados2.RData")

head(dados)

n = nrow(dados)

Q = 100000 # numero de iteracoes

# ~~~~~~~~~~~~~~~~
# valores iniciais
# ~~~~~~~~~~~~~~~~

G.max = 100 # maximo de componentes a ser verificado
G.samp = Gplus.samp = 10 # numero de componentes e clusters

lprio_G = logbnb(1:G.max) # valor da log priori de G para 1:G.max

# propostas geradas a partir das prioris
set.seed(2)

# para guardar amostras
prob.samp = mu.samp = sigma2.samp = matrix(NA, nrow = 2, ncol = G.max)
z.samp = matrix(NA, nrow = 2, ncol = n)

# hiperparametros
gammaProb = rf(1,6,3)
xi = rep(gammaProb/G.samp,G.samp) # prob
alpha = 0; beta = 10
r = 2.1; s = 1.1
a = 2.1; b = 1.1 

set.seed(1)
# valores iniciais
prob.samp[1,1:G.samp] = gtools::rdirichlet(1, alpha = xi) # pesos
eta.samp = rnorm(1, mean = alpha, sd = beta)
omega.samp = 1/sqrt(rgamma(1, shape = r, rate = s))
sigma2.samp[1,1:G.samp] = 1/rgamma(G.samp,a,b)
mu.samp[1,1:G.samp] = sapply(1:G.samp, function(a) rnorm(1, mean = eta.samp, sd = omega.samp))
z.samp[1,] = sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T)

# ~~~~~~~~~~~~~~~~~~~
# Barra de progresso
# ~~~~~~~~~~~~~~~~~~~
library(progress)
format = "(:spin) [:bar] :percent [Decorrido: :elapsedfull || Estimado: :eta] Taxa gammaP :prop"
pb = progress_bar$new(format, clear = FALSE, total = Q, complete = "=", incomplete = "-", width = 100)

# ~~~~~~~~~~~~~~~~~~~
# Config output
# ~~~~~~~~~~~~~~~~~~~
test = 2
cat(paste0("mu_",1:G.max),"\n", file = paste0("Outputs2/test",test,"/mu.txt"))
cat(paste0("sigma2_",1:G.max),"\n", file = paste0("Outputs2/test",test,"/sigma2.txt"))
cat(paste0("p_",1:G.max),"\n", file = paste0("Outputs2/test",test,"/prob.txt"))
cat("eta","omega","gammaP","Gplus","G","\n", file = paste0("Outputs2/test",test,"/param.txt"))
cat(paste0("z_",1:n),"\n", file = paste0("Outputs2/test",test,"/z.txt"))

cat(mu.samp[1,],"\n", file = paste0("Outputs2/test",test,"/mu.txt"), append=T)
cat(sigma2.samp[1,],"\n", file = paste0("Outputs2/test",test,"/sigma2.txt"), append=T)
cat(prob.samp[1,],"\n", file = paste0("Outputs2/test",test,"/prob.txt"), append=T)
cat(eta.samp, omega.samp, gammaProb[1],Gplus.samp[1],G.samp[1],"\n", file = paste0("Outputs2/test",test,"/param.txt"), append=T)
cat(z.samp[1,],"\n", file = paste0("Outputs2/test",test,"/z.txt"), append=T)

# -------------------------------------------------------------------------------
#                                    Amostrador
# -------------------------------------------------------------------------------

phigamma = 2.5 # "desvio padrao" da proposta de MH para gammaProb
contgamma = 0 # contador de aceites de MH para gammaProb

library(compiler)
enableJIT(3)
t_tmp = Sys.time()
set.seed(2)
for(i in 2:Q){
  # ===========================================================================
  #                                  Passo 1
  # ===========================================================================
  
  # ----------------------
  # a) Atualizando Z
  # ----------------------
  z.samp[2,] = full_Z.TS(G.samp[1],prob.samp[1,1:G.samp[1]],
                         mu.samp[1,1:G.samp[1]],
                         sigma2.samp[1,1:G.samp[1]], n, y)
  
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
    S1 = sum(y[z.samp[2,] == j])
    mu.samp[2,j] = full_mu.TS(j, M, S1, eta.samp[1], omega.samp[1], sigma2.samp[1,j])
    
    S2 = sum((y[z.samp[2,] == j]-mu.samp[2,j])^2) # sum(y[j]-mu[j])^2
    
    sigma2.samp[2,j] = full_sig2.TS(j, M, S2, a, b)
  }
  
  # ------------------------------
  # b) atualizando hiperparametros
  # ------------------------------
  
  eta.samp[2] = full_eta.TS(G.samp[1], alpha, beta, mu.samp[2,1:Gplus.samp[2]], omega.samp[1])
  omega.samp[2] = sqrt(full_omega2.TS(G.samp[1], r, s, mu.samp[2,1:Gplus.samp[2]], eta.samp[2]))
  
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
    
    sigma2.samp[2,mtyindex] = 1/rgamma(dif,a,b)
    mu.samp[2,mtyindex] = sapply(1:dif, 
                                 function(a) rnorm(1,eta.samp[2],omega.samp[2]))
  }
  
  # ----------------------
  # Atualizando prob
  # ----------------------
  compindex = 1:G.samp[2]
  prob.samp[2,compindex] = full_prob.TS(z.samp[2,],G.samp[2],gammaProb[2])
  
  # output
  cat(mu.samp[2,],"\n", file = paste0("Outputs2/test",test,"/mu.txt"), append=T)
  cat(sigma2.samp[2,],"\n", file = paste0("Outputs2/test",test,"/sigma2.txt"), append=T)
  cat(prob.samp[2,],"\n", file = paste0("Outputs2/test",test,"/prob.txt"), append=T)
  cat(eta.samp[2],omega.samp[2],gammaProb[2],Gplus.samp[2],G.samp[2],"\n", file = paste0("Outputs2/test",test,"/param.txt"), append=T)
  cat(z.samp[2,],"\n", file = paste0("Outputs2/test",test,"/z.txt"), append=T)
  
  # o novo se torna antigo na proxima iteracao
  z.samp[1,] = z.samp[2,]
  mu.samp[1,] = mu.samp[2,]
  sigma2.samp[1,] = sigma2.samp[2,]
  eta.samp[1] = eta.samp[2]
  omega.samp[1] = omega.samp[2]
  Gplus.samp[1] = Gplus.samp[2]
  G.samp[1] = G.samp[2]
  gammaProb[1] = gammaProb[2]
  prob.samp[1,] = prob.samp[2,]
  
  # limpando o novo
  z.samp[2,] = rep(NA,n)
  prob.samp[2,] = mu.samp[2,] = sigma2.samp[2,] = rep(NA, G.max)
  eta.samp[2] = omega.samp[2] = Gplus.samp[2] = G.samp[2] = gammaProb[2] = NA
  
  pb$tick(tokens = list(prop = paste0(formatC(contgamma/i * 100,2,format="f"),"%")))
};time_ok = Sys.time()-t_tmp


burn = 75000
thin = 5

id = seq(burn+1,Q, by = thin); length(id)

saida = read.table(paste0("Outputs2/test",test,"/param.txt"),h=T)
plot((saida$G), type="l")
plot(saida$Gplus, type="l")
plot(log(saida$gammaP)[id], type="l")

barplot(prop.table(table(saida$G[id[1:501]])),  ylim = c(0,.8))
barplot(prop.table(table(saida$G[id[(length(id)-500):length(id)]])), add = T, col =2)

barplot(prop.table(table(saida$Gplus[id[1:501]])),  ylim = c(0,1))
barplot(prop.table(table(saida$Gplus[id[(length(id)-500):length(id)]])), add = T, col =2)

barplot(prop.table(table(saida$G[id])), xlim = c(0,7), ylim = c(0,1))
barplot(prop.table(table(saida$Gplus[id])), xlim = c(0,7), ylim = c(0,1))

round(prop.table(table(saida$Gplus[id])),4)
round(prop.table(table(saida$G[id])),4)

source("Re-labeling File.R")

mu.samp = read.table(paste0("Outputs2/test",test,"/mu.txt"),h=T)[id,]
sigma2.samp = read.table(paste0("Outputs2/test",test,"/sigma2.txt"),h=T)[id,]
prob.samp = read.table(paste0("Outputs2/test",test,"/prob.txt"),h=T)[id,]
z.samp = read.table(paste0("Outputs2/test",test,"/z.txt"),h=T)[id,]
gammaP.samp = saida$gammaP[id]
Gplus.samp = saida$Gplus[id]
G.samp = saida$G[id]

coda::traceplot(coda::mcmc(log(gammaP.samp)), xlab = "Iterações", ylab = expression(log(gamma)))

barplot(prop.table(table(Gplus.samp[4001:5000])),  xlim = c(0,8),ylim = c(0,.7), xlab = expression(G["+"]), ylab = "Proporção")
barplot(prop.table(table(Gplus.samp[1:1000])), add = T, col =2)
legend("topright", legend = c("Mil primeiras", "Mil últimas"), fill = c("red","gray"))

barplot(prop.table(table(G.samp[4001:5000])),  xlim = c(0,11),ylim = c(0,.5), xlab = expression(G), ylab = "Proporção")
barplot(prop.table(table(G.samp[1:1000])), add = T, col =2)
legend("topright", legend = c("Mil primeiras", "Mil últimas"), fill = c("red","gray"))

relabel(mu.samp, sigma2.samp, prob.samp, z.samp, gammaP.samp,Gplus.samp, G.samp, try=500, test, saida = 2)

saida = read.table(paste0("Outputs2/test",test,"/corrigido/param.txt"),h=T)

mu.samp = read.table(paste0("Outputs2/test",test,"/corrigido/mu.txt"),h=T)
sigma2.samp = read.table(paste0("Outputs2/test",test,"/corrigido/sigma2.txt"),h=T)
prob.samp = read.table(paste0("Outputs2/test",test,"/corrigido/prob.txt"),h=T)
z.samp = read.table(paste0("Outputs2/test",test,"/corrigido/z.txt"),h=T)
gammaP.samp = saida$gammaP
G.samp = saida$G

aux = list(prob.samp, mu.samp, sigma2.samp)
nomes = c("p",expression(mu),expression(sigma^2))

range1 = c(0,1)
range2 = mu.samp

nperm = length(gammaP.samp)
GplusHat = unique(sort(Gplus.samp))[which.max(table(Gplus.samp ))]

par(mfrow = c(3,GplusHat), mar = c(2.1, 4.1, 2.1, 1.1))
for(i in 1:3){
  for(j in 1:GplusHat){
    if(i == 1){
      plot(1, type = "n", xlim = c(0,nperm), ylim = c(0,1), xlab = "",
           ylab = ifelse(j==1,nomes[i],""))
      lines(aux[[1]][,j])
    }else{
      if(i == 3){par(mar = c(5.1, 4.1, 2.1, 1.1))}
      plot(1, type = "n", xlim = c(0,nperm), ylim = range(aux[[i]][,j]), 
           xlab = ifelse(i==3, "Iterações",""),
           ylab = ifelse(j==1,nomes[i],""))
      lines(aux[[i]][,j])
    }
  }
}

aux2 = cbind(prob.samp, mu.samp, sigma2.samp)
round(colMeans(aux2),4)
round(coda::HPDinterval(coda::mcmc(aux2)),4)

mean(G.samp)
coda::HPDinterval(coda::mcmc(G.samp))

moda = function(x) unique(sort(x))[which.max(table(x))]
classifica = apply(z.samp, 2, moda)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 1.1))
plot(y, pch = 16, cex = .7, xlab = "Índice", ylab = "Amostra", col = z)
plot(y, pch = 16, cex = .7, xlab = "Índice", ylab = "Amostra", col = classifica)
