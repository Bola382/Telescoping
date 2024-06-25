# ==============================================================================
# Gerando amostras da posteriori utilizando Telescoping
# salva resultados em um .txt
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
invisible(sapply(list.files("Funcoes auxiliares",pattern="*.R$",full.names=TRUE, 
                            ignore.case=TRUE),source,.GlobalEnv))

load("~/Marcus/Misturas de normais/dados.RData")

head(dados)

n = nrow(dados)

Q = 50000 # numero de iteracoes

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
eta = 0; omega = 0.01 # media eta var b/((a-1)*omega)
a = 2.1; b = 1.1 # media b/(a-1) var b^2/((a-1)^2*(a-2))

set.seed(1)
# valores iniciais
prob.samp[1,1:G.samp] = gtools::rdirichlet(1, alpha = xi) # pesos
sigma2.samp[1,1:G.samp] = 1/rgamma(G.samp,a,b)
mu.samp[1,1:G.samp] = sapply(1:G.samp, function(a) rnorm(1, mean = eta, sd = sqrt(sigma2.samp[1,a]/omega)))
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
cat(paste0("mu_",1:G.max),"\n", file = "Outputs/test1/mu.txt")
cat(paste0("sigma2_",1:G.max),"\n", file = "Outputs/test1/sigma2.txt")
cat(paste0("p_",1:G.max),"\n", file = "Outputs/test1/prob.txt")
cat("gammaP","Gplus","G","\n", file = "Outputs/test1/param.txt")
cat(paste0("z_",1:n),"\n", file = "Outputs/test1/z.txt")

cat(mu.samp[1,],"\n", file = "Outputs/test1/mu.txt", append=T)
cat(sigma2.samp[1,],"\n", file = "Outputs/test1/sigma2.txt", append=T)
cat(prob.samp[1,],"\n", file = "Outputs/test1/prob.txt", append=T)
cat(gammaProb[1],Gplus.samp[1],G.samp[1],"\n", file = "Outputs/test1/param.txt", append=T)
cat(z.samp[1,],"\n", file = "Outputs/test1/z.txt", append=T)

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
  S2 = sum(y[z.samp[2,] == j]^2)
  
  thetaj = full_theta.TS(j, M, S1, S2, eta, omega, a, b)
  mu.samp[2,j] = thetaj[1]
  sigma2.samp[2,j] = thetaj[2]
 }
 
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
                               function(a) rnorm(1,eta,
                                                 sqrt(sigma2.samp[2,mtyindex[a]]/omega)))
 }
 
 # ----------------------
 # Atualizando prob
 # ----------------------
 compindex = 1:G.samp[2]
 prob.samp[2,compindex] = full_prob.TS(z.samp[2,],G.samp[2],gammaProb[2])
 
 # output
 cat(mu.samp[2,],"\n", file = "Outputs/test1/mu.txt", append=T)
 cat(sigma2.samp[2,],"\n", file = "Outputs/test1/sigma2.txt", append=T)
 cat(prob.samp[2,],"\n", file = "Outputs/test1/prob.txt", append=T)
 cat(gammaProb[2],Gplus.samp[2],G.samp[2],"\n", file = "Outputs/test1/param.txt", append=T)
 cat(z.samp[2,],"\n", file = "Outputs/test1/z.txt", append=T)
 
 # o novo se torna antigo na proxima iteracao
 z.samp[1,] = z.samp[2,]
 mu.samp[1,] = mu.samp[2,]
 sigma2.samp[1,] = sigma2.samp[2,]
 Gplus.samp[1] = Gplus.samp[2]
 G.samp[1] = G.samp[2]
 gammaProb[1] = gammaProb[2]
 prob.samp[1,] = prob.samp[2,]
 
 # limpando o novo
 z.samp[2,] = rep(NA,n)
 prob.samp[2,] = mu.samp[2,] = sigma2.samp[2,] = rep(NA, G.max)
 Gplus.samp[2] = G.samp[2] = gammaProb[2] = NA
 
 pb$tick(tokens = list(prop = paste0(formatC(contgamma/i * 100,2,format="f"),"%")))
};time_ok = Sys.time()-t_tmp

saida = read.table("Outputs/test1/param.txt",h=T)
plot((saida$G)[40001:Q], type="l")
plot(saida$Gplus, type="l")
plot(log(saida$gammaP), type="l")

barplot(prop.table(table(saida$G[(Q/2):(Q/2+1000)])),  ylim = c(0,.5))
barplot(prop.table(table(saida$G[(Q-1000):(Q)])), add = T, col =2)

round(prop.table(table(saida$Gplus[40001:Q])),4)
round(prop.table(table(saida$G[40001:Q])),4)

source("Re-labeling File.R")

burn = Q/2
thin = 5

id = seq(Q/2+1,Q, by = thin); length(id)

mu.samp = read.table("Outputs/test1/mu.txt",h=T)[id,]
sigma2.samp = read.table("Outputs/test1/sigma2.txt",h=T)[id,]
prob.samp = read.table("Outputs/test1/prob.txt",h=T)[id,]
z.samp = read.table("Outputs/test1/z.txt",h=T)[id,]
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

relabel(mu.samp, sigma2.samp, prob.samp, z.samp, gammaP.samp,Gplus.samp, G.samp, try=10)

saida = read.table("Outputs/test1/corrigido/param.txt",h=T)

mu.samp = read.table("Outputs/test1/corrigido/mu.txt",h=T)
sigma2.samp = read.table("Outputs/test1/corrigido/sigma2.txt",h=T)
prob.samp = read.table("Outputs/test1/corrigido/prob.txt",h=T)
z.samp = read.table("Outputs/test1/corrigido/z.txt",h=T)
gammaP.samp = saida$gammaP
G.samp = saida$G

aux = list(prob.samp, mu.samp, sigma2.samp)
nomes = c("p",expression(mu),expression(sigma^2))

range1 = c(0,1)
range2 = mu.samp

par(mfrow = c(3,3), mar = c(2.1, 4.1, 2.1, 1.1))
for(i in 1:3){
 for(j in 1:3){
  if(i == 1){
   plot(1, type = "n", xlim = c(0,3000), ylim = c(0,1), xlab = "",
        ylab = ifelse(j==1,nomes[i],""))
   lines(aux[[1]][,j])
  }else{
   if(i == 3){par(mar = c(5.1, 4.1, 2.1, 1.1))}
   plot(1, type = "n", xlim = c(0,3000), ylim = range(aux[[i]][,j]), 
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
