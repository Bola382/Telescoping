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
G.samp = Gplus.samp = 3
burn = 0 # burnin
thin = 1 # salto
Q = 100000

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

gammaProb = 1 #rf(1,6,3) # hiperparametro dos pesos das componentes
prob.samp[1,1:G.samp] = c(0.3,0.4,0.3) # gtools::rdirichlet(1, alpha = rep(gammaProb,G.samp)/G.samp) # pesos
beta.samp[1,,1:G.samp] = matrix(c(30,4,15,4,-30,-4),nrow=2,ncol=3) # matrix(rnorm(k*G.samp,sd=c),ncol=G.samp) # coef reg
tau2.samp[1,1:G.samp] = c(0.009900990099, 0.061538461538, 0.138461538462)# 1/rgamma(G.samp,r,s)# escala 
Delta.samp[1,1:G.samp] = c(0.9950371902, -1.9845557534, -2.9768336301) #rnorm(G.samp,eta,omega) # forma
alpha.samp = 0.1474441583# hiperparametro da priori de nu
nu.samp = 3 # 1 + rexp(1,rate=alpha.samp) # gl 

u.samp[1,] = c(0.46691579380, 0.22596012137, 0.97174887617, 1.68439693217, 1.24973703179, 
               1.10577492859, 0.47861808640, 0.24039646400, 0.31682958946, 0.48211153586, 
               0.63704743219, 0.20463979803, 1.32654817195, 1.12138734404, 1.42007148266, 
               1.29004706597, 0.71730331170, 1.66948797391, 1.58699542606, 0.28560030830, 
               0.38600454246, 0.51878855917, 0.59990240125, 0.30398780547, 1.30907795744, 
               0.41107151221, 1.60169288449, 1.27252336948, 0.56149560478, 0.39897957260, 
               1.08942111858, 0.28656027143, 0.03171434983, 1.27746242569, 0.59387280380, 
               0.85589642366, 0.43918234398, 0.91347356826, 0.12633476974, 0.84385599547, 
               1.14253299072, 0.55655081051, 1.10057649394, 0.57967022175, 0.65932334812, 
               1.20578297299, 0.68546545934, 0.26322703792, 0.05508530929, 2.00167841642, 
               0.77274998975, 2.90178135861, 1.02135790461, 0.27737303518, 0.58958769071, 
               0.09284385021, 0.76595936029, 0.59047525287, 0.71714865485, 0.33157519354, 
               1.76484281061, 1.28079538305, 0.52835410056, 1.12142613739, 0.90710957242, 
               1.56376343552, 0.47929869414, 1.16748088822, 1.81424361164, 0.21158971147, 
               1.66469078222, 1.21519229217, 2.14422903844, 1.09097546497, 0.08721980548, 
               0.82808805595, 1.49654591928, 1.28979008716, 0.24937923884, 0.19768194260, 
               0.12603814930, 2.36540889500, 1.23008325171, 1.41151900767, 0.94738997008, 
               2.25973678132, 0.31020250385, 0.68369596005, 0.68510937291, 0.01704718184, 
               0.43051105861, 0.90918462525, 1.03733431993, 0.76187823087, 0.58982390705, 
               1.86265143110, 0.53128778140, 0.13090203893, 1.22642641557, 0.61852568417)
  #rgamma(n, shape = nu.samp/2, rate = nu.samp/2) # t e u da representacao aumentada
t.samp[1,] = c(1.80430101265, 0.00783320333, 1.53348886437, 0.36652986279, 0.71375323959, 
               0.92624681880, 0.99645921185, 1.94948885098, 2.18823803712, 1.37812657728, 
               1.08974520628, 2.01312654329, 0.64360380868, 0.06469725608, 0.27167875993, 
               0.95659514180, 1.19953221744, 0.59422498506, 0.88883597010, 0.83862451060, 
               0.75928209243, 1.63895628005, 1.89824801658, 2.37855893957, 0.08436391155, 
               3.69605416185, 0.70372985925, 0.22355438936, 1.15538365562, 0.92232655666, 
               0.01200412591, 0.70025310388, 1.78502019323, 0.43247612414, 3.44997351136, 
               1.81622892638, 1.17636064308, 0.74625630828, 1.52736890719, 0.96425278680, 
               0.32612672252, 1.35123733558, 1.79507377851, 1.22014571289, 0.36231657772, 
               0.56002243253, 1.14391024066, 1.16746415849, 6.49168308479, 0.14573650279, 
               0.65330479312, 0.81608299945, 0.06967723378, 0.81813268335, 0.77128130770, 
               3.21991092707, 0.60833501604, 0.11771644153, 0.18479208305, 1.28044256954, 
               0.15155827385, 0.97389314235, 0.02304130397, 0.15277865065, 2.12590559910, 
               0.56272793854, 1.38779854858, 1.65708890712, 0.79006211391, 0.03834125045, 
               0.30220111304, 0.44525728272, 0.71413258919, 0.85803085625, 4.29819926358, 
               0.65257717721, 0.63403319076, 1.37129973497, 0.73171259591, 1.83654981568, 
               0.17079328194, 0.32599612370, 0.83497438237, 0.03109040480, 1.09540408853, 
               0.15862804458, 2.68462387598, 1.41760374481, 1.76112566996, 0.72803838098, 
               1.29190945605, 1.70356008562, 1.38298354523, 0.62067527363, 0.36284472144, 
               0.14212645240, 2.16239429319, 4.07830807124, 0.13057865599, 1.21201110066)
  #abs(rnorm(n))/sqrt(u.samp[1,])
z.samp[1,] = c(1, 2, 1, 2, 1, 3, 3, 1, 3, 2, 2, 3, 2, 3, 2, 3, 1, 3, 1, 3, 1, 3, 3, 2, 2, 
               3, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 3, 3, 3, 1, 2, 2, 1, 2, 2, 2, 
               3, 1, 2, 1, 3, 1, 1, 1, 3, 3, 3, 1, 3, 2, 3, 2, 3, 2, 3, 1, 2, 1, 2, 2, 2, 
               2, 3, 3, 2, 2, 2, 3, 2, 2, 3, 1, 1, 2, 3, 1, 2, 2, 3, 1, 3, 2, 1, 3, 1, 2)
  #data[,1]
  #sample(1:G.samp, n, replace = T)
  #sample(1:G.samp,n,prob=prob.samp[1,1:G.samp], replace = T) # latente que indica os grupos

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
  
  # ----------------------
  # Atualizando Z
  # ----------------------
  
  z.samp[i,] = full_Z.TS(G.samp[i-1], prob.samp[i-1,1:G.samp[i-1]], b, y, X, beta.samp[i-1,,1:G.samp[i-1]], tau2.samp[i-1,1:G.samp[i-1]], Delta.samp[i-1,1:G.samp[i-1]], nu.samp[i-1])
  # z.samp[i,z.samp[i,]==2] = 1
  # plot(X[,-1], y, col = z.samp[i,])
  
  # ------------------------------------
  # atualizando pametros por componente
  # ------------------------------------
  
  U = diag(u.samp[i-1,])
  m = contagem(z.samp[i,],G.samp[i-1]) # n de amostras por comp
  for(j in 1:G.samp[i-1]){
    if(m[j] != 0){
      index = which(z.samp[i,]==j)
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
    }
  }

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
  alpha.samp[i] = full_alpha(nu.samp[i-1])
  
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
  nu.samp[i] = full_nu(prob.samp[i,1:G.samp[i]], b, n, y, X, phi, beta.samp[i,,1:G.samp[i]], tau2.samp[i,1:G.samp[i]], Delta.samp[i,1:G.samp[i]], alpha.samp[i], nu.samp[i-1])
  cont = ifelse(nu.samp[i]==nu.samp[i-1],cont,cont + 1)

  # A INDICADORA DE NU > 1 FOI OMITIDA, ESTOU ASSUMINDO QUE A PROPOSTA TEM SUPORTE 
  # EM (1, Inf)
  
  if(nu.samp[i] < 1){print(i)}
  pb$tick(tokens = list(taxa = paste0(formatC(cont/i * 100,2,format="f"),"%"),
                        prop = paste0(formatC(contgamma/i * 100,2,format="f"),"%")))
}
plot.ts(Gplus.samp)

rm(list=ls())

# ========================
# G
# ========================

barplot(prop.table(table(G.samp[3000:4000])))
barplot(prop.table(table(G.samp[(Q-1000):Q])), col = rgb(1, 0, 0, alpha=0.5), add = T)

# ========================
# G+
# ========================

barplot(prop.table(table(Gplus.samp[3000:4000])))
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
plot.ts(beta.samp[,2,1]) # culpa da ordenacao

ii = which.max(beta.samp[,2,1])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

plot.ts(beta.samp[,1,2])
plot.ts(beta.samp[,2,2]) # culpa da ordem

ii = which.max(beta.samp[,2,2])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

plot.ts(beta.samp[,1,3])
plot.ts(beta.samp[,2,3]) # culpa da ordem

ii = which.min(beta.samp[,2,3])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

# ========================
# tau2
# ========================

plot.ts(cumsum(tau2.samp[,1])/1:Q)
plot.ts(tau2.samp[,1])

plot.ts(log(tau2.samp[,1]))

ii = which.max(tau2.samp[,1])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

plot.ts(cumsum(tau2.samp[,2])/1:Q)
plot.ts(tau2.samp[,2])

plot.ts(log(tau2.samp[,2]))

ii = which.max(tau2.samp[,2])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

plot.ts(cumsum(tau2.samp[,3])/1:Q)
plot.ts(tau2.samp[,3])

plot.ts(log(tau2.samp[,3]))

ii = which.max(tau2.samp[,2])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

# ========================
# Delta
# ========================

plot.ts(cumsum(Delta.samp[,1])/1:Q)
plot.ts(Delta.samp[,1])

ii = which.min(Delta.samp[,1])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

plot.ts(cumsum(Delta.samp[,2])/1:Q)
plot.ts(Delta.samp[,2])

ii = which.min(Delta.samp[,2])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

plot.ts(cumsum(Delta.samp[,3])/1:Q)
plot.ts(Delta.samp[,3])

ii = which.min(Delta.samp[,3])
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}

# ========================
# nu
# ========================

plot.ts(cumsum(nu.samp)/1:Q)
plot.ts(nu.samp)
plot.ts(log(nu.samp))
acf(nu.samp)

# ========================
# zmean
# ========================

plot.ts(rowMeans(z.samp))
plot.ts(cumsum(rowMeans(z.samp))/1:Q)

# ========================
# tmean
# ========================

plot.ts(rowMeans(t.samp))
plot.ts(cumsum(rowMeans(t.samp))/1:Q)

# ========================
# umean
# ========================

plot.ts(rowMeans(u.samp))
plot.ts(cumsum(rowMeans(u.samp))/1:Q)

z.samp[ii,]
zcopy.samp[ii,]

for(ii in 1:Q){
  plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
  for(jj in 1:G.samp[ii]){
    abline(beta.samp[ii,,jj], col = jj)
  }
  Sys.sleep(.3)
}

ii = which(abs(nu.samp-1)<.1)[2]-2
plot(X[,-1], y, col = z.samp[ii,], pch = 16, main = paste("iteracao",ii))
for(jj in 1:G.samp[ii]){
  abline(beta.samp[ii,,jj], col = jj)
}
