# ==================================================================
# Trata o problema de label switching (Sylvia Fruhwirth-Schnatter)
# as cadeias geradas nao estao disponiveis por conta de seu tamanho
# ===================================================================

# aplicando burn-in e thin
burn = Q/2; thin = 5
index = seq(burn+1,Q,by=5)

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
GplusHat = unique(Gplus.samp)[which.max(table(Gplus.samp))]

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

z_real = data$comp
# cortesia do chatgpt v
rm(list = setdiff(ls(),c(grep("_ok$",ls(),value = T),"z_real","X","y")))
# cortesia do chatgpt ^