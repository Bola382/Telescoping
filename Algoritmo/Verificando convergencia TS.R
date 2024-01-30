# ==================================================================
# Verifica a convergência de uma cadeia MCMC
# as cadeias geradas nao estao disponiveis por conta de seu tamanho
# ==================================================================
library(coda)
nsamp = length(G_ok)
GplusHat = ncol(prob_ok)

# --------------------------------------------
# G
# --------------------------------------------
barplot(prop.table(table(G_ok[1:floor(.1*nsamp)])),main="",col=2, xlab = "G", ylab = "Probabilidade",ylim = c(0,.3))
barplot(prop.table(table(G_ok[floor(.9*nsamp):nsamp])),add=T,xaxt="n")
legend("topright", legend=c("primeiros 10%", "ultimos 10%"), col = c(2,1), pch = 15)

# --------------------------------------------
# prob
# --------------------------------------------

par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
for (j in 1:GplusHat) {
 traceplot(mcmc(prob_ok[, j]), xlab = "Iterações",
           main = eval(expression(paste("p", j))), ylim = c(.3,.7))
}
for (j in 1:GplusHat) {
 acf(prob_ok[, j], xlab = "Iterações",
     main = eval(expression(paste("p", j))))
}

# --------------------------------------------
# gammaProb
# --------------------------------------------
traceplot(mcmc(gammaProb_ok), xlab = "Iterações",
          main = "gamaP")

acf(gammaProb_ok, xlab = "Iterações", main = "gamaP")

# --------------------------------------------
# Beta
# --------------------------------------------
p = ncol(X)

par(mfrow = c(p, GplusHat), mar = c(4, 4, 3, 1))
for (i in 1:p) {
 for (j in 1:GplusHat) {
  traceplot(mcmc(beta_ok[, i, j]), xlab = "Iterações",
            main = eval(expression(paste("beta", i, j))))
 }
}

for (i in 1:p) {
 for (j in 1:GplusHat) {
  acf(beta_ok[, i, j], xlab = "Iterações",
            main = eval(expression(paste("beta", i, j))))
 }
}

# --------------------------------------------
# tau2
# --------------------------------------------

par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
for (j in 1:GplusHat) {
 traceplot(mcmc(tau2_ok[, j]), xlab = "Iterações",
           main = eval(expression(paste("tau^2", j))))
}

for (j in 1:GplusHat) {
 acf(tau2_ok[, j], xlab = "Iterações",
     main = eval(expression(paste("tau^2", j))))
}

# --------------------------------------------
# Delta
# --------------------------------------------

par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
for (j in 1:GplusHat) {
 traceplot(mcmc(Delta_ok[, j]), xlab = "Iterações",
           main = eval(expression(paste("Delta", j))))
}

for (j in 1:GplusHat) {
 acf(Delta_ok[, j], xlab = "Iterações",
     main = eval(expression(paste("Delta", j))))
}

# --------------------------------------------
# nu
# --------------------------------------------

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
traceplot(mcmc(nu_ok), xlab = "Iterações",
          main = "nu")

acf(nu_ok, xlab = "Iterações", main = "nu")

# --------------------------------------------
# alpha
# --------------------------------------------

traceplot(mcmc(alpha_ok), xlab = "Iterações",
          main = "alpha")

acf(alpha_ok, xlab = "Iterações", main = "alpha")

# --------------------------------------------
# z
# --------------------------------------------

moda = function(vec){
 unique(vec)[which.max(table(vec))]
}

plot3d(x = X[,2],y = X[,3],z = y, col = apply(z_ok,2,moda), cex= 1.1)
