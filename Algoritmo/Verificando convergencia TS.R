# ==================================================================
# Verifica a convergência de uma cadeia MCMC
# ==================================================================

source("Funcoes auxiliares/contagem.R")
my_prop.table = function(vec,categ){
  out = contagem(vec,categ)/length(vec)
  names(out) = 1:categ
  return(out)
}

convergencia = function(file,repl){
  caminho = paste0("Outputs/Paralelo/",file,"/",repl,"/")
  caminho2 = paste0("Outputs/Paralelo/",file,"/",repl,"/corrigido/")
  
  # corrigidos
  beta.samp = as.matrix(read.table(paste0(caminho2,"beta.txt"), h =T))
  tau2.samp = read.table(paste0(caminho2,"tau2.txt"), h =T)
  Delta.samp = read.table(paste0(caminho2,"Delta.txt"), h =T)
  prob.samp = read.table(paste0(caminho2,"prob.txt"), h =T)
  zmean.samp = rowMeans(read.table(paste0(caminho2,"z.txt"), h =T))
  param = read.table(paste0(caminho,"param.txt"), h =T)
  nu.samp = param$nu; alpha.samp = param$alpha; rm("param")
  ut.samp = read.table(paste0(caminho2,"ut.txt"), h =T)
  umean.samp = ut.samp[,1]
  tmean.samp = ut.samp[,2]
  
  # sem correcao
  param = read.table(paste0(caminho,"param.txt"), h =T)
  Gplus.samp = param$Gplus; G.samp = param$G; gammaProb = param$gammaP; rm("param")
  
  # estimando G+
  GplusHat = unique(sort(Gplus.samp))[which.max(table(Gplus.samp))]
  Q = length(nu.samp)
  nsamp = nrow(tau2.samp)
  
  # --------------------------------------------
  # G
  # --------------------------------------------
  ymax1 = max(prop.table(table(G.samp[4001:5000]))) + .15
  ymax2 = max(prop.table(table(G.samp[1:1000]))) + .15
  ymax = max(ymax1,ymax2)
  barplot(my_prop.table(G.samp[4001:5000], max(G.samp)),  xlim = c(0,10.5),ylim = c(0,ymax), xlab = expression(G), ylab = "Proporção")
  barplot(my_prop.table(G.samp[1:1000], max(G.samp)), add = T, col = rgb(1, 0, 0, alpha=0.5), axes = F)
  legend("topright", legend = c("Mil primeiras", "Mil últimas"), fill = c(rgb(1, 0, 0, alpha=0.5),"gray"))
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # Gplus
  # --------------------------------------------
  ymax1 = max(prop.table(table(Gplus.samp[4001:5000]))) + .15
  ymax2 = max(prop.table(table(Gplus.samp[1:1000]))) + .15
  ymax = max(ymax1,ymax2)
  barplot(my_prop.table(Gplus.samp[4001:5000], max(G.samp)),  xlim = c(0,10.5),ylim = c(0,ymax), xlab = expression(G["+"]), ylab = "Proporção")
  barplot(my_prop.table(Gplus.samp[1:1000], max(G.samp)), add = T, col = rgb(1, 0, 0, alpha=0.5), axes = F)
  legend("topright", legend = c("Mil primeiras", "Mil últimas"), fill = c(rgb(1, 0, 0, alpha=0.5),"gray"))
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # gammaProb
  # --------------------------------------------
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
  coda::traceplot(coda::mcmc(log(gammaProb)), xlab = "Iterações",
            main = "traceplot de log(gama)")
  invisible(readline(prompt="Press [enter] to continue"))
  
  acf(gammaProb, xlab = "Iterações", main = "acf de gama")
  invisible(readline(prompt="Press [enter] to continue"))
  
  plot(cumsum(gammaProb)/1:Q, type = "l",
       xlab = expression(gama), 
       ylab = "media acumulada")
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # prob
  # --------------------------------------------
  
  par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
  for (j in 1:GplusHat) {
    ymin = min(prob.samp[, j])-.1
    ymax = max(prob.samp[, j])+.1
    coda::traceplot(coda::mcmc(prob.samp[, j]), xlab = "Iterações",
              main = eval(expression(paste("p", j))), ylim = c(ymin,ymax))
  }
  invisible(readline(prompt="Press [enter] to continue"))
  for (j in 1:GplusHat) {
    acf(prob.samp[, j], xlab = "Iterações",
        main = eval(expression(paste("p", j))))
  }
  
  for(j in 1:GplusHat){
    plot(cumsum(prob.samp[, j])/1:nsamp, type = "l",
         xlab = eval(expression(paste("p", j))), 
         ylab = "media acumulada")
  }
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # Beta
  # --------------------------------------------
  k = ncol(beta.samp)/GplusHat
  
  beta.aux = array(NA, dim = c(nsamp,k,GplusHat), dimnames = list(1:nsamp,1:k,1:GplusHat))
  enableJIT(3)
  for(i in 1:nsamp){
    qual = c(1:(k-1),0)
    for(j in 1:k){
      quais = which(1:ncol(beta.samp) %% k == qual[j])
      beta.aux[i,j,] = beta.samp[i,quais]
    }
  }
  
  par(mfrow = c(k, GplusHat), mar = c(4, 4, 3, 1))
  for (i in 1:k) {
    for (j in 1:GplusHat) {
      coda::traceplot(coda::mcmc(beta.aux[, i, j]), xlab = "Iterações",
          main = eval(expression(paste("beta", i, j))))
    }
  }
  invisible(readline(prompt="Press [enter] to continue"))
  
  for (i in 1:k) {
    for (j in 1:GplusHat) {
      plot(cumsum(beta.aux[, i, j])/1:nsamp, type = "l",
           xlab = eval(expression(paste("beta", i, j))), 
           ylab = "media acumulada")
    }
  }
  invisible(readline(prompt="Press [enter] to continue"))
  
  for (i in 1:k) {
    for (j in 1:GplusHat) {
      acf(beta.aux[, i, j], xlab = "Iterações",
          main = eval(expression(paste("beta", i, j))))
    }
  }
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # tau2
  # --------------------------------------------
  
  par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
  for (j in 1:GplusHat) {
    ymin = min(tau2.samp[, j])-.1
    ymax = max(tau2.samp[, j])+.1
    coda::traceplot(coda::mcmc(tau2.samp[, j]), xlab = "Iterações",
                    main = eval(expression(paste("tau^2", j))), ylim = c(ymin,ymax))
  }
  invisible(readline(prompt="Press [enter] to continue"))
  for (j in 1:GplusHat) {
    acf(tau2.samp[, j], xlab = "Iterações",
        main = eval(expression(paste("tau^2", j))))
  }
  
  for(j in 1:GplusHat){
    plot(cumsum(tau2.samp[, j])/1:nsamp, type = "l",
         xlab = eval(expression(paste("tau^2", j))), 
         ylab = "media acumulada")
  }
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # Delta
  # --------------------------------------------
  
  par(mfrow = c(1, GplusHat), mar = c(5.1, 4.1, 4.1, 2.1))
  for (j in 1:GplusHat) {
    ymin = min(Delta.samp[, j])-.1
    ymax = max(Delta.samp[, j])+.1
    coda::traceplot(coda::mcmc(Delta.samp[, j]), xlab = "Iterações",
                    main = eval(expression(paste("Delta", j))), ylim = c(ymin,ymax))
  }
  invisible(readline(prompt="Press [enter] to continue"))
  for (j in 1:GplusHat) {
    acf(Delta.samp[, j], xlab = "Iterações",
        main = eval(expression(paste("Delta", j))))
  }
  
  for(j in 1:GplusHat){
    plot(cumsum(Delta.samp[, j])/1:nsamp, type = "l",
         xlab = eval(expression(paste("Delta", j))), 
         ylab = "media acumulada")
  }
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # nu
  # --------------------------------------------
  
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
  ymin = min(nu.samp)-.1
  ymax = max(nu.samp)+.1
  coda::traceplot(coda::mcmc(nu.samp), xlab = "Iterações",
                  main = expression(nu), ylim = c(ymin,ymax))
  invisible(readline(prompt="Press [enter] to continue"))
  
  acf(nu.samp, xlab = "Iterações",
      main = expression(nu))

  plot(cumsum(nu.samp)/1:Q, type = "l",
       xlab = expression(nu), 
       ylab = "media acumulada")
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
  
  # --------------------------------------------
  # zmean
  # --------------------------------------------
  
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
  ymin = min(zmean.samp)-.1
  ymax = max(zmean.samp)+.1
  coda::traceplot(coda::mcmc(zmean.samp), xlab = "Iterações",
                  main = expression(bar(Z)), ylim = c(ymin,ymax))
  invisible(readline(prompt="Press [enter] to continue"))
  
  acf(zmean.samp, xlab = "Iterações",
      main = expression(bar(Z)))
  
  plot(cumsum(zmean.samp)/1:nsamp, type = "l",
       xlab = expression(bar(Z)), 
       ylab = "media acumulada")
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # umean
  # --------------------------------------------
  
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
  ymin = min(umean.samp)-.1
  ymax = max(umean.samp)+.1
  coda::traceplot(coda::mcmc(umean.samp), xlab = "Iterações",
                  main = expression(bar(U)), ylim = c(ymin,ymax))
  invisible(readline(prompt="Press [enter] to continue"))
  
  acf(umean.samp, xlab = "Iterações",
      main = expression(bar(U)))
  
  plot(cumsum(umean.samp)/1:nsamp, type = "l",
       xlab = expression(bar(U)), 
       ylab = "media acumulada")
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  # --------------------------------------------
  # tmean
  # --------------------------------------------
  
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
  ymin = min(tmean.samp)-.1
  ymax = max(tmean.samp)+.1
  coda::traceplot(coda::mcmc(tmean.samp), xlab = "Iterações",
                  main = expression(bar("T")), ylim = c(ymin,ymax))
  invisible(readline(prompt="Press [enter] to continue"))
  
  acf(tmean.samp, xlab = "Iterações",
      main = expression(bar("T")))
  
  plot(cumsum(tmean.samp)/1:nsamp, type = "l",
       xlab = expression(bar("T")), 
       ylab = "media acumulada")
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  par(mfrow = c(1,1), mar = c(4, 4, 3, 1))
}

convergencia("reescrita n100type1",1)

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
 unique(sort(vec))[which.max(table(vec))]
}

library(rgl)
plot3d(x = X[,2],y = X[,3],z = y, col = apply(z_ok,2,moda), cex= 1.1)
plot3d(x = X[,2],y = X[,3],z = y, col = z_real, cex= 1.1)

table(apply(z_ok,2,moda),
z_real)