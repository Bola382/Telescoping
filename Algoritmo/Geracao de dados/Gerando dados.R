#=======================================================================
# Gerando um conjunto de dados do modelo de misturas de distribuicao St
#=======================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

library("rgl")
source("rregmixst.R")

tmp_n = 1000
probs = c(.3,.4,.3)
beta = cbind(c(3,-10,2), c(-2,-4,8), c(3,5,7))
sigma = c(1,2,3) # 1, 4, 9
lambda = c(-10,-8,8)
nu = 3

set.seed(1)

tmp_x1 = rnorm(tmp_n, mean = 10, 2)
tmp_x2 = rgamma(tmp_n, shape = 10, rate = 2)

tmp_covs = cbind(1, tmp_x1, tmp_x2)

tmp_resul = rregmixst(tmp_covs, probs, beta, sigma, lambda, nu)

data = data.frame(group = tmp_resul$grupo, resp = tmp_resul$resposta)
data = cbind.data.frame(data,tmp_covs[,-1])
colnames(data) = c("comp", "resp", "x1", "x2")

plot3d(x = data$x1,y = data$x2,z = data$resp, col = data$comp, cex = 1.1)

plot3d(x = data$x1,y = data$x2,z = data$resp, cex = 1.1)

rm(list = ls(pattern="tmp_"))
rm(list = lsf.str())

save.image("dadosSeparados.RData")
