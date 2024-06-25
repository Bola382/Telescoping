n = 100
G = 3
p = c(.2,.3,.5)
mu = c(-6,-3,2)
sigma = c(1,.5,1)

set.seed(100)
z = sample(1:G, n, prob = p, replace = T)
y = sapply(1:n, function(a) rnorm(1, mean = mu[z[a]], sd = sigma[z[a]]))

plot(y, pch = 16, cex = .6, xlab = "Índice", ylab = "Amostra")
hist(y,breaks=16,freq = F, main = "", ylim = c(0,.25), xlab = "Amostra", ylab = "Densidade")
curve(p[1]*dnorm(x,mu[1],sigma[1])+p[2]*dnorm(x,mu[2],sigma[2])+p[3]*dnorm(x,mu[3],sigma[3]),add=T)

dados = data.frame(y = y, aloc = z)
