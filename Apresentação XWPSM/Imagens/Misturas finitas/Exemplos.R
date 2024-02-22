mixnorm = function(x,G){
 sum(sapply(1:G, function(a) prob[a] * dnorm(x,mu[a],sigma[a])))
}
mixnorm = Vectorize(mixnorm,"x")
# padrao
G = 3
prob = c(.3,.6,.1)
mu = c(-5,0,5)
sigma = rep(1,3)

curve(mixnorm(x,G), from=-8, to =8, ylab = "Densidade", xlab = "", lwd = 2)

# assimetria
G = 2
prob = c(.4,.6)
mu = c(0,2)
sigma = rep(1,2)

curve(mixnorm(x,G), from=-5, to =5, ylab = "Densidade", xlab = "", lwd = 2)

# curtose
G = 2
prob = c(.4,.6)
mu = c(0,0)
sigma = c(1,2)

curve(mixnorm(x,G), from=-5, to =5, ylim = c(0,.4),
      ylab = "Densidade", xlab = "", lwd = 2)
curve(dnorm(x,0,1),add=T)

# n comp e n cluster
set.seed(3)

prob = c(.475,.475,.05)
mu = c(-5,0,5)
sigma = rep(1,3)

curve(mixnorm(x,3), from = -10, to = 10, ylab = "Densidade", xlab = "", lwd = 2)

z = sample(1:3,50, replace = T, prob = prob)
table(z)

y = rnorm(50, mean = mu[z], sd = sigma[z])
hist(y, freq = F, breaks = 20, xlim = c(-8,2), main = "", ylab = "Densidade", xlab = "")
