# gera observacoes de uma normal padrao truncada em (a,inf)

rtruncnorm = function(n,a){
 alpha = (a+sqrt(a^2+4))/2 # alpha otimo segundo o artigo
 samps = NULL
 i=1
 cont = 0
 while(length(samps) < n){
  prop = a + rexp(1, rate = alpha)
  aceit = ifelse(a<alpha, exp(-(alpha-prop)^2/2), exp((a-alpha)^2/2-(alpha-prop)^2/2))
  
  if(aceit>1 | aceit < 0){stop("Problema na prob de aceitacao")}
  
  if(runif(1)<=aceit){
   samps[i] = prop
   i = i + 1
  }
  cont = cont + 1
 }
 return(list(amostra = samps, taxa = round(n/cont*100,2)))
}

# n = 1000
# a = 1
# 
# set.seed(1)
# xx = rtruncnorm(n,a)
# x = xx$amostra
# y = Runuran::urnorm(n = n, lb = a)
# 
# ks.test(x,y)
# 
# hist(x,freq=F, breaks = 30)
# hist(y,add=T,freq=F,col=2,breaks=30)
# 
# plot(sort(x),sort(y));abline(a=0,b=1)
