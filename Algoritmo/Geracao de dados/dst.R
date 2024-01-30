# densidade da st por azallini
dst = function(x,mu,sigma,lambda,nu, log = F){
 a = (x-mu)/sigma
 if(log==T){
  log(2) + dt(a, nu,log=T)-log(sigma) + pt(a*lambda*sqrt((nu+1)/(a^2+nu)), nu+1, log = T)
 }else{
  2*dt(a, nu)/sigma*pt(a*lambda*sqrt((nu+1)/(a^2+nu)), nu+1)
 }
}
