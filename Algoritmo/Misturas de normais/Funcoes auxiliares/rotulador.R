# vec: vetor que indica a que compontente uma observacao pertence
# rerotula vec por ordem de aparicao, ou seja, a primeira observacao esta na
# componente 1 juntamente com as demais em vec com o mesmo componente, 
# a primeira observacao que nao esta no componente 1 esta no 2 juntamente com as 
# demais em vec que estao no mesmo componente e etc
rotulador = function(vec){
 len = length(vec)
 aux = vector("numeric",length = len)
 
 for(i in 1:len){
  index1 = which(aux==0)[1]
  index2 = which(vec==vec[index1])
  aux[index2] = i
  if(all(aux!=0)){break}
 }
 return(aux)
}
