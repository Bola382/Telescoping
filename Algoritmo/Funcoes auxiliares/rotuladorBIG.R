# vec: vetor que indica a que compontente uma observacao pertence
# rerotula vec por tamanho da componente, ou seja, a primeira componente e aquela com o maior numero de membros
# a segunda a segunda maior e etc, em caso de empate, rotula por ordem de aparicao
rotuladorBIG = function(vec){
  len = length(vec)
  aux = vector("numeric",length = len)
  
  m = sort(table(vec),decreasing = T) # tamanhos de componentes 
  id = as.numeric(names(m)) # nomes das componentes antes da mudanca de rotulos
  Gplus = length(id) # numero de componentes com membros
  
  i = 1

  repeat{
    niguais = sum(m == m[i]) # numero de componentes com aquele tamanho
    if(niguais == 1){
      test = vec == id[i]
      index = which(test)
      aux[index] = i
      i = i + 1
    }else{
      idIguais = i:(i+niguais)
      index = rep(FALSE, len)
      
      for(j in 1:niguais){ # constroi um subconjunto de vec onde estao as componentes de mesmo tamanho
        index = index | vec == id[idIguais[j]]
      }
      
      aux[index] = rotulador(vec[index]) + i - 1 # rotula por ordem de aparicao o subconjunto
      i = i + length(table(aux[index])) # pula componentes ja preenchidas 
    }
    if(all(aux!=0)){break}
  }
  return(aux)
}
