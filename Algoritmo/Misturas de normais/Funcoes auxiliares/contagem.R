# dado um vetor e um numero de categorias, conta quantos elementos estao em cada 
# categoria
contagem = function(vec,categ){
 sapply(1:categ, function(a) sum(vec==a))
}
