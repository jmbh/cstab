
instabLookup = function(x,y){
  a = stabExp(x,lkup)
  b = stabExp(y,lkup)
  return(a * (1-b) + (1-a) * b)
  }

#
# rootComb = function(x,root = 1){
#   n  = sum(x)
#   co = 1
#   for(i in x){
#     #cat(n,' ',i,'\n')
#     co = co * rootChoose(n,i,root)
#     n = n - i
#   }
#   return(co)
# }
#
# stabilityExp = function(x,root = 1){
#   co  = 0
#   n  = sum(x)
#   totComb = rootComb(x,root)
#   for(i in 1:length(x)){
#     k = x[i]
#     if(k > 1){
#       rootCo     = rootChoose(n-2,k-2,root) * rootComb(x[-i],root)
#       #print(rootCo)
#       normRootCo = rootCo / totComb
#       normCo     = normRootCo ** root
#       co = co + normCo
#     }
#   }
#   return(co)
# }
#
# instab = function(x,y,root){
#   a = stabilityExp(x,root)
#   b = stabilityExp(y,root)
#   return(a * (1-b) + (1-a) * b)
# }
