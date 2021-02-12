# March 2017
grt=read.table("../greater_rank_fisher_exactJ.merged.20pop.polymorphic.txt", stringsAsFactors=FALSE, header=TRUE)
les=read.table("../lesser_rank_fisher_exactJ.merged.20pop.polymorphic.txt", stringsAsFactors=FALSE, header=TRUE)

#library(utils) #combn
# combolist = list()
# for (i in 1:24){
#   combolist[[i]] = c(1,2)
# }
# all.combos = expand.grid(combolist) # total of 16777216 combos
nperm = 200
perm.L = list()
perm.L.both = list()

for (j in 1:nperm){
  grt_new = grt
  les_new = les
  ids = sample(c(0,1), size=ncol(grt), replace=TRUE) ## 1 equals switch S/F
  my.switch = which(ids==1)
  grt_new[,my.switch] = les[,my.switch]
  les_new[,my.switch] = grt[,my.switch]
  
  greater = apply(grt_new, MARGIN=1, FUN= function(X) -sum(log(X)) )
  lesser = apply(les_new, MARGIN=1, FUN= function(X) -sum(log(X)) )
  perm.L[[j]] = list(greater, lesser)
  perm.L.both[[j]] = unlist(greater, lesser)
}

save(perm.L, file="pop20_polymorphic_permutations_June2020.permL_grt_les.Rdata")
save(perm.L.both, file="pop20_polymorphic_permutations_June2020.permL_both.Rdata")

obs.L.both = c(apply(grt, MARGIN=1, FUN= function(X) -sum(log(X)) ), apply(les, MARGIN=1, FUN= function(X) -sum(log(X)) ) )

save(obs.L.both, file="pop20_polymorphic_June2020.obsL_both.Rdata")
