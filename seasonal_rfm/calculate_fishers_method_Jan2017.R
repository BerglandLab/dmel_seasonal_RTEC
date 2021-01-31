## Jan 2017

## Performing fisher's method

########## Setup
# args: dpfile Nsites prop S Npops outfile
# args=c("mel176_082013_CO.dp.txt", 100000, 0.01, 0.0001, 25, "testselected.txt")
# eg, Rscript fishers_method.R greater_rank_fisher_exactJ_pop$i\_effectS$S.$outfile.txt lesser_rank_fisher_exactJ_pop$i\_effectS$S.$outfile.txt L_rank_fisher_exactJ_pop$i\_effectS$S.$outfile.txt

args = commandArgs(trailingOnly=TRUE)
require(pracma,lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")

grt = read.table(args[1], header=TRUE)
grtL = apply(grt, MARGIN=1, FUN=function(X) (-1)*sum(log(X)) )
les = read.table(args[2], header=TRUE)
lesL = apply(les, MARGIN=1, FUN=function(X) (-1)*sum(log(X)) )
df = data.frame(greater=grtL, lesser=lesL)
write.table(df, file=paste("results/", args[3], sep=""), quote=FALSE, col.names=TRUE, row.names=FALSE)


k=ncol(grt) ## number of populations
g1 = vector()
for (i in 0:(7*k)){
  g1[i+1] = suppressWarnings(gammainc(i,k))[2]
}
g2all = g1/(gamma(k))

updown = list(grtL, lesL)
enrichall = data.frame()
N=length(updown[[1]])

for (i in 1:2){
  p_up=updown[[i]]
  maxZ = round(max(p_up))
  nenrich = matrix(nrow=maxZ, ncol=4)
  colnames(nenrich) = c("Z","fdr","Nobs","Nexp")
  for (l in 1:maxZ){
    e1 = g2all[l+1]*N
    s1 = sum(p_up>l)
    fdr = e1/s1
    enrichall[l,1+ (i-1)*4] = l
    enrichall[l,2+ (i-1)*4] = fdr
    enrichall[l,3+ (i-1)*4] = s1
    enrichall[l,4+ (i-1)*4] = e1
  }
}
colnames(enrichall) = c("L.gr","fdr.gr","obs.gr","exp.gr","L.less","fdr.less","obs.less","exp.less")
enrichall$exp.tot = enrichall$exp.gr + enrichall$exp.less
enrichall$obs.tot = enrichall$obs.gr + enrichall$obs.less
enrichall$fdr.tot = enrichall$exp.tot/enrichall$obs.tot
enrichall$seas.tot = enrichall$obs.tot-enrichall$exp.tot
write.table(enrichall, file=paste("results/enrich_", args[3], sep=""), quote=FALSE, row.names=FALSE)

# greater=list()
# lesser=list()
# greater[[1]] = data2[,2+(ind*2)-1]
# lesser[[1]] = data2[,2+(ind*2)]
# greater[[2]] = apply(greater[[1]], MARGIN=1, FUN= function(X) 2*(-sum(log(X))) )
# lesser[[2]] = apply(lesser[[1]], MARGIN=1, FUN= function(X) 2*(-sum(log(X))) )
# 
