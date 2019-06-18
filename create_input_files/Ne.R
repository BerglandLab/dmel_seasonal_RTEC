### adjusting for N and read depth
## the below correction captures the majority of the extra sequencing error due to pooled sequencing

## Usage: Rscript Ne.R mel_freqdp_X_042016.Rdata mel_nescent2_samplesize.txt mel_freqdp_X_042016_Ne.Rdata x

args <- commandArgs(trailingOnly = TRUE)

filein = args[1]
popsizes = args[2]
fileout = args[3]
chrom = as.character(args[4])  ##

#filein = "mel_freqdp_X_042016.Rdata"
#popsizes = "mel_nescent2_samplesize.txt"
#fileout = "mel_freqdp_X_042016_Ne.Rdata"
#chrom = "x"  ##

#### 158B is mislabelled - should be 159B

load(filein)
popN = read.table(popsizes, stringsAsFactors=FALSE, header=TRUE)

Ne = function(N,R){
    1/ ( (1/N) + (1/R) )
}


if (chrom=="X"){
    Ns = as.numeric(popN[,2])
}
if (chrom!="X"){
    Ns = as.numeric(popN[,2])*2
}

dpnew = dp
for (i in 1:ncol(dp)){
    out1 = unlist(lapply(dp[,i], FUN=function(X) Ne(N=Ns[i], R=X) ) )
    dpnew[,i] = out1
}

dp = dpnew

save(freq, dp, info, popinfo, file=fileout)
