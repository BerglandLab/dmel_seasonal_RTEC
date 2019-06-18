## Mar 2018
# fisher's method using fisher's exact test
# including option for increasing read depth

#setwd("/Users/heathermachado/nescent_melCA/fishers_method_Oct2016")
#
#load("../../data/mel_freqdp_042016_2L_fixed.Rdata")
#focalPopYear = "co_13"
#outprefix = "co_13_fishersmethod_testout"
#fileout = paste(outprefix, ".txt", sep="")
#add_dp = as.numeric(10)
## Usage: spring_fall_pvalue.R data.Rdata focalPopulationYear outfile
#while read line; do
#  Rscript spring_fall_pvalue.R ../data/mel_freqdp_2L_042016_fixed.Rdata $line $line.testout
#done < paired_spring_fall_populations.txt
#args = c("../../data/mel_freqdp_042016_Ne_fixed.Rdata", "co_13", "spring_fall_pvalue_testout.co_13", 10)

#install.packages("data.table", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")

require(data.table,lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")

args = commandArgs(trailingOnly=TRUE)

load(args[1])
focalPopYear = args[2]
fileout = args[3]
add_dp = as.numeric(args[4])

popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

focalS = which( S=="s" & PY==focalPopYear )
focalF = which( S=="f" & PY==focalPopYear )

focalA1 = cbind(info[,1:2],freq[,c(focalS,focalF)],dp[,c(focalS,focalF)]) # 852832
filter = read.table("/scratch/users/hmachado/nescent_melCA/data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
focalA = na.omit(merge(focalA1, filter, by=c(1,2)))
Nsnps = nrow(focalA)
add_dp1 = rpois(Nsnps, add_dp)
add_dp2 = rpois(Nsnps, add_dp)
focalA[,5] = round(focalA[,5] + add_dp1)
focalA[,6] = round(focalA[,6] + add_dp2)

############
focalA$dptotal= focalA[,6] + focalA[,5]
focalA$altS = round(focalA[,3]*focalA[,5])
focalA$altF = round(focalA[,4]*focalA[,6])
focalA$alt.total = focalA$altF + focalA$altS
#focalB = focalA[ focalA[,5]>=10 & focalA[,6]>=10, ]  # remove anything with dp < 10
#focalC = na.omit(focalB) # 851653
colnames(focalA) = c("chrom","pos","sFreq","fFreq","sDP","fDP","dp.total","sAlt","fAlt", "alt.total")
#q1 = quantile(focalC$dp.total, probs=c(0.01,0.99))  # remove top and bottom 5% dp
focalA$sRef = focalA$sDP-focalA$sAlt
focalA$fRef = focalA$fDP-focalA$fAlt
focalA$sRef[focalA$sRef<0] = 0
focalA$fRef[focalA$fRef<0] = 0

#focal = na.omit(focalC[focalC$dp.total>q1[1] & focalC$dp.total<q1[2], ])
focalA$greater=NA
focalA$less=NA
for (i in 1:nrow(focalA)){
  focalline = focalA[i,]
  mat1 = matrix(as.numeric(focalline[c(8,9,11,12)]), ncol=2) # 1st row is spring, 2nd row is fall. col1: alt reads, col2: ref reads
  greater = fisher.test(mat1, alternative = "greater")$p.value
  less = fisher.test(mat1, alternative = "less")$p.value
  focalA$greater[i] = greater
  focalA$less[i] = less
}

write.table(focalA, file=fileout, row.names=FALSE, quote=FALSE, col.names=TRUE)

