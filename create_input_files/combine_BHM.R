### April 2016

## combining the wild and F1 flies of BHM

args <- commandArgs(trailingOnly = TRUE)

filein = args[1]  ## mel_freqdp_X_042016_Ne.Rdata
fileout = args[2] ## mel_freqdp_X_042016_Ne_fixed.Rdata
#filein = "mel_freqdp_X_042016_Ne.Rdata"
#fileout1 = "mel_freqdp_X_042016_Ne_fixed.Rdata"

load(filein)
#[51,] "mel14BHM11f1_FAT"     "BHMf1"   "14" "other" "f"  "BHMf1_14"   "TRUE"
#[52,] "mel14BHM11wild_FAT"   "BHMwild" "14" "other" "f"  "BHMwild_14" "TRUE"
#[53,] "mel14BHM7f1_SPT"      "BHMf1"   "14" "other" "s"  "BHMf1_14"   "TRUE"
#[54,] "mel14BHM7wild_SPT"    "BHMwild" "14" "other" "s"  "BHMwild_14" "TRUE"
###[1] "dp"       "filein"   "fileout1" "freq"     "info"     "popinfo"


## BHM FAT
## for each SNP
ind1 = 51
ind2 = 52
outcol = 76  ### normally 75 columns
dpnew = dp
freqnew = freq

# fall
ind1 = 51
ind2 = 52
dpnew$mel14BHM11_FAT = apply(dp, MARGIN=1, FUN=function(X) sum(X[ind1],X[ind2], na.rm=TRUE) )
newtab = cbind(dp[,ind1],dp[,ind2],freq[,ind1],freq[,ind2])
freqnew$mel14BHM11_FAT = apply(newtab, MARGIN=1, FUN= function(X) sum(c(X[3]*(X[1]/sum(X[1],X[2], na.rm=TRUE)), X[4]*(X[2]/sum(X[1],X[2], na.rm=TRUE))), na.rm=TRUE))

freqnew$mel14BHM11_FAT[dpnew$mel14BHM11_FAT==0] = NA
dpnew$mel14BHM11_FAT[dpnew$mel14BHM11_FAT==0] = NA

# spring
ind1 = 53
ind2 = 54
dpnew$mel14BHM7_SPT = apply(dp, MARGIN=1, FUN=function(X) sum(X[ind1],X[ind2], na.rm=TRUE))
newtab = cbind(dp[,ind1],dp[,ind2],freq[,ind1],freq[,ind2])
freqnew$mel14BHM7_SPT = apply(newtab, MARGIN=1, FUN= function(X) sum(c(X[3]*(X[1]/sum(X[1],X[2], na.rm=TRUE)), X[4]*(X[2]/sum(X[1],X[2], na.rm=TRUE))), na.rm=TRUE))

freqnew$mel14BHM7_SPT[dpnew$mel14BHM7_SPT==0] = NA
dpnew$mel14BHM7_SPT[dpnew$mel14BHM7_SPT==0] = NA

BHMf = c("mel14BHM11_FAT", "BHM", "14","other","f","BHM_14","TRUE","FALSE")
BHMs = c("mel14BHM7_SPT", "BHM", "14","other","s","BHM_14","TRUE","FALSE")
popinfonew = rbind(popinfo, BHMf, BHMs)

dptmp = dpnew[,c(-51, -52, -53, -54)]
freqtmp = freqnew[,c(-51, -52, -53, -54)]
popinfotmp = popinfonew[c(-51, -52, -53, -54),]

dp = dptmp
freq = freqtmp
popinfo = popinfotmp

save(dp, freq, info, popinfo, file=fileout)