#### Creating a file with the mean dp (Ne) and freq and a record of the
#     10 bin quantile of each (for each SNP).
#   Mean freq calculated only using populations WITHOUT simulans contamination > 5%
#   Calculate the freq bins for the average fall freq across the three population groups

args <- commandArgs(trailingOnly = TRUE)
filein=args[1]
fileout=args[2]
##filein="mel_freqdp_X_042016_Ne_fixed.Rdata"
#args = c("mel_freqdp_042016_Ne_fixed.Rdata", "","data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt")

load(filein)
popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

pops = read.table(args[3], stringsAsFactors=FALSE)[,1]
focal = which( PY %in% pops & (S=="s" | S=="f") )

#library(data.table, lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
P.Y = PY[focal]
S1 = S[focal]
PYu = unique(P.Y)

focalfreq = freq[, focal]
focaldp = dp[, focal]
focalinfo = info


freq2 = info[,1:2]
freq2$Smean = apply(focalfreq[, which(S1=="s")], MARGIN=1, FUN=mean, na.rm=TRUE)
freq2$Fmean = apply(focalfreq[, which(S1=="f")], MARGIN=1, FUN=mean, na.rm=TRUE)
meanNe = apply(focaldp, MARGIN=1, FUN=mean, na.rm=TRUE)
means_dpfreq = cbind(freq2, meanNe)
save(means_dpfreq, file=fileout)

