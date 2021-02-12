## Oct 2016
# make input files for fisher's exact test (jamie's)

args = commandArgs(trailingOnly=TRUE)

load(args[1])
focalPopYear = args[2]
fileout = args[3]

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
filter = read.table("../../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
focalA = merge(focalA1, filter, by=c(1,2))

focalA$dptotal= focalA[,6] + focalA[,5]
focalA$altS = round(focalA[,3]*focalA[,5])
focalA$altF = round(focalA[,4]*focalA[,6])
focalA$alt.total = focalA$altF + focalA$altS
focalB = focalA[ focalA[,5]>=10 & focalA[,6]>=10, ]  # remove anything with dp < 10
focalC = na.omit(focalB) # 851653
colnames(focalC) = c("chrom","pos","sFreq","fFreq","sDP","fDP","dp.total","sAlt","fAlt", "alt.total")
q1 = quantile(focalC$dp.total, probs=c(0.01,0.99))  # remove top and bottom 5% dp
focalC$sRef = focalC$sDP-focalC$sAlt
focalC$fRef = focalC$fDP-focalC$fAlt
focal = na.omit(focalC[focalC$dp.total>q1[1] & focalC$dp.total<q1[2], ])
out1 = focal[,c(1,2,8,9,5,6)]

write.table(out1, file=fileout, row.names=FALSE, quote=FALSE, col.names=FALSE)

