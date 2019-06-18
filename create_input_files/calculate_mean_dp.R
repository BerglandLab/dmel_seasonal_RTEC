## Jan 2015
# Script to run seasonal analysis on mel seasonal data

args <- commandArgs(trailingOnly = TRUE)

filein = args[1]  ## mel_freqdp_X_042016_Ne.Rdata
fileout1 = args[2] ## "mel_ca.seas_pop_year.f_s.glm"

load(filein)

popID = popinfo[,1]

meandp = apply(dp, MARGIN=2, FUN=mean, na.rm=TRUE)
meandptab = cbind(popinfo, meandp)

write.table(meandptab, file=fileout1, quote=FALSE, row.names=FALSE, col.names=FALSE)
