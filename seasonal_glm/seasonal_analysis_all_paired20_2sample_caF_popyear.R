## April 2017
# Script to run seasonal analysis on mel seasonal data

args <- commandArgs(trailingOnly = TRUE)

filein = args[1]  ## mel_freqdp_042016_Ne.Rdata
fileout1 = args[2] ## "mel_ca.seas_pop_year.f_s.glm"
pops24 = read.table("../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]

load(filein)
popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

####### Change this to target analysis ########################
focal = which( PY %in% pops24 & (S=="s" | S=="f") )
#co_13
#rd_13
#rd_12
#BA_12
#MA_12
#NY_12
#PA_10
#PA_11
#WI_12
#VI_12
#OUK_13
#AGA_14
#CUA_14
#CWI_14
#LMA_14
#SCPA_14
#TKA_14
#CUA_15
#SON_15
#BHM_14
################################################################

freqfocal = freq[, focal]
dpfocal = dp[, focal] ### only looking at fall to frost
bothA=cbind(info[,1:2], freqfocal, dpfocal)
P1 = P[focal]
Y1 = Y[focal]
R1 = R[focal]
S1 = S[focal]
PY1 = PY[focal]
nsamples = length(S1)
filter=read.table("../data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
both=merge(bothA, filter, by=c(1,2))

seas.glm = function(x, fileout){
  # x is a vector- a row of vcf table
  freq = as.numeric(x[3:(nsamples+2)])
  dp = as.numeric(x[(nsamples+3):(2*nsamples + 2)])
  chrompos = c(as.character(x[1]), as.numeric(x[2]))
  out1 = NA
  snp1 = data.frame( freq, dp, S1, PY1, P1, Y1)
  colnames(snp1) = c("freq","dp","seas","popyear","pop","year")
  naN = tapply(is.na(snp1$dp), INDEX=snp1$seas, FUN=sum)
  Ns = tapply(snp1$seas, INDEX=snp1$seas, FUN=length) - naN
  out1 = try( 
    summary(glm(freq~seas+popyear, data=snp1, weights=snp1$dp, family=binomial() )),
    silent=TRUE
  )
  if (is.na(out1[2])){
    out2 = c(NA, NA)
  } else {
    out2 = out1$coefficient[2,c(1,4)]
  }
  writeLines(paste(c(chrompos, out2, Ns), collapse=" "),  con=fileout)
}
my.seas.wrapper = function(data, fileout){
  writetofile = file(fileout, open="wt")
  writeLines("chrom pos seas.coef seas.p seas1.N seas2.N",  con=writetofile)
  apply(data, MARGIN=1, FUN=seas.glm, fileout=writetofile)
  close(writetofile)
}

#### Change target file
my.seas.wrapper(both, fileout=fileout1 )
