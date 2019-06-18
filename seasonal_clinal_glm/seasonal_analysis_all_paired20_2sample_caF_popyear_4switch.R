## April 2017
# Script to run seasonal analysis on mel seasonal data

args <- commandArgs(trailingOnly = TRUE)

filein = args[1]  ## "data/mel_freqdp_042016_Ne_fixed.Rdata"
#filein = "data/mel_freqdp_042016_Ne_fixed.Rdata"

fileout1 = args[2] ## "mel_ca.seas_pop_year.f_s.glm"
pops24 = read.table("/scratch/users/hmachado/nescent_melCA/data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]

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
#rd_12  ## switch
#BA_12
#MA_12
#NY_12
#PA_10 ## switch ## remove
#PA_11 ## remove
#WI_12
#VI_12
#OUK_13
#AGA_14
#CUA_14
#CWI_14
#LMA_14
#SCPA_14
#TKA_14 ## switch
#CUA_15
#SON_15
#BHM_14 ## switch
################################################################

freqfocal = freq[, focal]
dpfocal = dp[, focal] ### only looking at fall to frost
bothA=cbind(info[,1:2], freqfocal, dpfocal)
P1a = P[focal]
Y1a = Y[focal]
R1a = R[focal]
S1a = S[focal]
PY1a = PY[focal]
nsamples = length(S1a)
filter=read.table("/scratch/users/hmachado/nescent_melCA/data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
bothB=merge(bothA, filter, by=c(1,2))

# switch the 4 samples that look spring/fall switched:
########### esparto 2012, MI 2014, KA2014, PA10
S1b = S1a
S1b[3] = S1a[4] # rd_12
S1b[4] = S1a[3]

S1b[13] = S1a[14] # pa_10
S1b[14] = S1a[13]

S1b[33] = S1a[34] # ka_14
S1b[34] = S1a[33]

S1b[39] = S1a[40] # BHM_14
S1b[40] = S1a[39]

both = bothB

P1 = P1a
Y1 = Y1a
R1 = R1a
S1 = S1b
PY1 = PY1a
nsamples = length(S1)




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
