## Jan 2015
# Script to run seasonal analysis on mel seasonal data

args <- commandArgs(trailingOnly = TRUE)

filein = args[1]  ## data/mel_freqdp_X_042016_Ne_fixed.Rdata
fileout1 = args[2] ## "mel_ca.seas_pop_year.f_s.glm"

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
#focal = which( (S=="s" | S=="f") & P!="FL" & P!="GA" & P!="ME" & P!="SC" & P!="TWA" & P!="UIL" & P!="VA")

### Only north america
focal = which(popID=="melPA_72011_SPT" | popID=="melSC_072010_SPT" | popID=="melGA_072008_SPT" | popID=="melFL_072010_SPT")
lat = c(25.5, 30.99, 39.88, 33.39)

#"melFL_072010_SPT": 25.5
#"melGA_072008_SPT": 30.99
#"melPA_72011_SPT": 39.88
#"melSC_072010_SPT": 33.39

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
filter=read.table("/scratch/users/hmachado/nescent_melCA/data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
both=merge(bothA, filter, by=c(1,2))



seas.glm = function(x, fileout){
  # x is a vector- a row of vcf table
  freq = as.numeric(x[3:(nsamples+2)])
  dp = as.numeric(x[(nsamples+3):(2*nsamples + 2)])
  chrompos = c(as.character(x[1]), as.numeric(x[2]))
  out1 = NA
  snp1 = data.frame( freq, dp, S1, PY1, P1, Y1, lat)
  colnames(snp1) = c("freq","dp","seas","popyear","pop","year","lat")
  naN = sum(is.na(snp1$dp))
  Ns = nrow(snp1) - naN
  out1 = try( 
    summary(glm(freq~lat, data=snp1, weights=snp1$dp, family=binomial() )),
    silent=TRUE
  )
  if (is.na(out1[2])){
    out2 = c(NA, NA)
  } else if (nrow(out1$coefficient) != 2){
    out2 = c(NA, NA)
  } else {
    out2 = out1$coefficient[2,c(1,4)]
  }
  writeLines(paste(c(chrompos, out2, Ns), collapse=" "),  con=fileout)
}
my.seas.wrapper = function(data, fileout){
  writetofile = file(fileout, open="wt")
  writeLines("chrom pos clinal.coef clinal.p clinal.N",  con=writetofile)
  apply(data, MARGIN=1, FUN=seas.glm, fileout=writetofile)
  close(writetofile)
}

#### Change target file
my.seas.wrapper(both, fileout=fileout1 )
