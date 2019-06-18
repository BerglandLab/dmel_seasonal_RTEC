## March 2017

## making a simulated dataset for testing of power of different analyses using REAL parameters

# Real paramenters
# 1) # individuals sampled: per population
# 2) # reads sampled: per SNP, per season
# 3) spring freq (use as true spring freq): per SNP

## Setup
#setwd("/Users/heathermachado/nescent_melCA/simulated_data/neutral_from_realdata")

args = commandArgs(trailingOnly=TRUE)

repN=args[1]
pop=args[2]

## Read in real data
load("../../data/mel_freqdp_042016_fixed.Rdata")
popID = popinfo[,1]
S = popinfo[,5]
PY = popinfo[,6]
samplesizesA = read.table("../../data/mel_nescent2_samplesize.txt", stringsAsFactors=FALSE, header=TRUE)
samplesizesB = data.frame(sample=c("mel14BHM7_SPT","mel14BHM11_FAT"), individuals=c(46+60,18+16))
samplesizes = rbind(samplesizesA,samplesizesB)


## Divide samples into populations to use 
#pops = read.table("../../pproduct_seasonal/paired_spring_fall_populations_noPA12.txt", stringsAsFactors=FALSE)[,1]
# Filter 
filter = read.table("../../data/chrom_pos_medfreq01_RRgrt0.txt")
auto = filter[filter[,1]!="X", ]
#pops_list = list()
#N_mat = matrix(ncol=2, nrow=length(pops))
N_mat = vector()

#for (i in 1:length(pops)){

focalA = which( S=="s" & PY == pop )
focalB = which( S=="f" & PY == pop )
N_mat[1] = samplesizes[samplesizes[,1]== popID[focalA], 2]
N_mat[2] = samplesizes[samplesizes[,1]== popID[focalB], 2]
  focal = c(focalA, focalB)
  freqfocal = freq[, focal]
  dpfocal = dp[, focal] ### only looking at fall to frost
  both=cbind(info[,1:2], freqfocal, dpfocal)
  colnames(both) = c("chrom","pos","freq.s","freq.f","dp.s","dp.f")
  pops_list = merge(both, auto, by=c(1,2))
#}
#names(pops_list) = pops


## Sampling (one round at population, one at sequencer)
#for (i in 1:length(pops)){
  focal = na.omit(pops_list)
  sind = N_mat[1]
  find = N_mat[2]
  focalnew = focal
  focalnew$NeS = 1 / ( (1/focalnew[,5]) + (1/(sind*2)) )  # calculate Neff sp
  focalnew$NeF = 1 / ( (1/focalnew[,6]) + (1/(find*2)) )  # calculate Neff fall
  
  for (j in 1:nrow(focal)){    
    # sample at the level of the population
    s1 = rbinom(1, size=sind*2, prob=focal[j,3])  # sample from spring for both sp/fall
    f1 = rbinom(1, size=find*2, prob=focal[j,3])  # sample from spring for both sp/fall
    sfreq1 = s1/(sind*2)
    ffreq1 = f1/(find*2)
    # sample at the level of the sequencer
    focalnew[j,3] = rbinom(1,size=focal[j,5], prob=sfreq1)
    focalnew[j,4] = rbinom(1,size=focal[j,6], prob=ffreq1)  
  }
  
  focalnew$dp.total= focalnew[,6] + focalnew[,5]
  focalB = focalnew[ focalnew$NeS>=10 & focalnew$NeF>=10, ]  # remove anything with dp < 10
  focalC = na.omit(focalB) # 851653
  q1 = quantile(focalC$dp.total, probs=c(0.01,0.99))  # remove top and bottom 5% dp
  focalD = na.omit(focalC[focalC$dp.total>q1[1] & focalC$dp.total<q1[2], ])
  out1 = focalD[,c(1:6)]
  write.table(out1, file=paste("neutralsample.rep",repN,".",pop, ".txt",sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE)
  # chrom, pos, salt, falt, sdp, fdp)
#}
  





