### April 2016

## match list of SNPs for read depth, allele freq and chromosome
## Can match on recombination rate
## (only for mel, also matching on inversion status)


## USAGE: Rscript ../match_SNPs/match_snps_dp_ch_SNPbySNP_recomb.R means_dpfreq.2L.Rdata mel_2L.bootstrap_fmean_dp.txt 100
## changed to only take one file
# 1) number of matched sets to produce
# 2) the Rdata object with the mean and dps
# 2) the output filename
# 3) the file to match

### OUTPUT:
# 1) chrom (focal SNP)
# 2) pos (focal SNP)
# 3) pos of matched SNP (one column per matched set)



## must have run export R_LIBS=/hsgs/projects/petrov/hmachado/software/R_libs/
## install packages like this:
#install.packages("foreach", repos="http://cran.r-project.org", lib="/hsgs/projects/petrov/hmachado/software/R_libs/")

#export R_LIBS=/home/hmachado/R/x86_64-unknown-linux-gnu-library/
#install.packages("foreach", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
#install.packages("doParallel", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
#install.packages("data.table", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
#install.packages("doMC", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
library(foreach)
library(iterators)
library(data.table)
library(doMC)
library(doParallel)
cl <- makeCluster(1)
registerDoParallel(cl)


args <- commandArgs(trailingOnly = TRUE)
#args = c("means_dpfreq_otherFreqBin.X.Rdata", "bootstrap_otherfmean_dp.mel_X.txt", 100)

load(args[1])
snpinfo = means_dpfreq
n = as.numeric(args[3])
fileout = args[2]

snpinfo$RRbin = trunc(snpinfo$RR)  #### only if including recombination rate
snpinfo$RRbin[snpinfo$RRbin>5] = 6   #### 7 bins of recombination rate
snpinfo$chrom = as.character(snpinfo[,1])
snpinfo$pos = as.numeric(snpinfo[,2])

inv <- data.frame(#inv = c("In2Lt", "In2RNS", "In3RK", "In3RMo", "In3RP", "In3LP", "InXA", "InXBe"),  # combined the overlapping inversions on 3R and X
                    chr=c("2L", "2R", "3R",  "3L", "X" ),
                    prox = c(2225744,11278659,7576289,3173046,13519769),
                    dist = c(13154180,16163839,24857019,16301941,19487744))

if (length(args)==4){
    focalA = read.table(args[4], stringsAsFactors=FALSE, header=FALSE)   # read in file of SNPs that need to be matched to a control. File contains header
    focal = na.omit(focalA[,1:2])  ## there should be no NA's for chrom or pos
    focalinfo = merge(snpinfo, focal, by=c(1,2) )
} else if (length(args)==3){
    focalinfo = snpinfo
} else warning("input arguments not of length 3 or 4")

snpinfo = focalinfo

registerDoMC(12)
out1 = foreach(s=1:nrow(focalinfo), .combine=rbind) %dopar%  {
    # can only match major chromosomal arms, so exlude other arms
    focalsnp = focalinfo[s,]
    if ( (focalsnp[1] %in% inv$chr) == FALSE ){  # if the chromosome is not among 2L 2R 3L 3R or X, skip
        return()
    }
    # same chrom, dp quantile, freq quantile, and rec rate
    potentialmatchesa = snpinfo[ which(as.character(snpinfo$chrom)==as.character(focalsnp$chrom) & as.numeric(snpinfo$NeQuant10)==as.numeric(focalsnp$NeQuant10) & as.numeric(snpinfo$FfreqQuant10)==as.numeric(focalsnp$FfreqQuant10) & as.numeric(snpinfo$RRbin)==as.numeric(focalsnp$RRbin) ),]
    ###  exclude the focal snp
    potentialmatches = potentialmatchesa[ which(as.numeric(potentialmatchesa[,2]) != as.numeric(focalsnp[2]) ),]
    # extract the inversion coordinates
    invCh = inv[inv[,1]==focalsnp$chrom, ]
    start = as.numeric(invCh[2]) # inversion start
    end = as.numeric(invCh[3]) # inversion end
    # if focal SNP is in the inversion, use inversion SNPs, if not, use SNPs before or after inversion
    if (focalsnp$pos >= start & focalsnp$pos <= end){
      potentialmatches2 = potentialmatches[ which(potentialmatches$pos >= start & potentialmatches$pos <= end),  ]
    } else potentialmatches2 = potentialmatches[ which(potentialmatches$pos < start | potentialmatches$pos > end),  ]
    if (nrow(potentialmatches2) < n/10){   # there has to be at least 10 different matches if doing 100 sets
      bychromSNP = c(focalsnp[1], rep(NA, times=n) )   } else { 
          bychromSNP = c(focalsnp[1:2], sample(potentialmatches2[,2], size=n, replace=TRUE) )
      }
  unlist(bychromSNP)
}

write.table(out1, file=fileout, quote=FALSE, row.names=FALSE, col.names=FALSE)

#stopCluster(cl)
