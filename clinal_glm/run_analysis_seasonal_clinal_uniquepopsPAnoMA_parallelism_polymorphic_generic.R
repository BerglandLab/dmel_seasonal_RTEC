## April 2016
## testing parallelism in seasonal sites


############ Extended quantile range
#############
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(scales)
library(RColorBrewer)

# Command line args: 1) focal seasonal regression, 2) name of output file
args = commandArgs(trailingOnly=TRUE)
focalseasonalfile = args[1]
outfile = args[2]


quant = 10^(-seq(from=0, to=3, by=0.1))

##### By quantile
compare_glm_glm_quantile_sym = function(data, pops){
  #try(detach("package:data.table", unload=TRUE)) # this package interferes with the base function "order"
  ## requires a dataframe in the format
  #1) chrom
  #2) pos
  #3) glm coef
  #4) glm quant
  #5:(4+npops) pops glm coefs
  #(5+npops):4+2*npops) pops glm coefs
  quant = 10^(-seq(from=0, to=3, by=0.25))

  comparedev2 = matrix(ncol=1+2*length(pops), nrow=length(quant) )
  comparedev2[,1] = quant
  npops = length(pops)
  for (j in 1:npops){   # 3 pops
    for (i in 1:length(quant)){
      nsites = nrow(data)
      focal = data[ intersect(order(as.numeric(data$seas.p))[1:(quant[i]*nsites)], order(as.numeric(data[, 4+npops+j]))[1:(quant[i]*nsites)]), ]
      sameup = sum(  na.omit(as.numeric(focal[,j+4]) > 0 & as.numeric(focal$seas.coef) > 0) )
      samedown = sum( na.omit(as.numeric(focal[,j+4]) < 0 & as.numeric(focal$seas.coef) < 0 ) )
      zeros = trunc(sum( na.omit( as.numeric(focal[,j+4]) == 0 & !is.na(focal$seas.coef)) ) /2)
      non.na = sum( na.omit(!is.na(focal[,j+4]) & !is.na(focal$seas.coef) ))
      comparedev2[i,j+1] = non.na
      comparedev2[i,npops+j+1] = (zeros+sameup+samedown)/non.na
    }
  }
  comparedev2SD= comparedev2
  for (i in 1:nrow(comparedev2SD) ){
      for (j in 1:length(pops) ){
          comparedev2SD[i,npops+j+1] = sqrt( comparedev2[i,npops+j+1]*(1-comparedev2[i,npops+j+1])/ comparedev2[i,j+1])
      }}
  colnames(comparedev2) = c("quantile",pops,pops)
  colnames(comparedev2SD) = c("quantile",pops,pops)
  list(comparedev2, comparedev2SD)
}
compare_glm_glm = function(data, pops, p2=NULL){
  ## requires a dataframe in the format
  #1) chrom
  #2) pos
  #3) glm coef
  #4) glm pvalue
  #5:(4+npops) pops glm coefs
  #(5+npops):4+2*npops) pops glm coefs
  allrows=nrow(na.omit(data))
  pvalue = c(1,0.1,0.05, 0.01, 0.005, 0.001,0.0005, 0.0001)
  comparedev2 = matrix(ncol=1+3*length(pops), nrow=length(pvalue) )
  comparedev2[,1] = pvalue
  npops = length(pops)
  for (j in 1:npops){   # 3 pops
    for (i in 1:length(pvalue)){
      if (is.null(p2)){
        focal = data[data$seas.p <= pvalue[i], ]
      } else {
        focal = data[data$seas.p <= pvalue[i] & data[, 4+npops+j] < p2, ]
      }
      sameup = sum(  na.omit(focal[,j+4] > 0 & focal$seas.coef > 0) )
      samedown = sum( na.omit(focal[,j+4] < 0 & focal$seas.coef < 0 ) )
      zeros = trunc(sum( na.omit( focal[,j+4] == 0 & !is.na(focal$seas.coef)) ) /2)
      non.na = sum( na.omit(!is.na(focal[,j+4]) & !is.na(focal$seas.coef) ))
      comparedev2[i,j+1] = nrow(na.omit(focal))/allrows
      comparedev2[i,npops+j+1] = (zeros+sameup+samedown)/non.na
      comparedev2[i,j+3] = nrow(na.omit(focal))
    }
  }
  colnames(comparedev2) = c("pvalue","psigBoth","pconcord","NsigBoth")
  comparedev2
}



######### read in full files

# clinal regression
tempA = read.table("../glm/mel_clinal_uniquepops_springPA_noMA.glm.noheader", stringsAsFactors=FALSE) ## only spring pops
my.filter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
my.filterA = my.filter[my.filter[,1]!="X",]
tempB = merge(my.filterA, tempA, by=c(1,2))
temp = tempB[order(tempB[,1],tempB[,2]),]

## focal seasonal regression
all_fsA = read.table(focalseasonalfile, stringsAsFactors=FALSE, header=TRUE) # no year factor
all_fsB = merge(my.filter, all_fsA, by=c(1,2))
all_fsC = all_fsB[all_fsB[,1]!="X",]
all_fs = all_fsC[order(all_fsC[,1],all_fsC[,2]),]
colnames(all_fs) = c("chrom", "pos", "seas.coef", "seas.p", "seas1.N")
allcoef = na.omit(merge(all_fs[,1:4], temp[,1:4], by=c(1,2)) )
rm(tempA, my.filter, my.filterA, temp, tempB, all_fsA, all_fsB, all_fsC, all_fs)

# combine the control SNPs with the
controls = read.table("../bootstrap/bootstrap_fmean_dp.mel.medfreq01_RRgrt0.recRate.txt", stringsAsFactors=FALSE)
controls2 = merge(controls, allcoef[,1:2], by=c(1,2))
ncontrol = ncol(controls2) - 2 # first 2 columns are read chrom and pos
rm(controlbootlist,controls)


nescentClinal_Cquant_control = matrix(ncol=101, nrow=13)

for (i in 1:101){
  ## for each control, get the coef and pvalue for the NA other dataset
  focal_control = controls2[,c(1,i+1,1,2)] ## the control chrom, pos, then the real chrom and pos
  newNA = merge(focal_control, allcoef[,1:4], by=c(1,2) ) # fetch the coef and pvalue for the matched control
  # the 4 columns of the control and then real chrom/pos, then control coef/pvalue
  #newNA2 = newNA[order(newNA[,3], newNA[,4]),]
  allcoef3 = na.omit(merge(newNA[,c(3:6)], allcoef[,c(1,2,5:ncol(allcoef))], by=c(1,2)))
  rm(newNA,focal_control)

  nescentcompareQ = compare_glm_glm_quantile_sym(allcoef3, c("clinal") )
  nescentClinal_Cquant_control[,i] = nescentcompareQ[[1]][,3]
}

nescentClinal_Cquant_controlr = matrix(ncol=4, nrow=13)

for (i in 1:13){
  nescentClinal_Cquant_controlr[i,] = c(nescentClinal_Cquant_control[i,1],median(nescentClinal_Cquant_control[i,2:101], na.rm=TRUE), sort(nescentClinal_Cquant_control[i,2:101])[c(3,98)])
}

rm(tempA, my.filter, my.filterA, temp, tempB, all_fsA, all_fsB, all_fsC, all_fs, out2, result_list,pvalue,out1, NmatP, plistP, mycols,j,i,maxmin,ncontrol,chroms,clistP)
save.image(outfile)
