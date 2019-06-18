### Feb 2017

## Sample by genomic position,  1 per every X bp (ex. 10kb), creating a dataset
# of the same initial length

# USAGE: Rscript block_sampling.R inputfile.txt outputfile.txt bp_size
# ex. Rscript block_sampling.R inputfile.txt outputfile.txt 10000

## input file
# header=TRUE
# Col 1: chrom
# Col 2: pos
# Col ...: any

## testing
#input="/Users/heathermachado/nescent_melCA/fishers_method_Oct2016/jamie_fisher_method/rank_fisher_exactJ_allpops_100K.tmp"
#output="/Users/heathermachado/nescent_melCA/fishers_method_Oct2016/jamie_fisher_method/testout.txt"
#bp_size = 10000

## set-up
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]
bp_size = as.numeric(args[3])
library(data.table, lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")

## Read in input file
infileA = read.table(input, header=TRUE, stringsAsFactors=FALSE)
# remove "NA's", if they exist in the first two columns
nas = apply(infileA[,1:2], MARGIN=1, FUN=function(X) any(is.na(X)))
infile = infileA[!nas,]
chroms = unique(infile[,1])

block_sample = function(x){
    #focal = infile[infile[,1]==chroms[i],]
  focal = x
  focal_min = min(focal[,2])
  focal_max = max(focal[,2])
  bins = seq(from=focal_min, to=focal_max, by=bp_size)
  bins2 = bins+bp_size
  focal_bin = list()
  for (j in 1:length(bins)){
    focal_bin[[j]] = focal[as.numeric(focal[,2])>=bins[j] & as.numeric(focal[,2])<bins2[j], ]
  }
  bins1 = unlist(lapply(focal_bin, FUN=nrow))
  realbins = which(bins1>0)  # which bins have at least 1 entry
  nrealbins = length(realbins)
  #nsamples = nrow(focal)/nrealbins
  #nsamples_full = floor(nsamples)
  #nsamples_resid = nsamples - nsamples_full
  #additional_bins = round(nrealbins*nsamples_resid)
  focal_out_mat = matrix(nrow=nrealbins, ncol=ncol(focal))
  for (k in 1:nrealbins){  ## for each bin with at least 1 snp
    focal_bin_k = focal_bin[[ realbins[k] ]]
    #if (k <= additional_bins){  # for a subset we sample one extra time
    focal_out_mat[k,] = unlist(focal_bin_k[sample(1:nrow(focal_bin_k), size=1, replace=TRUE), ])
      #} else {  # otherwise, we sample to nsamples_full
      #focal_out_list[[k]] = focal_bin_k[sample(1:nrow(focal_bin_k), size=nsamples_full, replace=TRUE), ]
      #}
  }
  as.data.frame(focal_out_mat)
}

bychrom = list()
for (i in 1:length(chroms)){
  focal_chrom = infile[infile[,1]==chroms[i], ]
  bychrom[[i]] = block_sample(focal_chrom)
}

my_out = data.frame(rbindlist(bychrom))
my_out2 = my_out[ , 3:ncol(my_out)]
write.table(my_out2, file=output, col.names=TRUE, row.names=FALSE, quote=FALSE)

