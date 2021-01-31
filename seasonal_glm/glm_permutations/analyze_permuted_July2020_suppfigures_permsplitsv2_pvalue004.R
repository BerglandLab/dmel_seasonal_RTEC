## Nov 2020
#library(ggplot2)


## Re-worked the second part to calculate enrichment at obs quantiles (not p-values)
args = commandArgs(trailingOnly=TRUE)
permN = args[1]
#print(permN)

# delineation of major inversions (combined the three on 3R)
mychrom_inv = data.frame(chrom=c("2L","2R","3L","3R"), start=c(2225744,11278659,3173046,7576289), end=c(13154180,16163839,16301941,24857019) )

# quantiles to test
quant = c(10^(-seq(from=0, to=4, by=1/4)), 0.003894319)

# chromosomes
mychrom = c("2L","2R","3L","3R")


## for 500 permutations
#setwd("/lustre/scratch116/casm/cgp/users/hm8/nescent_melCA/glm_permutations")
perm = na.omit(read.table(paste("results/permute",permN,"_mel_all_paired20_2sample_caF_popyear.glm", sep=""), header=TRUE, stringsAsFactors=FALSE))
#glm_pvalue = read.table("perm500_pvalue.txt", stringsAsFactors=F, header=T, sep="\t")
myfilter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
myfilter$ID = paste(myfilter[,1], myfilter[,2], sep="_")
perm$ID =  paste(perm[,1], perm[,2], sep="_")
myinds = which(perm$ID %in% myfilter$ID)
perm_filter = perm[myinds,]
perm_filter_inv = subset(perm_filter, (chrom=="2L" & pos>2225744 & pos<13154180) |
    (chrom=="2R" & pos>11278659 & pos<16163839) |
    (chrom=="3L" & pos>3173046 & pos<16301941) |
    (chrom=="3R" & pos>7576289 & pos<24857019) )
perm_filter_noninv = subset(perm_filter, (chrom=="2L" & pos<2225744 | pos>13154180) |
    (chrom=="2R" & pos<11278659 | pos>16163839) |
    (chrom=="3L" & pos<3173046 | pos>16301941) |
    (chrom=="3R" & pos<7576289 | pos>24857019) )



## supp figure showing proportion of permutations over observed
## enrichment per chromosome and within and outside of inversions (total 12 regions)
# Inversion Chromosome Proximal   Distal
#1   In(2L)t         2L  2225744 13154180
#2  In(2R)NS         2R 11278659 16163839
#3   In(3L)P         3L  3173046 16301941
#4   In(3R)K         3R  7576289 21966092
#5  In(3R)Mo         3R 17232639 24857019
#6   In(3R)P         3R 12257931 20569732
# for each quantile, 
perm_filter_quant = quantile(perm_filter$seas.p, probs=quant, na.rm=TRUE)
perm_filter_quant_inv = quantile(perm_filter_inv$seas.p, probs=quant, na.rm=TRUE)
perm_filter_quant_noninv = quantile(perm_filter_noninv$seas.p, probs=quant, na.rm=TRUE)


## supp figure showing the enrichment in seasonals SNPS (#SNPS at a given pvalue) of observed over permutations
## enrichment per chromosome
# Fold enrichment of seasonal SNPs over permutation (obs/perm)
# Using p=0.004 as cutoff (obs quantile of 0.01)
## Using pvalues from observed quantiles
perm_filter_pvalueN = unlist(lapply(quant, FUN=function(X) mean(perm_filter$seas.p<X, na.rm=TRUE)))
perm_filterN = nrow(perm_filter)
perm_filter_pvalueN_inv = unlist(lapply(quant, FUN=function(X) mean(perm_filter_inv$seas.p<X, na.rm=TRUE)))
perm_filterN_inv = nrow(perm_filter_inv)
perm_filter_pvalueN_noninv = unlist(lapply(quant, FUN=function(X) mean(perm_filter_noninv$seas.p<X, na.rm=TRUE)))
perm_filterN_noninv = nrow(perm_filter_noninv)



##### per chromosome
perm_filter_quant_bychrom = perm_filter_quant_bychromN = list()
perm_filter_quant_bychrom_inv = perm_filter_quant_bychrom_invN = list()
perm_filter_quant_bychrom_noninv = perm_filter_quant_bychrom_noninvN = list()
perm_filter_focalN = perm_filter_pvalueN_focal = list()
perm_filter_focalN_inv = perm_filter_pvalueN_focal_inv = list()
perm_filter_focalN_noninv = perm_filter_pvalueN_focal_noninv = list()

for (i in 1:length(mychrom)){
    
    # whole chromsome
    focal = subset(perm_filter, chrom==mychrom[i])
    perm_filter_quant_bychromN[[i]] = nrow(focal)
    perm_filter_quant_bychrom[[i]] = quantile(focal$seas.p, probs=quant, na.rm=TRUE)
    # number of snps at observed p-value, for each permutation
    perm_filter_focalN[[i]] = nrow(focal)
    perm_filter_pvalueN_focal[[i]] = unlist(lapply(quant, FUN=function(X) mean(focal$seas.p<X, na.rm=TRUE)))

    # for the inversion
    focal = subset(perm_filter_inv, chrom==mychrom[i])
    perm_filter_quant_bychrom_invN[[i]] = nrow(focal)
    perm_filter_quant_bychrom_inv[[i]] = quantile(focal$seas.p, probs=quant, na.rm=TRUE)
    # number of snps at observed p-value, for each permutation
    perm_filter_focalN_inv[[i]] = nrow(focal)
    perm_filter_pvalueN_focal_inv[[i]] = unlist(lapply(quant, FUN=function(X) mean(focal$seas.p<X, na.rm=TRUE)))
    
    # outside of the inversions
    focal = subset(perm_filter_noninv, chrom==mychrom[i])
    perm_filter_quant_bychrom_noninvN[[i]] = nrow(focal)
    perm_filter_quant_bychrom_noninv[[i]] = quantile(focal$seas.p, probs=quant, na.rm=TRUE)
    # number of snps at observed p-value, for each permutation
    perm_filter_focalN_noninv[[i]] = nrow(focal)
    perm_filter_pvalueN_focal_noninv[[i]] = unlist(lapply(quant, FUN=function(X) mean(focal$seas.p<X, na.rm=TRUE)))
}
names(perm_filter_quant_bychrom) = names(perm_filter_quant_bychrom_inv) = names(perm_filter_quant_bychrom_noninv) = mychrom
names(perm_filter_quant_bychromN) = names(perm_filter_quant_bychrom_invN) = names(perm_filter_quant_bychrom_noninvN) = mychrom


save(perm_filter_quant, perm_filter_pvalueN, perm_filterN, perm_filter_focalN, perm_filter_pvalueN_focal, perm_filter_quant_bychrom, perm_filter_quant_bychrom_inv, perm_filter_quant_bychrom_noninv, perm_filter_quant_bychromN, perm_filter_quant_bychrom_invN, perm_filter_quant_bychrom_noninvN, perm_filter_focalN_inv, perm_filter_pvalueN_focal_inv, perm_filter_focalN_noninv, perm_filter_pvalueN_focal_noninv, perm_filter_quant_inv, perm_filter_quant_noninv, perm_filter_pvalueN_inv, perm_filterN_inv, perm_filter_pvalueN_noninv, perm_filterN_noninv, file=paste("results/perm", permN,"_objects_v2_pvalue004.Rdata", sep=""))

 

