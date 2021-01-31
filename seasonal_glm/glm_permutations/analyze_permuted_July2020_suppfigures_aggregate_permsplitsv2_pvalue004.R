## Nov 2020
#library(ggplot2)

#setwd("/lustre/scratch116/casm/cgp/users/hm8/nescent_melCA/glm_permutations")

# each perm*_objects.Rdata (one for each of 500 permutations) has the following objects:
# perm_filter_quant, perm_filter_pvalueN, perm_filterN, perm_filter_focalN, perm_filter_pvalueN_focal, perm_filter_quant_bychrom, perm_filter_quant_bychrom_inv, perm_filter_quant_bychrom_noninv, perm_filter_quant_bychromN, perm_filter_quant_bychrom_invN, perm_filter_quant_bychrom_noninvN, perm_filter_focalN_inv, perm_filter_pvalueN_focal_inv, perm_filter_focalN_noninv, perm_filter_pvalueN_focal_noninv

## for 500 permutations
totalperms = 500
#totalperms = 10

perm_filter_quant_list = perm_filter_pvalueN_focal_list = perm_filter_pvalueN_list = perm_filter_focalN_list = perm_filter_quant_bychrom_list = list()
perm_filter_quant_bychrom_invN_list = perm_filter_quant_bychrom_noninvN_list = perm_filter_quant_bychrom_inv_list = perm_filter_quant_bychrom_noninv_list = perm_filter_focalN_inv_list = perm_filter_pvalueN_focal_inv_list = perm_filter_focalN_noninv_list =  perm_filter_pvalueN_focal_noninv_list = list()
perm_filterN_vec = vector()

perm_filter_quant_inv_list = perm_filter_quant_noninv_list = perm_filter_pvalueN_inv_list = perm_filterN_inv_list = perm_filter_pvalueN_noninv_list = perm_filterN_noninv_list = list()

#for (i in 1:500){
for (i in 1:totalperms){
    permN = i
    load(file=paste("results/perm", permN,"_objects_v2_pvalue004.Rdata", sep=""))
    perm_filter_quant_list[[i]] = perm_filter_quant
    perm_filter_pvalueN_focal_list[[i]] = perm_filter_pvalueN_focal
    perm_filter_pvalueN_list[[i]] = perm_filter_pvalueN
    perm_filter_focalN_list[[i]] = perm_filter_focalN
    perm_filterN_vec[i] = perm_filterN
    perm_filter_quant_bychrom_list[[i]] = perm_filter_quant_bychrom
    
    perm_filter_quant_bychrom_invN_list[[i]] = perm_filter_quant_bychrom_invN
    perm_filter_quant_bychrom_noninvN_list[[i]] = perm_filter_quant_bychrom_noninvN
    perm_filter_quant_bychrom_inv_list[[i]] = perm_filter_quant_bychrom_inv
    perm_filter_quant_bychrom_noninv_list[[i]] = perm_filter_quant_bychrom_noninv
    
    perm_filter_focalN_inv_list[[i]] = perm_filter_focalN_inv
    perm_filter_pvalueN_focal_inv_list[[i]] = perm_filter_pvalueN_focal_inv
    perm_filter_focalN_noninv_list[[i]] = perm_filter_focalN_noninv
    perm_filter_pvalueN_focal_noninv_list[[i]] = perm_filter_pvalueN_focal_noninv
    perm_filter_quant_inv_list[[i]] = perm_filter_quant_inv
    perm_filter_quant_noninv_list[[i]] = perm_filter_quant_noninv
    perm_filter_pvalueN_inv_list[[i]] = perm_filter_pvalueN_inv
    perm_filterN_inv_list[[i]] = perm_filterN_inv
    perm_filter_pvalueN_noninv_list[[i]] = perm_filter_pvalueN_noninv
    perm_filterN_noninv_list[[i]] = perm_filterN_noninv
}
#myfilter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
#myfilter$ID = paste(myfilter[,1], myfilter[,2], sep="_")


# quantiles to test
quant = c(10^(-seq(from=0, to=4, by=1/4)), 0.003894319)

## focal snps
myfilter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")

# observed seasonal data
seasonal = read.table("../glm/mel_all_paired20_2sample_caF_popyear.f_s.glm", header=T, stringsAsFactors = F)
seasonal_filter = merge(na.omit(seasonal), myfilter, by=c(1,2))
seasonal_filter_inv = subset(seasonal_filter, (chrom=="2L" & pos>2225744 & pos<13154180) |
    (chrom=="2R" & pos>11278659 & pos<16163839) |
    (chrom=="3L" & pos>3173046 & pos<16301941) |
    (chrom=="3R" & pos>7576289 & pos<24857019) )
seasonal_filter_noninv = subset(seasonal_filter, (chrom=="2L" & pos<2225744 | pos>13154180) |
    (chrom=="2R" & pos<11278659 | pos>16163839) |
    (chrom=="3L" & pos<3173046 | pos>16301941) |
    (chrom=="3R" & pos<7576289 | pos>24857019) )

mychrom_inv = data.frame(chrom=c("2L","2R","3L","3R"), start=c(2225744,11278659,3173046,7576289), end=c(13154180, 16163839,16301941,24857019) )
mychrom = c("2L","2R","3L","3R")


## enrichment per chromosome
# Fold enrichment of seasonal SNPs over permutation (obs/perm)
# Using p=0.004 as cutoff (obs quantile of 0.01)
# Or using a range of p-values

seas_pvalueN = unlist(lapply(quant, FUN=function(X) mean(seasonal_filter$seas.p<X, na.rm=TRUE)))
glm_pvalue_filter_pvalueN = data.frame(do.call(cbind, perm_filter_pvalueN_list))
enrich=seas_pvalueN/glm_pvalue_filter_pvalueN
mean_enrich = apply(enrich, MARGIN=1, FUN=function(X) mean(X, na.rm=TRUE) )
median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
pvalueN_mean_obs_less_perm = data.frame(chrom="all_chroms", group="all_regions", pvalue=quant, obs=seas_pvalueN,  mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)

seas_pvalueN_inv = unlist(lapply(quant, FUN=function(X) mean(seasonal_filter_inv$seas.p<X, na.rm=TRUE)))
glm_pvalue_filter_pvalueN_inv = data.frame(do.call(cbind, perm_filter_pvalueN_inv_list))
enrich=seas_pvalueN_inv/glm_pvalue_filter_pvalueN_inv
mean_enrich = apply(enrich, MARGIN=1, FUN=function(X)mean(X, na.rm=TRUE) )
median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
pvalueN_mean_obs_less_perm_inv = data.frame(chrom="all_chroms", group="inv", pvalue=quant, obs=seas_pvalueN_inv,  mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)

seas_pvalueN_noninv = unlist(lapply(quant, FUN=function(X) mean(seasonal_filter_noninv$seas.p<X, na.rm=TRUE)))
glm_pvalue_filter_pvalueN_noninv = data.frame(do.call(cbind, perm_filter_pvalueN_noninv_list))
enrich=seas_pvalueN_noninv/glm_pvalue_filter_pvalueN_noninv
mean_enrich = apply(enrich, MARGIN=1, FUN=function(X) mean(X, na.rm=TRUE) )
median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
pvalueN_mean_obs_less_perm_noninv = data.frame(chrom="all_chroms", group="noninv", pvalue=quant, obs=seas_pvalueN_noninv,  mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)

mychrom = c("2L","2R","3L","3R")
mychrom_list = obs_perm500_enrich_chrom_v2_pvalue004_allchroms = obs_perm500_enrich_chrom_v2_pvalue004_noninv = obs_perm500_enrich_chrom_v2_pvalue004_inv = list()
for (i in 1:length(mychrom)){

    # full chrom
    list2 = list()
    for (j in 1:totalperms){
        list2[[j]] = perm_filter_pvalueN_focal_list[[j]][[i]]
    }
    
    seasonal_filter_focal = subset(seasonal_filter, chrom==mychrom[i])
    seas_pvalueN_focal = unlist(lapply(quant, FUN=function(X) mean(seasonal_filter_focal$seas.p<X, na.rm=TRUE)))
    
    perm_filter_pvalueN_focal_chrom = data.frame(do.call(cbind, list2))
    enrich=seas_pvalueN_focal/perm_filter_pvalueN_focal_chrom
    mean_enrich = apply(enrich, MARGIN=1, FUN=function(X) mean(X, na.rm=TRUE) )
    median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
    lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
    upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
    lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
    upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
    mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
    pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
    mychrom_list[[i]] = data.frame(chrom=mychrom[i], group="all_regions", pvalue=quant, obs=seas_pvalueN_focal,  mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)
    obs_perm500_enrich_chrom_v2_pvalue004_allchroms[[i]] = enrich
    
    # inv
    list2 = list()
    for (j in 1:totalperms){
        list2[[j]] = perm_filter_pvalueN_focal_inv_list[[j]][[i]]
    }
    seasonal_filter_inv_focal = subset(seasonal_filter_inv, chrom==mychrom[i])
    seas_pvalueN_inv_focal = unlist(lapply(quant, FUN=function(X) mean(seasonal_filter_inv_focal$seas.p<X, na.rm=TRUE)))
    perm_filter_pvalueN_focal_inv_chrom = data.frame(do.call(cbind, list2))
    enrich=seas_pvalueN_inv_focal/perm_filter_pvalueN_focal_inv_chrom
    mean_enrich = apply(enrich, MARGIN=1, FUN=function(X) mean(X, na.rm=TRUE) )
    median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
    lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
    upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
    lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
    upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
    mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
    pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
    mychrom_list[[i+4]] = data.frame(chrom=mychrom[i], group="inv", pvalue=quant, obs=seas_pvalueN_inv_focal,  mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)
    obs_perm500_enrich_chrom_v2_pvalue004_noninv[[i]] = enrich
    
    # noninv
    focal_chrom_noninv = mychrom_inv[i,]
    list2 = list()
    for (j in 1:totalperms){
        list2[[j]] = perm_filter_pvalueN_focal_noninv_list[[j]][[i]]
    }
    seasonal_filter_noninv_focal = subset(seasonal_filter_noninv, chrom==mychrom[i])
    seas_pvalueN_noninv_focal = unlist(lapply(quant, FUN=function(X) mean(seasonal_filter_noninv_focal$seas.p<X, na.rm=TRUE)))
  
    perm_filter_pvalueN_focal_noninv_chrom = data.frame(do.call(cbind, list2))
    enrich=seas_pvalueN_noninv_focal/perm_filter_pvalueN_focal_noninv_chrom
    mean_enrich = apply(enrich, MARGIN=1, FUN=function(X) mean(X, na.rm=TRUE) )
    median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
    lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
    upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
    lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
    upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
    mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
    pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
    mychrom_list[[i+8]] = data.frame(chrom=mychrom[i], group="noninv", pvalue=quant, obs=seas_pvalueN_noninv_focal,  mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)
    obs_perm500_enrich_chrom_v2_pvalue004_inv[[i]] = enrich

}

mychrom_list[[13]] = pvalueN_mean_obs_less_perm
mychrom_list[[14]] = pvalueN_mean_obs_less_perm_inv
mychrom_list[[15]] = pvalueN_mean_obs_less_perm_noninv
pvalueN_mean_obs_less_perm_chrom = do.call(rbind, mychrom_list)
write.table(pvalueN_mean_obs_less_perm_chrom, file="results/pvalueN_mean_obs_less_perm_chrom_v2_pvalue004.txt", quote=F, col.names=T, row.names=F)

save(obs_perm500_enrich_chrom_v2_pvalue004_allchroms, obs_perm500_enrich_chrom_v2_pvalue004_noninv, obs_perm500_enrich_chrom_v2_pvalue004_inv, file="results/obs_perm500_enrich_chrom_v2_pvalue004.Rdata")
