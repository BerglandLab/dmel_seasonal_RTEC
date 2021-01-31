## Nov 2020
#library(ggplot2)

#setwd("/lustre/scratch116/casm/cgp/users/hm8/nescent_melCA/glm_permutations")

# each perm*_objects.Rdata (one for each of 500 permutations) has the following objects:
# perm_filter_quant, perm_filter_pvalueN, perm_filterN, perm_filter_focalN, perm_filter_pvalueN_focal, perm_filter_quant_bychrom, perm_filter_quant_bychrom_inv, perm_filter_quant_bychrom_noninv, perm_filter_quant_bychromN, perm_filter_quant_bychrom_invN, perm_filter_quant_bychrom_noninvN, perm_filter_focalN_inv, perm_filter_pvalueN_focal_inv, perm_filter_focalN_noninv, perm_filter_pvalueN_focal_noninv
# observed seasonal data
myfilter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
seasonal = read.table("../results/mel_all_paired20_2sample_caF_popyear.f_s.glm", header=T, stringsAsFactors = F)
seasonal_filter = merge(na.omit(seasonal), myfilter, by=c(1,2))
max_chrom = tapply(seasonal_filter$pos, INDEX=seasonal_filter$chrom, FUN=max)
#2L       2R       3L       3R
#22899773 20999969 24085366 27399947
min_chrom = tapply(seasonal_filter$pos, INDEX=seasonal_filter$chrom, FUN=min)
#2L     2R     3L     3R
#300026 703047 203260 700262

## Write function per region
# Input:
# 1) seasonal (or permutation) regression results
# 2) dataframe or vector of region, with columns c("region","start","end")
N_perpvalue_perregion = function(glm, regionDF, quant = c(10^(-seq(from=0, to=4, by=1/4)), 0.003894319) ){
        #glm = perm_filter
        #regionDF = mychrom_inv
        #quant = c(10^(-seq(from=0, to=4, by=1/4)))
        if ( !(nrow(regionDF)>0) ) {warning("nrow is not greater than 0. regionDF requires a data.frame with nrow >= 1")}
        if ( any( !(c("region","chrom","start","end") %in% colnames(regionDF)) ) ) {warning("regionDF requires the column headers: region, chrom, start, end" )}
        subregions_list = list()
        for (i in 1:nrow(regionDF)){
            subregions_list[[i]] = subset(glm, chrom==regionDF$chrom[i] & pos>regionDF$start[i] & pos<regionDF$end[i])
        }
        allregions1 = data.frame(do.call(rbind, subregions_list))
        allregions = allregions1[!(duplicated(allregions1)),]
        N_perpvalue = unlist(lapply(quant, FUN=function(X) mean(allregions$seas.p<X, na.rm=TRUE)))
        N = nrow(allregions)
        myregion = regionDF
        myout = list(N_perpvalue=N_perpvalue, N=N, region=regionDF, pvalue=quant)
        return(myout)
}


# delineation of major inversions (combined the three on 3R)
#mychrom_inv = data.frame(region=c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(2225744,11278659,3173046,7576289,17232639,12257931), end=c(13154180,16163839,16301941,21966092,24857019,20569732) )

mychrom_inv_bkrpt1_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_bp1", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(2225744,11278659,3173046,7576289,17232639,12257931)-0.5*10^6, end=c(2225744,11278659,3173046,7576289,17232639,12257931)+0.5*10^6 )
mychrom_inv_bkrpt2_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_bp2", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(13154180,16163839,16301941,21966092,24857019,20569732)-0.5*10^6, end=c(13154180,16163839,16301941,21966092,24857019,20569732)+0.5*10^6 )

mychrom_inv_pre_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_prebp", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=min_chrom[c(1:4,4,4)], end=c(2225744,11278659,3173046,7576289,17232639,12257931)-0.5*10^6 )
mychrom_inv_post_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_postbp", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(13154180,16163839,16301941,21966092,24857019,20569732)+0.5*10^6, end=max_chrom[c(1:4,4,4)] )

mychrom_inv_interior_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_interior", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(2225744,11278659,3173046,7576289,17232639,12257931)+0.5*10^6, end=c(13154180,16163839,16301941,21966092,24857019,20569732)-0.5*10^6 )

## the three cases 3R inversions together
threeRbkpt = rbind(data.frame(region=paste(c("In(3R)K","In(3R)Mo", "In(3R)P"), "_bp1", sep=""), chrom=c("3R","3R","3R"), start=c(7576289,17232639,12257931)-0.5*10^6, end=c(7576289,17232639,12257931)+0.5*10^6 ),
    data.frame(region=paste(c("In(3R)K","In(3R)Mo", "In(3R)P"), "_bp2", sep=""), chrom=c("3R","3R","3R"), start=c(21966092,24857019,20569732)-0.5*10^6, end=c(21966092,24857019,20569732)+0.5*10^6 ) )

threeRext = rbind( data.frame(region=paste(c("In(3R)"), "_prebp", sep=""), chrom=c("3R"), start=min_chrom[c(4)], end=c(7576289)-0.5*10^6 ),
    data.frame(region=paste(c("In(3R)"), "_postbp", sep=""), chrom=c("3R"), start=c(24857019)+0.5*10^6, end=max_chrom[c(4)] ) )

threeRint = data.frame(region=paste(c("In(3R)","In(3R)", "In(3R)", "In(3R)", "In(3R)"), "_intbp", sep=""), chrom=c("3R","3R","3R","3R","3R"),
    start=c(7576289,12257931,17232639,20569732,21966092)+0.5*10^6, end=c(12257931,17232639,20569732,21966092,24857019)-0.5*10^6 )



## the three cases with all combined
allbkpt = rbind(mychrom_inv_bkrpt1_1MB, mychrom_inv_bkrpt2_1MB)

allext = rbind( data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)"), "_prebp", sep=""), chrom=c("2L","2R","3L","3R"), start=min_chrom, end=c(2225744,11278659,3173046,7576289)-0.5*10^6 ),
    data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)"), "_postbp", sep=""), chrom=c("2L","2R","3L","3R"), start=c(13154180,16163839,16301941,24857019)+0.5*10^6, end=max_chrom ) )

allint = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)","In(3R)", "In(3R)", "In(3R)", "In(3R)"), "_intbp", sep=""), chrom=c("2L","2R","3L","3R","3R","3R","3R","3R"),
    start=c(2225744,11278659,3173046,7576289,12257931,17232639,20569732,21966092)+0.5*10^6, end=c(13154180,16163839,16301941,12257931,17232639,20569732,21966092,24857019)-0.5*10^6 )

myregions = list(
    rbind(mychrom_inv_bkrpt1_1MB[1,], mychrom_inv_bkrpt2_1MB[1,]), rbind(mychrom_inv_bkrpt1_1MB[2,], mychrom_inv_bkrpt2_1MB[2,]), rbind(mychrom_inv_bkrpt1_1MB[3,], mychrom_inv_bkrpt2_1MB[3,]), rbind(mychrom_inv_bkrpt1_1MB[4,], mychrom_inv_bkrpt2_1MB[4,]), rbind(mychrom_inv_bkrpt1_1MB[5,], mychrom_inv_bkrpt2_1MB[5,]), rbind(mychrom_inv_bkrpt1_1MB[6,], mychrom_inv_bkrpt2_1MB[6,]),
    rbind(mychrom_inv_pre_1MB[1,], mychrom_inv_post_1MB[1,]), rbind(mychrom_inv_pre_1MB[2,], mychrom_inv_post_1MB[2,]), rbind(mychrom_inv_pre_1MB[3,], mychrom_inv_post_1MB[3,]), rbind(mychrom_inv_pre_1MB[4,], mychrom_inv_post_1MB[4,]), rbind(mychrom_inv_pre_1MB[5,], mychrom_inv_post_1MB[5,]), rbind(mychrom_inv_pre_1MB[6,], mychrom_inv_post_1MB[6,]),
    mychrom_inv_interior_1MB[1,], mychrom_inv_interior_1MB[2,], mychrom_inv_interior_1MB[3,], mychrom_inv_interior_1MB[4,], mychrom_inv_interior_1MB[5,], mychrom_inv_interior_1MB[6,],
    threeRbkpt, threeRext, threeRint,
    allbkpt, allext, allint
)

names(myregions) = c(
    paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "bkpts", sep="_"),
    paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "exterior", sep="_"),
    paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "interior", sep="_"),
    paste("In(3R)KMoP", c("bkpts", "exterior", "interior"), sep="_"),
    paste("all", c("bkpts", "exterior", "interior"), sep="_") )
myregions_sizes = data.frame(region=names(myregions), length_bp=unlist(lapply(myregions, FUN=function(X) sum(X$end-X$start) )) )
myinvs = unlist(lapply(myregions_sizes$region, FUN=function(X) unlist(strsplit(as.character(X), split="_", fixed=TRUE))[1]))
myclass = unlist(lapply(myregions_sizes$region, FUN=function(X) unlist(strsplit(as.character(X), split="_", fixed=TRUE))[2]))
myregions_sizes$inv = myinvs
myregions_sizes$class = myclass
write.table(myregions_sizes, file="myregions_sizes.txt", col.names=T, row.names=F, quote=F)

## Running the function to calculate the stats per region
# N_perpvalue_perregion(glm=perm_filter, regionDF=mychrom_inv)
obs_results_N_perpvalue_perregion = lapply(myregions, FUN=function(X) N_perpvalue_perregion(glm=seasonal_filter, regionDF=X) )
names(obs_results_N_perpvalue_perregion) = names(myregions)


## for 500 permutations
totalperms = 500
#totalperms = 10


perm_filter_pvalueN = list()
perm_filter_N = list()
for (j in 1:24){
    perm_filter_pvalueN[[j]] = list()
    perm_filter_N[[j]] = list()
}

#for (i in 1:500){
for (i in 1:totalperms){
    permN = i
    load(file=paste("results/perm", permN,"_pvalue004_inversionbrkpts_500Kb.Rdata", sep=""))
    for (j in 1:24){
        perm_filter_pvalueN[[j]][[i]] = results_N_perpvalue_perregion[[j]]$N_perpvalue
        perm_filter_N[[j]][[i]] = results_N_perpvalue_perregion[[j]]$N
    }
}
perm_filter_pvalueN_byregion = lapply(perm_filter_pvalueN, FUN=function(X) do.call(rbind, X) )
#perm_filter_N_byregion = lapply(perm_filter_N, FUN=function(X) do.call(rbind, X) )

# quantiles to test
quant = c(10^(-seq(from=0, to=4, by=1/4)), 0.003894319)

## comparing with observed
#myinvs = c("In(2L)t_bkpts",     "In(2R)NS_bkpts",    "In(3L)P_bkpts",
#"In(3R)K_bkpts" ,    "In(3R)Mo_bkpts" ,   "In(3R)P_bkpts",
#"In(2L)t_exterior" , "In(2R)NS_exterior" ,"In(3L)P_exterior",
#"In(3R)K_exterior" , "In(3R)Mo_exterior", "In(3R)P_exterior",
#"In(2L)t_interior" , "In(2R)NS_interior", "In(3L)P_interior",
#"In(3R)K_interior" , "In(3R)Mo_interior", "In(3R)P_interior",
#"In(3R)KMoP_bkpts"    "In(3R)KMoP_exterior" "In(3R)KMoP_interior",
#"all_bkpts"    ,     "all_exterior"  ,    "all_interior")

myinvs = unlist(lapply(names(obs_results_N_perpvalue_perregion), FUN=function(X) unlist(strsplit(X, split="_", fixed=TRUE))[1]))
myclass = unlist(lapply(names(obs_results_N_perpvalue_perregion), FUN=function(X) unlist(strsplit(X, split="_", fixed=TRUE))[2]))
obs_perm_filter_pvalueN_byregion_list = list()
for (j in 1:24){
    focal_obs = obs_results_N_perpvalue_perregion[[j]]$N_perpvalue
    focal_N = obs_results_N_perpvalue_perregion[[j]]$N
    focal_perm = perm_filter_pvalueN_byregion[[j]]
    enrich=focal_obs/t(focal_perm)
    mean_enrich = apply(enrich, MARGIN=1, FUN=function(X) mean(X, na.rm=TRUE) )
    median_enrich = apply(enrich, MARGIN=1, FUN=function(X) median(X, na.rm=TRUE) )
    lower95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.025) )
    upper95_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.975) )
    lower90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.05) )
    upper90_enrich = apply(enrich, MARGIN=1, FUN=function(X) quantile(X, probs=0.95) )
    mean_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) sum(X>1) )
    pvalue_obs_grt_perm = apply(enrich, MARGIN=1, FUN=function(X) (sum(X<=1)+1)/501 )
    obs_perm_filter_pvalueN_byregion_list[[j]] = data.frame(region=names(obs_results_N_perpvalue_perregion)[j], inv=myinvs[j] , class=myclass[j], pvalue=quant, obs=focal_obs, N=focal_N, mean_enrich, median_enrich, lower95_enrich, upper95_enrich, lower90_enrich, upper90_enrich, mean_obs_grt_perm, pvalue_obs_grt_perm)
}
obs_perm_filter_pvalueN_byregion = do.call(rbind, obs_perm_filter_pvalueN_byregion_list)
write.table(obs_perm_filter_pvalueN_byregion, file="results/obs_perm_filter_pvalue004N_byregion_inversionbrkpts_500Kb.txt", col.names=T, row.names=F, quote=F)


