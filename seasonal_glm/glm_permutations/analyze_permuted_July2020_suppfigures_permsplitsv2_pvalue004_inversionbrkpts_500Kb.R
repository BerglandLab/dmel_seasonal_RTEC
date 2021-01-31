## Nov 2020
#library(ggplot2)

## Regions to run this analysis on:
# Total (three regions): 1MB +/- each inversion breakpoint, interior to 1MB +/- breakpoints, exterior to 1MB +/- breakpoints
# Per inversion (3*6 = 18 regions): 1MB +/- each inversion breakpoint, interior to 1MB +/- breakpoints, exterior to 1MB +/- breakpoints
# Inversion Chromosome Proximal   Distal
# 1   In(2L)t         2L  2225744 13154180
# 2  In(2R)NS         2R 11278659 16163839
# 3   In(3L)P         3L  3173046 16301941
# 4   In(3R)K         3R  7576289 21966092
# 5  In(3R)Mo         3R 17232639 24857019
# 6   In(3R)P         3R 12257931 20569732

## Re-worked the second part to calculate enrichment at obs quantiles (not p-values)
args = commandArgs(trailingOnly=TRUE)
permN = args[1]
#print(permN)

## read in permutation and filter
perm = na.omit(read.table(paste("results/permute",permN,"_mel_all_paired20_2sample_caF_popyear.glm", sep=""), header=TRUE, stringsAsFactors=FALSE))
#glm_pvalue = read.table("perm500_pvalue.txt", stringsAsFactors=F, header=T, sep="\t")
myfilter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
myfilter$ID = paste(myfilter[,1], myfilter[,2], sep="_")
perm$ID =  paste(perm[,1], perm[,2], sep="_")
myinds = which(perm$ID %in% myfilter$ID)
perm_filter = perm[myinds,]

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

mychrom_inv_pre_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_prebp", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=rep(0, times=6), end=c(2225744,11278659,3173046,7576289,17232639,12257931)-0.5*10^6 )
mychrom_inv_post_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_postbp", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(13154180,16163839,16301941,21966092,24857019,20569732)+0.5*10^6, end=rep(10^9, times=6) )

mychrom_inv_interior_1MB = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_interior", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(2225744,11278659,3173046,7576289,17232639,12257931)+0.5*10^6, end=c(13154180,16163839,16301941,21966092,24857019,20569732)-0.5*10^6 )

list(rbind(mychrom_inv_bkrpt1_1MB[1,], mychrom_inv_bkrpt2_1MB[1,]), rbind(mychrom_inv_bkrpt1_1MB[2,], mychrom_inv_bkrpt2_1MB[2,]), rbind(mychrom_inv_bkrpt1_1MB[3,], mychrom_inv_bkrpt2_1MB[3,]), rbind(mychrom_inv_bkrpt1_1MB[4,], mychrom_inv_bkrpt2_1MB[4,]), rbind(mychrom_inv_bkrpt1_1MB[5,], mychrom_inv_bkrpt2_1MB[5,]), rbind(mychrom_inv_bkrpt1_1MB[6,], mychrom_inv_bkrpt2_1MB[6,]) )

list(rbind(mychrom_inv_pre_1MB[1,], mychrom_inv_post_1MB[1,]), rbind(mychrom_inv_pre_1MB[2,], mychrom_inv_post_1MB[2,]), rbind(mychrom_inv_pre_1MB[3,], mychrom_inv_post_1MB[3,]), rbind(mychrom_inv_pre_1MB[4,], mychrom_inv_post_1MB[4,]), rbind(mychrom_inv_pre_1MB[5,], mychrom_inv_post_1MB[5,]), rbind(mychrom_inv_pre_1MB[6,], mychrom_inv_post_1MB[6,]) )

list(mychrom_inv_interior_1MB[1,], mychrom_inv_interior_1MB[2,], mychrom_inv_interior_1MB[3,], mychrom_inv_interior_1MB[4,], mychrom_inv_interior_1MB[5,], mychrom_inv_interior_1MB[6,])



## the three cases 3R inversions together
threeRbkpt = rbind(data.frame(region=paste(c("In(3R)K","In(3R)Mo", "In(3R)P"), "_bp1", sep=""), chrom=c("3R","3R","3R"), start=c(7576289,17232639,12257931)-0.5*10^6, end=c(7576289,17232639,12257931)+0.5*10^6 ),
    data.frame(region=paste(c("In(3R)K","In(3R)Mo", "In(3R)P"), "_bp2", sep=""), chrom=c("3R","3R","3R"), start=c(21966092,24857019,20569732)-0.5*10^6, end=c(21966092,24857019,20569732)+0.5*10^6 ) )

threeRext = rbind( data.frame(region=paste(c("In(3R)"), "_prebp", sep=""), chrom=c("3R"), start=rep(0, times=1), end=c(7576289)-0.5*10^6 ),
    data.frame(region=paste(c("In(3R)"), "_postbp", sep=""), chrom=c("3R"), start=c(24857019)+0.5*10^6, end=rep(10^9, times=1) ) )

threeRint = data.frame(region=paste(c("In(3R)","In(3R)", "In(3R)", "In(3R)", "In(3R)"), "_intbp", sep=""), chrom=c("3R","3R","3R","3R","3R"),
    start=c(7576289,12257931,17232639,20569732,21966092)+0.5*10^6, end=c(12257931,17232639,20569732,21966092,24857019)-0.5*10^6 )



## the three cases with all combined
allbkpt = rbind(mychrom_inv_bkrpt1_1MB, mychrom_inv_bkrpt2_1MB)

allext = rbind( data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)"), "_prebp", sep=""), chrom=c("2L","2R","3L","3R"), start=rep(0, times=4), end=c(2225744,11278659,3173046,7576289)-0.5*10^6 ),
    data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)"), "_postbp", sep=""), chrom=c("2L","2R","3L","3R"), start=c(13154180,16163839,16301941,24857019)+0.5*10^6, end=rep(10^9, times=4) ) )

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
    


## Running the function to calculate the stats per region
# N_perpvalue_perregion(glm=perm_filter, regionDF=mychrom_inv)
results_N_perpvalue_perregion = lapply(myregions, FUN=function(X) N_perpvalue_perregion(glm=perm_filter, regionDF=X) )
names(results_N_perpvalue_perregion) = names(results_N_perpvalue_perregion)


save(results_N_perpvalue_perregion, file=paste("results/perm", permN,"_pvalue004_inversionbrkpts_500Kb.Rdata", sep=""))

 

