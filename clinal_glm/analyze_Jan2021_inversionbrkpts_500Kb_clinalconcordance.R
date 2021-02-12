## Nov 2020
#library(ggplot2)

library(data.table)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)

### load in orginal data 
cline = fread("../results/mel_clinal_uniquepops_springPA_noMA.glm.noheader")
colnames(cline) = c("chrom","pos","clinal.coef","clinal.p","N")
seas = fread("../results/mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm")
filter = fread("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
setnames(filter, names(filter), c("chrom", "pos"))
setkey(filter, chrom, pos)
setkey(cline, chrom, pos)
setkey(seas, chrom, pos)
seas_clinal = merge(cline, seas)
seas_clinal = merge(seas_clinal, filter)
seas_clinal[,clinal.q:=rank(clinal.p)/(length(clinal.p)+1)]
seas_clinal[,seas.q:=rank(seas.p)/(length(seas.p)+1)]
#setkey(seas_clinal, chrom, pos)

max_chrom = tapply(seas_clinal$pos, INDEX=seas_clinal$chrom, FUN=max)
#2L       2R       3L       3R
#22899773 20999969 24085366 27399947

min_chrom = tapply(seas_clinal$pos, INDEX=seas_clinal$chrom, FUN=min)
# 2L     2R     3L     3R 
# 300026 703047 203260 700262 


## Write function per region
# Input:
# 1) seasonal (or permutation) regression results
# 2) dataframe or vector of region, with columns c("region","start","end")
allchroms = c("2L","2R","3L","3R")
concordance_perregion = function(glm, regionDF, quant = 0.05){
        #glm = seas_clinal
        #regionDF = "2L"
        #quant = c(10^(-seq(from=0, to=4, by=1/4)))
        if ( regionDF=="all" ){
            allregions = glm
            } else
        
        if ( regionDF %in% allchroms ){
            allregions = subset(glm, chrom==regionDF)
            } else
        
        if ( nrow(regionDF)>0 ){
            subregions_list = list()
            for (i in 1:nrow(regionDF)){
                subregions_list[[i]] = subset(glm, chrom==regionDF$chrom[i] & pos>regionDF$start[i] & pos<regionDF$end[i])
            }
            allregions1 = data.frame(do.call(rbind, subregions_list))
            allregions = allregions1[!(duplicated(allregions1)),]
            }
        
        focalquant = subset(allregions, clinal.q <= quant &  seas.q <= quant)
        N = nrow(focalquant)
        concord = mean( (focalquant$clinal.coef > 0 & focalquant$seas.coef > 0) |
            (focalquant$clinal.coef <= 0 & focalquant$seas.coef <= 0) )
        myregion = regionDF
        myout = list(concord=concord, N=N, region=regionDF, jointquantile=quant)
        return(myout)
}




mychrom_inv_bkrpt1 = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_bp1", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(2225744,11278659,3173046,7576289,17232639,12257931), end=c(2225744,11278659,3173046,7576289,17232639,12257931))
mychrom_inv_bkrpt2 = data.frame(region=paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "_bp2", sep=""), chrom=c("2L","2R","3L","3R","3R","3R"), start=c(13154180,16163839,16301941,21966092,24857019,20569732), end=c(13154180,16163839,16301941,21966092,24857019,20569732))
data.frame(mychrom_inv_bkrpt1[,1:3],mychrom_inv_bkrpt2[,4])


# delineation of major inversions (combined the three on 3R)
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

myregions = c(list(
    rbind(mychrom_inv_bkrpt1_1MB[1,], mychrom_inv_bkrpt2_1MB[1,]), rbind(mychrom_inv_bkrpt1_1MB[2,], mychrom_inv_bkrpt2_1MB[2,]), rbind(mychrom_inv_bkrpt1_1MB[3,], mychrom_inv_bkrpt2_1MB[3,]), rbind(mychrom_inv_bkrpt1_1MB[4,], mychrom_inv_bkrpt2_1MB[4,]), rbind(mychrom_inv_bkrpt1_1MB[5,], mychrom_inv_bkrpt2_1MB[5,]), rbind(mychrom_inv_bkrpt1_1MB[6,], mychrom_inv_bkrpt2_1MB[6,]),
    rbind(mychrom_inv_pre_1MB[1,], mychrom_inv_post_1MB[1,]), rbind(mychrom_inv_pre_1MB[2,], mychrom_inv_post_1MB[2,]), rbind(mychrom_inv_pre_1MB[3,], mychrom_inv_post_1MB[3,]), rbind(mychrom_inv_pre_1MB[4,], mychrom_inv_post_1MB[4,]), rbind(mychrom_inv_pre_1MB[5,], mychrom_inv_post_1MB[5,]), rbind(mychrom_inv_pre_1MB[6,], mychrom_inv_post_1MB[6,]),
    mychrom_inv_interior_1MB[1,], mychrom_inv_interior_1MB[2,], mychrom_inv_interior_1MB[3,], mychrom_inv_interior_1MB[4,], mychrom_inv_interior_1MB[5,], mychrom_inv_interior_1MB[6,],
    threeRbkpt, threeRext, threeRint,
    allbkpt, allext, allint, "all"), list("2L","2R","3L","3R")
)

names(myregions) = c(
    paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "bkpts", sep="_"),
    paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "exterior", sep="_"),
    paste(c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K","In(3R)Mo", "In(3R)P"), "interior", sep="_"),
    paste("In(3R)KMoP", c("bkpts", "exterior", "interior"), sep="_"),
    paste("all", c("bkpts", "exterior", "interior"), sep="_"), "all_all", paste(allchroms, "all", sep="_") )

## Running the function to calculate the stats per region
# N_perpvalue_perregion(glm=perm_filter, regionDF=mychrom_inv)
focalquants = c(0.01,0.05,0.1)
concord_res_list = list()
for (i in 1:length(focalquants)){
    obs_results_N_perpvalue_perregion = lapply(myregions, FUN=function(X) concordance_perregion(glm=seas_clinal, regionDF=X, quant=focalquants[i]) )
    names(obs_results_N_perpvalue_perregion) = names(myregions)
    concord_res_list[[i]] = obs_results_N_perpvalue_perregion
}
names(concord_res_list) = focalquants
#save(concord_res_list, file="concord_results_list.Rdata")




######### Plotting
reslist = list()
for (i in 1:3){
    myN = myC = c()
    for (j in 1:length(concord_res_list[[i]])){
        myC[j] = concord_res_list[[i]][[j]]$concord 
        myN[j] = concord_res_list[[i]][[j]]$N
    }
    reslist[[i]] = data.frame(region=names(concord_res_list[[i]]), concord = myC, N=myN, jointquantile=names(concord_res_list)[i])
}
concordDF = do.call(rbind, reslist)
concordDF$se = sqrt(concordDF$concord*(1-concordDF$concord)/concordDF$N)
myinvs = unlist(lapply(concordDF$region, FUN=function(X) unlist(strsplit(as.character(X), split="_", fixed=TRUE))[1]))
myclass = unlist(lapply(concordDF$region, FUN=function(X) unlist(strsplit(as.character(X), split="_", fixed=TRUE))[2]))
mychrom = revalue(myinvs, c("In(2L)t" = "2L", "In(2R)NS" = "2R", "In(3L)P" = "3L", "In(3R)K" = "3R", "In(3R)KMoP" = "3R", "In(3R)Mo" = "3R", "In(3R)P" = "3R") )
concordDF$inv = myinvs
concordDF$class = myclass
concordDF$chrom = mychrom
concordDF$lower95 = concordDF$concord - 1.96*concordDF$se
concordDF$upper95 = concordDF$concord + 1.96*concordDF$se

5
x = concordDF[i,2]
n = concordDF[i,3]
tryCatch(binom.test(x=round(x*n), n=n, p = 0.5,
                    alternative = c("greater"),
                    conf.level = 0.95), error=function(err) NA)

binomFUN = function(x, n){
    myres = tryCatch(binom.test(x=round(x*n), n=n, p = 0.5,
                        alternative = c("greater"),
                        conf.level = 0.95), error=function(err) NA)
    if (is.na(myres[1])){
        return(rep(NA, times=3))
    } else {
        
    }
    return(unlist(myres[c("p.value","conf.int")]))
}
binomALL = t(apply(concordDF, MARGIN=1, FUN=function(X) binomFUN(as.numeric(X[2]), as.numeric(X[3])) ))
colnames(binomALL) = c("binom_p","binom_lowerCI","binom_upperCI")
concordDF = data.frame(concordDF, binomALL)
concordDF$pvaluelabel = ""
concordDF$pvaluelabel[concordDF$binom_p<0.05] = "*"
write.table(concordDF, file="concordDF_Jan2021.txt", col.names=T, row.names=F, quote=F, sep="\t")
