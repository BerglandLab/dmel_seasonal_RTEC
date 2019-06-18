## July 2017
# script to calculate the average fit of sampled (per genomic block) datasets to the simulated datasets.
# INPUT: "enrich" files for sampled datasets, simualted datasets
# OUTPUT: mean, median, and sd (of 100 samples) SSE fit to simulated datasets
# USAGE: Rscript average_samples_simulated_observed_realdata_20pop_dist2_July2017.R sample_prefix Nsamples simulated_dataset_dir simulated_dataset_list outprefix
# eg: Rscript average_samples_simulated_observed_realdata_20pop_dist2_July2017.R /Users/heathermachado/nescent_melCA/fishers_method_Oct2016/jamie_fisher_method/results/enrich_L_rank_fisher_exactJ.merged.20pop.sample1K.txt 2 inlist_rep1_pop20.txt 20pop_sample1K
# eg: Rscript ../../scripts/average_samples_simulated_observed_realdata_20pop_dist2_July2017.R results/enrich_L_rank_fisher_exactJ.merged.20pop.sample1K.txt 2 /Users/heathermachado/nescent_melCA/simulated_data/one_script_from_realdata inlist_rep1_pop20.txt 20pop_sample1K


#install.packages("gplots", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
#install.packages("RColorBrewer", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
library(gplots, lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
library(RColorBrewer, lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")

### real data
#setwd("/Users/heathermachado/nescent_melCA/simulated_data/one_script_from_realdata")
#args = c("/Users/heathermachado/nescent_melCA/fishers_method_Oct2016/jamie_fisher_method/results/enrich_L_rank_fisher_exactJ.merged.20pop.chrom_pos.sample1K.txt", 100, "/Users/heathermachado/nescent_melCA/simulated_data/one_script_from_realdata/", "inlist_rep1_pop20.txt", "20pop_sample1K")
args <- commandArgs(trailingOnly = TRUE)
samplefile = args[1]
sampleN = args[2]
simdir = args[3]
simlist = args[4]
outpre = args[5]

samplelist = list()

for (i in 1:sampleN){
  samplelist[[i]] = read.table(paste(samplefile, i, sep=""), stringsAsFactors=FALSE, header=TRUE)
}

## calculate deviation of each simulated dataset
ss2 = function(enrichall, slist, npops){
  ss_all = data.frame(inlist, ss=NA)
  obsP = enrichall[(npops+1):max(which(!is.na(enrichall[,12]))),12]/enrichall[1,10]
  
  for (i in 1:length(slist)){
    focal = slist[[i]]
    focalP = focal[(npops+1):max(which(!is.na(focal[,12]))),12]/focal[1,10]
    
    if (length(focalP)>length(obsP)){
      obsPnew = rep(0, times=length(focalP))
      obsPnew[1:length(obsP)] = obsP
      obsP = obsPnew
    }
    
    if (length(focalP)<length(obsP)){
      focalPnew = rep(0, times=length(obsP))
      focalPnew[1:length(focalP)] = focalP
      focalP = focalPnew
    }
    
    ss = sum((focalP-obsP)^2)
    ss_all[i,2] = ss
  }
  ss_all  
}

###### distance between simulated vs obs
#ls results/enrich_L_rank_fisher_exactJ.Nsites971520.prop0.0*Ppop1* > inlist_rep1.txt
inlist = read.table(paste(simdir, "/", simlist, sep=""), stringsAsFactors=FALSE)[,1]
slist = list()
for (i in 1:length(inlist)){
  slist[[i]] = read.table(paste(simdir, "/", inlist[i], sep=""), stringsAsFactors=FALSE, header=TRUE)
}

ss_all_list = list()
for (i in 1:sampleN){
  ss_all_list[[i]] = ss2(samplelist[[i]], slist, 20)
}

allSS = matrix(ncol=as.numeric(sampleN), nrow=length(slist))
for (i in 1:sampleN){
  allSS[,i] = ss_all_list[[i]][,2]
}

# Have to normalize each run separately, or can't compare.
# Z normalize?
my.scale = function(x){
  xnew1 = scale(x)
  xnew2 = xnew1 - min(xnew1)
  xnew2
}

allSSscale = apply(allSS, MARGIN=2, my.scale)


mean1 = apply(allSS, MARGIN=1, FUN=mean, na.rm=TRUE)
med1 = apply(allSS, MARGIN=1, FUN=median, na.rm=TRUE)
sd1 = apply(allSS, MARGIN=1, FUN=mean, na.rm=TRUE)
meanS = apply(allSSscale, MARGIN=1, FUN=mean, na.rm=TRUE)
medS = apply(allSSscale, MARGIN=1, FUN=median, na.rm=TRUE)
sdS = apply(allSSscale, MARGIN=1, FUN=mean, na.rm=TRUE)
ss_all = data.frame(sim=ss_all_list[[1]][,1], mean=mean1, median=med1, sd=sd1, meanS=meanS, medianS=medS, sdS=sdS)


Ppop = rep(c(1), times=77)
prop = c(rep(0.001, times=11),rep(0.005, times=11),rep(0.01, times=11),rep(0.02, times=11),rep(0.03, times=11),rep(0.04, times=11),rep(0.05, times=11))
s = rep( c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), times=7) 

ss_all$Ppop = Ppop
ss_all$prop = prop
ss_all$s = s
ss_all2 = ss_all[order(ss_all$medianS),]
write.table(ss_all2, file=paste("SSE_table_", outpre,".txt", sep=""), col.names=TRUE, row.names=FALSE, quote=FALSE )

m_sprop = matrix(as.numeric(ss_all$medianS), ncol=11, byrow=TRUE)
rownames(m_sprop) = c(0.001,0.005,0.01,0.02,0.03,0.04,0.05)
colnames(m_sprop) = c(0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

library(gplots)
library(RColorBrewer)
mycols = colorRampPalette(brewer.pal(6,"Blues"))(20)

pdf(paste("heatmap_s_prop_realdata_Pprop1_dist2_normalized_", outpre, ".pdf", sep=""), width=7, height=6)
h1 =heatmap.2(-log(m_sprop), Rowv=NA, labRow=rownames(m_sprop), Colv=NA, dendrogram="none", trace="none", labCol=colnames(m_sprop),xlab="Seasonal selection coefficient", ylab="Proportion of sites", main="", col=mycols, notecol="black",density.info="none",key.xlab="Dist from observed",cexRow=0.9, cexCol=0.9, scale="none") #,key=FALSE)
dev.off()

## Now separate by chromosome: NO DONE AT grt/les stage
# inversions = data.frame(
#   Inversion=c("In(2L)t","In(2R)NS","In(3L)P","In(3R)K,Mo,P"), 
#   Chromosome=c("2L","2R","3L","3R"),
#   Proximal=c(2225744,11278659,3173046,7576289),
#   Distal=c(13154180,16163839,16301941,24857019) ) 
# write.table(inversions, file="inverted_regions_dmel.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
