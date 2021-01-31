## Oct 2016

## the most significant SNPs from the glm
# plot genomic location
# check distance between sites

glm = read.table("results/mel_all_paired20_2sample_caF_popyear.f_s.glm", stringsAsFactors=FALSE, header=TRUE) # PRE-FILTERED
filter = read.table("data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
glm2A = na.omit(merge(glm, filter, by=c(1,2)))
glm2 = glm2A[glm2A[,4]<1 & glm2A[,4]>0  & glm2A[,1]!="X",]
# 774651      6


##################   SNPs by genomic location
library(caTools)
chroms=c("2L","2R","3L","3R")

############# checking to see how close (in bp) the nearest sig snps are
nearest = function(focal){
  require(data.table)
  chroms=unique(focal[,1])
  bychrom=list()
  for (i in 1:length(chroms)){
    focal2 = focal[focal[,1]==chroms[i],]
    focal2$nearest = NA
    for (j in 1:nrow(focal2)){
      nearest = min(abs(as.numeric(focal2[j,2])-as.numeric(focal2[-j,2])))
      focal2$nearest[j] = nearest
    }
    bychrom[[i]] = focal2
  }  
  data.frame(rbindlist(bychrom))
}
focal = glm2[ glm2[,4] < 0.00001, ]
allNearest = nearest(focal)
sum(allNearest$nearest<1000)
dim(glm2)
min(allNearest$nearest)
allNearest[allNearest$nearest<5,]


############### histogram of P-value by bin (compared to permutations)
head(glm2)
hist(glm2[,4])
make_bins = function(x, size){
  x = na.omit(x)
  my_seq = seq(from=0, to=1, by=size)
  my_sum=vector()
  for (i in 1:(length(my_seq)-1) ){
    my_sum[i] = sum(x>=my_seq[i] & x<my_seq[i+1])
  }
  list(my_sum, my_seq)
}
  
mybins = make_bins(glm2[,4], 0.001)
trimmed = mybins[[1]][1:(length(mybins[[1]])-1)]

perms = read.table("results/permutation_prop_by_bin001_polymorphic.txt", header=TRUE)
perms_trim = perms[,1:(ncol(perms)-1)]
perms_trim2 = t(apply(perms_trim, MARGIN=1, FUN=function(X) X/sum(X)))
perms_trim2_med = apply(perms_trim2, MARGIN=2, FUN=median)
q01 = quantile(glm2[,4], probs=0.01)


trimmed_prop = trimmed/sum(trimmed)
obs_propP01 = sum(trimmed_prop[1:11])
perm_propP01 = apply(perms_trim2[,1:11], MARGIN=1, FUN=function(X) sum(X) )
mean(obs_propP01<perm_propP01) # 0.03

obs_propP001 = sum(trimmed_prop[1:21])
perm_propP001 = apply(perms_trim2[,1:21], MARGIN=1, FUN=function(X) sum(X) )
mean(obs_propP001<perm_propP001) # 0.02


pdf("glm_paired20_vs_100perm_hist_polymorphic.pdf", width=4,height=4)
#par(mar=c(5,5,2,2))
plot(mybins[[2]][2:(length(mybins[[2]])-1)], trimmed/sum(trimmed), type="l", log="x", xlab="P-value bin", ylab="Proportion SNPs", ylim=c(0.0008,0.0035), cex.axis=0.8)
for (i in 1:100){
  points(mybins[[2]][2:(length(mybins[[2]])-1)], perms_trim2[i,], type="l", col=grey(0.5, 0.5), log="x")
}
legend("topright", legend=c("Observed","Permuted","Permuted (med)"), col=c("red","grey","black"), bty="n", lty=1)
points(mybins[[2]][2:(length(mybins[[2]])-1)], trimmed/sum(trimmed), type="l", log="x", col="red", lwd=2)
points(mybins[[2]][2:(length(mybins[[2]])-1)], perms_trim2_med, type="l", col="black", log="x", lwd=2)
text(x=0.015, y=0.0009, labels="Quantile=0.01", cex=0.8)
abline(v=q01, lty=2, col="black")
dev.off()

glm_df = data.frame(pvalue=mybins[[2]][2:(length(mybins[[2]])-1)], prop=trimmed/sum(trimmed))
perm_df = data.frame(pvalue=mybins[[2]][2:(length(mybins[[2]])-1)], prop=t(perms_trim2) )
perm_df2 = melt(perm_df, id="pvalue")
colnames(perm_df2) = c("pvalue","perm","prop")
glm_df_orig = glm_df

pdf("glm_paired20_vs_100perm_hist_polymorphic_ggplot.pdf", width=5,height=4)
glm_hist =
ggplot(perm_df2, aes(pvalue, prop, group=perm)) +
  geom_line(aes(color="Permuted",size="Permuted"), alpha=0.1) +
  geom_line(data=glm_df, aes(pvalue, prop, group=NULL, color="Observed", size="Observed")) +
  scale_color_manual('',values=c('Permuted'='grey','Observed'='red')) + 
  scale_size_manual('',values=c('Permuted'=0.5,'Observed'=1)) + 
  scale_x_continuous(trans="log10", breaks=c(10^c(-3, -2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(0.001, .01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept=0.001, lty=2) +
  xlab("Seasonal P-value") + 
  ylab("Proportion of SNPs") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.5), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.title.align=.5,
        axis.text=element_text(size=8),
        axis.title=element_text(size=10))                     
dev.off()



####### spring_fall mean allele frequency change for all data compared with Q<0.01
freqdiff = read.table("../data/mel_meanfreqdiff_paired20_042016.txt", header=TRUE, stringsAsFactors=FALSE)
glm2F = merge(glm2, freqdiff, by=c(1,2))
Q01 = quantile(glm2F$seas.p, probs=0.01, na.rm=TRUE)

# h1 = hist(abs(glm2F$meanf_s[glm2F$seas.p>Q01]), plot=FALSE, breaks=50)
# h1$density = h1$counts/sum(h1$counts)
# h2 = hist(abs(glm2F$meanf_s[glm2F$seas.p<Q01]), plot=FALSE, breaks=h1$breaks)
# h2$density = h2$counts/sum(h2$counts)
# 
# pdf("glm_seasonalQ01_meanspAF_hist_polymorphic.pdf", width=4, height=4)
# plot(h1, freq=FALSE, col="grey", xlim=c(0,0.15), main="", xlab=expression(paste("Mean spring/fall ",Delta,"AF")), ylab="Proportion of SNPs", ylim=c(0,0.45))
# plot(h2, freq=FALSE, add=TRUE, col=rgb(1,0,0,0.2))
# legend("topright",bty="n",legend=c("Seasonal (Q<0.01)","Non-seasonal"), fill=c(rgb(1,0,0,0.2), "grey"), cex=1)
# dev.off()

glm2F$seas = NA
glm2F$seas[glm2F$seas.p>Q01] = "Non-seasonal"
glm2F$seas[glm2F$seas.p<=Q01] = "Seasonal (Quant=0.01)"
glm2F_df = na.omit(data.frame(AF=abs(glm2F$meanf_s), seas=as.factor(glm2F$seas) ))
glm2F_df$seas = factor(glm2F_df$seas, levels=rev(levels(glm2F_df$seas)))


pdf("glm_seasonalQ01_meanspAF_hist_polymorphic_ggplot.pdf", width=4, height=3)
#glm_AF_hist = 
ggplot(glm2F_df, aes(x=AF, fill=seas)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.001,
                 position="identity", alpha=0.6) + 
  scale_fill_manual(name='', values=c("Non-seasonal"="grey","Seasonal (Quant=0.01)"=rgb(1,0,0,0.2)) ) +
  xlab(expression(paste("Mean spring/fall ",Delta,"AF"))) + 
  ylab("Density of SNPs") + 
  xlim(c(0,0.15)) +
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.5), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.title.align=.5,
        axis.text=element_text(size=8),
        axis.title=element_text(size=10))     
dev.off()

save(glm_hist, perm_df2, glm_df, glm_AF_hist, glm2F, file="figure_objects_sig_snps_Jan2018.Rdata")

save(glm_df_orig,glm2FM_samp_orig, file="/Users/hm8/nescent_melCA/significant_snps/ggplot_glmsig_deltaAF.Rdata")






#################### mean spring allele frequency for all data compared with Q<0.01
#freqdiff = read.table("../data/mel_meanfreqdiff_paired20_042016.txt", header=TRUE, stringsAsFactors=FALSE)
load("../data/means_dpfreq.paired20.allmean.Rdata")
# "filter"       "glm"          "glm2"         "glm2A"        "means_dpfreq" "pvalues" 
meansall = means_dpfreq
glm2FM = merge(glm2F, meansall[,c(1:5)], by=c(1,2))
Q01 = quantile(glm2FM$seas.p, probs=0.01, na.rm=TRUE)
top1500 = glm2FM[order(glm2FM$seas.p)[1:1500],]
summary(abs(top1500$meanf_s))
glm2FM$top1500 = "no"
glm2FM$top1500[order(glm2FM$seas.p)[1:1500]] = "top1500"
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03099 0.06646 0.07609 0.07635 0.08580 0.15384 

summary(abs(glm2FM$meanf_s[glm2FM$seas== "(Quant=0.01)"]))


h1 = hist(abs(glm2FM$Smean[glm2FM$seas.p>Q01]), plot=FALSE, breaks=50)
h1$density = h1$counts/sum(h1$counts)
h2 = hist(abs(glm2FM$Smean[glm2FM$seas.p<Q01]), plot=FALSE, breaks=h1$breaks)
h2$density = h2$counts/sum(h2$counts)

pdf("glm_seasonalQ01_Smean_hist_polymorphic.pdf", width=4, height=4)
plot(h1, freq=FALSE, col="grey", xlim=c(0,1), main="", xlab=expression(paste(" Spring AF (mean)")), ylab="Proportion of SNPs")
plot(h2, freq=FALSE, add=TRUE, col=rgb(1,0,0,0.2))
legend("topright",bty="n",legend=c("Seasonal (Q<0.01)","Non-seasonal"), fill=c(rgb(1,0,0,0.2), "grey"), cex=1)
dev.off()

h1 = hist(abs(glm2$Fmean[glm2$seas.p>Q01]), plot=FALSE, breaks=50)
h1$density = h1$counts/sum(h1$counts)
h2 = hist(abs(glm2$Fmean[glm2$seas.p<Q01]), plot=FALSE, breaks=h1$breaks)
h2$density = h2$counts/sum(h2$counts)

pdf("glm_seasonalQ01_Fmean_hist_polymorphic.pdf", width=4, height=4)
plot(h1, freq=FALSE, col="grey", xlim=c(0,1), main="", xlab=expression(paste(" Fall AF (mean)")), ylab="Proportion of SNPs")
plot(h2, freq=FALSE, add=TRUE, col=rgb(1,0,0,0.2))
legend("topright",bty="n",legend=c("Seasonal (Q<0.01)","Non-seasonal"), fill=c(rgb(1,0,0,0.2), "grey"), cex=1)
dev.off()


# S/F change as a function of spring allele frequency 
fold = function(x){
  xnew = x
  if(x > 0.5){
    xnew = 1-x
  }
  xnew
}

glm2FM$Smean_fold = unlist(lapply(glm2FM$Smean, FUN=fold))
glm2FM$Smean_bin = round(glm2FM$Smean_fold, digits=2)
glm2FM$seas = factor(glm2FM$seas, levels=c("Seasonal (Quant=0.01)","Non-seasonal") )

meanfs_by_Safbin = tapply(abs(glm2FM$meanf_s), INDEX=glm2FM$Smean_bin, FUN=mean)
plot(names(meanfs_by_Safbin), meanfs_by_Safbin)
#meanfs_by_Safbin_df = data.frame(names(meanfs_by_Safbin), meanfs_by_Safbin)
#colnames(meanfs_by_Safbin_df) = c()
glm2FMsig = glm2FM[glm2FM$seas=="Seasonal (Quant=0.01)", ]
meanfs_by_Safbin_sig = tapply(abs(glm2FMsig$meanf_s), INDEX=glm2FMsig$Smean_bin, FUN=mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.03733 0.05485 0.06506 0.06223 0.07203 0.07576       1 

glm2FM1500 = glm2FM[glm2FM$top1500=="top1500", ]
meanfs_by_Safbin_1500 = tapply(abs(glm2FM1500$meanf_s), INDEX=glm2FM1500$Smean_bin, FUN=mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03882 0.06294 0.07593 0.07197 0.08346 0.08970 
seas_inds = which(glm2FM$seas=="Seasonal (Quant=0.01)")
glm2FM_samp = glm2FM[c(sample(which(glm2FM$seas=="Non-seasonal"), size=length(seas_inds)),seas_inds), ]


## Somehow lost the code for Figure 2B, so need to recreate here
#glm_AF_scatter = ### this is what is was called- this currently isn't a perfect replication
ggplot(glm2FM_samp) + 
  geom_point(aes(Smean_fold, abs(meanf_s), color=seas), pch=20, alpha=0.1) +
  geom_smooth(aes(Smean_fold, abs(meanf_s), color=seas)) + 
  #geom_smooth(mapping=aes(Smean_fold, abs(meanf_s), color="seas") ) + 
  scale_color_manual("", values=c("Non-seasonal"="black","Seasonal (Quant=0.01)"="red") ) +
  ylim(c(0,0.15)) + 
  ylab(expression(paste("Spring/fall ",Delta,"AF (abs)"))) +
  xlab("Spring AF") +
  theme(legend.direction="vertical", 
        legend.justification=c("left"), 
        legend.position=c(0.15,0.9), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.title.align=0,
        legend.margin=margin(0,0,0,0),
        #      legend.box.just = "left",
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title=element_text(size=10))  

glm2FM_samp_orig = glm2FM_samp

# spring_fall mean allele frequency change for all data compared with Q<0.01. FISHERS
fisher = read.table("/Users/hm8/nescent_melCA/fishers_method_Oct2016/jamie_fisher_method/fisher_exactJ_pop20_polymorphic/results/L_rank_fisher_exactJ.merged.20pop.polymorphic_chrom_pos.txt", header=TRUE, stringsAsFactors=FALSE)
Lhigh = apply(fisher[,3:4], MARGIN=1, FUN=max, na.rm=TRUE)
fisher$Lhigh = Lhigh
coef = apply(fisher[,3:4], MARGIN=1, FUN=function(X) X[1]>X[2] )
coef[coef==TRUE] = 1
coef[coef==FALSE] = -1
fisher$coef = coef
rankF = rank(1/fisher$Lhigh)
fisher$Q = rankF / nrow(fisher)
glm2 = merge(fisher[,c(1,2,6,7)], freqdiff, by=c(1,2))
colnames(glm2) = c("chrom","pos","seas.coef","seas.p","meanf_s")
Q01 = quantile(glm2$seas.p, probs=0.01, na.rm=TRUE)

h1 = hist(abs(glm2$meanf_s[glm2$seas.p>Q01]), plot=FALSE, breaks=10)
h1$density = h1$counts/sum(h1$counts)
h2 = hist(abs(glm2$meanf_s[glm2$seas.p<Q01]), plot=FALSE, breaks=seq(from=0, to = 0.130, by=0.01) )
h2$density = h2$counts/sum(h2$counts)

# h1 = hist(abs(glm2$meanf_s[glm2$seas.p>Q01]), breaks=20)
# h1$density = h1$counts/sum(h1$counts)
# h2 = hist(abs(glm2$meanf_s[glm2$seas.p<Q01]), breaks=18 )
# h2$density = h2$counts/sum(h2$counts)


pdf("fisher_seasonalQ01_meanspAF_hist_polymorphic.pdf", width=4, height=4)
plot(h1, freq=FALSE, col="grey", xlim=c(0,0.15), main="", xlab=expression(paste("Mean spring/fall ",Delta,"AF")), ylab="Proportion of SNPs", ylim=c(0,0.4))
plot(h2, freq=FALSE, add=TRUE, col=rgb(1,0,0,0.2))
legend("topright",bty="n",legend=c("Seasonal (Q<0.01)","Non-seasonal"), fill=c(rgb(1,0,0,0.2), "grey"), cex=1)
dev.off()

core20 = c(3,4,6:10,12:15,17:19,21:25) # not including WI_13, PA_9, PA_14, PA_15




################# By population
# S/F change as a function of spring allele frequency 
fold = function(x){
  xnew = x
  if(x > 0.5){
    xnew = 1-x
  }
  xnew
}
freqdiffpop = read.table("../data/mel_freqdiff_042016_paired20.txt", header=T, stringsAsFactors = F)
meanf_s = apply(freqdiffpop[, 3:ncol(freqdiffpop)], MARGIN=1, FUN=mean, na.rm=TRUE)
freqdiffpop2 = data.frame(freqdiffpop[,1:2], meanf_s, freqdiffpop[, 3:ncol(freqdiffpop)])
load("../data/means_dpfreq.paired20.allmean.Rdata")
# "filter"       "glm"          "glm2"         "glm2A"        "means_dpfreq" "pvalues" 
meansall = means_dpfreq
glm2M = merge(glm2, meansall[,c(1:5)], by=c(1,2))
Q01 = quantile(glm2M$seas.p, probs=0.01, na.rm=TRUE)
glm2M$Smean_fold = unlist(lapply(glm2M$Smean, FUN=fold))
glm2M$Smean_bin = round(glm2M$Smean_fold, digits=2)
glm2M$seas = "Non-seasonal"
glm2M$seas[glm2M$seas.p<Q01] = "Seasonal (Quant=0.01)"
glm2M$seas = factor(glm2M$seas, levels=c("Seasonal (Quant=0.01)","Non-seasonal") )
glm2FM = merge(glm2M, freqdiffpop2, by=c(1,2))
# top1500 = glm2FM[order(glm2FM$seas.p)[1:1500],]
# summary(abs(top1500$freqmean))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03092 0.06441 0.07357 0.07412 0.08366 0.15363
# glm2FM$top1500 = "no"
# glm2FM$top1500[order(glm2FM$seas.p)[1:1500]] = "top1500"
# meanfs_by_Safbin = tapply(abs(glm2FM$meanf_s), INDEX=glm2FM$Smean_bin, FUN=mean)
# plot(names(meanfs_by_Safbin), meanfs_by_Safbin)
glm2FMsig = glm2FM[glm2FM$seas=="Seasonal (Quant=0.01)", ]
meanfs_by_Safbin_sig = tapply(abs(glm2FMsig$meanf_s), INDEX=glm2FMsig$Smean_bin, FUN=mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.03733 0.05485 0.06506 0.06223 0.07203 0.07576       1 
#glm2FM1500 = glm2FM[glm2FM$top1500=="top1500", ]
#meanfs_by_Safbin_1500 = tapply(abs(glm2FM1500$meanf_s), INDEX=glm2FM1500$Smean_bin, FUN=mean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03882 0.06294 0.07593 0.07197 0.08346 0.08970 
seas_inds = which(glm2FM$seas=="Seasonal (Quant=0.01)")
glm2FM_samp = glm2FM[c(sample(which(glm2FM$seas=="Non-seasonal"), size=length(seas_inds)),seas_inds), ]

meanfs_by_Safbin_sigPmat = matrix(ncol=20,nrow=47)
meanfs_by_Safbin_sampPmat = matrix(ncol=20,nrow=47)
for (i in 1:20){
  focal = glm2FMsig[,c(1:13,i+13)]
  focal2 = glm2FM_samp[,c(1:13,i+13)]
  colnames(focal)[14] = "popf_s"
  colnames(focal2)[14] = "popf_s"
  meanfs_by_Safbin_sigPmat[,i] = unlist(tapply(abs(focal$popf_s), INDEX=focal$Smean_bin, FUN=mean, na.rm=TRUE))
  meanfs_by_Safbin_sampPmat[,i] = unlist(tapply(abs(focal2$popf_s), INDEX=focal2$Smean_bin, FUN=mean, na.rm=TRUE))
}
Smean_bin = names(unlist(tapply(abs(focal$popf_s), INDEX=focal$Smean_bin, FUN=mean, na.rm=TRUE)))
colnames(meanfs_by_Safbin_sigPmat) = colnames(glm2FMsig[,14:33])
meanfs_by_Safbin_sigPmat2 = data.frame(Smean_bin, meanfs_by_Safbin_sigPmat)
meanfs_by_Safbin_sigPmat3 = melt(meanfs_by_Safbin_sigPmat2)
colnames(meanfs_by_Safbin_sampPmat) = colnames(glm2FMsig[,14:33])
meanfs_by_Safbin_sampPmat2 = data.frame(Smean_bin, meanfs_by_Safbin_sampPmat)
meanfs_by_Safbin_sampPmat3 = melt(meanfs_by_Safbin_sampPmat2)
meanfs_by_Safbin_sigPmat3$sig = "sig"
meanfs_by_Safbin_sampPmat3$sig = "nonsig"

meanfs_by_Safbin_Pmat3 = rbind(meanfs_by_Safbin_sigPmat3, meanfs_by_Safbin_sampPmat3)
meanfs_by_Safbin_Pmat3$group2 = paste(meanfs_by_Safbin_Pmat3$variable, meanfs_by_Safbin_Pmat3$sig)

ggplot(meanfs_by_Safbin_Pmat3, aes(Smean_bin, value, group=group2, col=sig))+
  geom_line() + 
  scale_color_manual(values=c("sig"=rgb(1,0,0,0.3), "nonsig"=rgb(0,0,0,0.3)) )
  
  geom_line(meanfs_by_Safbin_sampPmat3, aes(Smean_bin, value, group=variable), col="black")

#ggplot(meanfs_by_Safbin_sampPmat3, aes(Smean_bin, value, group=variable))+
  geom_line(col=rgb(0,0,0,alpha=0.5))

ggplot(meanfs_by_Safbin_sigPmat3, aes(Smean_bin, value, group=variable))+
    geom_line(col=rgb(1,0,0,alpha=0.5)) 
    

glm2FM_sampP = merge(glm2FM_samp, core20, by=c(1,2))
glm2FM_sampPO = glm2FM_sampP[order(Smean_fold$Smean), ]

library(caTools)
k1 = runmean(glm2FM_sampPO[,14+i], k=1000)

plot(glm2FM_sampP$Smean_fold, abs(glm2FM_sampP[,14+1]), col=rgb(1,0,0,alpha=0.1), pch=20)
plot(glm2FM_sampP$Smean_fold, abs(glm2FM_sampP[,14+1]), col=rgb(1,0,0,alpha=0.1), pch=20)
  
#glm_AF_scatter = ### this is what is was called- this currently isn't a perfect replication
ggplot(glm2FM_samp) + 
  geom_point(aes(Smean_fold, abs(meanf_s), color=seas), pch=20, alpha=0.1) +
  geom_smooth(aes(Smean_fold, abs(meanf_s), color=seas)) + 
  #geom_smooth(mapping=aes(Smean_fold, abs(meanf_s), color="seas") ) + 
  scale_color_manual("", values=c("Non-seasonal"="black","Seasonal (Quant=0.01)"="red") ) +
  ylim(c(0,0.15)) + 
  ylab(expression(paste("Spring/fall ",Delta,"AF (abs)"))) +
  xlab("Spring AF") +
  theme(legend.direction="vertical", 
        legend.justification=c("left"), 
        legend.position=c(0.15,0.9), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.title.align=0,
        legend.margin=margin(0,0,0,0),
        #      legend.box.just = "left",
        axis.text=element_text(size=8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10)),
        axis.title=element_text(size=10))  
