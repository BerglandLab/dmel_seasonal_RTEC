## Bayenv analysis of clinal variation
# Feb 2019

setwd("/Users/hm8/stanford/nescent_melCA/clinal/bayenv")
library(ggplot2)
library(reshape2)

# Read in files
snps = read.table("clinal_allsnps_non0_unformatted_Feb2019.txt", stringsAsFactors = F)
clinal = read.table("environ_corr.clinal_environfile_Feb2019.txt", stringsAsFactors = F)
ID = unlist(lapply(clinal[,1], FUN=function(X) unlist(strsplit(unlist(strsplit(X, split="clinal_snpsfile_Feb2019_snp"))[2], split=".txt", fixed=T))[1]   ))
clinal$ID = ID
snps$ID = 1:nrow(snps)
clinal_snps = merge(clinal, snps, by="ID")
clinal_snps$FL = clinal_snps$V3/clinal_snps$V7
clinal_snps$GA = clinal_snps$V4/clinal_snps$V8
clinal_snps$PA = clinal_snps$V5/clinal_snps$V9
clinal_snps$SC = clinal_snps$V6/clinal_snps$V10
clinal_snps$chr_pos = paste(clinal_snps$V1.y, clinal_snps$V2.y, sep="_")
colnames(clinal_snps)[4] = "chrom"
colnames(clinal_snps)[5] = "pos"

### Polymorphic sites only
my.filter = read.table("../../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
colnames(my.filter) = c("chrom","pos")
filter = my.filter[my.filter[,1]!="X",]
clinal_snps2 = merge(clinal_snps, filter, by=c("chrom","pos") )  ## 781110

# Latitude
# lat = c(25.5, 30.99, 39.88, 33.39)
#"melFL_072010_SPT": 25.5
#"melGA_072008_SPT": 30.99
#"melPA_72011_SPT": 39.88
#"melSC_072010_SPT": 33.39

# Comparison with old analysis
oldclinal = read.table("../../glm/mel_clinal_uniquepops_springPA_noMA.glm.noheader",stringsAsFactors = F)
#oldclinal$chr_pos = paste(oldclinal$V1, oldclinal$V2, sep="_")
colnames(oldclinal) = c("chrom","pos","lm_beta","lm_pvalue","N")
clinal_snps3 = merge(oldclinal, clinal_snps2, by=c("chrom","pos") )
colnames(clinal_snps3)[colnames(clinal_snps3)=="V2.x"] = "Z"
summary(lm(clinal_snps3$lm_pvalue~clinal_snps3$Z))
# (Intercept)     0.4757446  0.0005918   803.9   <2e-16 ***
#   clinal_snps3$Z -0.9252675  0.0041199  -224.6   <2e-16 ***
# Residual standard error: 0.2993 on 781108 degrees of freedom
# Multiple R-squared:  0.06066,	Adjusted R-squared:  0.06065 
# F-statistic: 5.044e+04 on 1 and 781108 DF,  p-value: < 2.2e-16
## -0.6832792 strong negative correlation (the lm is pvalue, the bayenv is 0-0.5)

summary(lm(clinal_snps3$lm_beta~clinal_snps3$Z))
# (Intercept)     5.134e-02  9.537e-05   538.3   <2e-16 ***
#   clinal_snps3$Z -4.715e-01  6.639e-04  -710.1   <2e-16 ***
# Residual standard error: 0.04823 on 781108 degrees of freedom
# Multiple R-squared:  0.3923,	Adjusted R-squared:  0.3923 
# F-statistic: 5.043e+05 on 1 and 781108 DF,  p-value: < 2.2e-16
## Better correlation with the lm beta



# Read in seasonal data
all = read.table("../../glm/mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm", stringsAsFactors = F, header=T)
ca = read.table("../../glm/mel_ca_popyear_paired.f_s.glm", stringsAsFactors = F, header=T)
eur = read.table("../../glm/mel_eur.seas_pop_year.f_s.glm.noheader", stringsAsFactors = F)
colnames(eur) = c("chrom","pos","seas.coef",  "seas.p", "seas1.N", "seas2.N")

clinal_all = na.omit(merge(all[,c("seas.coef","seas.p","chrom","pos")], clinal_snps3[,c("lm_beta","lm_pvalue","Z","chrom","pos")], by= c("chrom","pos") ) )
myN = nrow(clinal_all)
clinal_all$seas.q = rank(clinal_all$seas.p)/myN ## 
clinal_all$bayclinal.q = (myN-rank(abs(clinal_all$Z))) /myN ## 0.5 is the greatest correlation
clinal_all$lmclinal.q = rank(clinal_all$lm_pvalue)/myN ## 
# summary(lm(seas.q~lmclinal.q))
# summary(lm(seas.q~bayclinal.q))

clinal_ca = na.omit(merge(ca[,c("seas.coef","seas.p","chrom","pos")], clinal_snps3[,c("lm_beta","lm_pvalue","Z","chrom","pos")], by=c("chrom","pos")))
myN = nrow(clinal_ca)
clinal_ca$seas.q = rank(clinal_ca$seas.p)/myN ## 
clinal_ca$bayclinal.q = (myN-rank(abs(clinal_ca$Z)) )/myN ## 0.5 is the greatest correlation
clinal_ca$lmclinal.q = rank(clinal_ca$lm_pvalue)/myN ## 

clinal_eur = na.omit(merge(eur[,c("seas.coef","seas.p","chrom","pos")], clinal_snps3[,c("lm_beta","lm_pvalue","Z","chrom","pos")], by=c("chrom","pos")))
myN = nrow(clinal_eur)
clinal_eur$seas.q = rank(clinal_eur$seas.p)/myN ## 
clinal_eur$bayclinal.q = (myN-rank(abs(clinal_eur$Z)) )/myN ## 0.5 is the greatest correlation
clinal_eur$lmclinal.q = rank(clinal_eur$lm_pvalue)/myN ## 


## For each quantile, measure concordance (using glm beta)
quant = 10^(-seq(from=0, to=3, by=0.1))
inds = 1:28
concord_lm = vector()
concord_bay = vector()
concord_lm_ca = vector()
concord_bay_ca = vector()
concord_lm_eur = vector()
concord_bay_eur = vector()

sd_lm = vector()
sd_bay = vector()
sd_lm_ca = vector()
sd_bay_ca = vector()
sd_lm_eur = vector()
sd_bay_eur = vector()

# For cumulative
for (i in 1:length(inds)){
  focalquant = quant[i]
  focallm = clinal_all[clinal_all$seas.q<focalquant & clinal_all$lmclinal.q<focalquant, ]
  focalbay = clinal_all[clinal_all$seas.q<focalquant & clinal_all$bayclinal.q<focalquant, ]
  concord_lm[i] = mean( (focallm$seas.coef<=0 & focallm$lm_beta<=0) | (focallm$seas.coef>0 & focallm$lm_beta>0) )
  concord_bay[i] = mean( (focalbay$seas.coef<=0 & focalbay$lm_beta<=0) | (focalbay$seas.coef>0 & focalbay$lm_beta>0) )
  sd_lm[i] = sqrt( (1-concord_lm[i])*concord_lm[i]/nrow(focallm) )
  sd_bay[i] = sqrt( (1-concord_bay[i])*concord_bay[i]/nrow(focalbay) )
  
  
  focallm = clinal_ca[clinal_ca$seas.q<focalquant & clinal_ca$lmclinal.q<focalquant, ]
  focalbay = clinal_ca[clinal_ca$seas.q<focalquant & clinal_ca$bayclinal.q<focalquant, ]
  concord_lm_ca[i] = mean( (focallm$seas.coef<=0 & focallm$lm_beta<=0) | (focallm$seas.coef>0 & focallm$lm_beta>0) )
  concord_bay_ca[i] = mean( (focalbay$seas.coef<=0 & focalbay$lm_beta<=0) | (focalbay$seas.coef>0 & focalbay$lm_beta>0) )
  sd_lm_ca[i] = sqrt( (1-concord_lm_ca[i])*concord_lm_ca[i]/nrow(focallm) )
  sd_bay_ca[i] = sqrt( (1-concord_bay_ca[i])*concord_bay_ca[i]/nrow(focalbay) )
  
  
  focallm = clinal_eur[clinal_eur$seas.q<focalquant & clinal_eur$lmclinal.q<focalquant, ]
  focalbay = clinal_eur[clinal_eur$seas.q<focalquant & clinal_eur$bayclinal.q<focalquant, ]
  concord_lm_eur[i] = mean( (focallm$seas.coef<=0 & focallm$lm_beta<=0) | (focallm$seas.coef>0 & focallm$lm_beta>0) )
  concord_bay_eur[i] = mean( (focalbay$seas.coef<=0 & focalbay$lm_beta<=0) | (focalbay$seas.coef>0 & focalbay$lm_beta>0) )
  sd_lm_eur[i] =sqrt( (1-concord_lm_eur[i])*concord_lm_eur[i]/nrow(focallm) )
  sd_bay_eur[i] = sqrt( (1-concord_bay_eur[i])*concord_bay_eur[i]/nrow(focalbay) )
}

concordDF =  data.frame(quantile = quant[inds], bayenv=concord_bay, lm=concord_lm, bayenv_ca=concord_bay_ca, lm_ca=concord_lm_ca, bayenv_eur=concord_bay_eur, lm_eur=concord_lm_eur)
concordDF_sd =  data.frame(quantile = quant[inds], bayenv=sd_bay, lm=sd_lm, bayenv_ca=sd_bay_ca, lm_ca=sd_lm_ca, bayenv_eur=sd_bay_eur, lm_eur=sd_lm_eur)

melt_concordDF = melt(concordDF, id="quantile")
melt_concordDF$model = "bayenv"
melt_concordDF$model[melt_concordDF$variable=="lm" | melt_concordDF$variable=="lm_ca" | melt_concordDF$variable=="lm_eur"] = "lm"
melt_concordDF$dataset = "all"
melt_concordDF$dataset[melt_concordDF$variable=="bayenv_ca" | melt_concordDF$variable=="lm_ca"] = "ca"
melt_concordDF$dataset[melt_concordDF$variable=="bayenv_eur" | melt_concordDF$variable=="lm_eur"] = "eur"
melt_concordDF_sd = melt(concordDF_sd, id="quantile")
melt_concordDF$sd = melt_concordDF_sd$value

## preparing for plotting
melt_concordDF$model = revalue(melt_concordDF$model, c("lm"="GLM"))
melt_concordDF$model = revalue(melt_concordDF$model, c("bayenv"="Bayenv"))

# polygon dataframes
mysubset_bayenv = subset(melt_concordDF, quantile>=0.01 & dataset=="all" & model=="Bayenv")
mysubset_lm = subset(melt_concordDF, quantile>=0.01 & dataset=="all" & model=="GLM")
sd_mysubset_bayenv = data.frame(x=c(rev(mysubset_bayenv$quantile),mysubset_bayenv$quantile), y=c(rev(mysubset_bayenv$value + 2*mysubset_bayenv$sd),mysubset_bayenv$value - 2*mysubset_bayenv$sd) )
sd_mysubset_lm= data.frame(x=c(rev(mysubset_lm$quantile),mysubset_lm$quantile), y=c(rev(mysubset_lm$value + 2*mysubset_lm$sd),mysubset_lm$value - 2*mysubset_lm$sd) )


clinal_bayenv_thresh = ggplot(subset(melt_concordDF, quantile>=0.01 & dataset=="all"), aes(x=quantile, y=value, group=model))+
  geom_line(aes(color=model))+
  #geom_line(col="black")+
  #geom_point(cex=0.8)+
  #facet_grid(.~dataset)+
  scale_x_log10()+
  coord_cartesian(ylim=c(0.25, 1))+
  geom_polygon(data=sd_mysubset_bayenv, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="#d95f02", alpha=0.2) +
  geom_polygon(data=sd_mysubset_lm, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="#1b9e77", alpha=0.2) +
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2) +
  scale_color_manual("", values=c("Bayenv"="#d95f02","GLM"="#1b9e77","eur"="#7570b3"))+
  ylab("Concordance")+
  xlab("Joint quantile threshold")+
  geom_hline(yintercept=0.5, linetype=2) + #, color, size)
  theme_light()
ggsave("bayenv_vs_lm_clinalconcordance_poly_original_errorbars_Nov2020.pdf", width=6,height=4)
ggsave("bayenv_vs_lm_clinalconcordance_poly_original_errorbars_Nov2020.png", width=6,height=4)
# ggplot(melt_concordDF, aes(quantile, value, group=dataset, color=dataset))+
#   geom_line()+
#   geom_point()+
#   facet_wrap(.~model)+
#   scale_x_log10()+
#   ylab("Concordance")+
#   xlab("Joint quantile")+
#   geom_hline(yintercept=0.5, linetype=2) + #, color, size)
#   theme_light()
# ggsave("bayenv_vs_lm_clinalconcordance_poly.pdf", width=6,height=3)




# For non-cumulative
for (i in 1:(length(inds)-1)){
  focalquant = quant[i]
  focalquantend = quant[i+1]
  focallm = clinal_all[clinal_all$seas.q<focalquant & clinal_all$lmclinal.q<focalquant & clinal_all$seas.q>focalquantend & clinal_all$lmclinal.q>focalquantend, ]
  focalbay = clinal_all[clinal_all$seas.q<focalquant & clinal_all$bayclinal.q<focalquant & clinal_all$seas.q>focalquantend & clinal_all$bayclinal.q>focalquantend, ]
  concord_lm[i] = mean( (focallm$seas.coef<=0 & focallm$lm_beta<=0) | (focallm$seas.coef>0 & focallm$lm_beta>0) )
  concord_bay[i] = mean( (focalbay$seas.coef<=0 & focalbay$lm_beta<=0) | (focalbay$seas.coef>0 & focalbay$lm_beta>0) )
  sd_lm[i] = sqrt( (1-concord_lm[i])*concord_lm[i]/nrow(focallm) )
  sd_bay[i] = sqrt( (1-concord_bay[i])*concord_bay[i]/nrow(focalbay) )
  
  # CA
  focallm = clinal_ca[clinal_ca$seas.q<focalquant & clinal_ca$lmclinal.q<focalquant & clinal_ca$seas.q>focalquantend & clinal_ca$lmclinal.q>focalquantend, ]
  focalbay = clinal_ca[clinal_ca$seas.q<focalquant & clinal_ca$bayclinal.q<focalquant & clinal_ca$seas.q>focalquantend & clinal_ca$bayclinal.q>focalquantend, ]
  concord_lm_ca[i] = mean( (focallm$seas.coef<=0 & focallm$lm_beta<=0) | (focallm$seas.coef>0 & focallm$lm_beta>0) )
  concord_bay_ca[i] = mean( (focalbay$seas.coef<=0 & focalbay$lm_beta<=0) | (focalbay$seas.coef>0 & focalbay$lm_beta>0) )
  sd_lm_ca[i] = sqrt( (1-concord_lm_ca[i])*concord_lm_ca[i]/nrow(focallm) )
  sd_bay_ca[i] = sqrt( (1-concord_bay_ca[i])*concord_bay_ca[i]/nrow(focalbay) )
  
  # EUR
  focallm = clinal_eur[clinal_eur$seas.q<focalquant & clinal_eur$lmclinal.q<focalquant & clinal_eur$seas.q>focalquantend & clinal_eur$lmclinal.q>focalquantend, ]
  focalbay = clinal_eur[clinal_eur$seas.q<focalquant & clinal_eur$bayclinal.q<focalquant & clinal_eur$seas.q>focalquantend & clinal_eur$bayclinal.q>focalquantend, ]
  
  concord_lm_eur[i] = mean( (focallm$seas.coef<=0 & focallm$lm_beta<=0) | (focallm$seas.coef>0 & focallm$lm_beta>0) )
  concord_bay_eur[i] = mean( (focalbay$seas.coef<=0 & focalbay$lm_beta<=0) | (focalbay$seas.coef>0 & focalbay$lm_beta>0) )
  sd_lm_eur[i] =sqrt( (1-concord_lm_eur[i])*concord_lm_eur[i]/nrow(focallm) )
  sd_bay_eur[i] = sqrt( (1-concord_bay_eur[i])*concord_bay_eur[i]/nrow(focalbay) )
  
}

concordDFnoncum =  data.frame(quantile = quant[inds], bayenv=concord_bay, lm=concord_lm, bayenv_ca=concord_bay_ca, lm_ca=concord_lm_ca, bayenv_eur=concord_bay_eur, lm_eur=concord_lm_eur)
concordDFnoncum_sd =  data.frame(quantile = quant[inds], bayenv=sd_bay, lm=sd_lm, bayenv_ca=sd_bay_ca, lm_ca=sd_lm_ca, bayenv_eur=sd_bay_eur, lm_eur=sd_lm_eur)

melt_concordDFnoncum = melt(concordDFnoncum, id="quantile")
melt_concordDFnoncum$model = "bayenv"
melt_concordDFnoncum$model[melt_concordDFnoncum$variable=="lm" | melt_concordDFnoncum$variable=="lm_ca" | melt_concordDFnoncum$variable=="lm_eur"] = "lm"
melt_concordDFnoncum$dataset = "all"
melt_concordDFnoncum$dataset[melt_concordDFnoncum$variable=="bayenv_ca" | melt_concordDFnoncum$variable=="lm_ca"] = "ca"
melt_concordDFnoncum$dataset[melt_concordDFnoncum$variable=="bayenv_eur" | melt_concordDFnoncum$variable=="lm_eur"] = "eur"
melt_concordDFnoncum_sd = melt(concordDFnoncum_sd, id="quantile")
melt_concordDFnoncum$sd = melt_concordDFnoncum_sd$value

######### Nov 2020- supplementary figure with error bars
library(plyr)
melt_concordDFnoncum$model = revalue(melt_concordDFnoncum$model, c("lm"="GLM"))
melt_concordDFnoncum$model = revalue(melt_concordDFnoncum$model, c("bayenv"="Bayenv"))

nCA = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )
n2 = data.frame(quant, nescentClinal_Cquant_controlrEU)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nEU = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )

  geom_polygon(data=nCA, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="purple", alpha=0.8) +
  geom_polygon(data=nEU, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="light green", alpha=0.5)

mysubset_bayenv = subset(melt_concordDFnoncum, quantile>=0.01 & dataset=="all" & model=="Bayenv")
mysubset_lm = subset(melt_concordDFnoncum, quantile>=0.01 & dataset=="all" & model=="GLM")
sd_mysubset_bayenv = data.frame(x=c(rev(mysubset_bayenv$quantile),mysubset_bayenv$quantile), y=c(rev(mysubset_bayenv$value + 2*mysubset_bayenv$sd),mysubset_bayenv$value - 2*mysubset_bayenv$sd) )
sd_mysubset_lm= data.frame(x=c(rev(mysubset_lm$quantile),mysubset_lm$quantile), y=c(rev(mysubset_lm$value + 2*mysubset_lm$sd),mysubset_lm$value - 2*mysubset_lm$sd) )


clinal_bayenv_bin = ggplot(subset(melt_concordDFnoncum, quantile>=0.01 & dataset=="all"), aes(x=quantile, y=value, group=model))+
  geom_line(aes(color=model))+
  #geom_line(col="black")+
  #geom_point(cex=0.8)+
  #facet_grid(.~dataset)+
  scale_x_log10()+
  coord_cartesian(ylim=c(0.25, 1))+
  geom_polygon(data=sd_mysubset_bayenv, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="#d95f02", alpha=0.2) +
  geom_polygon(data=sd_mysubset_lm, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="#1b9e77", alpha=0.2) +
  #geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2) +
  scale_color_manual("", values=c("Bayenv"="#d95f02","GLM"="#1b9e77","eur"="#7570b3"))+
  ylab("Concordance")+
  xlab("Joint quantile bin")+
  geom_hline(yintercept=0.5, linetype=2) + #, color, size)
  theme_light()
ggsave("bayenv_vs_lm_clinalconcordance_poly_original_noncum_errorbars_Nov2020.pdf", width=6,height=4)
ggsave("bayenv_vs_lm_clinalconcordance_poly_original_noncum_errorbars_Nov2020.png", width=6,height=4)

library(cowplot)
save(clinal_bayenv_bin, clinal_bayenv_thresh, file="clinal_bayenv_plot_objects.Rdata")
clinal_bayenv_thresh2 = clinal_bayenv_thresh + theme(legend.position = "none")
clinal_bayenv_bin2 = clinal_bayenv_bin + ylab("")
p1 = plot_grid(clinal_bayenv_thresh2, clinal_bayenv_bin2, rel_widths = c(1,1.3), labels=c("A","B"))
save_plot(p1, file="bayenv_vs_lm_clinalconcordance_poly_original_cum_noncum_errorbars_Nov2020.pdf", base_width=8,base_height=3.2)
save_plot(p1, file="bayenv_vs_lm_clinalconcordance_poly_original_cum_noncum_errorbars_Nov2020.png", base_width=8,base_height=3.2)
