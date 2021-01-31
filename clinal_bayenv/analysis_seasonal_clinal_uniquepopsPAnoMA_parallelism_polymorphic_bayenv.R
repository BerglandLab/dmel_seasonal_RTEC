## April 2016
## testing parallelism in seasonal sites
 
setwd("/Users/hm8/stanford/nescent_melCA/clinal/bayenv")

############# Small quantile range
# load("results/clinal_uniquepopsPAnoMA_parallelism_vs_controlsDec2017_quant_polymorphic.Rdata")

############ Extended quantile range
#############
load("clinal_uniquepopsPAnoMA_parallelism_vs_controlsMar2019_bayenv_polymorphic.Rdata")
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(scales)
library(RColorBrewer)

quant = 10^(-seq(from=0, to=3, by=0.1))
inds = 1:28
n2 = data.frame(quant, nescentClinal_Cquant_controlr)
colnames(n2) = c("quant","obs","med","lower96","upper96")
n3 <- melt(n2[inds,], id=c("quant"))
#n4 = data.frame(x=c(n2$quant[inds],rev(n2$quant[inds])), y=c(n2$X3[inds], rev(n2$X4[inds])) )
n4 = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )

pdf("figures/clinal_uniquepopsPAnoMA_seasonalparallel_polymorphic_Dec2017_96CI_ext_ggplot.pdf", width=3.65,height=4)
ggplot(data=n2[inds,], aes(quant, obs) ) +
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_polygon(n4, mapping=aes(x=x, y=y), fill="grey") +
  #geom_line(aes(quant, value, group=variable)) + 
  geom_line() + 
  geom_point(color="black", fill="red", shape=21, size=1) + 
  scale_colour_manual(name='', values=c("Observed"="red") ) +
  #guides(color = guide_legend(override.aes = list(linetype=c(1,0), shape=c(NA, 16)))) +
  ylab("Concordance") + 
  xlab("Joint quantile") +
  geom_hline(yintercept = 0.5, lty=2) +
  theme(legend.direction="vertical", 
        legend.justification=c(1,0), 
        legend.position=c(1,0), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.title.align=.5,
        axis.text=element_text(size=8),
        axis.title=element_text(size=10))
dev.off()
#scale_colour_manual(name="Error Bars",values=cols, guide = guide_legend(fill = NULL,colour = NULL)) + 

######### Heatmap

# Each list item is a different latitudinal quantile. 
# Each column is a seasonal quantile
# Creating a matrix for which the columns are latitudinal and rows are seasonal.
enrichmatP = matrix(nrow=length(quant), ncol=length(quant))
colnames(enrichmatP) = c(quant)
rownames(enrichmatP) = c(quant)
for (i in 1:length(quant)){
  enrichmatP[,i] = ((nescentClinal_Npvalue_controlr[[i]][,1] - nescentClinal_Npvalue_controlr[[i]][,2])/nescentClinal_Npvalue_controlr[[i]][,2])
}

# 21  = 1%
#enrichmatP[21,21]

enrichmatP96 = matrix(nrow=length(quant), ncol=length(quant))
colnames(enrichmatP96) = c(quant)
rownames(enrichmatP96) = c(quant)
for (i in 1:length(quant)){
  focal = nescentClinal_Npvalue_controlr[[i]][,1]
  naind = (nescentClinal_Npvalue_controlr[[i]][,1] > nescentClinal_Npvalue_controlr[[i]][,3] & 
             nescentClinal_Npvalue_controlr[[i]][,1] < nescentClinal_Npvalue_controlr[[i]][,4]) 
  focal[naind] = NA
  enrichmatP96[,i] = ((focal - nescentClinal_Npvalue_controlr[[i]][,2])/nescentClinal_Npvalue_controlr[[i]][,2])
}

enrichmatC = matrix(nrow=length(quant), ncol=length(quant))
colnames(enrichmatC) = c(quant)
rownames(enrichmatC) = c(quant)
for (i in 1:length(quant)){
  enrichmatC[,i] = ((nescentClinal_Cpvalue_controlr[[i]][,1] - nescentClinal_Cpvalue_controlr[[i]][,2])/nescentClinal_Cpvalue_controlr[[i]][,2])
}

matC = matrix(nrow=length(quant), ncol=length(quant))
colnames(matC) = c(quant)
rownames(matC) = c(quant)
for (i in 1:length(quant)){
  matC[,i] = nescentClinal_Cpvalue_controlr[[i]][,1]
}

enrichmatC96 = matrix(nrow=length(quant), ncol=length(quant))
colnames(enrichmatC96) = c(quant)
rownames(enrichmatC96) = c(quant)
for (i in 1:length(quant)){
  focal = nescentClinal_Cpvalue_controlr[[i]][,1]
  naind = (nescentClinal_Cpvalue_controlr[[i]][,1] > nescentClinal_Cpvalue_controlr[[i]][,3] & 
      nescentClinal_Cpvalue_controlr[[i]][,1] < nescentClinal_Cpvalue_controlr[[i]][,4]) 
  focal[naind] = NA
  enrichmatC96[,i] = ((focal - nescentClinal_Cpvalue_controlr[[i]][,2])/nescentClinal_Cpvalue_controlr[[i]][,2])
}

NmatP = matrix(nrow=length(quant), ncol=length(quant))
colnames(NmatP) = c(quant)
rownames(NmatP) = c(quant)
for (i in 1:length(quant)){
  NmatP[,i] = nescentClinal_Npvalue_controlr[[i]][,1]
}
NmatP_Nsnps = round(NmatP*773862)  ### same as the switched?
# Columns are latitudinal and rows are seasonal.
#write.table(NmatP_Nsnps, file="NmatP_Nsnps_seasonal_clinal.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

NmatC = matrix(nrow=length(quant), ncol=length(quant))
colnames(NmatC) = c(quant)
rownames(NmatC) = c(quant)
for (i in 1:length(quant)){
  NmatC[,i] = nescentClinal_Cpvalue_controlr[[i]][,1]
}
#write.table(NmatC, file="NmatC_concord_seasonal_clinal.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)


# Var1 = latitudinal
# Var2 = seasonal
myinds = 1:24
eN = melt(enrichmatP96[myinds,myinds])
eN$Var1 = as.numeric(eN$Var1)
eN$Var2 = as.numeric(eN$Var2)
eN$value = as.numeric(eN$value)*100
colnames(eN) = c("Latitudinal", "Seasonal","Enrichment")

pdf("figures/clinal_uniquepopsPAnoMA_seas_signficant_heatmap_all_caF_nonclinal_polymorphic_Dec2017_ggplot.pdf", width=5.5,height=4)
#heat_sig_poly = 
ggplot(eN, aes(Latitudinal, Seasonal)) + 
  geom_tile(aes(fill = Enrichment)) + 
  scale_fill_gradient2(expression(paste("Percent\nEnrichment")), low = muted("blue"), mid = "white", high = "black", midpoint = 0, space = "Lab", na.value = "white", guide = "colourbar") + 
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  scale_y_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  xlab("Latitudinal Quantile") + 
  ylab("Seasonal Quantile") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        #      legend.position=c(1.1,0.5), 
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
dev.off()

pdf("figures/clinal_uniquepopsPAnoMA_seas_signficant_heatmap_all_caF_nonclinal_polymorphic_Dec2017_ggplot_red.pdf", width=5.5,height=4)
#heat_sig_poly_red = 
ggplot(eN, aes(Latitudinal, Seasonal)) + 
  geom_tile(aes(fill = Enrichment)) + 
  scale_fill_gradientn(expression(paste("Percent\nEnrichment")), colors=rev(myPalette(13)),na.value = "white") + 
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  scale_y_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  xlab("Latitudinal Quantile") + 
  ylab("Seasonal Quantile") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        #      legend.position=c(1.1,0.5), 
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
dev.off()

myinds = 1:24
eN = melt(enrichmatP[myinds,myinds])
eN$Var1 = as.numeric(eN$Var1)
eN$Var2 = as.numeric(eN$Var2)
eN$value = as.numeric(eN$value)*100
colnames(eN) = c("Latitudinal", "Seasonal","Enrichment")
myPalette <- colorRampPalette(rev(brewer.pal(9, "Reds")))

heat_sig_poly_red_all = 
ggplot(eN, aes(Latitudinal, Seasonal)) + 
  geom_tile(aes(fill = Enrichment)) + 
  scale_fill_gradientn(expression(paste("Percent\nEnrichment")), colors=rev(myPalette(13)),na.value = "white") + 
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  scale_y_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  xlab("Latitudinal Quantile") + 
  ylab("Seasonal Quantile") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        #      legend.position=c(1.1,0.5), 
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
save(heat_sig_poly_red_all, file="/Users/hm8/nescent_melCA/clinal/heat_sig_poly_red_all.Rdata")

# Var1 = latitudinal
# Var2 = seasonal
myinds = 1:24
eC = melt(enrichmatC96[myinds,myinds])
eC$Var1 = as.numeric(eC$Var1)
eC$Var2 = as.numeric(eC$Var2)
eC$value = as.numeric(eC$value)
colnames(eC) = c("Latitudinal", "Seasonal","Enrichment")

pdf("figures/clinal_uniquepopsPAnoMA_enrich_concordance_heatmap_all_caF_nonclinal_polymorphic_Dec2017_ggplot.pdf", width=5.5,height=4)
#heat_poly = 
ggplot(eC, aes(Latitudinal, Seasonal)) + 
  geom_tile(aes(fill = Enrichment)) + 
  scale_fill_gradient2(expression(paste("Concordance\n Enrichment")), low = muted("red"), mid = "white", high = "blue", midpoint = 0, space = "Lab", na.value = "white", guide = "colourbar") + 
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  scale_y_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  xlab("Latitudinal Quantile") + 
  ylab("Seasonal Quantile") + 
  theme(legend.direction="vertical", 
      legend.justification=c("right"), 
#      legend.position=c(1.1,0.5), 
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
dev.off()

cC = melt(matC[myinds,myinds])
cC$Var1 = as.numeric(cC$Var1)
cC$Var2 = as.numeric(cC$Var2)
cC$value = as.numeric(cC$value)
colnames(cC) = c("Latitudinal", "Seasonal","Enrichment")

pdf("figures/clinal_uniquepopsPAnoMA_seas_concordance_heatmap_all_caF_nonclinal_polymorphic_Dec2017_ggplot.pdf", width=5.5,height=4)
#heat_poly_con = 
ggplot(cC, aes(Latitudinal, Seasonal)) + 
  geom_tile(aes(fill = Enrichment)) + 
  scale_fill_gradient2("Concordance", low = muted("red"), mid = "white", high = "black", midpoint = 0.52, space = "Lab", na.value = "white", guide = "colourbar") + 
  scale_x_continuous(trans="log10", breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  scale_y_continuous(trans="log10", breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) +
  xlab("Latitudinal Quantile") + 
  ylab("Seasonal Quantile") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        #      legend.position=c(1.1,0.5), 
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
dev.off()




myPalette <- colorRampPalette(rev(brewer.pal(9, "Reds")))

pdf("figures/clinal_uniquepopsPAnoMA_seas_concordance_heatmap_all_caF_nonclinal_polymorphic_Dec2017_ggplot_red.pdf", width=5.5,height=4)
heat_poly_con_red = 
ggplot(cC, aes(Latitudinal, Seasonal)) + 
  geom_tile(aes(fill = Enrichment)) + 
  scale_fill_gradientn("Concord.", colors=rev(myPalette(13)),na.value = "white") + 
  scale_x_continuous(trans="log10", breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  scale_y_continuous(trans="log10", breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) +
  xlab("Latitudinal Quantile") + 
  ylab("Seasonal Quantile") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        #      legend.position=c(1.1,0.5), 
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
dev.off()

######## CA and EU with clinal
load("results/clinal_uniquepopsPAnoMA_CAparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata")
nescentClinal_Cquant_controlrCA = nescentClinal_Cquant_controlr
nescentClinal_Cquant_controlCA = nescentClinal_Cquant_control
load("results/clinal_uniquepopsPAnoMA_EUparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata")
nescentClinal_Cquant_controlrEU = nescentClinal_Cquant_controlr
nescentClinal_Cquant_controlEU = nescentClinal_Cquant_control
load("results/clinal_uniquepopsPAnoMA_parallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata")
inds=1:24
df1 = data.frame(quant[inds], nescentClinal_Cquant_controlrCA[inds,1], nescentClinal_Cquant_controlrEU[inds,1], nescentClinal_Cquant_controlr[inds,1])
colnames(df1) = c("quant","CA","EU","ALL")
df2 = melt(df1, id=quant)
colnames(df2) = c("quant","pops","concord")


n2 = data.frame(quant, nescentClinal_Cquant_controlrCA)
colnames(n2) = c("quant","obs","med","lower96","upper96")
#n4 = data.frame(x=c(n2$quant[inds],rev(n2$quant[inds])), y=c(n2$X3[inds], rev(n2$X4[inds])) )
nCA = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )
n2 = data.frame(quant, nescentClinal_Cquant_controlrEU)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nEU = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )


pdf("figures/clinal_uniquepopsPAnoMA_CA_EUseasonalparallel_polymorphic_Dec2017_96CI_ext_ggplot.pdf", width=4,height=4)
#g_lat_line =
ggplot(data=df2, aes(quant, concord, group=pops) ) +
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept = 0.5, lty=2) +
  geom_polygon(data=nCA, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="purple", alpha=0.8) +
  geom_polygon(data=nEU, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="light green", alpha=0.5) +
  geom_line(aes(quant, concord)) + 
  #scale_color_manual(name='', values=c("ALL"="blue","CA"="black","EU"="black") ) +
  geom_point(aes(quant, concord, fill=pops, size=pops), color="black", shape=21) + 
  scale_fill_manual(name='', values=c("ALL"="black","CA"="purple","EU"="light green") ) +
  scale_alpha_manual(name='', values=c("ALL"=1,"CA"=0.7,"EU"=0.2) ) +
  scale_size_manual(name='', values=c("ALL"=2,"CA"=2,"EU"=2) ) +
  ylab("Concordance") + 
  xlab("Seasonal/latitudinal quantile") +
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.25), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        legend.title.align=.5,  
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10) )
#        axis.title=element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10) ) )
dev.off()


df2$pops <- factor(df2$pops, levels = c("ALL","CA","EU"))

pdf("figures/clinal_uniquepopsPAnoMA_CA_EUseasonalparallel_polymorphic_Dec2017_96CI_ext_ggplot_red.pdf", width=4,height=4)
g_lat_line_red =
ggplot(data=df2, aes(quant, concord, group=pops) ) +
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept = 0.5, lty=2) +
  geom_polygon(data=nCA, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="purple", alpha=0.8) +
  geom_polygon(data=nEU, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="light green", alpha=0.5) +
  geom_line(aes(quant, concord)) + 
  #scale_color_manual(name='', values=c("ALL"="blue","CA"="black","EU"="black") ) +
  geom_point(aes(quant, concord, fill=pops), color="black", shape=21, size=1.5) + 
  scale_fill_manual(name='', values=c("ALL"="red","CA"="purple","EU"="light green") ) +
  scale_alpha_manual(name='', values=c("ALL"=1,"CA"=0.7,"EU"=0.2) ) +
  scale_size_manual(name='', values=c("ALL"=2,"CA"=2,"EU"=2) ) +
  ylab("Concordance") + 
  xlab("Seasonal/latitudinal quantile") +
  ylim(0.2,1) +
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.2), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=0),
        legend.title.align=.5,  
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10) )
#        axis.title=element_text(size=10, margin = margin(t = 0, r = 0, b = 0, l = 10) ) )
dev.off()

#save(cC, eC, heat_poly, heat_sig_poly, heat_poly_con, g_lat_line, df2, nCA, nEU, heat_poly_con_red, heat_sig_poly_red, g_lat_line_red, file="figure_objects_Jan2018.Rdata")
#load("figure_objects_Jan2018.Rdata")
save(heat_poly_con_red, heat_sig_poly_red, g_lat_line_red, file="clinal_objects_red_May2018.Rdata")


####### Reversed 
library("scales")
library("ggpubr")

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

g_lat_line_red_rev =
  ggplot(data=df2, aes(quant, concord, group=pops) ) +
  scale_x_continuous(trans=reverselog_trans(10), 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept = 0.5, lty=2) +
  geom_polygon(data=nCA, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="purple", alpha=0.8) +
  geom_polygon(data=nEU, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="light green", alpha=0.5) +
  geom_line(aes(quant, concord)) + 
  #scale_color_manual(name='', values=c("ALL"="blue","CA"="black","EU"="black") ) +
  geom_point(aes(quant, concord, fill=pops), color="black", shape=21, size=1.5) + 
  scale_fill_manual(name='', values=c("ALL"="red","CA"="purple","EU"="light green") ) +
  scale_alpha_manual(name='', values=c("ALL"=1,"CA"=0.7,"EU"=0.2) ) +
  scale_size_manual(name='', values=c("ALL"=2,"CA"=2,"EU"=2) ) +
  ylab("Concordance") + 
  xlab("Seasonal/latitudinal quantile") +
  ylim(0.2,1) +
  grids(axis = "y", linetype = "solid") + 
  ggtitle("Original") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(0.1,0.2), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=0),
        legend.title.align=.5,  
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10) )

save(nCA, nEU, df2, g_lat_line_red_rev, file="clinal_objects_red_rev_Aug2018.Rdata")




############ Testing difference between original and switched clinal concordance
## Performing a Fishers exact test with multiple testing correction (bonferroni)
bayenv = read.table("environ_corr.clinal_environfile_Feb2019.txt", stringsAsFactors = F)
snps = read.table("clinal_allsnps_non0_unformatted_Feb2019.txt", stringsAsFactors = F)
ID = unlist(lapply(bayenv[,1], FUN=function(X) unlist(strsplit(unlist(strsplit(X, split="clinal_snpsfile_Feb2019_snp"))[2], split=".txt", fixed=T))[1]   ))
bayenv$ID = ID
snps$ID = 1:nrow(snps)
clinal_snps = merge(bayenv, snps, by="ID")
clinal_snps$FL = clinal_snps$V3/clinal_snps$V7
clinal_snps$GA = clinal_snps$V4/clinal_snps$V8
clinal_snps$PA = clinal_snps$V5/clinal_snps$V9
clinal_snps$SC = clinal_snps$V6/clinal_snps$V10
clinal_snps$chr_pos = paste(clinal_snps$V1.y, clinal_snps$V2.y, sep="_")

# merge old lm to get beta 
tempA = read.table("../../glm/mel_clinal_uniquepops_springPA_noMA.glm.noheader", stringsAsFactors=FALSE) ## only spring pops
my.filter = read.table("../../data/chrom_pos_medfreq01_RRgrt0.txt")
my.filterA = my.filter[my.filter[,1]!="X",]
tempB = merge(my.filterA, tempA, by=c(1,2))
temp = merge(tempB[,c(1:3)], clinal_snps[,c(4,5,3)], by=c(1,2))
colnames(temp) = c("chrom", "pos", "clinal.coef", "clinal.p")

# read in the seasonal glms
clinal = temp
allS = read.table("/Users/hm8/stanford/nescent_melCA/glm/mel_all_caF_switch_nonclinal_uniquepopsSpringPA_noMA.glm.noheader", stringsAsFactors = F)
allO = read.table("/Users/hm8/stanford/nescent_melCA/glm/mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm", header=T, stringsAsFactors = F)
caS = read.table("/Users/hm8/stanford/nescent_melCA/glm/mel_ca_popyear_paired_switch.f_s.glm", stringsAsFactors = F, header=TRUE)
caO = read.table("/Users/hm8/stanford/nescent_melCA/glm/mel_ca_popyear_paired.f_s.glm", stringsAsFactors = F, header=TRUE)


Nquant_thresh = function(x, c, quants){
  #c = clinal # columns in order, chrom, pos, coef, pval
  #x = allO # columns in order, chrom, pos, coef, pval
  #quants = quant[inds]
  # Output:
  # List: N concordant and N sites per quantile
  cx = merge(c[,1:4], x[,1:4], by=c(1,2)) # col 4: pval
  q4 = quantile(cx[,4], probs=quants, na.rm=T)
  q6 = quantile(cx[,6], probs=quants, na.rm=T)
  
  Nquant = vector()
  concord = vector()
  for (i in 1:length(quants)){
    focal = cx[cx[,4]<q4[i] & cx[,6]<q6[i], ] # pval under focal quantile
    Nquant[i] = nrow(focal)
    concord[i] = sum( (focal[,3]>=0 & focal[,5]>=0) | (focal[,3]<0 & focal[,5]<0), na.rm=T  ) # N concordant coef
  }
  list(concord, Nquant)
}

Nquant_bins = function(x, c, quants){
  #c = clinal # columns in order, chrom, pos, coef, pval
  #x = allO # columns in order, chrom, pos, coef, pval
  #quants = quant[inds2]
  # Output:
  # List: N concordant and N sites per quantile
  cx = na.omit(merge(c[,1:4], x[,1:4], by=c(1,2))) # col 4: pval
  q4 = quantile(cx[,4], probs=quants, na.rm=T)
  q6 = quantile(cx[,6], probs=quants, na.rm=T)
  
  Nquant = vector()
  concord = vector()
  for (i in 1:(length(quants)-1) ){
    focal = cx[cx[,4]<q4[i] & cx[,4]>q4[i+1] & cx[,6]<q6[i] & cx[,6]>q6[i+1], ] # pval under focal quantile
    Nquant[i] = nrow(focal)
    concord[i] = sum( (focal[,3]>=0 & focal[,5]>=0) | (focal[,3]<0 & focal[,5]<0), na.rm=T  ) # N concordant coef
  }
  list(concord, Nquant)
}

concord_difference_thresh = function(dataset1, dataset2, clinal, quants){
  #clinal = clinal # columns in order, chrom, pos, coef, pval
  #dataset1 = allO # columns in order, chrom, pos, coef, pval
  #dataset2 = allS # columns in order, chrom, pos, coef, pval
  #quants = quant[inds]
  ## Output:
  # Table, rows = quantiles, 
  #   col 1: quantiles
  #   col 2: pvalue
  #   col 3: pvalue bonferroni corrected
  #   col 4: odds ratio
  #   col 5: lower 95% CI
  #   col 6: upper 95% CI
  outmat = matrix(ncol=6, nrow=length(quants))
  colnames(outmat) = c("quantile","pval","pval.Bcor", "oddsratio", "l95CI","u95CI")
  D1quant = Nquant_thresh(dataset1, clinal, quants)
  D2quant = Nquant_thresh(dataset2, clinal, quants)
  Nq = length(quants)
  
  for (i in 1:length(quants)){
    f1 = fisher.test(matrix( c(D1quant[[1]][i], D1quant[[2]][i], D2quant[[1]][i], D2quant[[2]][i]), ncol=2) )
    outmat[i,] = c(quants[i], f1$p.value, p.adjust(f1$p.value, method = "bonferroni", n = Nq), f1$estimate, f1$conf.int[1], f1$conf.int[2])
  }
  data.frame(outmat)
}

concord_difference_bins = function(dataset1, dataset2, clinal, quants){
  #clinal = clinal # columns in order, chrom, pos, coef, pval
  #dataset1 = allO # columns in order, chrom, pos, coef, pval
  #dataset2 = allS # columns in order, chrom, pos, coef, pval
  #quants = quant[inds2]
  ## Output:
  # Table, rows = quantiles, 
  #   col 1: quantiles
  #   col 2: pvalue
  #   col 3: pvalue bonferroni corrected
  #   col 4: odds ratio
  #   col 5: lower 95% CI
  #   col 6: upper 95% CI
  outmat = matrix(ncol=6, nrow=length(quants)-1)
  colnames(outmat) = c("quantile","pval","pval.Bcor", "oddsratio", "l95CI","u95CI")
  D1quant = Nquant_bins(dataset1, clinal, quants)
  D2quant = Nquant_bins(dataset2, clinal, quants)
  Nq = length(quants)
  
  for (i in 1:(length(quants)-1) ){
    f1 = fisher.test(matrix( c(D1quant[[1]][i], D1quant[[2]][i], D2quant[[1]][i], D2quant[[2]][i]), ncol=2) )
    outmat[i,] = c(quants[i], f1$p.value, p.adjust(f1$p.value, method = "bonferroni", n = Nq), f1$estimate, f1$conf.int[1], f1$conf.int[2])
  }
  data.frame(cbind(outmat, D1quant[[1]], D1quant[[2]], D1quant[[1]]/D1quant[[2]], D2quant[[1]], D2quant[[2]], D2quant[[1]]/D2quant[[2]]))
}
concord_diff_all = concord_difference_thresh(allO, allS, clinal, quant[inds])
concord_diff_ca = concord_difference_thresh(caO, caS, clinal, quant[inds])
concord_diff_all$sig05 =  "notsig"
concord_diff_all$sig05[concord_diff_all$pval < 0.05] = "sig05"
concord_diff_ca$sig05 = "notsig"
concord_diff_ca$sig05[concord_diff_ca$pval < 0.05] = "sig05"




###############################################################
####### Figure like final figure in paper
###############################################################
quant = 10^(-seq(from=0, to=3, by=0.1))
inds=1:24
############ Plotting clouds: CA
load("/Users/hm8/stanford/nescent_melCA/clinal/bayenv/clinal_uniquepopsPAnoMA_CAparallelism_vs_controlsMar2019_bayenv_polymorphic.Rdata")
nescentClinal_Cquant_controlrCA = nescentClinal_Cquant_controlr
nescentClinal_Cquant_controlCA = nescentClinal_Cquant_control
load("/Users/hm8/stanford/nescent_melCA/clinal/bayenv/clinal_uniquepopsPAnoMA_CAparallelism_vs_controlsMar2019_bayenv_polymorphic_switch.Rdata")
nescentClinal_Cquant_controlrCAs = nescentClinal_Cquant_controlr
nescentClinal_Cquant_controlCAs = nescentClinal_Cquant_control

df1 = data.frame(quant[inds], nescentClinal_Cquant_controlrCA[inds,1], nescentClinal_Cquant_controlrCAs[inds,1])
colnames(df1) = c("quant","Original","Flipped")
df2 = data.frame(melt(df1, id=quant), sig05=c(rep("notsig",times=24), concord_diff_ca$sig05) )
colnames(df2) = c("quant","pops","concord", "sig05")

n2 = data.frame(quant, nescentClinal_Cquant_controlrCA)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nCA = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )
n2 = data.frame(quant, nescentClinal_Cquant_controlrCAs)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nCAs = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )

#ca_plot = 
  ggplot(data=df2, aes(quant, concord, group=pops) ) +
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept = 0.5, lty=2) +
  geom_polygon(data=nCA, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="red", alpha=0.5) +
  geom_polygon(data=nCAs, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="#4292c6", alpha=0.5) +
  geom_line(aes(quant, concord)) + 
  geom_point(aes(quant, concord, fill=pops), color="black", shape=21, size=1.2) + 
  #geom_point(aes(quant, concord, fill=pops, size=sig05), color="black", shape=21) + 
  #scale_size_manual(name='', values=c("notsig"=1.2,"sig05"=2) ) +
  scale_fill_manual(name='', values=c("Original"="red","Flipped"="#4292c6") ) +
  geom_text(aes(x=quant, y=concord+0.02, label="*"), data=df2[df2$sig05=="sig05",], size=2) + 		
  xlab("Seasonal / latitudinal quantile threshold") +
  ylim(0.2,1) +
  ggtitle("California") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.2), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=0),
        legend.title.align=.5,  
        #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10) ,
        plot.title = element_text(size=10),
        plot.margin = margin(5, 0, 0, 5, "pt"),
        axis.title.y = element_blank())


############ Plotting clouds: EU (no switch)
load("/Users/hm8/stanford/nescent_melCA/clinal/results/clinal_uniquepopsPAnoMA_EUparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata")
nescentClinal_Cquant_controlrEU = nescentClinal_Cquant_controlr
nescentClinal_Cquant_controlEU = nescentClinal_Cquant_control

df1 = data.frame(quant[inds], nescentClinal_Cquant_controlrEU[inds,1])
colnames(df1) = c("quant","Original")
df2 = melt(df1, id=quant)
colnames(df2) = c("quant","pops","concord")

n2 = data.frame(quant, nescentClinal_Cquant_controlrEU)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nEU = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )

#eu_plot = 
  ggplot(data=df2, aes(quant, concord, group=pops) ) +
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept = 0.5, lty=2) +
  geom_polygon(data=nEU, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="red", alpha=0.5) +
  geom_line(aes(quant, concord)) + 
  geom_point(aes(quant, concord, fill=pops), color="black", shape=21, size=1.2) + 
  scale_fill_manual(name='', values=c("Original"="red") ) +
  #ylab("") + 
  xlab("") +
  ylim(0.2,1) +
  ggtitle("Europe") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.2), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=0),
        legend.title.align=.5,  
        #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        plot.title = element_text(size=10),
        plot.margin = margin(5, 0, 0, 5, "pt"),
        axis.title.y = element_blank()
  )



############ Plotting clouds: ALL
load("/Users/hm8/stanford/nescent_melCA/clinal/results/clinal_uniquepopsPAnoMA_parallelism_vs_controlsDec2017_quant_polymorphic_ext_switch.Rdata")
nescentClinal_Cquant_controlrs = nescentClinal_Cquant_controlr
nescentClinal_Cquant_controls = nescentClinal_Cquant_control
load("/Users/hm8/stanford/nescent_melCA/clinal/results/clinal_uniquepopsPAnoMA_parallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata")

df1 = data.frame(quant[inds], nescentClinal_Cquant_controlr[inds,1], nescentClinal_Cquant_controlrs[inds,1])
colnames(df1) = c("quant","Original","Flipped")
#df2 = melt(df1, id=quant)
#colnames(df2) = c("quant","pops","concord")
df2 = data.frame(melt(df1, id=quant), sig05=c(rep("notsig",times=24), concord_diff_all$sig05) )
colnames(df2) = c("quant","pops","concord", "sig05")

n2 = data.frame(quant, nescentClinal_Cquant_controlr)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nALL = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )
n2 = data.frame(quant, nescentClinal_Cquant_controlrs)
colnames(n2) = c("quant","obs","med","lower96","upper96")
nALLs = data.frame(x=c(rev(n2$quant[inds]),n2$quant[inds]), y=c(rev(n2$upper96[inds]), n2$lower96[inds]) )

#all_plot = 
  ggplot(data=df2, aes(quant, concord, group=pops) ) +
  scale_x_continuous(trans="log10", 
                     breaks=c(10^c(-2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), 
                     labels=c(.01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept = 0.5, lty=2) +
  geom_polygon(data=nCA, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="red", alpha=0.5) +
  geom_polygon(data=nCAs, mapping=aes(x=x, y=y, group=NULL, fill=NULL), fill="#4292c6", alpha=0.5) +
  geom_line(aes(quant, concord)) + 
  geom_point(aes(quant, concord, fill=pops), color="black", shape=21, size=1.2) + 
  scale_fill_manual(name='', values=c("Original"="red","Flipped"="#4292c6") ) +
  geom_text(aes(x=quant, y=concord+0.02, label="*"), data=df2[df2$sig05=="sig05",], size=2) + 		
  ylab("Concordance") + 
  xlab("") +
  ylim(0.2,1) +
  ggtitle("All Populations") + 
  theme(legend.direction="vertical", 
        legend.justification=c("right"), 
        legend.position=c(1,0.2), 
        legend.key.size = unit(0.35, "cm"),
        legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
        legend.text=element_text(size=8),
        legend.title=element_text(size=0),
        legend.title.align=.5,  
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
        axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        plot.title = element_text(size=10),
        plot.margin = margin(5, 0, 0, 0, "pt") )







############## PERFORMED ON SHERLOCK
################ saved as files:
# clinal_uniquepopsPA_parallelism_vs_controlsNov2017.Rdata  
# clinal_uniquepopsPA_CAparallelism_vs_controlsNov2017.Rdata              
#clinal_uniquepopsPA_EUparallelism_vs_controlsNov2017.Rdata 


#allcoef2 = allcoef[allcoef[,4]<1 & allcoef[,6]<1,] ## are no sites with Pvalue==0 (or==1)
# rm(tempA, my.filter, my.filterA, temp, tempB, all_fsA, all_fsB, all_fsC, all_fs)
# 
# # combine the control SNPs with the 
# #controls = data.frame(rbindlist(controlbootlist))
# controls = read.table("../bootstrap/bootstrap_fmean_dp.mel.medfreq01_RRgrt0.recRate.txt", stringsAsFactors=FALSE)
# controls2 = merge(controls, allcoef[,1:2], by=c(1,2))
# ncontrol = ncol(controls2) - 2 # first 2 columns are read chrom and pos
# rm(controls)
# 
# nescentClinal_Npvalue_control = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control1 = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control05 = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control01 = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control005 = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control001 = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control0005 = matrix(ncol=101, nrow=8)
# nescentClinal_Npvalue_control0001 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control1 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control05 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control01 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control005 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control001 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control0005 = matrix(ncol=101, nrow=8)
# nescentClinal_Cpvalue_control0001 = matrix(ncol=101, nrow=8)
# nescentClinal_Cquant_control = matrix(ncol=101, nrow=13)
# 
# for (i in 1:101){
#   ## for each control, get the coef and pvalue for the NA other dataset
#   focal_control = controls2[,c(1,i+1,1,2)] ## the control chrom, pos, then the real chrom and pos
#   newNA = merge(focal_control, allcoef[,1:4], by=c(1,2) ) # fetch the coef and pvalue for the matched control
#   # the 4 columns of the control and then real chrom/pos, then control coef/pvalue
#   #newNA2 = newNA[order(newNA[,3], newNA[,4]),]
#   allcoef3 = na.omit(merge(newNA[,c(3:6)], allcoef[,c(1,2,5:ncol(allcoef))], by=c(1,2)))
#   rm(newNA,focal_control)
#   
#   nescentcompare = compare_glm_glm(allcoef3, c("clinal") )
#   nescentcompareP1 = compare_glm_glm(allcoef3, c("clinal") , p2=0.1)
#   nescentcompareP05 = compare_glm_glm(allcoef3, c("clinal") , p2=0.05)
#   nescentcompareP01 = compare_glm_glm(allcoef3, c("clinal") , p2=0.01)
#   nescentcompareP005 = compare_glm_glm(allcoef3,c("clinal") , p2=0.005)
#   nescentcompareP001 = compare_glm_glm(allcoef3,c("clinal") , p2=0.001)
#   nescentcompareP0005 = compare_glm_glm(allcoef3,c("clinal") , p2=0.0005)
#   nescentcompareP0001 = compare_glm_glm(allcoef3,c("clinal") , p2=0.0001)
#   
#   nescentClinal_Npvalue_control[,i] = nescentcompare[,2]
#   nescentClinal_Npvalue_control1[,i] = nescentcompareP1[,2]
#   nescentClinal_Npvalue_control05[,i] = nescentcompareP05[,2]
#   nescentClinal_Npvalue_control01[,i] = nescentcompareP01[,2]
#   nescentClinal_Npvalue_control005[,i] = nescentcompareP005[,2]
#   nescentClinal_Npvalue_control001[,i] = nescentcompareP001[,2]
#   nescentClinal_Npvalue_control0005[,i] = nescentcompareP0005[,2]
#   nescentClinal_Npvalue_control0001[,i] = nescentcompareP0001[,2]
#   
#   nescentClinal_Cpvalue_control[,i] = nescentcompare[,3]
#   nescentClinal_Cpvalue_control1[,i] = nescentcompareP1[,3]
#   nescentClinal_Cpvalue_control05[,i] = nescentcompareP05[,3]
#   nescentClinal_Cpvalue_control01[,i] = nescentcompareP01[,3]
#   nescentClinal_Cpvalue_control005[,i] = nescentcompareP005[,3]
#   nescentClinal_Cpvalue_control001[,i] = nescentcompareP001[,3]
#   nescentClinal_Cpvalue_control0005[,i] = nescentcompareP0005[,3]
#   nescentClinal_Cpvalue_control0001[,i] = nescentcompareP0001[,3]
#   
#   nescentcompareQ = compare_glm_glm_quantile_sym(allcoef3, c("clinal") )
#   nescentClinal_Cquant_control[,i] = nescentcompareQ[[1]][,3]
# }
# 
# nescentClinal_Npvalue_controlr = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control1r = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control05r = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control01r = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control00r = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control001r = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control0005r = matrix(ncol=4, nrow=8)
# nescentClinal_Npvalue_control0001r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_controlr = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control1r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control05r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control01r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control005r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control001r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control0005r = matrix(ncol=4, nrow=8)
# nescentClinal_Cpvalue_control0001r = matrix(ncol=4, nrow=8)
# nescentClinal_Cquant_controlr = matrix(ncol=4, nrow=13)
# 
# for (i in 1:8){
#   nescentClinal_Npvalue_controlr[i,] = c(nescentClinal_Npvalue_control[i,1],median(nescentClinal_Npvalue_control[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control1r[i,] = c(nescentClinal_Npvalue_control1[i,1],median(nescentClinal_Npvalue_control1[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control1[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control05r[i,] = c(nescentClinal_Npvalue_control05[i,1],median(nescentClinal_Npvalue_control05[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control05[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control01r[i,] = c(nescentClinal_Npvalue_control01[i,1],median(nescentClinal_Npvalue_control01[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control01[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control00r[i,] = c(nescentClinal_Npvalue_control005[i,1],median(nescentClinal_Npvalue_control005[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control005[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control001r[i,] = c(nescentClinal_Npvalue_control001[i,1],median(nescentClinal_Npvalue_control001[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control001[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control0005r[i,] = c(nescentClinal_Npvalue_control0005[i,1],median(nescentClinal_Npvalue_control0005[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control0005[i,2:101])[c(3,98)])
#   nescentClinal_Npvalue_control0001r[i,] = c(nescentClinal_Npvalue_control0001[i,1],median(nescentClinal_Npvalue_control0001[i,2:101], na.rm=TRUE), sort(nescentClinal_Npvalue_control0001[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_controlr[i,] = c(nescentClinal_Cpvalue_control[i,1],median(nescentClinal_Cpvalue_control[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control1r[i,] = c(nescentClinal_Cpvalue_control1[i,1],median(nescentClinal_Cpvalue_control1[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control1[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control05r[i,] = c(nescentClinal_Cpvalue_control05[i,1],median(nescentClinal_Cpvalue_control05[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control05[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control01r[i,] = c(nescentClinal_Cpvalue_control01[i,1],median(nescentClinal_Cpvalue_control01[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control01[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control005r[i,] = c(nescentClinal_Cpvalue_control005[i,1],median(nescentClinal_Cpvalue_control005[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control005[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control001r[i,] = c(nescentClinal_Cpvalue_control001[i,1],median(nescentClinal_Cpvalue_control001[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control001[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control0005r[i,] = c(nescentClinal_Cpvalue_control0005[i,1],median(nescentClinal_Cpvalue_control0005[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control0005[i,2:101])[c(3,98)])
#   nescentClinal_Cpvalue_control0001r[i,] = c(nescentClinal_Cpvalue_control0001[i,1],median(nescentClinal_Cpvalue_control0001[i,2:101], na.rm=TRUE), sort(nescentClinal_Cpvalue_control0001[i,2:101])[c(3,98)])
# }
# 
# for (i in 1:13){
#   nescentClinal_Cquant_controlr[i,] = c(nescentClinal_Cquant_control[i,1],median(nescentClinal_Cquant_control[i,2:101], na.rm=TRUE), sort(nescentClinal_Cquant_control[i,2:101])[c(3,98)])
# }
# 
# rm(tempA, my.filter, my.filterA, temp, tempB, all_fsA, all_fsB, all_fsC, all_fs, out2, result_list,pvalue,out1, NmatP, plistP, mycols,j,i,maxmin,ncontrol,chroms,clistP)
# save.image("clinal_uniquepopsPA_parallelism_vs_controlsNov2017.Rdata")



################# For only the SNPS with allele frequency 0<freq>1
