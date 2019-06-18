######### March 2018
### New PCA (using prcomp)

library(ggplot2)
library(RColorBrewer)
library(cowplot)

pops = read.table("data/paired_spring_fall_populations_noPA12.txt", stringsAsFactors=FALSE)[,1]
pops = c(pops, "PA_12")
load("data/mel_freqdp_042016_fixed.Rdata")

popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

####### Change this to target analysis ########################
focalA = which( (S=="s" | S=="f") & PY %in% pops)  # 25 paired samples (including PA_12)
focalB = which(popID=="melSC_072010_SPT" | popID=="melGA_072008_SPT" | popID=="melFL_072010_SPT")  # clinal (popID=="melPA_72011_SPT" already in seasonal)
################################################################
focal = c(focalA, focalB)

freqfocal = freq[, focal]
both=cbind(info[,1:2], freqfocal)
P1 = P[focal]
Y1 = Y[focal]
R1 = R[focal]
S1 = S[focal]
PY1 = PY[focal]
PYS1 = vector()
for (i in 1:length(PY1)){
  PYS1[i] = paste(PY1[i],S1[i], sep="")
}

lats = c(38.0, 38.7, 38.7, 38.7, 38.7, 38.0, 41.4, 41.4, 42.5, 42.5, 39.9, 42.4, 42.4,
         39.9, 39.9, 39.9, 39.9, 39.9, 39.9, 39.9, 43.1, 46.5, 46.5, 43.1, 48.2, 43.1, 43.1, 48.2, 
         34.0, 34.0, 38.0, 38.0, 43.1, 43.1, 42.5, 42.5, 39.9, 39.9, 40.8, 40.8, 39.2,
         39.2, 38.0, 38.0, 39.9, 39.9, 46.3, 46.3, 42.1, 42.1, 25.5, 31.0, 33.4)

mypop = P1
mypop[P1=="OUK"] = "UA_od"
mypop[P1=="VI"] = "AT_gr"
mypop[P1=="BA"] = "ES_ba"
mypop[P1=="rd"] = "CA_es"
mypop[P1=="co"] = "CA_tu"
mypop[P1=="SON"] = "ON_su"
mypop[P1=="WI"] = "WI_cp"
mypop[P1=="CWI"] = "WI_cp"
mypop[P1=="MA"] = "MA_la"
mypop[P1=="LMA"] = "MA_la"
mypop[P1=="NY"] = "NY_it"
mypop[P1=="BHM"] = "MI_bh"
mypop[P1=="SCPA"] = "PA_sc"
mypop[P1=="PA"] = "PA_li"
mypop[P1=="TKA"] = "KA_to"
mypop[P1=="CUA"] = "VA_ch"
mypop[P1=="AGA"] = "GA_at"
mypop[P1=="SC"] = "SC_eu"
mypop[P1=="GA"] = "GA_ha"
mypop[P1=="FL"] = "FL_ho"

mycol = brewer.pal(n=9, "Spectral")
mycol2 = colorRampPalette(mycol)(13)
mycolEUR = c("#cccccc","#969696","#525252")
mycolCA = c("#9e9ac8","#6a51a3")
mycol3 = c(mycolEUR, mycolCA, mycol2)
f1 = c("UA_od",  "AT_gr", "ES_ba",   "CA_es", "CA_tu", "ON_su",  "WI_cp",   "MA_la",   "NY_it",   "MI_bh", "PA_sc", "PA_li",   "KA_to",  "VA_ch",    "GA_at",   "SC_eu",   "GA_ha",  "FL_ho")  

nsamples = length(S1) # 51
filter = read.table("data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
filtered = merge(both, filter, by=c(1,2))
filtered2 = na.omit(filtered[filtered[,1] != "X", 3:ncol(filtered)])
colnames(filtered2) = PYS1
filtered2zeros = apply(filtered2, MARGIN=1, FUN=function(X) any(X==0 | X==1))
filtered3 = filtered2[!filtered2zeros, ]

S2 = S1
S2[S1=="s"] = "spring"
S2[S1=="f"] = "fall"

# SNPs on columns, populations on rows
prcomp1 = prcomp(t(filtered3), scale=TRUE)
df1 =data.frame(pc1=prcomp1$x[,1], pc2=prcomp1$x[,2], prcomp1$x[,3:51], sample= rownames(prcomp1$x), Season=S2, popyear=PY1, lat=lats, pop=mypop)
df2 = df1[order(df1$lat),]
df2$pop = factor(df2$pop, levels=f1)
save(prcomp1, df2, file="prcomp_df2_objects_PCA_Mar2018_correct.Rdata")
summary(prcomp1)
#                           PC1       PC2      PC3      PC4       PC5
# Standard deviation     257.4572 189.05360 157.1605 153.7745 143.94782
# Proportion of Variance   0.1197   0.06454   0.0446   0.0427   0.03741
# Cumulative Proportion    0.1197   0.18422   0.2288   0.2715   0.30893

pca_plot_prcomp =
ggplot(data=df2, aes(pc1, pc2)) + 
  geom_point(aes(pc1, pc2, color=pop), size=1.5) +
  scale_color_manual('Location', values = mycol3) + 
  #guides(fill=guide_legend(override.aes=list(colour=c(a="green",b="red"))))
  geom_point(aes(fill=pop, shape=Season), size=2) +
  scale_fill_manual('Location', values = mycol3, guide=FALSE) + 
  scale_shape_manual('Season', values = c(fall=21,spring=24)) + 
  ylab("PC 2") + 
  xlab("PC 1") +
  theme(#legend.direction="vertical", 
    #legend.justification=c("left"), 
    #legend.position=c(0.01,0.9), 
    legend.key.size = unit(0.35, "cm"),
    legend.background=element_rect(fill="transparent", colour="white", linetype="solid", size=.25),
    legend.text=element_text(size=8),
    legend.title=element_text(size=10),
    legend.title.align=.5,  
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 12)),
    axis.text=element_text(size=8),
    axis.title=element_text(size=10) )
save_plot(pca_plot_prcomp, file="PCA_allsamples_ggplot_prcomp_Mar2018.pdf", base_aspect_ratio = 1.6, base_height=3)
