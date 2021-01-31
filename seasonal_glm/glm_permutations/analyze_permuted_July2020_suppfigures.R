## Nov 2020

## for 500 permutations
## First part done on farm
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)

## Does the analyses per chromosome and interior/exterior to inversions
# Inversion Chromosome Proximal   Distal
# 1   In(2L)t         2L  2225744 13154180
# 2  In(2R)NS         2R 11278659 16163839
# 3   In(3L)P         3L  3173046 16301941
# 4   In(3R)K         3R  7576289 21966092
# 5  In(3R)Mo         3R 17232639 24857019
# 6   In(3R)P         3R 12257931 20569732
# chrom    start      end
# 1    2L  2225744 13154180
# 2    2R 11278659 16163839
# 3    3L  3173046 16301941
# 4    3R  7576289 24857019


###### Enrichment, by by pvalues not quant
pvalueN_mean_obs_less_perm_chromv2 = read.table("results/pvalueN_mean_obs_less_perm_chrom_v2_pvalue004.txt", header=T, stringsAsFactors = F)


##### Getting the number of SNPs per region
seasonal = read.table("../results/mel_all_paired20_2sample_caF_popyear.f_s.glm", header=T, stringsAsFactors = F)
myfilter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
seasonal_filter = merge(na.omit(seasonal), myfilter, by=c(1,2))
seasonal_filter_inv = subset(seasonal_filter, (chrom=="2L" & pos>2225744 & pos<13154180) |
                               (chrom=="2R" & pos>11278659 & pos<16163839) |
                               (chrom=="3L" & pos>3173046 & pos<16301941) |
                               (chrom=="3R" & pos>7576289 & pos<24857019) )
seasonal_filter_noninv = subset(seasonal_filter, (chrom=="2L" & pos<2225744 | pos>13154180) |
                                  (chrom=="2R" & pos<11278659 | pos>16163839) |
                                  (chrom=="3L" & pos<3173046 | pos>16301941) |
                                  (chrom=="3R" & pos<7576289 | pos>24857019) )



###### Enrichment, by by pvalues not quant, per inversion, and isolating the inversion breakpoints, 500Kb around each breakpoint
pvalueN_mean_obs_less_perm_chromv2_inv = read.table("results/obs_perm_filter_pvalue004N_byregion_inversionbrkpts_500Kb.txt", header=T, stringsAsFactors = F)
pvalueN_mean_obs_less_perm_chromv2_inv$inv = factor(pvalueN_mean_obs_less_perm_chromv2_inv$inv, levels=unique(pvalueN_mean_obs_less_perm_chromv2_inv$inv))
pvalueN_mean_obs_less_perm_chromv2_inv$class = factor(pvalueN_mean_obs_less_perm_chromv2_inv$class, levels=c("bkpts", "interior",   "exterior" ))
myregions_sizes = read.table(file="myregions_sizes.txt", header=T, stringsAsFactors = F)
myregions_sizes$length_Mb = round(myregions_sizes$length_bp/10^6)



######## Inv PLUS full chrom
#pvalueN_mean_obs_less_perm_chromv2_inv$class = factor(pvalueN_mean_obs_less_perm_chromv2_inv$class, levels=c("bkpts", "interior",   "exterior" ))
chronly = subset(pvalueN_mean_obs_less_perm_chromv2, group=="all_regions")
chronly$class = "chr"
chronly$inv = "In(2L)t"
chronly$inv[chronly$chrom=="2R"] = "In(2R)NS"
chronly$inv[chronly$chrom=="3R"] = "In(3R)KMoP"
chronly$inv[chronly$chrom=="3L"] = "In(3L)P"
chronly$inv[chronly$chrom=="all_chroms"] = "all"
chr_inv = rbind(pvalueN_mean_obs_less_perm_chromv2_inv[,c( "inv" , "class" ,"pvalue", "obs" ,"mean_enrich" ,   "median_enrich" , "lower95_enrich", "upper95_enrich" ,"lower90_enrich", "upper90_enrich","mean_obs_grt_perm", "pvalue_obs_grt_perm")], chronly[,c( "inv" , "class" ,"pvalue", "obs" ,"mean_enrich" ,   "median_enrich" , "lower95_enrich", "upper95_enrich" ,"lower90_enrich", "upper90_enrich","mean_obs_grt_perm", "pvalue_obs_grt_perm")])
# myregions_sizes = read.table(file="/Users/hm8/volumes/hm8_network/nescent_melCA/glm_permutations/myregions_sizes.txt", header=T, stringsAsFactors = F)
# myregions_sizes$length_Mb = round(myregions_sizes$length_bp/10^6)
chr_inv$chrom = revalue(chr_inv$inv, c("In(2L)t" = "2L", "In(2R)NS" = "2R", "In(3L)P" = "3L", "In(3R)K" = "3R", "In(3R)KMoP" = "3R", "In(3R)Mo" = "3R", "In(3R)P" = "3R") )
chr_inv$class = revalue(chr_inv$class, c("chr"= "all") )
sub1 = subset(chr_inv, lower90_enrich>1 & pvalue>=0.0001 & pvalue<1 & !(inv %in% c("In(3R)K",    "In(3R)Mo",    "In(3R)P")) )

ggplot(subset(chr_inv, pvalue>=0.0001 & pvalue<1 & !(inv %in% c("In(3R)K",    "In(3R)Mo",    "In(3R)P")) ), aes(pvalue, median_enrich))+
  geom_line()+
  geom_errorbar(aes(pvalue, ymin=lower95_enrich, ymax=upper95_enrich), col="grey", size=0.2)+ #width=0.5
  #geom_text(size = 2, data = subset(text2, !(inv %in% c("In(3R)K",  "In(3R)Mo",    "In(3R)P")) ), mapping = aes(x = Inf, y = Inf, label = paste(SNPs_K, "K SNPs, ", length_Mb, "Mb", sep="")), hjust   = 1.05, vjust   = 1.5)+
  geom_hline(yintercept=1, linetype="dashed")+
  facet_grid(class~inv)+
  scale_x_continuous(trans='log10', limits=c(0.00005,1))+
  #scale_y_continuous(limits=c(0.5,4))+
  theme_light()+
  ylab("Seasonal enrichment (# SNPs obs/perm)")+
  xlab("p-value")+
  geom_text(aes(pvalue, upper95_enrich+0.1), data = subset(sub1,!(inv %in% c("In(3R)K",  "In(3R)Mo",    "In(3R)P")) ), label = "*")+
  #ylim(c(0,2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )
ggsave("seasonalenrich_bychrom_v2_pvalue_inversionbrkpts_500Kb_fullchr.pdf", width=8, height=6)
ggsave("seasonalenrich_bychrom_v2_pvalue_inversionbrkpts_500Kb_fullchr.png", width=8, height=6)


sub2 = subset(chr_inv, lower90_enrich>1 & pvalue==0.001 & !(inv %in% c("In(3R)K",    "In(3R)Mo",    "In(3R)P")) )
chr_inv$sig_lower90 = ""
chr_inv$sig_lower90[chr_inv$lower90_enrich>1] = "*"
chr_inv$pvaluelabel = ""
chr_inv$pvaluelabel[chr_inv$pvalue_obs_grt_perm<0.1] = "."
chr_inv$pvaluelabel[chr_inv$pvalue_obs_grt_perm<0.05] = "*"
chr_inv$pvaluelabel[chr_inv$pvalue_obs_grt_perm<0.01] = "**"


ggplot(subset(chr_inv, pvalue==0.001 & !(inv %in% c("In(3R)K",    "In(3R)Mo",    "In(3R)P")) ), aes(x=inv, weight=median_enrich, ymin=lower95_enrich, ymax=upper95_enrich, fill=class)) +
  geom_bar      (position=position_dodge(), aes(y=median_enrich), stat="identity") +
  geom_errorbar (position=position_dodge(width=0.9), colour="grey", width=0.5) +
  geom_text    (position=position_dodge(width=0.9), aes(y=0.05, label=sig_lower90), size=8)+
  #geom_text    (position=position_dodge(width=0.9), aes(y=median_enrich+0.05, label=sig_lower90), size=8)+
  geom_hline(yintercept=1, linetype="dashed")+
  scale_fill_brewer("", palette = "Set2")+
  theme_light()+
  ylab("Seasonal enrichment (# SNPs obs/perm)")+
  xlab("Region")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )
ggsave("seasonalenrich_bychrom_v2_pvalue_inversionbrkpts_500Kb_fullchr_barplotp001.pdf", width=6, height=5)
ggsave("seasonalenrich_bychrom_v2_pvalue_inversionbrkpts_500Kb_fullchr_barplotp001.png", width=8, height=6)


ggplot(subset(chr_inv, pvalue==0.003894319 & !(inv %in% c("In(3R)K",    "In(3R)Mo",    "In(3R)P")) ), aes(x=inv, weight=median_enrich, ymin=lower95_enrich, ymax=upper95_enrich, fill=class)) +
  geom_bar      (position=position_dodge(), aes(y=median_enrich), stat="identity") +
  geom_errorbar (position=position_dodge(width=0.9), colour="grey", width=0.5) +
  geom_text    (position=position_dodge(width=0.9), aes(y=0.05, label=sig_lower90), size=8)+
  #geom_text    (position=position_dodge(width=0.9), aes(y=median_enrich+0.05, label=sig_lower90), size=8)+
  geom_hline(yintercept=1, linetype="dashed")+
  scale_fill_brewer("", palette = "Set2")+
  theme_light()+
  ylab("Seasonal enrichment (# SNPs obs/perm)")+
  xlab("Region")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() )
write.table(chr_inv, file="chr_inv_enrichment_pvalues_Jan2021.txt", col.names=T, row.names=F, quote=F, sep="\t")
ggsave("seasonalenrich_bychrom_v2_pvalue_inversionbrkpts_500Kb_fullchr_barplotp004.pdf", width=6, height=5)
ggsave("seasonalenrich_bychrom_v2_pvalue_inversionbrkpts_500Kb_fullchr_barplotp004.png", width=8, height=6)


# subset(chr_inv, pvalue==0.003894319 & !(inv %in% c("In(3R)K",    "In(3R)Mo",    "In(3R)P")) & pvalue_obs_grt_perm<=0.1)
# inv    class      pvalue        obs mean_enrich median_enrich lower95_enrich upper95_enrich
# 18      In(2L)t    bkpts 0.003894319 0.01412543    1.811800      1.849154      0.9859881       2.495575
# 54      In(3L)P    bkpts 0.003894319 0.01084980    1.315378      1.302702      1.0420830       1.651689
# 234     In(2L)t interior 0.003894319 0.01133474    1.349652      1.384193      0.9009021       1.654088
# 342  In(3R)KMoP    bkpts 0.003894319 0.01161282    1.240338      1.237460      0.8724551       1.585714
# 378  In(3R)KMoP interior 0.003894319 0.01153358    1.243103      1.258382      0.9583038       1.457294
# 396         all    bkpts 0.003894319 0.01158368    1.286676      1.293536      1.0286083       1.520139
# 432         all interior 0.003894319 0.01025779    1.161460      1.166804      0.9748726       1.294318
# 1810    In(2L)t      all 0.003894319 0.01116054    1.300860      1.322546      0.8899287       1.549032
# 721  In(3R)KMoP      all 0.003894319 0.01093348    1.177174      1.197526      0.9203795       1.363479
# 2341        all      all 0.003894319 0.01000058    1.126827      1.131445      0.9801304       1.230325
# lower90_enrich upper90_enrich mean_obs_grt_perm pvalue_obs_grt_perm chrom pvaluelabel
# 18        1.1363631       2.369629               486          0.02994012    2L           *
#   54        1.0752106       1.608516               495          0.01197605    3L           *
#   234       0.9569711       1.607914               466          0.06986028    2L           .
# 342       0.9453250       1.520757               464          0.07385230    3R           .
# 378       1.0140660       1.433116               476          0.04990020    3R           *
#   396       1.0843214       1.469120               494          0.01397206   all           *
#   432       1.0250926       1.280738               480          0.04191617   all           *
#   1810      0.9755944       1.512996               470          0.06187625    2L           .
# 721       0.9700327       1.335089               461          0.07984032    3R           .
# 2341      1.0128924       1.217957               484          0.03393214   all           *


#### testing for difference in enrichment among chromosomes
library(data.table)
seas = read.table("../results/mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm", stringsAsFactors = F, header=T)
filter = read.table("../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", stringsAsFactors = F, header=F)
seas_filter = merge(seas, filter, by=c(1,2))
snpschrom = data.frame(enrich=subset(chr_inv, class=="all" & pvalue==0.003894319)$mean_enrich[1:4], inv=
             subset(chr_inv, class=="all" & pvalue==0.003894319)$inv[1:4], N=table(seas_filter$chrom))
#           inv class      pvalue         obs mean_enrich median_enrich lower95_enrich upper95_enrich
# 1810    In(2L)t   all 0.003894319 0.011160535   1.3008598     1.3225457      0.8899287       1.549032
# 3610   In(2R)NS   all 0.003894319 0.009012914   0.9824823     0.9953951      0.8044230       1.105250
# 541     In(3L)P   all 0.003894319 0.008679926   1.0427718     1.0436152      0.9563403       1.126967
# 721  In(3R)KMoP   all 0.003894319 0.010933485   1.1771742     1.1975259      0.9203795       1.363479
# 2341        all   all 0.003894319 0.010000581   1.1268273     1.1314452      0.9801304       1.230325
# lower90_enrich upper90_enrich mean_obs_grt_perm pvalue_obs_grt_perm chrom pvaluelabel
# 1810      0.9755944       1.512996               470          0.06187625    2L           .
# 3610      0.8381220       1.082881               232          0.53692615    2R            
# 541       0.9668557       1.118099               403          0.19560878    3L            
# 721       0.9700327       1.335089               461          0.07984032    3R           .
# 2341      1.0128924       1.217957               484          0.03393214   all           *

# obs_perm500_enrich_chrom_v2_pvalue004_allchroms
load("results/obs_perm500_enrich_chrom_v2_pvalue004.Rdata")
quant = c(10^(-seq(from=0, to=4, by=1/4)), 0.003894319)
names(obs_perm500_enrich_chrom_v2_pvalue004_allchroms) = c("2L","2R","3L","3R")
rownames(obs_perm500_enrich_chrom_v2_pvalue004_allchroms) = quant
perm500_allchrom = cbind(c("2L","2R","3L","3R"), rbind(obs_perm500_enrich_chrom_v2_pvalue004_allchroms[[1]][18,],
                                         obs_perm500_enrich_chrom_v2_pvalue004_allchroms[[2]][18,],
                                         obs_perm500_enrich_chrom_v2_pvalue004_allchroms[[3]][18,],
                                         obs_perm500_enrich_chrom_v2_pvalue004_allchroms[[4]][18,]) )
colnames(perm500_allchrom)[1] = "chrom"
perm500_allchrom_melt = melt(perm500_allchrom, id.vars="chrom")
res.aov = aov(value~chrom , data=perm500_allchrom_melt)
summary(res.aov)
TukeyHSD(res.aov)
# chrom          3  30.36  10.120   888.9 <2e-16 ***
#   Residuals   1996  22.72   0.011                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## the reported p-value is 0
# $chrom
#             diff         lwr         upr p adj
# 2R-2L -0.3183775 -0.33572839 -0.30102657     0
# 3L-2L -0.2580880 -0.27543888 -0.24073706     0
# 3R-2L -0.1236856 -0.14103652 -0.10633470     0
# 3L-2R  0.0602895  0.04293859  0.07764041     0
# 3R-2R  0.1946919  0.17734096  0.21204278     0
# 3R-3L  0.1344024  0.11705145  0.15175327     0
