## Mar 2018

library(dplyr)
library(ggplot2)


L_list = list()
glm_list = list()
dpplus = c(0,10,100)

for (i in 1:length(dpplus)){
  L_list[[i]] = read.table(paste("results/L_rank_fisher_exact_dpplus",dpplus[i],".txt", sep=""), stringsAsFactors=FALSE, header=TRUE)
  glm_list[[i]] = read.table(paste("results/mel_all_paired20_2sample_caF_popyear_dpplus",dpplus[i],".f_s.glm", sep=""), stringsAsFactors=FALSE, header=TRUE)
}

#filter = read.table("../../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
old = read.table("/Users/hm8/stanford/nescent_melCA/fishers_method_Oct2016/jamie_fisher_method/fisher_exactJ_pop20_polymorphic/results/L_rank_fisher_exactJ.merged.20pop.polymorphic.txt", stringsAsFactors=FALSE, header=TRUE)

#enrichall_list = list()
# for (i in 1:length(dpplus)){
#   enrichall_list[[i]] = read.table(paste("results/enrich_L_rank_fisher_exact_dpplus",dpplus[i],".txt", sep=""), stringsAsFactors=FALSE, header=TRUE)
# }

quants = seq(from=0, to=1, by=0.01)
qchisq_quants = unlist(lapply(quants, FUN=qchisq, df=40))
prop_qchisq_quants = list()
for (i in 1:length(dpplus)){
  outvec = vector()
  focal = 2*unlist(L_list[[i]])
  for (q in 1:(length(quants)-1)){
    outvec[[q]] = mean(focal>=qchisq_quants[q] & focal<qchisq_quants[q+1], na.rm=TRUE)
  }
  prop_qchisq_quants[[i]] = outvec
}

oldvec = vector()
old3 = 2*unlist(old[,1])
for (q in 1:(length(quants)-1)){
  oldvec[[q]] = mean(old3>=qchisq_quants[q] & old3<qchisq_quants[q+1], na.rm=TRUE)
}

prop_glm_quants = list()
for (i in 1:length(dpplus)){
  outvec = vector()
  focal = glm_list[[i]]$seas.p
  for (q in 1:(length(quants)-1)){
    outvec[[q]] = mean(focal>=quants[q] & focal<quants[q+1], na.rm=TRUE)
  }
  prop_glm_quants[[i]] = outvec
}



# plot(quants[1:(length(quants)-1)], prop_glm_quants[[1]], log="x", type="l", ylim=c(0.008,0.02))
# plot(quants[1:(length(quants)-1)], prop_glm_quants[[1]], type="l", ylim=c(0.008,0.02))
# points(quants[1:(length(quants)-1)], prop_qchisq_quants[[1]], type="l", col="blue")
# points(quants[1:(length(quants)-1)], prop_qchisq_quants[[2]], type="l", col="red")
# points(quants[1:(length(quants)-1)], prop_glm_quants[[2]], type="l", col="grey")
# points(quants[1:(length(quants)-1)], prop_glm_quants[[3]], type="l", col="grey")
# points(quants[1:(length(quants)-1)], oldvec, type="l", lty=2)
# abline(h=0.01)


glm_df = data.frame(pvalue=quants[1:(length(quants)-1)]+0.01, "0"=prop_glm_quants[[1]], "10"=prop_glm_quants[[2]], "100"=prop_glm_quants[[3]])
fisher_df = data.frame(pvalue=quants[1:(length(quants)-1)]+0.01, "0"=prop_qchisq_quants[[1]], "10"=prop_qchisq_quants[[2]], "100"=prop_qchisq_quants[[3]])
glm_df2 = melt(glm_df[1:(nrow(glm_df)-1),], id="pvalue")
fisher_df2 = melt(fisher_df[1:(nrow(fisher_df)-1),], id="pvalue")

glm_df2 = melt(glm_df, id="pvalue")
fisher_df2 = melt(fisher_df, id="pvalue")


#pdf("glm_paired20_vs_100perm_hist_polymorphic_ggplot.pdf", width=5,height=4)
glm_dptest_plot = 
  ggplot(glm_df2, aes(pvalue, value, color=variable)) +
  geom_line(alpha=1) +
  scale_color_manual('DP+',values=c('black','darkgrey','lightgrey'), labels=c(0,10,100)) + 
  #scale_size_manual('',values=c('Permuted'=0.5,'Observed'=1)) + 
  #scale_x_log10() + 
  scale_x_continuous(trans="log10", breaks=c(10^c(-3, -2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(0.001, .01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
  geom_hline(yintercept=0.01, lty=2) +
  xlab("GLM p-value bin") + 
  ylab("Proportion of SNPs") + 
  ylim(c(0.009,0.21)) +
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
#dev.off()
  
fisher_dptest_plot = 
ggplot(fisher_df2, aes(pvalue, value, color=variable)) +
    geom_line(alpha=1) +
    scale_color_manual('DP+',values=c('black','darkgrey','lightgrey'), labels=c(0,10,100)) + 
    #scale_size_manual('',values=c('Permuted'=0.5,'Observed'=1)) + 
    #scale_x_log10() + 
    scale_x_continuous(trans="log10", breaks=c(10^c(-3, -2, -1, 0), expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1))) [,1]), labels=c(0.001, .01, .1, 1, rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])) ) ) + 
    geom_hline(yintercept=0.01, lty=2) +
    xlab(expression(paste(Chi^2, " expected quantile bin"))) + 
    ylab("Proportion of SNPs") + 
    ylim(c(0.0095,0.012)) +
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
  

save(fisher_dptest_plot, glm_dptest_plot,  file="test_samplesize_obs_exp_ggplot.Rdata")
multiPlot = plot_grid(glm_dptest_plot, fisher_dptest_plot, nrow=1, ncol=2, labels = c('A', 'B'), vjust=1.2, align='hv')
save_plot(multiPlot, file="test_samplesize_obs_exp_ggplot_Mar2018.png", base_aspect_ratio = 2.2, base_height=3)
save_plot(multiPlot, file="test_samplesize_obs_exp_ggplot_Mar2018.pdf", base_aspect_ratio = 2.2, base_height=3)

