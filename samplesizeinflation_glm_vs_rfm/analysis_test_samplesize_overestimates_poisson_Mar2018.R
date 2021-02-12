## Mar 2018
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

L_list = list()
glm_list = list()
dpplus = c(0,10,100)

for (i in 1:length(dpplus)){
  L_list[[i]] = read.table(paste("results/L_rank_fisher_exact_dpplus",dpplus[i],"_poisson.txt", sep=""), stringsAsFactors=FALSE, header=TRUE)
  glm_list[[i]] = read.table(paste("results/mel_all_paired20_2sample_caF_popyear_dpplus",dpplus[i],"_poisson.f_s.glm", sep=""), stringsAsFactors=FALSE, header=TRUE)
}

old = read.table("../results/L_rank_fisher_exactJ.merged.20pop.polymorphic.txt", stringsAsFactors=FALSE, header=TRUE)

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


glm_df = data.frame(pvalue=quants[1:(length(quants)-1)]+0.01, "0"=prop_glm_quants[[1]], "10"=prop_glm_quants[[2]], "100"=prop_glm_quants[[3]])
fisher_df = data.frame(pvalue=quants[1:(length(quants)-1)]+0.01, "0"=prop_qchisq_quants[[1]], "10"=prop_qchisq_quants[[2]], "100"=prop_qchisq_quants[[3]])
glm_df2 = melt(glm_df[1:(nrow(glm_df)-1),], id="pvalue")
fisher_df2 = melt(fisher_df[1:(nrow(fisher_df)-1),], id="pvalue")

glm_df2 = melt(glm_df, id="pvalue")
fisher_df2 = melt(fisher_df, id="pvalue")

glm_dptest_plot_pois =
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
  
fisher_dptest_plot_pois = 
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
save(fisher_dptest_plot_pois, glm_dptest_plot_pois,  file="test_samplesize_obs_exp_poisson_ggplot.Rdata")

#fisher_dptest_plot, glm_dptest_plot,  
load(file="test_samplesize_obs_exp_ggplot.Rdata")
multiPlot = plot_grid(fisher_dptest_plot+theme_cowplot()+, glm_dptest_plot+theme_cowplot(), fisher_dptest_plot_pois+theme_cowplot(), glm_dptest_plot_pois+theme_cowplot(), nrow=2, ncol=2, labels = c('A', 'B','C','D'), vjust=1.1, align='hv')
save_plot(multiPlot, file="test_samplesize_set_poisson_obs_exp_ggplot_4panel_Mar2018.pdf", base_aspect_ratio = 1.4, base_height=8)
