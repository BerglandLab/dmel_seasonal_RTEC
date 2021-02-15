### libraries
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(foreach)
  library(patchwork)

### set working directory
  setwd("/Users/alanbergland/Documents/work/Projects/2005-2019/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/actualMS/Alan's data for the MS")

### set daysPrior (only options are 14 and 21)
  daysPrior.use <- 21

##############
### Fig 3A ###
##############
  ### days prior irrelevant
  source("/Users/alanbergland/Documents/GitHub/drosRTEC_revisions/figure3/Fig3a_use.R")
  fig3a

##############
### Fig 3B ###
##############
  ### days prior relevant
  source("/Users/alanbergland/Documents/GitHub/drosRTEC_revisions/figure3/Fig3b_use.R")
  fig3b

##############
### Fig 3C ###
##############
  source("/Users/alanbergland/Documents/GitHub/drosRTEC_revisions/figure3/Fig3c_use.R")
  fig3c

##############
### Fig 3D ###
##############
  source("/Users/alanbergland/Documents/GitHub/drosRTEC_revisions/figure3/Fig3d_use.R")
  #fig3d.raw
  fig3d.loess

###############
#### Fig 3E ###
###############
#  source("/Users/alanbergland/Documents/GitHub/drosRTEC_revisions/figure3/Fig3e_use.R")
#  fig3e


####################
### Merge & plot ###
####################


#### version 2
  layout <- "
  AACCCC
  BBBDDD
  "

  fig3 <-
  fig3a + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
  fig3b + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
  fig3c + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
  fig3d.loess + theme(axis.text=element_text(size=16), axis.title=element_text(size=18)) +
  plot_layout(design=layout) +
  plot_annotation(tag_levels="A")
  #fig3

  ggsave(fig3, file=paste("~/Fig3_alt.fix", daysPrior.use, "daysPrior.pdf", sep=""), h=8.5, w=12)








  layout <- "
  AAACCCDDD
  AAACCCDDD
  AAACCCDDD
  BBBCCCEEE
  BBBCCCEEE
  BBBCCCEEE
  "

  fig3 <- fig3a + fig3b + fig3c + fig3d.loess + fig3e +
  plot_layout(design=layout) +
  plot_annotation(tag_levels="A")


  ggsave(fig3, file="~/Fig3.pdf", h=8.5, w=12)










        fig3AB <- plot_grid(fig3a,
                            fig3b,
                            nrow=2, labels=c("A", "B"),
                            align="hv", axis="tblr")

        fig3ABC <- plot_grid(fig3AB,
                            fig3c,
                            ncol=2,
                            rel_widths=c(1,1.3),
                            labels=c("", "C"))


        fig3DE <- plot_grid(fig3d.loess,
                            fig3e,
                            nrow=2,
                            labels=c("D", "E"),
                            align="hv", axis="tblr")

        fig3 <- plot_grid(fig3ABC, fig3DE,
                ncol=2,
                rel_widths=c(2.5,1),
                align="hv", axis="tblr")

        fig3
