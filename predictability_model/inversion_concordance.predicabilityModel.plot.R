scp aob2x@rivanna.hpc.virginia.edu:~/loo_inv.Rdata ~/.

### libraries
	library(ggplot2)
	library(data.table)
	library(cowplot); theme_set(theme_cowplot())

### load
	load("~/loo_inv.Rdata")
	o.ag <- o[win==500000][n>25][,list(pr=mean(frac>.5)), list(target, inv, win, log.sp.q.th)]
	o.ag[,y:=1]

### set working directory
  setwd("/Users/alanbergland/Documents/work/Projects/2005-2019/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/actualMS/Alan's data for the MS")

### get whole genome
  load(file="Fig3/m.fix.Rdata")
  pops <- unique(m[season.y=="spring"]$pop)

  #load(file="/mnt/pricey_1/dropPop/19_drop_1_concordance_data.Rdata") ### o19
  load("19_drop_1_concordance_data.Rdata")
    o19[,type := "Core 20 LOOCV"]


  #load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.Rdata")  ### o20
  load("20_add_1_concordance_data.Rdata")

  o20[pop%in%c("PA_12", "WI_13", "PA_14", "PA_9"), type := "New set CV"]
  o20[is.na(type) & pop!="PA_2015" , type := "frost"]

  concordDat <- rbind(o19, o20[p<.05][n>=50], fill=T)


	concordDat.ag <- concordDat[,list(pr=mean(frac>.5), y=1.1), list(log.sp.q.th)]



### plot
	ggplot() +
	geom_vline(xintercept=-1, linetype="dashed") +
	geom_line(data=concordDat[type=="Core 20 LOOCV"][pop%in%pops],
						 aes(y=frac, x=log.sp.q.th, group=pop, color=beta)) +
	geom_line(data=o[win==500000][n>25], aes(x=log.sp.q.th, y=frac, group=pop)) +
	geom_label(data=o.ag[log.sp.q.th==-1], aes(x=log.sp.q.th, y=y, label=round(pr, 2))) +
	geom_label(data=concordDat.ag[log.sp.q.th==-1], aes(x=log.sp.q.th, y=y, label=round(pr, 2), fill="orange")) +
	scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
							space = "Lab", na.value = "grey50", guide = "colourbar",
							limits=range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta),
							name=NULL) +
	facet_grid(target~inv)
