### libraries
	library(data.table)
	library(ggplot2)
  library(cowplot)
  library(doMC)
  registerDoMC(20)
  library(speedglm)


### dummy plot
	spring <- ggplot(data=weather.pop[pop=="MA_12"][season=="spring"]) +
						geom_line(aes(x=-1*daysPrior, y=tmax)) +
						geom_point(aes(x=-1*daysPrior, y=tmax, color=tmax>330)) +
						geom_hline(yintercept=330)

	fall <- ggplot(data=weather.pop[pop=="MA_12"][season=="fall"]) +
						geom_line(aes(x=-1*daysPrior, y=tmin)) +
						geom_point(aes(x=-1*daysPrior, y=tmin, color=tmin<10)) +
						geom_hline(yintercept=10)

	ggsave("~/temperature.pdf", plot_grid(spring, fall, nrow=1, labels=c("spring", "fall")))

### load data
  load(file="/mnt/pricey_1/dropPop/altModels_out.Rdata")


	### 14 day data
		load(file="/mnt/pricey_1/dropPop/thermalLimitModels_univariate_out.14.Rdata")
  	load(file="/mnt/pricey_1/dropPop/thermalLimitModels_bivariate_out.14.Rdata")

		o.tlm.multi.ag[,days:=14]
		o.tlm.uni.ag[,days:=14]

		o.tlm.multi.ag.14 <- o.tlm.multi.ag
		o.tlm.uni.ag.14 <- o.tlm.uni.ag

	### 21 day data
		load(file="/mnt/pricey_1/dropPop/thermalLimitModels_univariate_out.Rdata")
		load(file="/mnt/pricey_1/dropPop/thermalLimitModels_bivariate_out.Rdata")

		o.tlm.multi.ag[,days:=21]
		o.tlm.uni.ag[,days:=21]

		o.tlm.multi.ag <- rbind(o.tlm.multi.ag.14, o.tlm.multi.ag)
		o.tlm.uni.ag <- rbind(o.tlm.uni.ag.14, o.tlm.uni.ag)


		o.tlm.multi.ag[,modelSet:="tlm_bivariate"]
		o.tlm.uni.ag[,modelSet:="tlm_univariate"]


### save
	save(o.tlm.multi.ag, o.tlm.uni.ag, o.altModels, file="~/mods.Rdata")
	scp bergland@bergland-lab.bio.virginia.edu:~/mods.Rdata ~/.

	load("~/mods.Rdata")





















########
### contort data into same forms
	str(o.tlm.uni.ag)
	str(o.tlm.multi.ag)

	m <- rbind(o.tlm.uni.ag[,c("boot.num", "set", "x", "r2", "modelSet", "days"), with=F],
						 o.tlm.multi.ag[,c("boot.num", "set", "x", "r2", "modelSet", "days"), with=F] )

	m <- rbind(m, o.altModels[,c("boot.num", "set", "x", "r2", "modelSet"), with=F], fill=T)
	m[modelSet=="aveThermalLimit.bivariate", modelSet:="AveMinMax_bi"]
	m[modelSet=="aveThermalLimit.univariate", modelSet:="AveMinMax_uni"]
	m[modelSet=="tlm_bivariate", modelSet:="tlm_bi"]
	m[modelSet=="tlm_univariate", modelSet:="tlm_uni"]
	m[modelSet=="cumulDegreeDays", modelSet:="cdd"]

	m[modelSet=="tlm_bi" & x=="m1", x:="S-Min;F-Min"]
	m[modelSet=="tlm_bi" & x=="m2", x:="S-Max;F-Max"]
	m[modelSet=="tlm_bi" & x=="m3", x:="S-Min;F-Max"]
	m[modelSet=="tlm_bi" & x=="m4", x:="S-Max;F-Min"]

	m[modelSet=="tlm_uni" & x=="m1", x:="S-Min"]
	m[modelSet=="tlm_uni" & x=="m2", x:="S-Max"]
	m[modelSet=="tlm_uni" & x=="m3", x:="F-Min"]
	m[modelSet=="tlm_uni" & x=="m4", x:="F-Mmax"]


	m[,modelSet:=factor(modelSet, levels=c("tlm_bi", "tlm_uni",
																				 "aveTemp",
																				 "AveMinMax_bi", "AveMinMax_uni",
																				 "cdd",
																				 "geography",
																				 "method",
																				 "simulans",
																				 "substrate"))]


	ggplot() +
	geom_boxplot(data=m[set=="Perm"], aes(y=r2, x=x, fill=as.factor(days)), size=.5) +
	geom_point(data=m[set=="Obs"], aes(y=r2, x=x, fill=as.factor(days)),
						color="red", size=3, position=position_dodge(width=0.75)) +
	scale_fill_manual(values=c("white", "lightgrey")) + xlab("") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  +
	ylim(0,1) +
	facet_grid(.~modelSet, scales="free_x", space="free_x")






### function
  makePlot <- function(dat, mod) {
    tmp <- dat[modelSet==mod]

    pr <- tmp[,list(pr=mean(r2[set=="Obs"]>r2[set=="Perm"]), set="Obs"), list(x)]
    setkey(tmp, set, x)
    setkey(pr, set, x)
    tmp <- merge(tmp, pr, all.x=T)

    ggplot() +
    geom_density(data=tmp[set=="Perm"], aes(r2), fill="lightblue", n=10000) +
    geom_vline(data=tmp[set=="Obs"], aes(xintercept=r2), color="red", size=1) +
    geom_text(data=tmp[set=="Obs"], aes(x=r2, y=3, label=pr), hjust=-.5) +
    facet_grid(modelSet~x) +
    xlim(-.2,1.2)
  }

	makePlot2 <- function(dat, mod, xlabs) {
		#dat <- o.tlm.multi.ag; mod="tlm_bivariate"; xlabs=c("S-Min;F-Min", "S-Max;F-Max", "S-Min;F-Max", "S-Max;F-Min")
		tmp <- dat[modelSet==mod]

		pr <- tmp[,list(pr=mean(r2[set=="Obs"]>r2[set=="Perm"]), set="Obs"), list(x)]
		setkey(tmp, set, x)
		setkey(pr, set, x)
		tmp <- merge(tmp, pr, all.x=T)

		ggplot() +
		geom_boxplot(data=tmp[set=="Perm"], aes(y=r2, x=x, fill=as.factor(days)), size=.5) +
		geom_point(data=tmp[set=="Obs"], aes(y=r2, x=x, fill=as.factor(days)),
							color="red", size=5, position=position_dodge(width=0.75)) +
		scale_fill_manual(values=c("white", "lightgrey")) +
		scale_x_discrete(labels=xlabs) + xlab("") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))  + ylim(0,1)


	}


	makePlot3 <- function(dat, mod, xlabs) {
			#dat <- o.altModels; mod="tlm_bivariate"; xlabs=c("S-Min;F-Min", "S-Max;F-Max", "S-Min;F-Max", "S-Max;F-Min")
			tmp <- dat		#[modelSet==mod]

			pr <- tmp[,list(pr=mean(r2[set=="Obs"]>r2[set=="Perm"]), set="Obs"), list(x)]
			setkey(tmp, set, x)
			setkey(pr, set, x)
			tmp <- merge(tmp, pr, all.x=T)


			ggplot() +
			geom_boxplot(data=tmp[set=="Perm"], aes(y=r2, x=x ), size=.5) +
			geom_point(data=tmp[set=="Obs"], aes(y=r2, x=x),
								color="red", size=5, position=position_dodge(width=0.75)) +
			xlab("") +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
			facet_grid(.~modelSet, scales="free_x") + ylim(0,1)

		}


### plots
	#thermal limit models
	  tlm.bi <- makePlot2(dat=o.tlm.multi.ag, mod="tlm_bivariate", xlabs=c("S-Min;F-Min", "S-Max;F-Max", "S-Min;F-Max", "S-Max;F-Min"))
		tlm.uni <- makePlot2(dat=o.tlm.uni.ag, mod="tlm_univariate", xlabs=c("S-Min", "S-Max", "F-Min", "F-Max"))

	### ecology plot
		alt.plot <- makePlot3(dat=o.altModels)

		plot_grid(alt.plot, tlm.uni, tlm.bi, nrow=1)















	  geography.plot <- makePlot(dat=o.altModels, mod="geography")
    cumulDegreeDays.plot <- makePlot(dat=o.altModels, mod="cumulDegreeDays")
    aveTemp.plot <- makePlot(dat=o.altModels, mod="aveTemp")
    aveThermalLimit.univariate.plot <- makePlot(dat=o.altModels, mod="aveThermalLimit.univariate")
    aveThermalLimit.bivariate.plot <- makePlot(dat=o.altModels, mod="aveThermalLimit.bivariate")

    models.ecol.plot <- plot_grid(geography.plot,
              cumulDegreeDays.plot,
              aveTemp.plot,
              aveThermalLimit.univariate.plot,
              aveThermalLimit.bivariate.plot,
              ncol=1)
    ggsave(models.ecol.plot, file="~/models_ecol_plot.pdf", h=10.5, w=8)

  ### experimental attributes plot
    substrate.plot <- makePlot(dat=o.altModels, mod="substrate")
    simulans.plot <- makePlot(dat=o.altModels, mod="simulans")
    method.plot <- makePlot(dat=o.altModels, mod="method")

    models.exp.plot <- plot_grid(substrate.plot,
              simulans.plot,
              method.plot,
              ncol=1)
    ggsave(models.exp.plot, file="~/models_exp_plot.pdf", h=10.5, w=8)


scp bergland@bergland-lab.bio.virginia.edu:~/models_exp_plot.pdf .
scp bergland@bergland-lab.bio.virginia.edu:~/models_ecol_plot.pdf .
scp bergland@bergland-lab.bio.virginia.edu:~/tlm_mods.pdf .
scp bergland@bergland-lab.bio.virginia.edu:~/temperature.pdf .
scp bergland@bergland-lab.bio.virginia.edu:~/tlm_mods.14.pdf .
