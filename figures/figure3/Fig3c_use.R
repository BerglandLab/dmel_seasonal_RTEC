  ### scp aob2x@rivanna.hpc.virginia.edu:/home/aob2x/mods.perm.fix.Rdata ~/.


  ###

  ### new with Europe swap & lots of perms
    load("Fig3_data/mods.perm.fix.Rdata")

  	str(o.tlm.uni.ag.perm)
  	str(o.tlm.multi.ag.perm)

  	m <- rbind(o.tlm.uni.ag.perm[,c("bootSet", "sub.boot.num", "set", "x", "r2", "modelSet", "daysPrior"), with=F],
  						 o.tlm.multi.ag.perm[,c("bootSet", "sub.boot.num", "set", "x", "r2", "modelSet", "daysPrior"), with=F] )

  	m <- rbind(m, o.altModels.perm[,c("bootSet", "sub.boot.num", "set", "x", "r2", "modelSet", "daysPrior"), with=F], fill=T)
    m.orig <-m
    #m[is.na(daysPrior), daysPrior:=0]
    m <- m[is.na(daysPrior) | daysPrior==daysPrior.use]

  ### renames

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
  	m[modelSet=="tlm_uni" & x=="m4", x:="F-Max"]


  	m[,modelSet:=factor(modelSet, levels=c("tlm_bi", "tlm_uni",
  																				 "aveTemp",
  																				 "AveMinMax_bi", "AveMinMax_uni",
  																				 "cdd",
  																				 "geography",
  																				 "method",
  																				 "simulans",
  																				 "substrate"))]

    setnames(m, "daysPrior", "days")

    m[,perm:=bootSet*sub.boot.num]


    m.ag <- m[,list(nTests=length(r2[set!="Obs"]),
                    r2.max=round(r2[set=="Obs"], 3),
                    pr=mean(r2[set=="Obs"] >= r2[set!="Obs"])),
              list(x, modelSet, days)]
    m.ag[,prA:= 1 - p.adjust(1-pr)]
    m.ag[order(pr)]

    ### compre new & old

    #m.ag.perm <- m.ag
    #m.ag.perm[,version:="perm"]

      #m.ag.orig <- m.ag
      #m.ag.orig[,version:="orig"]
#
      #m.ag <- rbind(m.ag.perm, m.ag.orig)
      m.ag[x=="S-Max;F-Min"][modelSet=="tlm_bi"]

    ###### number of tests is 42. Bonferroni 0.05 = 0.05/42


    setkey(m, x, modelSet, days)
    setkey(m.ag, x, modelSet, days)

    m <- merge(m, m.ag)


    m[grepl("PercentSimulans", x), x:=gsub("PercentSimulans", "sim", x)]
    m[,modelSet.orig:=modelSet]

    m[modelSet=="tlm_bi", modelSet:="Thermal Limit Model Bi"]
    m[modelSet=="tlm_uni", modelSet:="Thermal Limit Model Uni"]

    m[modelSet=="aveTemp", modelSet:="Ave Temp"]
    m[modelSet=="AveMinMax_bi", modelSet:="Ave MinMax Bi"]
    m[modelSet=="AveMinMax_uni", modelSet:="Ave MinMax Uni"]
    m[modelSet=="cdd", modelSet:="Cumul Deg Day"]
    m[modelSet=="geography", modelSet:="Geo"]
    m[modelSet=="method", modelSet:="Collection Method"]
    m[modelSet=="simulans", modelSet:="Percent Sim"]
    m[modelSet=="substrate", modelSet:="Collection Substr"]

    m[,modelSet:=factor(modelSet,
      levels=c("Thermal Limit Model Uni",
      "Thermal Limit Model Bi",
      "Ave Temp",
      "Ave MinMax Uni",
      "Ave MinMax Bi",
      "Cumul Deg Day",
      "Geo",
      "Collection Method",
      "Percent Sim",
      "Collection Substr"))]



    m <- foreach(ms.i=levels(m$modelSet), .combine="rbind")%do%{
      print(ms.i)
      tmp <- m[modelSet==ms.i]
      tmp[,x2:=as.numeric(as.factor(as.character(x)))]
      tmp[,x2:=factor(x2, rev(unique(x2)))]
      tmp
    }
    m[,x:=gsub(".ave", "", x)]
    m[,x:=gsub("spring", "sp", x)]

    m[,x:=gsub("sp", "S", x)]
    m[,x:=gsub("fall", "F", x)]

    #m[,x:=gsub("tmin.", "", x)]
  #  m[,x:=gsub("tmax.", "", x)]

    m[,x:=gsub("max.cdd.", "", x)]
    m[,x:=gsub("max.cdd.", "", x)]

    m[,x:=gsub("longitude", "Long", x)]
    m[,x:=gsub("latitude", "Lat", x)]

    m[,x:=gsub("cm_", "", x)]
    m[,x:=gsub("cm_", "", x)]

    m[,x:=gsub("sim_", "", x)]
    m[,x:=gsub("sim_", "", x)]

    m[,x:=gsub("cs_", "", x)]
    m[,x:=gsub("cs_", "", x)]


    m[,x:=gsub("tmin.S", "S-Min", x)]
    m[,x:=gsub("tmin.F", "F-Min", x)]
    m[,x:=gsub("tmax.S", "S-Max", x)]
    m[,x:=gsub("tmax.F", "F-Max", x)]

    m[,x:=gsub("_", ";", x)]

    #m[,x.orig:=x]

    m[,x:=factor(x,
     levels=rev(c("S", "F", "S;F",
              "S-Min", "F-Min", "S-Max", "F-Max",
              "S-Min;F-Min", "S-Min;F-Max", "S-Max;F-Min", "S-Max;F-Max",
              "Lat", "Long", "Lat;Long")))]


    m[,days:=factor(days, c(14, 21))]

  #### Version 1
  #	fig3c1 <- ggplot() +
  #	geom_violin(data=m[set=="Perm"][grepl("Thermal", modelSet)], aes(y=r2, x=x2, fill=as.factor(days)), size=.5) +
  #  geom_point(data=m[set=="Obs"][grepl("Thermal", modelSet)][prA>.9], aes(y=r2, x=x2, fill=as.factor(days), color=as.factor(prA>.9)), size=5, shape=23, position=position_dodge(width=1.0)) +
  #	geom_point(data=m[set=="Obs"][grepl("Thermal", modelSet)], aes(y=r2, x=x2, fill=as.factor(days), color=as.factor(prA>.9)), size=4, shape=23, position=position_dodge(width=1.0)) +
  #	scale_fill_manual(values=c("white", "lightgrey")) +
  #  scale_color_manual(values=c("black", "red")) +
  #	theme(legend.position = "none",
  #        strip.text.y = element_text(angle=0, size=8)) +
  #	ylim(0,1) + xlab("") +
  #	facet_grid(modelSet~., scales="free_y", space="free_y", labeller = labeller(modelSet = label_wrap_gen(4))) + coord_flip()

  #  fig3c2 <- ggplot() +
  #  geom_violin(data=m[set=="Perm"][!grepl("Thermal", modelSet)], aes(y=r2, x=x2, fill=as.factor(days)), size=.5) +
  #  geom_point(data=m[set=="Obs"][!grepl("Thermal", modelSet)], aes(y=r2, x=x2, fill="white", color=as.factor(prA>.9)), size=4, shape=23, position=position_dodge(width=0.75)) +
  #  scale_fill_manual(values=c("white", "lightgrey")) +
  #  scale_color_manual(values=c("black", "red")) +
  #  theme(legend.position = "none", strip.text.y = element_text(angle=0, size=8))  +
  #  ylim(0,1) + xlab("") +
  #  facet_grid(modelSet~., scales="free_y", space="free_y", labeller = labeller(modelSet = label_wrap_gen(4))) + coord_flip()


  #  fig3c_1 <- plot_grid(fig3c1, fig3c2, nrow=2, rel_heights=c(1,1.5), align="v", axis = "b")
  #  fig3c_1
  #  #element_text(angle = 45, hjust = 1)

    ### Version 2

      pa.thresh <- .2
      fig3cA <- ggplot() +
      geom_violin(data=m[set=="Perm"][!grepl("Thermal", modelSet)], aes(y=r2, x=x), size=.5) +
      geom_point(data=m[set=="Obs"][!grepl("Thermal", modelSet)], aes(y=r2, x=x, fill="white", color=as.factor(prA>pa.thresh), fill=as.factor(prA>pa.thresh)), size=4, shape=23, position=position_dodge(width=0.75)) +
      scale_fill_manual(values=c("white", "lightgrey")) +
      scale_color_manual(values=c("black", "red")) +
      theme_cowplot() +
      theme(legend.position = "none",
            strip.text.y = element_text(angle=0, size=8),
            axis.text = element_text(angle=0, size=8),
            #axis.text.x = element_blank(),
            #axis.title.x=element_blank(),
            #axis.ticks.x=element_blank(),
            #axis.line.x=element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm"))  +
      ylim(0,1) + xlab("") +
      facet_grid(modelSet~., scales="free_y", space="free_y", labeller = labeller(modelSet = label_wrap_gen(4))) +
      coord_flip()

      fig3cB <- ggplot() +
      geom_violin(data=m[set=="Perm"][grepl("Thermal", modelSet)], aes(y=r2, x=x), size=.5) +
      geom_point(data=m[set=="Obs"][grepl("Thermal", modelSet)][prA>pa.thresh], aes(y=r2, x=x, color=as.factor(prA>pa.thresh), fill=as.factor(prA>pa.thresh)), size=5, shape=23, position=position_dodge(width=.9)) +
      geom_point(data=m[set=="Obs"][grepl("Thermal", modelSet)], aes(y=r2, x=x, color=as.factor(prA>pa.thresh), fill=as.factor(prA>pa.thresh)), size=4, shape=23, position=position_dodge(width=.9)) +
      geom_text(data=m[set=="Obs"][grepl("Thermal", modelSet)][prA>pa.thresh], aes(y=r2, x=x, label=paste(round(pr, 3), round(prA, 3), sep=";")), hjust="right", nudge_y = 0.05) +
      scale_fill_manual(values=c("white", "white")) +
      scale_color_manual(values=c("black", "red")) +
      theme_cowplot() +
      theme(legend.direction="vertical",
            legend.position = c(.95, 2),
            legend.justification = c("right", "top"),
            legend.text = element_text(size=8),
            legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
            strip.text.y = element_text(angle=0, size=8),
            axis.text = element_text(angle=0, size=8),
            plot.margin = unit(c(0,0,0,0), "cm"),
            axis.title=element_text(size=10)) +
      ylim(0,1) + xlab("") +
      facet_grid(modelSet~., scales="free_y", space="free_y", labeller = labeller(modelSet = label_wrap_gen(8)))   +
      coord_flip()


#fig3c_2 <-
      #fig3c <- plot_grid(fig3cA,
      #          fig3cB,
      #          nrow=2, rel_heights=c(1,.6), align="hv", axis="tblr")
#
      #          fig3c

    #fig3c <- fig3cA / fig3cB +
    #        plot_layout(heights=c(1,.7))
#
    fig3c <- fig3cA | fig3cB
    fig3c
