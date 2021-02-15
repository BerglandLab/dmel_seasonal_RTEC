
  load(file="Fig3_data/m.fix.Rdata")
  pops <- unique(m[season.y=="spring"]$pop)

  #load(file="/mnt/pricey_1/dropPop/19_drop_1_concordance_data.Rdata") ### o19
  load("Fig3_data/19_drop_1_concordance_data.Rdata")
    o19[,type := "Core 20 LOOCV"]


  #load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.Rdata")  ### o20
  load("20_add_1_concordance_data.Rdata")

    o20[pop%in%c("PA_12", "WI_13", "PA_14", "PA_9"), type := "New set CV"]
    o20[is.na(type) & pop!="PA_2015" , type := "frost"]




    concordDat <- rbind(o19, o20[p<.05][n>=50], fill=T)


    fig3a <- ggplot(data=concordDat[type=="Core 20 LOOCV"][pop%in%pops],
                 aes(y=frac, x=log.sp.q.th, group=pop, color=beta)) +
            geom_line() +
            scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
                        space = "Lab", na.value = "grey50", guide = "colourbar",
                        limits=range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta),
                        name=NULL) +
            geom_hline(yintercept=.5, linetype="dashed") +
            scale_x_continuous(breaks=c(c(-3, -2, -1, 0),
                          log10(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])),
                      labels=c(.0001, .001, .01, 1,
                          rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])))) +
            ylim(.27, .75) + ylab("Fraction concordant") + xlab("Joint significance quantile") +
            theme_cowplot() +
            theme(legend.direction="vertical",
            legend.position = c(.95, .05),
              legend.justification = c("right", "bottom"),
              legend.key.size = unit(0.35, "cm"),
              legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
              legend.text=element_text(size=8),
              legend.title=element_text(size=10),
              legend.title.align=.5,
              axis.text=element_text(size=8),
              axis.title=element_text(size=10),
              plot.background=element_blank())
