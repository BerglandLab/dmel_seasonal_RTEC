
swaps <- list(c("rd_12", "CA_es_12"),
     c("BHM_14", "MI_bh_14"),
     c("MA_12", "MA_la_12"),
     c("TKA_14", "KA_to_14"),
     c("CWI_14", "WI_cp_14"),
     c("AGA_14", "GA_at_14"),
     c("NY_12", "NY_it_12"),
     c("SON_15", "ON_su_15"),
     c("rd_13", "CA_es_13"),
     c("WI_12", "WI_cp_12"),
     c("LMA_12", "MA_la_12"),
     c("SCPA_14", "PA_st_14"),
     c("VI_12", "ES_ba_12"),
     c("BA_12", "AT_gr_12"),
     c("co_13", "CA_tu_13"),
     c("CUA_15", "VA_ch_15"),
     c("CUA_14", "VA_ch_14"))



swaps <- do.call("rbind", swaps)

swaps <- swaps[swaps[,2]!="ES_ba_12",]


### Fig 3b


  ### old figure

#sp.plot <- ggplot(data=m[season.y=="spring"], aes(x=pop.f, y=tmax, fill=beta)) + geom_boxplot() +
#            scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
#                        space = "Lab", na.value = "grey50", guide = "colourbar",
#                        limits=range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta),
#                        name=NULL) +
#            geom_hline(yintercept=320/10) +
#            theme_cowplot() +
#            theme(axis.text.x = element_blank(),
#                  legend.position="none",
#                  plot.margin = unit(c(.1,.1,.1,1), "cm")) +
#            ylab(expression("Spring Max (°C)\n14 days prior")) + xlab("")

#fall.plot <- ggplot(data=m[season.y=="fall"], aes(x=pop.f, y=tmin, fill=beta)) + geom_boxplot() +
#            scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
#                        space = "Lab", na.value = "grey50", guide = "colourbar",
#                        limits=range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta),
#                        name=NULL) +
#            geom_hline(yintercept=50/10) +
#            theme_cowplot() +
#            theme(axis.text.x = element_text(angle = 45, hjust = 1),
#                  legend.position="none",
#                  plot.margin = unit(c(.1,.1,.1,1), "cm"),
#                  plot.background = element_blank()) +
#            ylab(expression("Fall Min (°C)\n14 days prior")) + xlab("") +
#            scale_x_discrete(labels= swaps[,2])



#fig3b <- plot_grid(sp.plot, fall.plot, nrow=2, rel_heights=c(1, 1.35))
##fig3b




load(file="Fig3_data/m.fix.Rdata")
pops <- unique(m[season.y=="spring"]$pop)

m <- m[daysPrior<=daysPrior.use]

### new
  m[season.y=="spring", y:=tmax]
  m[season.y=="fall", y:=tmin]

  m[season.y=="spring", fl:="Spring Tmax (S-Max)"]
  m[season.y=="fall", fl:="Fall Tmin (F-Min)"]
  m[,fl:=factor(fl, levels=c("Spring Tmax (S-Max)", "Fall Tmin (F-Min)"))]

  m <- m[pop!="VI_12"]

  #m[season.y=="spring", fl:="Spring Tmax"]
  #m[season.y=="fall", fl:="Fall Tmin"]

  fig3b <- ggplot(data=m, aes(x=pop.f, y=y/10, fill=beta)) +
  geom_boxplot() +
  scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
              space = "Lab", na.value = "grey50", guide = "colourbar",
              limits=range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta),
              name=NULL) +
  facet_wrap(fl~., scales="free_x") + coord_flip() +
  geom_hline(data=m[fl=="Fall Tmin (F-Min)"], aes(yintercept=c(50/10))) +
  geom_hline(data=m[fl=="Spring Tmax (S-Max)"], aes(yintercept=c(320/10))) +
  theme(legend.position="none",
        plot.background = element_blank()) +
  ylab("°C") + xlab("") +
  scale_x_discrete(labels= swaps[,2]) +
  ggtitle(max(m$daysPrior))

  fig3b
    #expression("Fall Min (°C)\n14 days prior")
