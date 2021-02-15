
### fig 3E
load("Fig3/popStats.ag.Rdata")
fig3e <- ggplot() +
    geom_abline(slope=1, intercept=0) +
    geom_segment(data=popStats.ag[set=="Core20"],
          aes(x=pred.mu-pred.sd, xend=pred.mu+pred.sd, y=beta.obs, yend=beta.obs), alpha=.75) +
  #  geom_segment(data=popStats.ag[set=="NewSet"],
  #        aes(x=pred.mu-pred.sd, xend=pred.mu+pred.sd, y=beta.obs, yend=beta.obs), alpha=.75) +
    geom_point(data=popStats.ag[set=="Core20"], aes(x=pred.mu, y=beta.obs, fill=beta.obs),
          size=4.5, color="black", shape=21, alpha=.95) +
  # geom_point(data=popStats.ag[set=="NewSet"], aes(x=pred.mu, y=beta.obs, fill=beta.obs),
  #       size=6, color="black", shape=23) +


    scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
                          space = "Lab", na.value = "grey50", guide = "colourbar",
                          name=NULL) +
    ylim(-.25, .25) + xlim(-.25, .25) +
    theme_cowplot() +
    theme(legend.position="none",
        axis.text=element_text(size=8),
        axis.title=element_text(size=10)) +
    xlab("Predicted") + ylab("Observed") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_vline(xintercept=0, lty="dashed")
