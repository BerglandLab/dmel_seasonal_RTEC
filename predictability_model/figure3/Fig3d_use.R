

### Fig 3D
  load("Fig3_data/paramContourData.fix.Rdata")

  thresh.out.obs <- thresh.out.obs[daysPrior==daysPrior.use]

  setnames(thresh.out.obs, "pr", "p")

	fig3d.raw <- ggplot() +
						geom_tile(data=thresh.out.obs, aes(x=I(sp.th/10), y=I(fall.th/10), fill=r2)) +
						geom_point(data=thresh.out.obs[p>.95], aes(x=I(sp.th/10), y=I(fall.th/10)), size=.05) +
						scale_fill_distiller(palette="Spectral",
											breaks=seq(0, 1, by=0.15),
											minor_breaks=seq(0,1, by=0.1),
											name=expression(r^2),
											limits=c(.0, .85)) +
						ylab("Fall lower threshold, °C") +
						xlab("Spring upper threshold, °C") +
						theme(legend.direction="vertical",
                legend.position=c(0.05,.95),
								legend.justification=c(0,1),
								legend.background=element_rect(fill="white", colour="black", linetype="solid", size=.25),
								legend.key.size = unit(0.35, "cm"),
								legend.text=element_text(size=8),
								legend.title=element_text(size=10),
								axis.text=element_text(size=8),
								axis.title=element_text(size=10))


t1 <- loess(r2~I(sp.th/10) + I(fall.th/10), thresh.out.obs, span=.03)
thresh.out.obs[,loess:=t1$fitted]
nBins <- 10


library(metR)
pal <- NULL
if(is.null(pal)) pal <- "BrBG"

fig3d.loess <- ggplot() +
geom_contour_fill(data=thresh.out.obs, aes(x=I(sp.th/10), y=I(fall.th/10), z=loess), bins=nBins, alpha=.75, na.fill=TRUE) +
geom_tile(data=thresh.out.obs[p>.95], aes(x=I(sp.th/10), y=I(fall.th/10), fill=loess), size=.5, color="black") +
geom_contour(data=thresh.out.obs, aes(x=I(sp.th/10), y=I(fall.th/10), z=loess), bins=nBins,  color="darkgrey", size=1, na.fill=TRUE) +
scale_fill_distiller(palette=pal)+
ylab("Fall lower threshold, °C") +
xlab("Spring upper threshold, °C") +
theme_cowplot() +
theme(legend.direction="vertical", legend.position=c(0.05,.95),
   legend.justification=c(0,1),
   legend.background=element_rect(fill="white", colour="black", linetype="solid", size=.25),
   legend.key.size = unit(0.5, "cm"),
   legend.text=element_text(size=8),
   legend.title=element_text(size=10),
   axis.text=element_text(size=8),
   axis.title=element_text(size=10),
   plot.background=element_rect(fill="white"))



fig3d.loess



#rawData.topo <- ggplot() +
#geom_contour_fill(data=thresh.out.obs, aes(x=I(sp.th/10), y=I(fall.th/10), z=r2), bins=nBins, na.fill=T, alpha=.75) +
#geom_point(data=thresh.out.obs[p>.95], aes(x=I(sp.th/10), y=I(fall.th/10)), size=1, color="black") +
#geom_contour(data=thresh.out.obs, aes(x=I(sp.th/10), y=I(fall.th/10), z=r2), bins=nBins, na.fill=T, color="black") +
#scale_fill_distiller(palette=pal)  + ggtitle(paste("raw data:", pal))
#
#
#
#fig3d <- eyeGlasses.loess
#fig3DE <- plot_grid(fig3d, fig3e, nrow=2, labels=c("D", "E"))
#
#ggsave(plot_grid(fig3AB, fig3c_2, fig3DE, ncol=3, labels=c("", "C", "")), file="~/Fig3_B.pdf", h=8.5, w=12)
#
