##### are inversion frequencies estimated from pooled data accurate?

### libraries
	library(data.table)
	library(tidyr)
	library(ggplot2)

### load 'inversion-SNPs'
	fixedInv <- read.delim("/mnt/icy_1/inbredLines/metaData/inv_fixed_coord.txt", header=F, sep=" ", as.is=T)
	fixedInv <- as.data.frame(t(fixedInv))
	fixedInv$chr <- substr(fixedInv[,1], 4, 5)
	fixedInv$pos <- as.numeric(as.character(fixedInv[,2]))

	fixedInv <- as.data.table(fixedInv)

### load NESCent data
	load(file="~/pigmentation/inputData/dat.Rdata")

### extract
	setkey(fixedInv, chr, pos)
	setkey(dat, chr, pos)
	freqs <- merge(fixedInv, dat)
	freqs.ag <- freqs[,list(f=mean(af, na.rm=T),
							rd=median(dp, na.rm=T),
							nSites=length(af)), list(inversion=V1, pop=pop)]
	freqs.ag[pop=="melNC_2003"]

### tack in info
	popInfo <- fread("/mnt/icy_3/nescent/data/all_popinfo.csv")
	setnames(popInfo, "name", "pop")
	setkey(popInfo, pop)

	popInfo[,locale:=do.call("rbind", strsplit(pop_name, "_"))[,1]]

	setkey(freqs.ag, pop)
	setkey(popInfo, pop)

	freqs.ag <- merge(freqs.ag, popInfo)
	save(freqs.ag, file="~/nescent_inversions.Rdata")

### plot
	library(data.table)
	library(ggplot2)
	library(cowplot)

	load("/Users/alanbergland/Documents/work/Projects/2005-2019/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/InversionFreqs/nescent_inversions.Rdata")

	freqs.ag[,season:=factor(freqs.ag$season, c("spring", "fall", "frost"))]

	state.ag <- na.omit(freqs.ag[,list(lat=mean(as.numeric(latitude), na.rm=T)),
								 list(state)])
	state.ag <- state.ag[order(lat)]


	freqs.ag[,state.f:=factor(state, state.ag$state)]

	ggplot(data=freqs.ag[order(season)][!is.na(season)][season!="frost"],
			aes(x=season, y=f, color=state.f, group=pop_name)) +
	geom_point() +
	geom_line() +
	ylim(0,1) +
	facet_wrap(~inversion, scales="fixed")

	freqs.ag.wide <- dcast(freqs.ag[,c("locale", "year", "inversion", "f", "season"),with=F], year + locale + inversion ~ season, value.var="f")
	freqs.ag.wide.ag <- freqs.ag.wide[,list(pr=mean(fall>spring, na.rm=T)), list(inversion)]

	ggplot(data=freqs.ag.wide, aes(x=spring, y=fall)) + geom_point() + facet_wrap(~inversion) + geom_abline(intercept=0, slope=1) +
	geom_text(data=freqs.ag.wide.ag, aes(x=.4, y=.1, label=pr)) + geom_cowplot()
