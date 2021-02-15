##############################################################
### Analyses associated with seasonal predictability plots ###
############## Alan O. Bergland. April 2018 ##################
##############################################################




##############################################################
### Calculate per population environmental data statistics ###
##############################################################
	### libraries
		library(data.table)

	### load in data.
		### script that generates this file is here:
		load(file="/mnt/spicy_1/pigmentation/inputData/dat.Rdata")
		setkey(dat, chr, pos)

	### load in popInfo
		popInfo <- fread("/mnt/icy_3/nescent/data/all_popinfo.csv")

	### load SNPs to use
		snps <- fread("/mnt/pricey_1/dropPop/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
		setnames(snps, names(snps), c("chr", "pos"))
		setkey(snps, chr, pos)

	### summary stats
		ss <- dat[J(snps),list(median.rd=median(dp, na.rm=T)), list(chr, pop)]
		setkey(ss, pop)

		ss <- merge(ss, popInfo)

	### save
		save(ss, file="/mnt/pricey_1/dropPop/popSS.Rdata")

	### load in translation table
		load(file="/mnt/pricey_1/dropPop/popSS.Rdata")
		setnames(ss, "pop", "popOrig")
		setkey(ss, pop_name)

		translate <- fread("/mnt/pricey_1/dropPop/popTranslate.delim", header=T)
		setkey(translate, pop_name)

		ss <- merge(translate, ss)

	### save
		save(ss, file="/mnt/pricey_1/dropPop/popSS.translate.Rdata")


########################################
### load raw data & save merged data ###
########################################
	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(20)
		library(feather)
		library(ggplot2)
		library(RColorBrewer)

	### load SNPs to use
		snps <- fread("/mnt/pricey_1/dropPop/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
		setnames(snps, names(snps), c("chr", "pos"))
		setkey(snps, chr, pos)

		snps[,filter:="pass"]

	### single population
		singlePop.fl <- system("ls /mnt/pricey_1/dropPop/fisher_*", inter=T)

		singlePop <- foreach(i = singlePop.fl)%do%{
			print(i)

			dat <- fread(i)

			dat[,pop:=gsub("fisher_exactJ_", "", strsplit(rev(strsplit(i, "/")[[1]])[1], "\\.")[[1]][1])]

			setnames(dat, c("chrom", "coef", "minp2"), c("chr", "sp.coef", "sp.p"))

			setkey(dat, chr, pos)

			dat <- merge(dat, snps, all.y=T, all.x=T)

			dat[, sp.q.all := rank(sp.p)/(length(sp.p)+1)]
			dat[filter=="pass", sp.q.pass := rank(sp.p)/(length(sp.p)+1)]

			dat

		}
		singlePop <- rbindlist(singlePop)

	### drop population
		dropPop.fl <- system("ls /mnt/pricey_1/dropPop/mel_all*", inter=T)

		dropPop <- foreach(i = dropPop.fl)%do%{
			print(i)

			dat <- fread(i)

			dat[,pop:=gsub("mel_all_popyear_all20_drop", "", strsplit(rev(strsplit(i, "/")[[1]])[1], "\\.")[[1]][1])]

			setnames(dat, c("chrom", "seas.coef", "seas.p"), c("chr", "dp.coef", "dp.p"))

			#dat[, dp.q := rank(dp.p)/(length(dp.p)+1)]

			setkey(dat, chr, pos)


			dat <- merge(dat, snps, all.y=T, all.x=T)[,c("chr", "pos", "dp.coef", "dp.p", "pop", "filter"), with=F]

			dat[, dp.q.all := rank(dp.p)/(length(dp.p)+1)]
			dat[filter=="pass", dp.q.pass := rank(dp.p)/(length(dp.p)+1)]

			dat[,-"filter",with=F]

		}
		dropPop <- rbindlist(dropPop)

	### clinal analysis
		cline <- fread("/mnt/pricey_1/dropPop/mel_clinal_uniquepops.glm")
		setnames(cline, "chrom", "chr")
		setkey(cline, chr, pos)

		cline <- merge(cline, snps, all.y=T, all.x=T)
		cline[, clinal.q.all := rank(clinal.p)/(length(clinal.p)+1)]
		cline[filter=="pass", clinal.q.pass := rank(clinal.p)/(length(clinal.p)+1)]


	### merge
		setkey(singlePop, chr, pos, pop)
		setkey(dropPop, chr, pos, pop)

		m <- merge(dropPop, singlePop)

		setkey(m, chr, pos)
		m <- merge(m, cline[,-"filter", with=F])
		m[is.na(filter), filter:="fail"]

	### clean up
		rm(singlePop, dropPop, dat, snps)

	### save
		write_feather(m, path="/mnt/pricey_1/dropPop/dropPop.feather")

#############################
### analysis: concordance ###
#############################

	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(10)
		library(feather)
		library(ggplot2)
		library(RColorBrewer)

	### load data
		m <- as.data.table(read_feather(path="/mnt/pricey_1/dropPop/dropPop.feather"))

	### a few transformations
		m[,log.sp.q.all := log10(sp.q.all)]
		m[,log.dp.q.all := log10(dp.q.all)]
		m[,log.clinal.q.all := log10(clinal.q.all)]

		m[,log.sp.q.pass := log10(sp.q.pass)]
		m[,log.dp.q.pass := log10(dp.q.pass)]
		m[,log.clinal.q.pass := log10(clinal.q.pass)]

	### only work with SNPs that are in PASS

	### summarize along sp & dp q, conditional on cline.
		### for all sites
			### generate binned quantile values
			#	m[,log.sp.q.all.bin := round(log.sp.q.all, 2)]
			#	m[,log.dp.q.all.bin := round(log.dp.q.all, 2)]

			#	wins <- unique(c(m$log.sp.q.all.bin, m$log.dp.q.all.bin))
			#	wins <- wins[order(wins)]

			#	setkey(m, log.sp.q.all.bin)

			### calculate scores
			#	registerDoMC(5)
			#	o.all <- foreach(i = 1:(length(wins)-50), .combine="rbind", .export="m")%dopar%{
					print(i)

					m.temp <- m[J(wins[1:i]), nomatch=0]
					setkey(m.temp, log.dp.q.all.bin)
					m.temp <- m.temp[J(wins[1:i]), nomatch=0]

					m.ag <- m.temp[,
									list(TT=sum(sign(dp.coef)==sign(sp.coef)),
										n=length(dp.coef),
										log.sp.q.th=wins[i],
										log.dp.q.th=wins[i],
										class="all"),
								  list(pop)]

					m.ag
				}

		### for filtered sites
			### generate binned quantile values
				m[,log.sp.q.pass.bin := round(log.sp.q.pass, 2)]
				m[,log.dp.q.pass.bin := round(log.dp.q.pass, 2)]

				wins <- unique(c(m$log.sp.q.pass.bin, m$log.dp.q.pass.bin))
				wins <- wins[order(wins)]

				setkey(m, log.sp.q.pass.bin)

			### calculate scores
				registerDoMC(5)
				o.pass <- foreach(i = 1:(length(wins)), .combine="rbind", .export="m")%dopar%{
					print(i)

					m.temp <- m[J(wins[1:i]), nomatch=0]
					setkey(m.temp, log.dp.q.pass.bin)
					m.temp <- m.temp[J(wins[1:i]), nomatch=0]

					m.ag <- m.temp[,
									list(TT=sum(sign(dp.coef)==sign(sp.coef)),
										n=length(dp.coef),
										log.sp.q.th=wins[i],
										log.dp.q.th=wins[i],
										class="pass"),
								  list(pop)]
					m.ag

				}


		### merge them
			# o <- rbind(o.all, o.pass)
			o <- o.pass

		### summarize stats
			o[,p := dbinom(TT, n, .5)]
			o[,frac := TT/n ]
			o[,q:= p.adjust(p, "fdr")]


	### save object
		save(o, file="/mnt/pricey_1/dropPop/concordance.Rdata")

	### plot
		library(ggplot2)
		library(data.table)
		library(foreach)

		load(file="/mnt/pricey_1/dropPop/concordance.Rdata")


		### model fit to get slopes
			small.o <- o[p<.05]

			pred.beta <- foreach(i=unique(small.o$pop), .combine="rbind")%do%{
				foreach(filter.set=c("pass"), .combine="rbind")%do%{
					temp <- small.o[pop==i][class==filter.set]

				### playing
					#i <- "CUA_14"
					#temp <- small.o[pop==i][log.sp.q.th < 0]

					#t1 <- glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=small.o[pop==i]$n)
					#temp[,pred:=predict(t1)]
					#plot(frac~I(-1*log.sp.q.th), data=temp, type="l")
					#points(plogis(pred)~I(-1*log.sp.q.th), data=temp, type="l", col="red")



					glm.out <- summary(glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=temp$n))

					data.table(pop=i,
								beta=glm.out$coef[2,1],
								p=glm.out$coef[2,4],
								class=filter.set)
				}
			}

		### tack in slopes to o[] to color by slope
			setkey(o, pop)
			setkey(pred.beta, pop)

			o <- merge(o, pred.beta[,c("pop", "beta"), with=F])

			o19 <- o[p<.05][class=="pass"][n>35]

			save(o19, file="/mnt/pricey_1/dropPop/19_drop_1_concordance_data.Rdata")


######################################################################################
### analysis: concordance as a function of joint-q, per chromosome and genome-wide ###
######################################################################################

	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(10)
		library(feather)
		library(ggplot2)
		library(RColorBrewer)

	### load data
		m <- as.data.table(read_feather(path="/mnt/pricey_1/dropPop/dropPop.feather"))

	######################################################
	#### this verison uses the genome-ranked q-values. ###
	######################################################

	### a few transformations
		m[,log.sp.q.all := log10(sp.q.all)]
		m[,log.dp.q.all := log10(dp.q.all)]
		m[,log.clinal.q.all := log10(clinal.q.all)]

		m[,log.sp.q.pass := log10(sp.q.pass)]
		m[,log.dp.q.pass := log10(dp.q.pass)]
		m[,log.clinal.q.pass := log10(clinal.q.pass)]

	### only work with SNPs that are in PASS

		### for filtered sites
			### generate binned quantile values
				m[,log.sp.q.pass.bin := round(log.sp.q.pass, 2)]
				m[,log.dp.q.pass.bin := round(log.dp.q.pass, 2)]

				wins <- unique(c(m$log.sp.q.pass.bin, m$log.dp.q.pass.bin))
				wins <- wins[order(wins)]

				setkey(m, log.sp.q.pass.bin)

			### calculate scores

				### highly parallelized smaller subsets
					registerDoMC(18)
					m.ag <- foreach(chrs=list(c("2L", "2R", "3L", "3R"), "2L", "2R", "3L", "3R"))%do%{
						print(chrs)

						foreach(i = 1:(length(wins)-50), .combine="rbind", .export="m")%dopar%{
							print(i)
							m.temp <- m[J(wins[1:i]), nomatch=0]

							m.temp <- m[J(wins[1:i]), nomatch=0]
							setkey(m.temp, log.dp.q.pass.bin)
							m.temp <- m.temp[J(wins[1:i]), nomatch=0]
							setkey(m.temp, chr)

							m.temp <- m.temp[J(chrs)]

							m.ag <- m.temp[,
											list(TT=sum(sign(dp.coef)==sign(sp.coef)),
												n=length(dp.coef),
												log.sp.q.th=wins[i],
												log.dp.q.th=wins[i],
												class="all",
												chr=paste(chrs, collapse=";")),
										  list(pop)]
							m.ag[]
						}

					}
					o.one <- rbindlist(m.ag)

				### medium parallelized medium subsets
					registerDoMC(10)
					m.ag <- foreach(chrs=list(c("2L", "2R", "3L", "3R"), "2L", "2R", "3L", "3R"))%do%{
						print(chrs)

						foreach(i = ((length(wins)-50)+1):(length(wins)-20), .combine="rbind", .export="m")%dopar%{
							print(i)
							m.temp <- m[J(wins[1:i]), nomatch=0]

							m.temp <- m[J(wins[1:i]), nomatch=0]
							setkey(m.temp, log.dp.q.pass.bin)
							m.temp <- m.temp[J(wins[1:i]), nomatch=0]
							setkey(m.temp, chr)

							m.temp <- m.temp[J(chrs)]

							m.ag <- m.temp[,
											list(TT=sum(sign(dp.coef)==sign(sp.coef)),
												n=length(dp.coef),
												log.sp.q.th=wins[i],
												log.dp.q.th=wins[i],
												class="all",
												chr=paste(chrs, collapse=";")),
										  list(pop)]
							m.ag[]
						}

					}
					o.two <- rbindlist(m.ag)

				### low paralllized for large
					registerDoMC(5)
					m.ag <- foreach(chrs=list(c("2L", "2R", "3L", "3R"), "2L", "2R", "3L", "3R"))%do%{
						print(chrs)

						foreach(i = ((length(wins)-20)+1):(length(wins)-0), .combine="rbind", .export="m")%dopar%{
							print(i)
							m.temp <- m[J(wins[1:i]), nomatch=0]

							m.temp <- m[J(wins[1:i]), nomatch=0]
							setkey(m.temp, log.dp.q.pass.bin)
							m.temp <- m.temp[J(wins[1:i]), nomatch=0]
							setkey(m.temp, chr)

							m.temp <- m.temp[J(chrs)]

							m.ag <- m.temp[,
											list(TT=sum(sign(dp.coef)==sign(sp.coef)),
												n=length(dp.coef),
												log.sp.q.th=wins[i],
												log.dp.q.th=wins[i],
												class="all",
												chr=paste(chrs, collapse=";")),
										  list(pop)]
							m.ag[]
						}

					}
					o.three <- rbindlist(m.ag)

				### merge
					o <- rbindlist(list(o.one, o.two, o.three))

				### p-values

					o[,p := dbinom(TT, n, .5)]
					o[,frac := TT/n ]
					o[,q:= p.adjust(p, "fdr")]


			### model fit to get slopes: slopes here are only calcualted for the full genome model...
				small.o <- o[p<.05][n>50]

				pred.beta <- foreach(i=unique(small.o$pop), .combine="rbind")%do%{
					foreach(filter.set=c("all"), .combine="rbind")%do%{
						temp <- small.o[pop==i][class==filter.set][chr=="2L;2R;3L;3R"]

					### playing
						#i <- "CUA_14"
						#temp <- small.o[pop==i][log.sp.q.th < 0]

						#t1 <- glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=small.o[pop==i]$n)
						#temp[,pred:=predict(t1)]
						#plot(frac~I(-1*log.sp.q.th), data=temp, type="l")
						#points(plogis(pred)~I(-1*log.sp.q.th), data=temp, type="l", col="red")



						glm.out <- summary(glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=temp$n))

						data.table(pop=i,
									beta=glm.out$coef[2,1],
									p=glm.out$coef[2,4],
									class=filter.set)
					}
				}

			### tack in slopes to o[] to color by slope
				setkey(small.o, pop)
				setkey(pred.beta, pop)

				small.o <- merge(small.o, pred.beta[,c("pop", "beta"), with=F])

				o19.chr <- small.o[p<.05][class=="all"][n>35]

			### plot
				o19.chr[,chr := factor(chr, levels=c("2L;2R;3L;3R", "2L", "2R", "3L", "3R"))]

				setnames(o19.chr, "beta.x", "beta")

				ggplot(data=o19.chr,
						   aes(y=frac, x=log.sp.q.th, group=pop, color=beta)) +
					geom_line() +
					facet_grid(~chr) +
					scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
										  space = "Lab", na.value = "grey50", guide = "colourbar",
										  limits=range(o19.chr$beta),
										  name=NULL) +
					geom_hline(yintercept=.5, linetype="dashed") +
					scale_x_continuous(breaks=c(c(-3, -2, -1, 0),
												log10(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])),
										labels=c(.0001, .001, .01, 1,
												rep("", length(expand.grid(sapply(c(1:9), function(x) x*10^c(-3, -2, -1)))[,1])))) +
					ylim(.2, .8) + ylab("Fraction concordant") + xlab("Joint significance quantile") +
					theme(legend.direction="vertical",
						legend.justification=c(1,0),
						legend.position=c(1,0),
						legend.key.size = unit(0.35, "cm"),
						legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
						legend.text=element_text(size=8),
						legend.title=element_text(size=10),
						legend.title.align=.5,
						axis.text=element_text(size=8),
						axis.title=element_text(size=10))

		### save
			save(o19.chr, file="/mnt/pricey_1/dropPop/concordance_o19_chr.Rdata")




###########################################
### analysis: concordance  & enrichment ###
###########################################

	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(10)
		library(feather)
		library(ggplot2)
		library(RColorBrewer)

	### load data
		m <- as.data.table(read_feather(path="/mnt/pricey_1/dropPop/dropPop.feather"))

	### a few transformations
		m[,log.sp.q.all := log10(sp.q.all)]
		m[,log.dp.q.all := log10(dp.q.all)]
		m[,log.clinal.q.all := log10(clinal.q.all)]

		m[,log.sp.q.pass := log10(sp.q.pass)]
		m[,log.dp.q.pass := log10(dp.q.pass)]
		m[,log.clinal.q.pass := log10(clinal.q.pass)]

	### only work with SNPs that are in PASS

	### summarize along sp & dp q, conditional on cline.
		### for all sites
			### generate binned quantile values
			#	m[,log.sp.q.all.bin := round(log.sp.q.all, 2)]
			#	m[,log.dp.q.all.bin := round(log.dp.q.all, 2)]

			#	wins <- unique(c(m$log.sp.q.all.bin, m$log.dp.q.all.bin))
			#	wins <- wins[order(wins)]

			#	setkey(m, log.sp.q.all.bin)

			### calculate scores
			#	registerDoMC(5)
			#	o.all <- foreach(i = 1:(length(wins)-50), .combine="rbind", .export="m")%dopar%{
					print(i)

					m.temp <- m[J(wins[1:i]), nomatch=0]
					setkey(m.temp, log.dp.q.all.bin)
					m.temp <- m.temp[J(wins[1:i]), nomatch=0]

					m.ag <- m.temp[,
									list(TT=sum(sign(dp.coef)==sign(sp.coef)),
										n=length(dp.coef),
										log.sp.q.th=wins[i],
										log.dp.q.th=wins[i],
										class="all"),
								  list(pop)]

					m.ag
				}

		### for filtered sites
			### generate binned quantile values
				m[,log.sp.q.pass.bin := round(log.sp.q.pass, 2)]
				m[,log.dp.q.pass.bin := round(log.dp.q.pass, 2)]

				wins <- unique(c(m$log.sp.q.pass.bin, m$log.dp.q.pass.bin))
				wins <- wins[order(wins)]


				m.n <- m[,list(nTot=length(pos)), pop]
				setkey(m.n, pop)

				m2 <- m
				setkey(m2, chr, pos)
				m2 <- m2[!duplicated(m2)]
				m2 <- m2[!is.na(log.dp.q.pass.bin)]
				setkey(m2, log.dp.q.pass.bin)

				setkey(m, log.sp.q.pass.bin)


			### calculate scores
				registerDoMC(5)
				o.enrich.pass <- foreach(i = 1:(length(wins)-50), .combine="rbind", .export="m")%dopar%{
					print(i)

					m.dp.temp <- m2[J(wins[1:i]), nomatch=0]
					m.sp.temp <- m[J(wins[1:i]), nomatch=0]

					dp.frac <- dim(m.dp.temp)[1]/dim(m2)[1]

					sp.frac <- m.sp.temp[,list(n.sp=length(pos)), pop]
					setkey(sp.frac, pop)
					sp.frac <- merge(sp.frac, m.n)
					sp.frac[,dp.frac:=dp.frac]
					sp.frac[,sp.frac:=n.sp/nTot]


					setkey(m.sp.temp, log.dp.q.pass.bin)
					m.temp <- m.sp.temp[J(wins[1:i]), nomatch=0]



					m.ag <- m.temp[,
									list(n.sp.dp=length(pos),
										log.sp.q.th=wins[i],
										log.dp.q.th=wins[i],
										class="pass"),
								  list(pop)]
					setkey(m.ag, pop)

					m.ag <- merge(m.ag, sp.frac)

					m.ag[,exp.n.sp.dp := dp.frac*sp.frac*nTot]

				}


		### merge them
			# o <- rbind(o.all, o.pass)
			#o <- o.pass

		### summarize stats

			#o.ag <- o.enrich.pass[n.sp.dp>10,list(en=log2(n.sp.dp/exp.n.sp.dp)), list(pop, log.sp.q.th)]
			o.ag <- o.enrich.pass[n.sp.dp>10, list(en=(n.sp.dp-exp.n.sp.dp)/exp.n.sp.dp,
													en2=log2(n.sp.dp/exp.n.sp.dp)),
												list(pop, log.sp.q.th)]


			setkey(o.ag, pop, log.sp.q.th)

		### load in predictability scores
			load(file="/mnt/pricey_1/dropPop/19_drop_1_concordance_data.Rdata")
			setkey(o19,  pop, log.sp.q.th)

			foo <- merge(o.ag, o19)

		### save joint object
			save(foo, file="/mnt/pricey_1/dropPop/19_drop_1_enrichment.Rdata")




#########################################################################################
### analysis: model fit and correlation with experimental and environmental variables ###
#########################################################################################

	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(10)
		library(feather)
		library(ggplot2)
		library(RColorBrewer)
		library(cowplot)

	### load concordance object
		load(file="/mnt/pricey_1/dropPop/concordance.Rdata")

	### model fit
		small.o <- o[p<.05]


		pred.beta <- foreach(i=unique(small.o$pop), .combine="rbind")%do%{
			foreach(filter.set=c("pass"), .combine="rbind")%do%{
				temp <- small.o[pop==i][class==filter.set]

			### playing
				#i <- "CUA_14"
				#temp <- small.o[pop==i][log.sp.q.th < 0]

				#t1 <- glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=small.o[pop==i]$n)
				#temp[,pred:=predict(t1)]
				#plot(frac~I(-1*log.sp.q.th), data=temp, type="l")
				#points(plogis(pred)~I(-1*log.sp.q.th), data=temp, type="l", col="red")



				glm.out <- summary(glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=temp$n))

				data.table(pop=i,
							beta=glm.out$coef[2,1],
							p=glm.out$coef[2,4],
							class=filter.set)
			}
		}

	### merge with summary stats
		load(file="/mnt/pricey_1/dropPop/popSS.translate.Rdata")

		### clean a few doobies
			ss[pop=="OUK_13" & month==10, col_date:="10/15/2013"]
			ss[pop=="OUK_13" & month==5, col_date:="5/15/2013"]

			ss[pop=="TKA_14", latitude:="39.202831"]
			ss[,latitude:=as.numeric(as.character(latitude))]

			ss[pop=="AGA_14" & season=="fall", col_date:="10/15/2014"]

		### do merge
			setkey(ss, pop)

			ss.ag <- ss[season!="frost",list(median.rd=mean(median.rd),
							  latitude=mean(as.numeric(as.character(latitude)), na.rm=T),
							  longitude=mean(longitude, na.rm=T),
							  col_date=unique(col_date[!is.na(col_date)]),
							  dayTrue=!is.na(day)[1]),
						list(pop, season)]

			setkey(ss.ag, pop, season)
			setkey(pred.beta, pop)

			pbs <- merge(ss.ag, pred.beta, all.y=T, allow.cartesian=T)

		### clean up a few doobies
			pbs[,col_date:=gsub("2014", "14", col_date)]
			pbs[,col_date:=gsub("2013", "13", col_date)]

		### format date stamp
			pbs[,col_date := as.Date(col_date, format="%m/%d/%y")]
			setnames(pbs, "col_date", "date")

		### set key
			setkey(pbs, pop, date)

		### save pbs objct
			save(pbs, file="/mnt/pricey_1/dropPop/pbs.Rdata")

	### load climate data
		load("/mnt/pricey_1/dropPop/ghcnd_records.RData")

		weather <- as.data.table(ghcnd_records)
		setnames(weather, "pop_name", "pop")

		setkey(weather, pop, date)


	### Calculate days above and below thresholds & include permutations. Takes a while to run so load if possible
		daysPrior <- 21

		### make general spring & fall datasets
			weather.pop.spring <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{

				wt <- weather[J(data.table(pop=pbs[i]$pop,
					 												 date=seq(from=pbs[i]$date-(daysPrior), to=pbs[i]$date, by=1),
					 											 	 key="pop,date"))]


				wt[,daysPrior := pbs[i]$date - date]
				wt[,season := "spring"]

				wt
			}
			weather.pop.fall <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{

				wt <- weather[J(data.table(pop=pbs[i]$pop,
					 date=seq(from=pbs[i]$date-(daysPrior), to=pbs[i]$date, by=1),
					 key="pop,date"))]


				wt[,daysPrior := pbs[i]$date - date]
				wt[,season := "fall"]

				wt
			}

			weather.pop <- rbind(weather.pop.spring, weather.pop.fall)

		### clean up pbs
			pbs.small <- pbs[class=="pass"][season=="spring"]
			setkey(pbs.small, pop)

		### save weather.pop & pbs.small
			save(weather.pop, pbs.small, file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata")



		### iterate through temperature sets & do permutation test
			nPerms <- 10000
			registerDoMC(20)

			thresh.out <- foreach(sp.th=seq(from=-50, to=400, by=5), .combine="rbind", .errorhandling="remove")%dopar%{
				foreach(fall.th=seq(from=-50, to=400, by=5), .combine="rbind", .errorhandling="remove")%do%{

					print(paste(sp.th, fall.th, sep= " / "))

					weather.pop.ag <- weather.pop[daysPrior<=daysPrior,list(tmax=mean(tmax[season=="spring"] >= sp.th, na.rm=T),
									 				   																			tmin=mean(tmin[season=="fall"] <= fall.th, na.rm=T)),
								 list(pop)]

					setkey(weather.pop.ag, pop)

					m <- merge(pbs.small, weather.pop.ag)[pop!="OUK_13"][dayTrue==T]

					### make matrix of betas for observed & permuted
						set.seed(1234)
						beta.mat <- cbind(m$beta, replicate(nPerms, sample(m$beta, replace=F)))

					### run models
						#m4 <- lm(beta.mat ~ tmax + tmin, m)
						#p4 <- predict(m4)

						t4 <- summary(lm(beta.mat ~ tmax + tmin, m))

					### extract out data
						foreach(i=1:length(t4), .combine="rbind")%do%{

							data.table(sp.th=sp.th, fall.th=fall.th,
										boot.num =i,
										r2=(t4[[i]])$r.squared,
										int=coef(t4[[i]])[1,1],
										tmax=coef(t4[[i]])[2,1],
										tmin=coef(t4[[i]])[3,1])
						}

					}
				}


			save(thresh.out, file="/mnt/pricey_1/dropPop/threshDays_withPerms_10K.Rdata")



	### start here for revisions
		### generate permutation based p-values and plot
			library(data.table)
			library(ggplot2)
			load(file="/mnt/pricey_1/dropPop/threshDays_withPerms_10K.Rdata")

			thresh.out.ag <- thresh.out[,list(maxR2=max(r2)), list(boot.num)]
			thresh.out.obs <- thresh.out[boot.num==1,
										 list(r2=r2,
										 	  p=mean(r2 >= thresh.out.ag[boot.num!=1]$maxR2),
										 	  int=int, tmax=tmax, tmin=tmin),
										 list(sp.th=sp.th,
										 	  fall.th=fall.th)]

			save(thresh.out.obs, file="/mnt/pricey_1/dropPop/paramContourData.Rdata")


		### best model plot
			library(data.table)
			load(file="/mnt/pricey_1/dropPop/threshDays_withPerms_10K.Rdata")
			load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata")

			thresh.out.ag <- thresh.out[,list(maxR2=max(r2)), list(boot.num)]
			thresh.out.obs <- thresh.out[boot.num==1,
										 list(r2=r2,
										 	  p=mean(r2 >= thresh.out.ag[boot.num!=1]$maxR2),
										 	  int=int, tmax=tmax, tmin=tmin),
										 list(sp.th=sp.th,
										 	  fall.th=fall.th)]



			bestModel <- thresh.out.obs[which.max(r2)]

			weather.pop.ag <- weather.pop[,list(tmax.max=mean(tmax[season=="spring"] >= bestModel$sp.th, na.rm=T),
									 			tmin.min=mean(tmin[season=="fall"] <= bestModel$fall.th, na.rm=T),
									 			tmax.mean=mean(tmax[season=="spring"], na.rm=T),
									 			tmin.mean=mean(tmin[season=="fall"], na.rm=T),
									 			cdd.spring=cdd[season=="spring"][daysPrior==0],
									 			cdd.fall=cdd[season=="fall"][daysPrior==0]),
								 list(pop)]

			setkey(weather.pop.ag, pop)

			m <- merge(pbs.small, weather.pop.ag)[pop!="OUK_13"][dayTrue==T]

			t4 <- lm(beta ~ tmax.max + tmin.min, m)


			pred.full <- data.table(tmax.max = rep(seq(from=0, to=.4, by=.001), each=length(seq(from=0, to=.3, by=.001))),
								   tmin.min = rep(seq(from=0, to=.3, by=.001), length(seq(from=0, to=.4, by=.001))))
			pred.full[,pred := predict(t4, newdata=pred.full)]

		### save best fit model objects
			save(m, pred.full, file="/mnt/pricey_1/dropPop/pred_pred.plotData.Rdata")



#####################################################
### 20 pred 1 analysis for additional populations ###
#####################################################

	### save and export data
		### libraries
			library(data.table)
			library(foreach)

		### calculate Z scores for new population pairs
			### load newPopPairs file
				newPopPairs <- fread("/mnt/pricey_1/dropPop/newPopPairs.csv")
				setnames(newPopPairs, names(newPopPairs)[c(1,2)], c("pop.t1", "pop.t2"))

			### load allele frequency data
				load(file="/mnt/spicy_1/pigmentation/inputData/dat.Rdata") ### script that makes this is sftp://bergland:@bergland-lab.bio.virginia.edu//mnt/spicy_1/pigmentation/scripts/Sept2017_v1/pigmentationOrganizeData_generalized.R
				dat[,pop := as.character(pop)]
				setkey(dat, pop)

			### get set of SNPs that lines up with other single population sets
				singlePopSNPs <- fread(system("ls /mnt/pricey_1/dropPop/fisher_*", inter=T)[1])[,c("chrom", "pos"), with=F]
				setnames(singlePopSNPs, "chrom", "chr")
				setkey(singlePopSNPs, chr, pos)

			### calculate Z/FET scores and write to file
				foreach(i=1:dim(newPopPairs)[1])%do%{
					print(paste(i, dim(newPopPairs)[1], sep=" / "))

					temp.t0 <- dat[J(newPopPairs[i]$pop.t1)]
					temp.t1 <- dat[J(newPopPairs[i]$pop.t2)]

					setkey(temp.t0, chr, pos)
					setkey(temp.t1, chr, pos)

					temp.m <- merge(temp.t0[,-c("ref", "alt"), with=F],
									temp.t1[,-c("ref", "alt"), with=F])


					temp.m[,d := af.x - af.y]
					temp.m[,p0 := af.x/2 + af.y/2]
					temp.m[,scaleFactor := 1/dp.x + 1/dp.y]
					temp.m[,s:=sqrt(p0*(1-p0)*scaleFactor)]
					temp.m[,Z:=d/s]

					temp.m[,sp.p := 2*apply(cbind(pnorm(Z, 0, 1), 1-pnorm(Z, 0, 1)), 1, min)]
					temp.m[,sp.coef := sign(Z)]
					temp.m[,pop:=newPopPairs[i]$pop_name]


					setkey(temp.m, chr, pos)

					temp.o <- merge(temp.m, singlePopSNPs, all.y=T)[,c("chr", "pos", "sp.coef", "sp.p", "pop"), with=F]

					write.csv(temp.o, file=paste("/mnt/pricey_1/dropPop/Z_", newPopPairs[i]$pop_name, ".coef_minp2.txt", sep=""),
								quote=F, row.names=F)

				}

	### now, replicate analysis of "19 drop 1" using 20 drop 1 with new frost & s/f populations
		### libraries
			library(data.table)
			library(foreach)
			library(doMC)
			registerDoMC(20)
			library(ggplot2)

		### load SNPs to use
			snps <- fread("/mnt/pricey_1/dropPop/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
			setnames(snps, names(snps), c("chr", "pos"))
			setkey(snps, chr, pos)

			snps[,filter:="pass"]

		### functionalize it

			analysis_20_1 <- function(filename, stat_class) {

				### load full 20 dataset
					full20 <- fread(filename)
					#full20 <- fread("/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm")
					#full20 <- fread("/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")

					setnames(full20, "chrom", "chr")
					setkey(full20, chr, pos)

					full20 <- merge(full20, snps, all.x=T)
					full20[,full.q.all := rank(seas.p)/(length(seas.p)+1)]
					full20[filter=="pass",full.q.pass := rank(seas.p)/(length(seas.p)+1)]
					setkey(full20, chr, pos)

				### single population
					#singlePop.new.fl <- system("ls /mnt/pricey_1/dropPop/Z_*", inter=T)
					singlePop.new.fl <-c("/mnt/pricey_1/dropPop/fisher_exactJ_PA_12.coef_minp2.txt",
										 "/mnt/pricey_1/dropPop/fisher_exactJ_PA_9.coef_minp2.txt",
										 "/mnt/pricey_1/dropPop/fisher_exactJ_PA_14.coef_minp2.txt",
										 "/mnt/pricey_1/dropPop/fisher_exactJ_WI_13.coef_minp2.txt",
										 "/mnt/pricey_1/dropPop/fisher_exactJ_PA_15.coef_minp2.txt")



					singlePop <- foreach(i = singlePop.new.fl)%do%{
						print(i)

						dat <- fread(i)

						setnames(dat, c("chrom", "minp", "coef"), c("chr", "sp.p", "sp.coef")) #### specific to FET analysis
						dat[,pop := gsub(".coef_minp2.txt", "", gsub("/mnt/pricey_1/dropPop/fisher_exactJ_", "", i))]


						setkey(dat, chr, pos)

						dat <- merge(dat, snps, all.y=T, all.x=T)

						dat[, sp.q.all := rank(sp.p)/(length(sp.p)+1)]
						dat[filter=="pass", sp.q.pass := rank(sp.p)/(length(sp.p)+1)]

						dat

					}
					singlePop <- rbindlist(singlePop)
					setkey(singlePop, chr, pos)

				### merge singlePop & full20
					m <- merge(full20, singlePop)

				### a few data transformations
					m[,log.sp.q.all := log10(sp.q.all)]
					m[,log.full.q.all := log10(full.q.all)]

					m[,log.sp.q.pass := log10(sp.q.pass)]
					m[,log.full.q.pass := log10(full.q.pass)]

				### summarize along sp & dp q, conditional on cline. for filterted sites
					### generate binned quantile values
						m[,log.sp.q.pass.bin := round(log.sp.q.pass, 2)]
						m[,log.full.q.pass.bin := round(log.full.q.pass, 2)]

						wins <- unique(c(m$log.sp.q.pass.bin, m$log.full.q.pass.bin))
						wins <- wins[order(wins)]

						setkey(m, log.sp.q.pass.bin)

				### calculate scores
					registerDoMC(5)
					o <- foreach(i = 1:(length(wins)), .combine="rbind", .export="m")%dopar%{
						print(i)

						m.temp <- m[J(wins[1:i]), nomatch=0]
						setkey(m.temp, log.full.q.pass.bin)
						m.temp <- m.temp[J(wins[1:i]), nomatch=0]

						m.ag <- m.temp[,
										list(TT=sum(sign(seas.coef)==sign(sp.coef)),
											n=length(seas.coef),
											log.sp.q.th=wins[i],
											log.full.q.th=wins[i],
											class="pass"),
									  list(pop)]
						m.ag

					}

				### summarize stats
					o[,p := dbinom(TT, n, .5)]
					o[,frac := TT/n ]
					o[,q:= p.adjust(p, "fdr")]

				### tack in metadata and return
					o[,class:=stat_class]
					return(o)
			}

		### run both flipped & original models

			o <- list(analysis_20_1(filename="/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm", stat_class="orig"),
					  analysis_20_1(filename="/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm", stat_class="swap"))

			o <- rbindlist(o)


			### save
				save(o, file="/mnt/pricey_1/dropPop/concordance_full20_add1.orig.swap.FET.Rdata")

				#save(o, file="/mnt/pricey_1/dropPop/concordance_full20_add1.orig.swap.Rdata")
				#save(o, file="/mnt/pricey_1/dropPop/concordance_full20_add1.Rdata")

		### make plot frame
			load(file="/mnt/pricey_1/dropPop/concordance_full20_add1.orig.swap.FET.Rdata")

			#load(file="/mnt/pricey_1/dropPop/concordance_full20_add1.orig.swap.Rdata")
			#load(file="/mnt/pricey_1/dropPop/concordance_full20_add1.orig.Rdata")
			#load(file="/mnt/pricey_1/dropPop/concordance_full20_add1.Rdata")

		### functionalize it
			o20 <- foreach(k=unique(o$class))%do%{

				small.o <- o[class==k][p<.05]

				pred.beta <- foreach(i=unique(small.o$pop), .combine="rbind")%do%{
					temp <- small.o[pop==i]

					### playing
						#i <- "CUA_14"
						#temp <- small.o[pop==i][log.sp.q.th < 0]

						#t1 <- glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=small.o[pop==i]$n)
						#temp[,pred:=predict(t1)]
						#plot(frac~I(-1*log.sp.q.th), data=temp, type="l")
						#points(plogis(pred)~I(-1*log.sp.q.th), data=temp, type="l", col="red")



					glm.out <- summary(glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=temp$n))

					data.table(pop=i,
								beta=glm.out$coef[2,1],
								p=glm.out$coef[2,4],
								class=k)
				}

				### tack in slopes to o[] to color by slope
					setkey(small.o, pop, class)
					setkey(pred.beta, pop, class)

					small.o <- merge(small.o, pred.beta[,c("pop", "beta"), with=F])

				### return
					small.o
			}
			o20 <- rbindlist(o20)

		### add a bit of annotation
			o20[pop%in%c("PA_12", "WI_13", "PA_14", "PA_9"), type := "New set CV"]
			o20[is.na(type) & pop!="PA_2015" , type := "frost"]

		### simple plot
			ggplot(data=o20[type=="New set CV"][p<.05][n>=50], aes(x=log.sp.q.th, y=frac, group=pop, color=beta)) +
			geom_line() +
			facet_wrap(~class)

		### save
			save(o20, file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.orig.swap.FET.Rdata")

#############################
##### Weather 20 choose 1 ###
#############################

	### weather model fit using parameters from previous analysis
		### libraries
			library(data.table)
			library(foreach)
			library(doMC)
			registerDoMC(20)
			library(ggplot2)
			library(epiR)

		### load predictability data
			load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.orig.swap.FET.Rdata")

#			load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.orig.swap.Rdata")
#			load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.orig.Rdata")
#			load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.Rdata")

		### prediction model fits
			pred.beta <- o20[,list(beta=mean(beta), p=mean(p)), list(pop, class)]
			setkey(pred.beta, pop)

		### tack in translation table inorder to merge with popInfo
			translate <- fread("/mnt/pricey_1/dropPop/popTranslate.delim", header=T)
			setkey(translate, pop)

			pred.beta <- merge(pred.beta, translate)
			setkey(pred.beta, pop_name)

		### tack in collection date from popInfo file; requires sample name conversion
			popInfo <- fread("/mnt/icy_3/nescent/data/all_popinfo.csv")
			setkey(popInfo, pop_name)

			pred.beta <- merge(pred.beta, popInfo[,c("pop_name", "col_date", "season"), with=F])

		### convert col_date to Datetime class
			pred.beta[, col_date := as.Date(col_date, format="%m/%d/%y")]

		### load in weather data & gdd data
			load("/mnt/pricey_1/dropPop/ghcnd_records.RData")

			weather <- as.data.table(ghcnd_records)
			setnames(weather, "pop_name", "pop")

			setkey(weather, pop, date)

		### load in previously calculated threshold data from 19:1 analysis
			#load(file="/mnt/pricey_1/dropPop/threshDays_remade.Rdata")

		### calculate weather statistics for new populations
			### make general spring & fall datasets
				daysPrior <- 21
				frost_set <- unique(pred.beta[season=="frost"]$pop)
				sf_set <- unique(pred.beta[season=="spring"]$pop)

				weather.pop.spring.sf <- foreach(i=c(1:dim(pred.beta)[1])[pred.beta$pop%in%sf_set & pred.beta$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{

					wt <- weather[J(data.table(pop=pred.beta[i]$pop,
												 date=seq(from=pred.beta[i]$col_date-(daysPrior), to=pred.beta[i]$col_date, by=1),
						 key="pop,date"))]


					wt[,daysPrior := pred.beta[i]$col_date - date]
					wt[,season := "spring"]
					wt[,set := "spring_fall"]

					wt

				}

				weather.pop.fall.sf <- foreach(i=c(1:dim(pred.beta)[1])[pred.beta$pop%in%sf_set & pred.beta$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{

					wt <- weather[J(data.table(pop=pred.beta[i]$pop,
												 date=seq(from=pred.beta[i]$col_date-(daysPrior), to=pred.beta[i]$col_date, by=1),
						 key="pop,date"))]


					wt[,daysPrior := pred.beta[i]$col_date - date]
					wt[,season := "fall"]
					wt[,set := "spring_fall"]

					wt

				}

				weather.pop.pre.frost <- foreach(i=c(1:dim(pred.beta)[1])[pred.beta$pop%in%frost_set & pred.beta$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{

					wt <- weather[J(data.table(pop=pred.beta[i]$pop,
												 date=seq(from=pred.beta[i]$col_date-(daysPrior), to=pred.beta[i]$col_date, by=1),
						 key="pop,date"))]


					wt[,daysPrior := pred.beta[i]$col_date - date]
					wt[,season := "pre"]
					wt[,set := "frost"]

					wt

				}

				weather.pop.post.frost <- foreach(i=c(1:dim(pred.beta)[1])[pred.beta$pop%in%frost_set & pred.beta$season=="frost"], .combine="rbind", .errorhandling="remove")%do%{

					wt <- weather[J(data.table(pop=pred.beta[i]$pop,
												 date=seq(from=pred.beta[i]$col_date-(daysPrior), to=pred.beta[i]$col_date, by=1),
						 key="pop,date"))]


					wt[,daysPrior := pred.beta[i]$col_date - date]
					wt[,season := "post"]
					wt[,set := "frost"]

					wt

				}


				weather.pop <- rbind(weather.pop.spring.sf, weather.pop.fall.sf,
									weather.pop.pre.frost, weather.pop.post.frost)

		### save
			pred.beta.new <- pred.beta
			weather.pop.new <- weather.pop[set=="spring_fall"]
			save(pred.beta.new, weather.pop.new, file="/mnt/pricey_1/dropPop/weather_pred.new.orig.swap.FET.Rdata")

		### use best model(s) from 19:1 to predict 20:1 & also use best permuted model as null distribution for r
			library(data.table)
			library(foreach)
			library(doMC)
			registerDoMC(20)
			library(epiR)

			load(file="/mnt/pricey_1/dropPop/weather_pred.new.orig.swap.FET.Rdata")
			load(file="/mnt/pricey_1/dropPop/threshDays_withPerms_10K.Rdata")
				weather.pop.new <- weather.pop[set=="spring_fall"]
				weather.pop <- weather.pop.new
				pred.beta <- pred.beta.new

			thresh.out.null <- thresh.out[boot.num!=1,list(maxR2=max(r2),
											  int=int[which.max(r2)],
											  tmax=tmax[which.max(r2)],
											  tmin=tmin[which.max(r2)],
											  sp.th=sp.th[which.max(r2)],
											  fall.th=fall.th[which.max(r2)],
											  set="null"),
										list(boot.num)]

			thresh.out.obs <- thresh.out[boot.num==1,
										 list(r2=r2,
										 	  p=mean(r2 >= thresh.out.null$maxR2),
										 	  int=int, tmax=tmax, tmin=tmin, set="obs", boot.num=1),
										 list(sp.th=sp.th,
										 	  fall.th=fall.th)]

			thresholds <- rbind(thresh.out.obs[p>.95], thresh.out.null, fill=T)

			o <- foreach(i=1:dim(thresholds)[1], .combine="rbind")%dopar%{
				print(paste(i, 	dim(thresholds)[1], sep=" / "))

				### calculate threshold days for ith model
					weather.pop.ag <- weather.pop[daysPrior<=daysPrior, list(tmax.max=mean(tmax[season%in%c("spring", "pre")] >= thresholds[i]$sp.th, na.rm=T),
																			 tmin.min=mean(tmin[season%in%c("fall", "post")] <= thresholds[i]$fall.th, na.rm=T)),
										 list(pop, set)]

				### merge
					setkey(weather.pop.ag, pop)
					setkey(pred.beta, pop)
					m <- merge(pred.beta, weather.pop.ag)[!is.na(col_date)][set=="spring_fall"][season=="spring"]
					#m <- merge(pred.beta, weather.pop.ag)[!is.na(col_date)][set=="frost"][season=="fall"]

				### model it!
					# mod <- summary(lm(beta~tmax.max + tmin.min, m))

				### predict based on 19:1 model parameters
					m[,pred := thresholds[i]$int +
								thresholds[i]$tmax * m$tmax.max +
								thresholds[i]$tmin * m$tmin.min]

				### output

					cbind(m[,list(r=c(cor(pred, beta),
									cor(pred, beta, method="kendall"),
									cor(pred, beta, method="spearman"),
									epi.ccc(beta, pred)$rho.c$est),
									method=c("pearson", "kendall", "spearman", "ccc")),
							list(class)],
						thresholds[i])

			}

			o.ccc <- o

			save(o.ccc, file="/mnt/pricey_1/dropPop/20_weather_ccc.orig.swap.FET.Rdata")

		### summarize
			load(file="/mnt/pricey_1/dropPop/20_weather_ccc.orig.swap.FET.Rdata")

			o.ccc[,list(pr.mu=mean(median(r[set=="obs"], na.rm=T) > r[set=="null"], na.rm=T),
						pr.min=mean(min(r[set=="obs"], na.rm=T) > r[set=="null"], na.rm=T),
						pr.max=mean(max(r[set=="obs"], na.rm=T) > r[set=="null"], na.rm=T),
						mu=median(r[set=="obs"], na.rm=T)),
					list(method, class)]

		o.ccc[,list(pr.mu=mean(median(r[set=="obs"], na.rm=T) > r[set=="null"], na.rm=T),
						pr.min=mean(min(r[set=="obs"], na.rm=T) > r[set=="null"], na.rm=T),
						pr.max=mean(max(r[set=="obs"], na.rm=T) > r[set=="null"], na.rm=T),
						mu=median(r[set=="obs"], na.rm=T),
						max=max(r[set=="obs"], na.rm=T),
						min=min(r[set=="obs"], na.rm=T)),
					list(method, class)]





		### plot
			load(file="/mnt/pricey_1/dropPop/20_weather_ccc.orig.swap.FET.Rdata")

			o.ccc.plot <- o.ccc[method=="pearson"]
			o.ccc.plot[set=="null", x:="null"]
			o.ccc.plot[set=="obs" & class=="orig", x:="orig"]
			o.ccc.plot[set=="obs" & class=="swap", x:="swap"]



			ggplot(data=o.ccc.plot, aes(x=x, y=r, color=x)) + geom_boxplot()




##################################################################################
### model averaged prediction plot: averaged across set of significant models ####
### has 19+1 best models, plus new set   									  ####
##################################################################################

	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(20)

	### load threshold model data
		load(file="/mnt/pricey_1/dropPop/threshDays_withPerms_10K.Rdata")

		thresh.out.ag <- thresh.out[,list(maxR2=max(r2)), list(boot.num)]
		thresh.out.obs <- thresh.out[boot.num==1,
									 list(r2=r2,
										  p=mean(r2 >= thresh.out.ag[boot.num!=1]$maxR2),
										  int=int, tmax=tmax, tmin=tmin),
									 list(sp.th=sp.th,
										  fall.th=fall.th)]

	### load climate & pbs data (pbs data is from 20 pop set)
		load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata") ### weather.pop, pbs.small

	### load in new populations
		#load( file="/mnt/pricey_1/dropPop/weather_pred.new.Rdata") ### pred.beta.new, weather.pop.new
		#load( file="/mnt/pricey_1/dropPop/weather_pred.new.orig.Rdata") ### pred.beta.new, weather.pop.new
		load(file="/mnt/pricey_1/dropPop/weather_pred.new.orig.swap.FET.Rdata")

	### merge together
		pbs.small[,set:="Core20"]
		pred.beta.new[,set:="NewSet"]

		pop.pred <- rbind(pbs.small[pop!="OUK_13"][dayTrue==T][,c("pop", "beta", "set"), with=F],
						  pred.beta.new[pop!="PA_15"][season=="spring"][,c("pop", "beta", "set", "class"), with=F], fill=T)
		setkey(pop.pred, pop)

		weather.pop <- rbind(weather.pop, weather.pop.new, fill=T)

	### average over models with p>.95
		pred.full <- foreach(i=which(thresh.out.obs$p>=.95))%dopar%{
			print(i)

			weather.pop.ag <- weather.pop[,list(tmax.max=mean(tmax[season=="spring"] >= thresh.out.obs[i]$sp.th, na.rm=T),
											tmin.min=mean(tmin[season=="fall"] <= thresh.out.obs[i]$fall.th, na.rm=T),
											tmax.n=length(unique(date[season=="spring"][!is.na(tmax)]))-1,
											tmin.n=length(unique(date[season=="fall"][!is.na(tmin)]))-1),
							 list(pop)]

			setkey(weather.pop.ag, pop)

			m <- merge(pop.pred, weather.pop.ag)[pop!="OUK_13"]
			m[,modelN:=i]

			m[,pred := thresh.out.obs[i]$int +
						thresh.out.obs[i]$tmax * tmax.max +
						thresh.out.obs[i]$tmin * tmin.min]

		#	t4 <- lm(beta ~ tmax.max + tmin.min, m)
		#	m[,pred:=predict(t4)]
		#

		#	pred.full <- data.table(tmax.max = rep(seq(from=0, to=.5, by=.001), each=length(seq(from=0, to=.5, by=.001))),
		#						   tmin.min = rep(seq(from=0, to=.5, by=.001), length(seq(from=0, to=.5, by=.001))))
		#	pred.full[,pred := predict(t4, newdata=pred.full)]
		#	pred.full[,modelN := i]

		#	list(m, cbind(pred.full, thresh.out.obs[i]))

		}
		pred.full <- rbindlist(pred.full)


	### make supplemental table 3
		pred.full[,sp.th := thresh.out.obs[pred.full$modelN]$sp.th]
		pred.full[,fall.th := thresh.out.obs[pred.full$modelN]$fall.th]

		write.csv(pred.full, file="~/supplementalTable3.csv")
	### extract components of pred.full
		#is.even <- function(x) x %% 2 == 0

		#preds <- rbindlist(lapply(c(1:length(pred.full)), function(x) pred.full[[x]][[2]]))
		#popStats <- rbindlist(lapply(c(1:length(pred.full)), function(x) pred.full[[x]][[1]]))

	### average
		#preds.ag <- preds[,list(pred.mu=mean(pred), pred.sd=sd(pred)),
		#					list(tmax.max, tmin.min)]



		popStats.ag <- pred.full[,list(pred.mu=median(pred), beta.obs=mean(beta),
									  pred.sd=sd(pred),
										tmax.max=mean(tmax.max), tmin.min=mean(tmin.min)),
								 list(pop, set, class)]

		save(popStats.ag, file="/mnt/pricey_1/dropPop/popStats.ag.swap.orig.FET.Rdata")



		ggplot(data=popStats.ag[is.na(class) | class=="orig"], aes(x=pred.mu, y=beta.obs, color=set)) + geom_point()
		ggplot(data=popStats.ag[is.na(class) | class=="swap"], aes(x=pred.mu, y=beta.obs, color=set)) + geom_point()


##################
### joint plot ###
##################

	### libraries
		library(cowplot)
		library(data.table)

	### load in summary data
		### concordance plots
			load(file="/mnt/pricey_1/dropPop/19_drop_1_concordance_data.Rdata") ### o19
			#load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/19_drop_1_concordance_data.Rdata")
				o19[,type := "Core 20 LOOCV"]


			load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_data.Rdata")  ### o20
			#load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/20_add_1_concordance_data.Rdata")
			#load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/20_add_1_concordance_data.orig.Rdata")

				o20[pop%in%c("PA_12", "WI_13", "PA_14", "PA_9"), type := "New set CV"]
				o20[is.na(type) & pop!="PA_2015" , type := "frost"]

				o20[type!="frost",list(mb=beta[which.min(log.sp.q.th)]), list(pop)]



			concordDat <- rbind(o19, o20[p<.05][n>=50], fill=T)

		### weather parameter fits: LOOCV (19)
			load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/paramContourData.Rdata")

		### weather prediction data: 19 population
			load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/pred_pred.plotData.old.Rdata")

		### New population prediction plot
			#load(file="/mnt/pricey_1/dropPop/20_weather_ccc.Rdata") ### o.ccc
			load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/20_weather_ccc.Rdata")

		### drop N CV predictability
			load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/dropN.predBeta.Rdata")

		### model averaged predictions & observed
			load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/popStats.ag.orig.Rdata")
			#load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/popStats.ag.Rdata")

		### popInfo

			locs <- fread("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/MapFigure/all_popinfo_formap_set.csv")

			locs[,jDay := as.POSIXlt(col_date, format = "%m/%d/%y")$yday]
			locs[,lat := as.numeric(as.character(tstrsplit(latitude, "\\ ")[[1]]))]

			locs[set=="Cline" & pop_name=="PA_9" ,jDay:=jDay+1]
			locs[name=="melPA_112011_WIN_FRT", set:="extra"]

			locs[,season:=factor(season, levels=c("spring", "fall"))]
			setnames(locs, "pop_name", "pop")


	### make plot components
		### concordance plot LOOCV
			concord.19.plot <- ggplot(data=concordDat[type=="Core 20 LOOCV"],
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
							theme(legend.direction="vertical",
								legend.justification=c(1,0),
								legend.position=c(1,0),
								legend.key.size = unit(0.35, "cm"),
								legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
								legend.text=element_text(size=8),
								legend.title=element_text(size=10),
								legend.title.align=.5,
								axis.text=element_text(size=8),
								axis.title=element_text(size=10))


		### concordance plot New Set CV
			concord.20.plot <- ggplot(data=concordDat[type=="New set CV"],
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
							ylim(.27, .75) +
							ylab("Fraction concordant") + xlab("Joint significance quantile")  +
							theme(legend.direction="vertical",
								legend.justification=c(1,0),
								legend.position=c(1,0),
								legend.key.size = unit(0.35, "cm"),
								legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
								legend.text=element_text(size=8),
								legend.title=element_text(size=10),
								axis.text=element_text(size=8),
								axis.title=element_text(size=10))

		### 19 choose 1 parameter surface & prediction plot
			load("~/predSurface.plot.Rdata")
			prediction19.plot <- ggplot() +
							geom_tile(data=pred.full,
									 aes(x=I(21*tmax.max), y=I(21*tmin.min), fill=pred), alpha=.85) +

							geom_point(data=m[tmax.max>0 | tmin.min>0][beta>0],
										aes(x=I(21*tmax.max), y=I(21*tmin.min), fill=beta), colour="black", pch=21, size=3,
										position=position_jitter(width=.5, height=.0)) +

							geom_point(data=m[tmax.max>0 | tmin.min>0][beta<0],
										aes(x=I(21*tmax.max), y=I(21*tmin.min), fill=beta), colour="white", pch=21, size=3) +

							geom_point(data=m[tmax.max==0 & tmin.min==0],
									aes(x=I(21*tmax.max), y=I(21*tmin.min), fill=beta), colour="black", pch=21, size=3,
									position=position_jitter(width=.5, height=.5)) +

							scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
																  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
																  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +
							scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
																  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
																  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +

							geom_line(data=pred.full[, list(tmin.min=tmin.min[which.min(abs(pred-0))]), tmax.max][tmin.min!=0],
										aes(x=I(21*tmax.max), y=I(21*tmin.min)), linetype="dashed") +

							ylab("No. days below fall threshold") + xlab("No. days above spring threshold") +
							theme(legend.direction="vertical",
								legend.justification=c(1,1),
								legend.position=c(1,1),
								legend.key.size = unit(0.35, "cm"),
								legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
								legend.text=element_text(size=8),
								legend.title=element_text(size=10),
								axis.text=element_text(size=8),
								axis.title=element_text(size=10))




		### 19 choose 1 parameter esimate plot

			contour19_1.plot <- ggplot() +
								geom_tile(data=thresh.out.obs, aes(x=I(sp.th/10), y=I(fall.th/10), fill=r2)) +
								geom_point(data=thresh.out.obs[p>.95], aes(x=I(sp.th/10), y=I(fall.th/10)), size=.05) +
								geom_point(data=thresh.out.obs[which.max(p)], aes(x=I(sp.th/10), y=I(fall.th/10)), size=.15, color="blue") +

								scale_fill_distiller(palette="Spectral",
													breaks=seq(0, 1, by=0.15),
													minor_breaks=seq(0,1, by=0.1),
													name=expression(r^2),
													limits=c(.0, .85)) +
								ylab("Fall lower threshold, °C") +
								xlab("Spring upper threshold, °C") +
								theme(legend.direction="vertical", legend.position=c(0.05,.95),
										legend.justification=c(0,1),
										legend.background=element_rect(fill="white", colour="black", linetype="solid", size=.25),
										legend.key.size = unit(0.35, "cm"),
										legend.text=element_text(size=8),
										legend.title=element_text(size=10),
										axis.text=element_text(size=8),
										axis.title=element_text(size=10))

		### New population CCC plot
			pred20.boxplot <- ggplot(data=o.ccc, aes(y=ccc.r, x=set, color=set)) +
								geom_boxplot() +
								ylab("CCC") + xlab(NULL) +
								theme(legend.position="none",
								axis.text=element_text(size=8),
								axis.title=element_text(size=10))

		### drop N pred
			dropN.pred <- ggplot(pred.beta.np, aes(y = abs(beta), x = nPops.f, color=nPops.f)) +
							geom_boxplot() +
							theme(legend.position="none")

		### model averaged prediction

			ma.pred <- ggplot() +
						geom_abline(slope=1, intercept=0) +
						geom_segment(data=popStats.ag[set=="Core20"],
									aes(x=pred.mu-pred.sd, xend=pred.mu+pred.sd, y=beta.obs, yend=beta.obs), alpha=.75) +
						geom_segment(data=popStats.ag[set=="NewSet"],
									aes(x=pred.mu-pred.sd, xend=pred.mu+pred.sd, y=beta.obs, yend=beta.obs), alpha=.75) +

						geom_point(data=popStats.ag[set=="Core20"], aes(x=pred.mu, y=beta.obs, fill=beta.obs),
									size=2.5, color="black", shape=21, alpha=.95) +
						geom_point(data=popStats.ag[set=="NewSet"], aes(x=pred.mu, y=beta.obs, fill=beta.obs),
									size=4, color="black", shape=23) +


						scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
																  space = "Lab", na.value = "grey50", guide = "colourbar",
																  name=NULL) +
						ylim(-.35, .25) + xlim(-.35, .25) +
						theme(legend.justification=c(1,0),
							  legend.position=c(1,0),
							  legend.key.size = unit(0.35, "cm"),
							  axis.text=element_text(size=8),
							  axis.title=element_text(size=10)) +
						xlab("Predicted") + ylab("Observed") +
						geom_hline(yintercept=0, lty="dashed") +
						geom_vline(xintercept=0, lty="dashed")


	### pred score & sampling params
		setkey(popStats.ag, pop)
		setkey(locs, pop)
		m <- merge(popStats.ag, locs)



		ggplot(data=m[set.y=="Core20"], aes(x=jDay, y=beta.obs)) + geom_point() + facet_grid(~season)


	### 19 choose 1 enrichment
			load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/19_drop_1_enrichment.Rdata")

		pred.beta.en <- foreach(i=unique(foo$pop), .combine="rbind")%do%{
			temp <- foo[pop==i]

		### playing
			#i <- "CUA_14"
			#temp <- small.o[pop==i][log.sp.q.th < 0]

			#t1 <- glm(frac ~ I(-1*log.sp.q.th), data=temp, family=binomial(), weights=small.o[pop==i]$n)
			#temp[,pred:=predict(t1)]
			#plot(frac~I(-1*log.sp.q.th), data=temp, type="l")
			#points(plogis(pred)~I(-1*log.sp.q.th), data=temp, type="l", col="red")



			glm.out <- summary(glm(en ~ I(-1*log.sp.q.th), data=temp, , weights=temp$n))

			data.table(pop=i,
						beta.en=glm.out$coef[2,1])
		}

		setkey(foo, pop)
		setkey(pred.beta.en, pop)

		foo <- merge(foo, pred.beta.en)


			ggplot(foo, aes(x=log.sp.q.th, y=en, group=pop, color=beta)) +
			geom_line() +

			scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
												  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
												  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +
			scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
												  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
												  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +

			theme(legend.direction="vertical",
			legend.justification=c(1,0),
			legend.position=c(1,0),
			legend.key.size = unit(0.35, "cm"),
			legend.background=element_rect(fill="white", colour="white", linetype="solid", size=.25),
			legend.text=element_text(size=8),
			legend.title=element_text(size=10),
			legend.title.align=.5,
			axis.text=element_text(size=8),
			axis.title=element_text(size=10))




	#### collection times figure w/ pred beta ###

		load("/Users/alanbergland/Documents/work/Projects/2011_6Dimensions/Paper_3_NESCent_seasonal/MS/PredictabilityFigure/pred_pred.plotData.Rdata")
		m[,season:=factor(season, levels=c("spring", "fall"))]
		m[,jDay := as.POSIXlt(date, format = "%m/%d/%y")$yday]

		setkey(m, season, pop)
		m[,cdd := c(cdd.spring[season=="spring"], cdd.fall[season=="fall"])]

		### version 1

		lat <-	ggplot(data=m[order(jDay)][dayTrue==T], aes(x=jDay, y=latitude, color=beta)) +
			scale_x_reverse() +
			facet_grid(.~season) + coord_flip() + coord_flip() + theme(legend.position="none") +
			geom_vline(xintercept=as.POSIXlt("6/20/2010", format = "%m/%d/%y")$yday, lty="dashed", lwd=.5) +
			geom_vline(xintercept=as.POSIXlt("9/22/2010", format = "%m/%d/%y")$yday, lty="dashed", lwd=.5) +
			geom_vline(xintercept=as.POSIXlt("12/21/2010", format = "%m/%d/%y")$yday, lty="dashed", lwd=.5) +
			geom_point() +
			scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
								  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
								  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +
			scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
								  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
								  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL)

		### version 1
		cdd <- 	ggplot(data=m[order(jDay)][dayTrue==T], aes(x=jDay, y=cdd, color=beta)) +
			scale_x_reverse() +
			facet_grid(.~season) + coord_flip() + coord_flip() + theme(legend.position="none", axis.text = element_text(angle=90)) +
			geom_vline(xintercept=as.POSIXlt("6/20/2010", format = "%m/%d/%y")$yday, lty="dashed", lwd=.5) +
			geom_vline(xintercept=as.POSIXlt("9/22/2010", format = "%m/%d/%y")$yday, lty="dashed", lwd=.5) +
			geom_vline(xintercept=as.POSIXlt("12/21/2010", format = "%m/%d/%y")$yday, lty="dashed", lwd=.5) +
			geom_point() +
			scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
								  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
								  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +
			scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
								  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
								  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL)


		plot_grid(lat, cdd, labels=c('A', 'B'))


#### Supplemental table 2
#### model parameters for best-fit models

	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(20)

	### load threshold model data
		load(file="/mnt/pricey_1/dropPop/threshDays_withPerms_10K.Rdata")

		thresh.out.ag <- thresh.out[,list(maxR2=max(r2)), list(boot.num)]
		thresh.out.obs <- thresh.out[boot.num==1,
									 list(r2=r2,
										  p=mean(r2 >= thresh.out.ag[boot.num!=1]$maxR2),
										  int=int, tmax=tmax, tmin=tmin),
									 list(sp.th=sp.th,
										  fall.th=fall.th)]

	### write out Supplemental Table 2
		write.csv(thresh.out.obs, file="~/SupplementalTable2.csv", quote=F, row.names=F)


	### load climate & pbs data (pbs data is from 20 pop set)
		load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata") ### weather.pop, pbs.small

	### load in new populations
		#load( file="/mnt/pricey_1/dropPop/weather_pred.new.Rdata") ### pred.beta.new, weather.pop.new
		#load( file="/mnt/pricey_1/dropPop/weather_pred.new.orig.Rdata") ### pred.beta.new, weather.pop.new
		load(file="/mnt/pricey_1/dropPop/weather_pred.new.orig.swap.FET.Rdata")

	### merge together
		pbs.small[,set:="Core20"]
		pred.beta.new[,set:="NewSet"]

		pop.pred <- rbind(pbs.small[pop!="OUK_13"][dayTrue==T][,c("pop", "beta", "set"), with=F],
						  pred.beta.new[pop!="PA_15"][season=="spring"][,c("pop", "beta", "set", "class"), with=F], fill=T)
		setkey(pop.pred, pop)

		weather.pop <- rbind(weather.pop, weather.pop.new, fill=T)

	### iterate through to generate Supplemental Table 3
		o <- foreach(i=which(thresh.out.obs$p>.95), .combine="rbind")%do%{
			weather.pop[daysPrior<=21,list(nSpringDays=sum(tmax>=thresh.out.obs$sp.th[i], na.rm=T), nFallDays=sum(tmin<=thresh.out.obs$fall.th[i], na.rm=T),
								springThresh=thresh.out.obs$sp.th[i], fallThresh=thresh.out.obs$fall.th[i]), list(pop)]
		}





















	### average over models with p>.95
		pred.full <- foreach(i=which(thresh.out.obs$p>=.95))%dopar%{
			print(i)

			weather.pop.ag <- weather.pop[,list(tmax.max=sum(tmax[season=="spring"] >= thresh.out.obs[i]$sp.th, na.rm=T),
											    tmin.min=sum(tmin[season=="fall"] <= thresh.out.obs[i]$fall.th, na.rm=T),
											    spring.thresh=thresh.out.obs[i]$sp.th,
											    fall.thresh=thresh.out.obs[i]$fall.th),
							 list(pop)]

			setkey(weather.pop.ag, pop)

			m <- merge(pop.pred, weather.pop.ag)[pop!="OUK_13"]
			m[,modelN:=i]

			m[,pred := thresh.out.obs[i]$int +
						thresh.out.obs[i]$tmax * tmax.max +
						thresh.out.obs[i]$tmin * tmin.min]

		#	t4 <- lm(beta ~ tmax.max + tmin.min, m)
		#	m[,pred:=predict(t4)]
		#

		#	pred.full <- data.table(tmax.max = rep(seq(from=0, to=.5, by=.001), each=length(seq(from=0, to=.5, by=.001))),
		#						   tmin.min = rep(seq(from=0, to=.5, by=.001), length(seq(from=0, to=.5, by=.001))))
		#	pred.full[,pred := predict(t4, newdata=pred.full)]
		#	pred.full[,modelN := i]

		#	list(m, cbind(pred.full, thresh.out.obs[i]))

		}
		pred.full <- rbindlist(pred.full)

	### extract components of pred.full
		#is.even <- function(x) x %% 2 == 0

		#preds <- rbindlist(lapply(c(1:length(pred.full)), function(x) pred.full[[x]][[2]]))
		#popStats <- rbindlist(lapply(c(1:length(pred.full)), function(x) pred.full[[x]][[1]]))

	### average
		#preds.ag <- preds[,list(pred.mu=mean(pred), pred.sd=sd(pred)),
		#					list(tmax.max, tmin.min)]



		popStats.ag <- pred.full[,list(pred.mu=median(pred), beta.obs=mean(beta),
									  pred.sd=sd(pred),
										tmax.max=mean(tmax.max), tmin.min=mean(tmin.min)),
								 list(pop, set, class)]

		save(popStats.ag, file="/mnt/pricey_1/dropPop/popStats.ag.swap.orig.FET.Rdata")



		ggplot(data=popStats.ag[is.na(class) | class=="orig"], aes(x=pred.mu, y=beta.obs, color=set)) + geom_point()
		ggplot(data=popStats.ag[is.na(class) | class=="swap"], aes(x=pred.mu, y=beta.obs, color=set)) + geom_point()




























		### version 1
			setkey(m, season, pop)

			m[,cdd := c(cdd
			spring.plot <- ggplot(data=m[order(cdd.spring)][season=="spring"], aes(x=cdd.spring, y=jDay, color=beta)) +
							scale_x_reverse() +
							facet_grid(.~season) + coord_flip() + coord_flip() + theme(legend.position="none")  +
							geom_point() +
							scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
												  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
												  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +
							scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
												  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
												  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL)

			fall.plot <- ggplot(data=m[order(cdd.fall)][season=="fall"], aes(x=cdd.fall, y=jDay, color=beta)) +
						scale_x_reverse() +
						facet_grid(.~season) + coord_flip() + coord_flip() + theme(legend.position="none")  +
						geom_point()  +
						scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
											  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
											  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL) +
						scale_fill_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
											  space = "Lab", na.value = rgb(red=1, blue=1, green=1, alpha=1), guide = "colourbar",
											  limits=1*range(concordDat[type%in%c("Core 20 LOOCV", "New set CV")]$beta), name=NULL)


		plot_grid(spring.plot, fall.plot, ncol=2)




		### correlations
			cor.test(m[season=="spring"]$cdd, m[season=="spring"]$beta)
			cor.test(m[season=="fall"]$cdd, m[season=="fall"]$beta)

			cor.test(m[season=="fall"]$jDay, m[season=="fall"]$beta)

			cor.test(m[season=="fall"]$latitude, m[season=="fall"]$beta)

			cor.test(m[season=="fall"]$cdd - m[season=="spring"]$cdd , m[season=="fall"]$beta)

			summary(lm(

	### composite plot
		### 4 panel
			multiPlot4 <- plot_grid(concord.19.plot, contour19_1.plot, concord.20.plot, ma.pred,
									labels = c('A', 'B',
												'C', 'D'),
									nrow=2, ncol=2,
									vjust=.5)

			save_plot(multiPlot4, file="~/plot4.orig.pdf", base_height=6)


		### 6 panel
			multiPlot <- plot_grid(concord.19.plot, dropN.pred,
					contour19_1.plot, prediction19.plot,
					concord.20.plot, pred20.boxplot,
					labels = c('A', 'B',
								'C', 'D',
								'E', 'F'),
					nrow=3, ncol=2,
					vjust=.5)

			save_plot(multiPlot, file="~/plot.pdf", base_aspect_ratio = 1/1.5, base_height=8)



				plot_grid(contour19_1.plot, prediction19.plot, nrow=1)




		### doodles

			### how many sites are we talking about here

			plot(log10(abs(frac-.5)*n)~log.sp.q.th, concordDat, pch=19, cex=.5)









		bottomRow <- plot_grid(concord.20.plot, pred20.boxplot,
							labels=c('C', 'D'), nrow=1, ncol=2)

		plot_grid(topRow, bottomRow, nrow=2)


	prediction19.plot


	plot_grid(concord.19.plot, contour19_1.plot,

					concord.20.plot, pred20.boxplot ,
					labels = c('A', 'B', 'C', 'D'), nrow=2, ncol=2)



	ggplot() +
	geom_point(data=thresh.out.obs[p>.95], aes(x=tmax, y=tmin)) +
	ylim(-1.5, 1.5) + xlim(-1.5, 1.5) +
	geom_hline(yintercept=0, linetype="dashed") +
	geom_vline(xintercept=0, linetype="dashed")











	dat[seas.p<1e-5, list(]











qplot(10^(1:6), 10^(1:6)) + scale_x_log10(breaks = 10^(1:6),
											minor_breaks = log(c(sapply(10^(1:6), function(x) seq(0, x, x/10))), 10))






		ggplot(data=dat[seas.p<.5], aes(x=pos, y=-log10(seas.p), color=chrom)) + geom_point() + facet_grid(~chrom)
			###	make plot
				concord20.plot <-
				ggplot(data=o[p<.05][pop%in%c("PA_12", "PA_15", "WI_13", "PA_14", "PA_9")],
													aes(y=frac, x=log.sp.q.th, group=pop,
															color=beta)) +
												geom_line() +
				scale_colour_gradient2(low = "blue", mid="orange", high = "red",
				  space = "Lab", na.value = "grey50", guide = "colourbar")



			### save plot
				save(concord20.plot, file="/mnt/pricey_1/dropPop/20_add_1_concordance_plot.Rdata")











		###	make plot

			concord.plot <- ggplot(data=,
								aes(y=frac, x=log.sp.q.th, group=pop,
										color=beta)) +
							geom_line() + scale_colour_gradientn(colours=c("navy", "blue", "orange", "orangered", "red"),
													  space = "Lab", na.value = "grey50", guide = "colourbar") +
							geom_hline(yintercept=.5, linetype="dashed")



		### show plot
			concord.plot

		### save plot
			save(concord.plot, file="/mnt/pricey_1/dropPop/19_drop_1_concordance_plot.Rdata")

scale_colour_distiller(palette = "Spectral") +
							scale_fill_distiller(palette="Spectral", limits=c(0, 1),
												breaks=seq(0,1, by=0.2),
												minor_breaks=seq(0,1, by=0.1))

scale_colour_gradient2(low = "blue", mid="orange", high = "red",
													  space = "Lab", na.value = "grey50", guide = "colourbar")





	load(file="/mnt/pricey_1/dropPop/19_drop_1_concordance_plot.Rdata")	### concord.plot
	load(file="/mnt/pricey_1/dropPop/contour19_1.plot.Rdata")	### contour19_1.plot
	load(file="/mnt/pricey_1/dropPop/20_add_1_concordance_plot.Rdata") ### concord20.plot

	plot_grid(concord.plot, contour19_1.plot, concord20.plot)




















			ggplot(data=o, aes(y=sign.test, x=set, color=set)) + geom_boxplot()



			pred20.plot <- ggplot() +
			geom_tile(data=o, aes(x=sp.th, fall.th, fill=newPred.r)) +
			scale_colour_distiller(palette = "Spectral") +
			scale_fill_distiller(palette="Spectral", limits=c(-1, 1),
								breaks=seq(-1,1, by=0.2),
								minor_breaks=seq(-1,1, by=0.1)) +
			geom_point(data=thresh.out.obs[p>.95], aes(x=sp.th, y=fall.th), size=.5)

			save(pred20.plot, file="/mnt/pricey_1/dropPop/pred20.plot.Rdata")

			o <- rbindlist(o)































		### best model
			bestModel <- thresh.out[which.max(r2)]

			weather.pop.ag <- weather.pop[daysPrior<=daysPrior, list(tmax.max=mean(tmax[season%in%c("spring", "pre")] >= bestModel$sp.th, na.rm=T),
									 				    	 		 tmin.min=mean(tmin[season%in%c("fall", "post")] <= bestModel$fall.th, na.rm=T),
									 				    	 tmax.mean=mean(tmax[season=="spring"], na.rm=T),
									 				    	 tmin.mean=mean(tmin[season=="fall"], na.rm=T)),
								 list(pop, set)]

			setkey(weather.pop.ag, pop)
			setkey(pred.beta, pop)
			m <- merge(pred.beta, weather.pop.ag)[!is.na(col_date)]

			t4 <- lm(beta~tmax.max + tmin.min, m[set=="frost"])


			pred.full <- data.table(tmax.max = rep(seq(from=0, to=.9, by=.01), each=length(seq(from=0, to=.9, by=.01))),
								   tmin.min = rep(seq(from=0, to=.9, by=.01), length(seq(from=0, to=.9, by=.01))))
			pred.full[,pred := predict(t4, newdata=pred.full)]

		### plot: add it to old prediction plot
			load(file="/mnt/pricey_1/dropPop/pred_pred.plot.Rdata")

			pred_pred.plot +
			geom_point(data=m, aes(x=tmax.max, y=tmin.min, fill=beta), colour="black", pch=23, size=5)


			pred_pred.plot <- ggplot() +
							geom_tile(data=pred.full,
									 aes(x=tmax.max, y=tmin.min, fill=pred)) +
							geom_point(data=m, aes(x=tmax.max, y=tmin.min, fill=beta), colour="black",pch=21, size=5) +
							scale_colour_distiller(palette = "Spectral") +
							scale_fill_distiller(palette="Spectral")



			save(pred_pred.plot, file="/mnt/pricey_1/dropPop/pred_pred.plot.Rdata")


































	### calculate # days above and below threshold
		registerDoMC(20)
		threshDays.coarse <- foreach(sp.th=seq(from=0, to=400, by=5), .combine="rbind")%dopar%{
			foreach(fall.th=seq(from=0, to=400, by=5), .combine="rbind")%do%{
				print(paste(sp.th, fall.th, sep= " / "))

				spring.tMax <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{


						weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
											 date=seq(from=pbs[i]$date-(15), to=pbs[i]$date, by=1),
											 key="pop,date"))]

						weather.pop[,list(spring.tMax=sum(tmax>sp.th, na.rm=T),
										  date=pbs[i]$date),
									list(pop)]

				}

				fall.tMin <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{


						weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
											 date=seq(from=pbs[i]$date-(15), to=pbs[i]$date, by=1),
											 key="pop,date"))]

						weather.pop[,list(fall.tMin=sum(tmin<fall.th, na.rm=T),
										  date=pbs[i]$date),
									list(pop)]

				}

				setkey(spring.tMax, pop)
				setkey(fall.tMin, pop)

				tMax.tMin <- merge(spring.tMax, fall.tMin)

				setkey(tMax.tMin, pop)

				m <- merge(pbs[season=="spring"], tMax.tMin, allow.cartesian=T)[pop!="OUK_13"][dayTrue==T]

				foreach(filter.set=unique(m$class), .errorhandling="remove", .combine="rbind")%do%{

					tm <- lm(beta~spring.tMax + fall.tMin , m[class==filter.set])
					stm <- summary(tm)

					data.table(sp.th=sp.th, fall.th=fall.th,
							r2=stm$r.squared,
							class=filter.set)
				}

			}
		}

		threshDays.fine <- foreach(sp.th=seq(from=200, to=400, by=5), .combine="rbind")%dopar%{
			foreach(fall.th=seq(from=-50, to=150, by=5), .combine="rbind")%do%{
				print(paste(sp.th, fall.th, sep= " / "))

				spring.tMax <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{


						weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
											 date=seq(from=pbs[i]$date-(30), to=pbs[i]$date, by=1),
											 key="pop,date"))]

						weather.pop[,list(spring.tMax=sum(tmax>sp.th, na.rm=T),
										  date=pbs[i]$date),
									list(pop)]

				}

				fall.tMin <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{


						weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
											 date=seq(from=pbs[i]$date-(30), to=pbs[i]$date, by=1),
											 key="pop,date"))]

						weather.pop[,list(fall.tMin=sum(tmin<fall.th, na.rm=T),
										  date=pbs[i]$date),
									list(pop)]

				}

				setkey(spring.tMax, pop)
				setkey(fall.tMin, pop)

				tMax.tMin <- merge(spring.tMax, fall.tMin)

				setkey(tMax.tMin, pop)

				m <- merge(pbs[season=="spring"], tMax.tMin, allow.cartesian=T)[pop!="OUK_13"][dayTrue==T]
				setkey(m, pop)

				foreach(filter.set=unique(m$class), .errorhandling="remove", .combine="rbind")%do%{

					tm <- lm(beta~spring.tMax + fall.tMin , m[class==filter.set][!duplicated(m)])
					stm <- summary(tm)

					data.table(sp.th=sp.th, fall.th=fall.th,
							r2=stm$r.squared,
							class=filter.set)
				}

			}
		}

		theshDays <- rbind(threshDays.coarse, threshDays.fine)

		threshDays.model <- loess(r2~fall.th * sp.th, threshDays[class=="pass"], degree=2, span=.05)
		threshDays[class=="pass", pred := predict(threshDays.model)]

		ggplot() +
		geom_tile(data=threshDays[class=="pass"], aes(x=sp.th, y=fall.th, fill=r2)) +
		scale_colour_distiller(palette = "Spectral") +
		scale_fill_distiller(palette="Spectral", limits=c(0, 1),
							breaks=seq(0,1, by=0.2),
							minor_breaks=seq(0,1, by=0.1)) +
		geom_contour(data=threshDays[class=="pass"], aes(x=sp.th, y=fall.th, z=pred), binwidth=.05)

		save(threshDays, file="/mnt/pricey_1/dropPop/threshDays.Rdata")


		ggplot() +
		geom_tile(data=threshDays.fine, aes(x=sp.th, y=fall.th, fill=r2)) +
		scale_colour_distiller(palette = "Spectral") +
		scale_fill_distiller(palette="Spectral", limits=c(0, 1),
							breaks=seq(0,1, by=0.2),
							minor_breaks=seq(0,1, by=0.1)) +
		facet_wrap(~class)






				weather.pop.ag <- weather.pop[,list(spring.tMax = mean(tmax[season=="spring"] >= 350, na.rm=T),
													fall.tmin = mean(tmin[season=="fall"] <= 90, na.rm=T)),
												list(pop)]

				setkey(pred.beta, pop)
				setkey(weather.pop.ag, pop)

				pbw <- merge(pred.beta[season=="spring"], weather.pop.ag)























			getWeatherData <- function(wd
				### debug
					wd <- weather
					popName <- "PA_9"
					col_date=as.Date(c("7/15/09", "11/15/09"), format="%m/%d/%y")
					spring_fall <- T
					daysPrior <- 15
					tmax.thresh <- 	360
					tmin.thresh <- 90

			if(spring_fall==T) {
				### first, get out weather data for days prior to spring
					spring.weather.pop <- weather[J(data.table(pop=popName,
								 						date=seq(from=min(col_date)-(daysPrior), to=min(col_date), by=1),
								 						key="pop,date"))]

					spring.tMax <- spring.weather.pop[,list(spring.tMax=sum(tmax>tmax.thresh, na.rm=T),
										 			 date=col_date),
												list(pop)]


					fall.weather.pop <- weather[J(data.table(pop=popName,
								 						date=seq(from=min(col_date)-(daysPrior), to=min(col_date), by=1),
								 						key="pop,date"))]

					fall.tMin <- fall.weather.pop[,list(spring.tMax=sum(tmin<tmin.thresh, na.rm=T),
										 			 date=col_date),
												list(pop)]




			}





			spring.tMax <- foreach(i=c(1:dim(pred.beta)[1])[sea
				weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
									 date=seq(from=pbs[i]$date-(daysPrior), to=pbs[i]$date, by=1),
									 key="pop,date"))]

					weather.pop[,list(spring.tMax=sum(tmax>bestModel$sp.th, na.rm=T),
									  date=pbs[i]$date),
								list(pop)]

				}










































		### load in metadata

		daysPrior <- 15

		### calculate number of days above and below thresholds in 15 days prior to sampling


		### tack into predictability scores
			setkey(spring.tMax, pop)
			setkey(fall.tMin, pop)

			tMax.tMin <- merge(spring.tMax, fall.tMin)

			setkey(tMax.tMin, pop)

			m <- merge(pbs[season=="spring"], tMax.tMin, allow.cartesian=T)[pop!="OUK_13"][dayTrue==T]

			m <- m[class=="pass"]
			setkey(m, pop)
			m <- m[!duplicated(pop)]

		### make model & full predictions
			t4 <- lm(beta~spring.tMax + fall.tMin, m)












	### get out best model (sp.th=360; fall.th=90) == totally sensible
		### extract out best model from above
			#bestModel <- threshDays[class=="pass"][which.max(r2)]
			bestModel <- threshDays.fine[class=="pass"][which.max(r2)]

		### days prior
			daysPrior <- 15

		### calculate number of days above and below thresholds in 15 days prior to sampling

			spring.tMax <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{


					weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
										 date=seq(from=pbs[i]$date-(daysPrior), to=pbs[i]$date, by=1),
										 key="pop,date"))]

					weather.pop[,list(spring.tMax=sum(tmax>bestModel$sp.th, na.rm=T),
									  date=pbs[i]$date),
								list(pop)]

			}

			fall.tMin <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{


					weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
										 date=seq(from=pbs[i]$date-(daysPrior), to=pbs[i]$date, by=1),
										 key="pop,date"))]

					weather.pop[,list(fall.tMin=sum(tmin<bestModel$fall.th, na.rm=T),
									  date=pbs[i]$date),
								list(pop)]

			}

		### tack into predictability scores
			setkey(spring.tMax, pop)
			setkey(fall.tMin, pop)

			tMax.tMin <- merge(spring.tMax, fall.tMin)

			setkey(tMax.tMin, pop)

			m <- merge(pbs[season=="spring"], tMax.tMin, allow.cartesian=T)[pop!="OUK_13"][dayTrue==T]

			m <- m[class=="pass"]
			setkey(m, pop)
			m <- m[!duplicated(pop)]

		### make model & full predictions
			t4 <- lm(beta~spring.tMax + fall.tMin, m)


			pred.full <- data.table(spring.tMax = rep(c(10:15), each=length(0:1)),
														fall.tMin = rep(c(0:1), length(0:15)))
			pred.full[,pred := predict(t4, newdata=pred.full)]

		### plot
			ggplot() +
			geom_tile(data=pred.full,
					 aes(x=spring.tMax, y=fall.tMin, fill=pred)) +
			geom_point(data=m, aes(x=spring.tMax, y=fall.tMin, fill=beta), colour="black",pch=21, size=5, position="jitter") +
			scale_colour_distiller(palette = "Spectral") +
			scale_fill_distiller(palette="Spectral")

















	### for all sites
				### generate binned quantile values
					m[,log.sp.q.all.bin := round(log.sp.q.all, 2)]
					m[,log.full.q.all.bin := round(log.full.q.all, 2)]

					wins <- unique(c(m$log.sp.q.all.bin, m$log.full.q.all.bin))
					wins <- wins[order(wins)]

					setkey(m, log.sp.q.all.bin)

				### calculate scores
					registerDoMC(5)
					o.all <- foreach(i = 1:(length(wins)-50), .combine="rbind", .export="m")%dopar%{
						print(i)

						m.temp <- m[J(wins[1:i]), nomatch=0]
						setkey(m.temp, log.full.q.all.bin)
						m.temp <- m.temp[J(wins[1:i]), nomatch=0]

						m.ag <- m.temp[,
										list(TT=sum(sign(seas.coef)==sign(sp.coef)),
											n=length(seas.coef),
											log.sp.q.th=wins[i],
											log.full.q.th=wins[i],
											class="all"),
									  list(pop)]

						m.ag
					}





	### attemp #2
	### 5 day average at various points in time in the past. Gets to ~75% r2; convoluted
		spring.tMax <- foreach(j=0:90, .combine="rbind")%dopar%{
			foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{
				print(i)

				weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
									 date=seq(from=pbs[i]$date-(j), to=pbs[i]$date-(j-5), by=1),
									 key="pop,date"))]

				weather.pop[,list(spring.tMax=mean(tmax, na.rm=T),
									spring.tMin=mean(tmin, na.rm=T),
								  date=pbs[i]$date,
								  win=j),
							list(pop)]

			}
		}

		fall.tMin <- foreach(j=0:90, .combine="rbind")%dopar%{
			foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{

				weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
									 date=seq(from=pbs[i]$date-(j), to=pbs[i]$date-(j-5), by=1),
									 key="pop,date"))]

				weather.pop[,list(fall.tMin=mean(tmin, na.rm=T),
								  fall.tMax=mean(tmax, na.rm=T),
								  date=pbs[i]$date,
								  win=j),
							list(pop)]

			}
		}





		wwn <- foreach(i=c("obs", "shuf1", "shuf2", "shuf3", "shuf4"), .combine="rbind")%do%{

			pbs.use <- pbs[season=="spring"][pop!="OUK_13"][dayTrue==T][class=="pass"]
			#pbs.use <- pbs[season=="spring"][pop!="OUK_13"]

			if(i=="obs") pbs.use[,beta := beta]
			if(grepl("shuf", i)) pbs.use[,beta := sample(beta)]


			o <- foreach(j=1:90, .combine="rbind")%dopar%{
				foreach(k=1:90, .combine="rbind")%do%{
					spring.temp <- spring.tMax[win==j]
					fall.temp <- fall.tMin[win==k]

					setkey(spring.temp, pop)
					setkey(fall.temp, pop)

					sf.temp <- merge(spring.temp, fall.temp)

					setkey(sf.temp, pop)
					setkey(pbs.use, pop)

					temp <- merge(pbs.use, sf.temp)

					#temp <- temp[beta>0]

					t1 <- lm(beta~spring.tMax + fall.tMin, temp)
					t2 <- lm(beta~spring.tMin + fall.tMax, temp)

					t3 <- lm(beta~spring.tMax + fall.tMax, temp)
					t4 <- lm(beta~spring.tMin + fall.tMin, temp)

					t5 <- lm(beta~spring.tMin + fall.tMin + spring.tMax + fall.tMax, temp)

					temp[,pred:=predict(t1)]

					data.table(r2 = c(summary(t1)$r.squared, summary(t2)$r.squared,
										summary(t3)$r.squared, summary(t4)$r.squared,
										summary(t5)$r.squared),

								springCoef=c(coef(t1)[2], coef(t2)[2],
											coef(t3)[2], coef(t4)[2],
											coef(t5)[2]),

								fallCoef=c(coef(t1)[3], coef(t2)[3],
											coef(t3)[3], coef(t4)[3],
											coef(t5)[3]),

								springCoef.p=c(summary(t1)$coef[2,4], summary(t2)$coef[2,4],
											summary(t3)$coef[2,4], summary(t4)$coef[2,4],
											summary(t5)$coef[2,4]),

								fallCoef.p=c(summary(t1)$coef[3,4], summary(t2)$coef[3,4],
											summary(t3)$coef[3,4], summary(t4)$coef[3,4],
											summary(t5)$coef[3,4]),

								model.p = c(1 - pf(summary(t1)$fstatistic[1], summary(t1)$fstatistic[2], summary(t1)$fstatistic[3]),
											1 - pf(summary(t2)$fstatistic[1], summary(t2)$fstatistic[2], summary(t2)$fstatistic[3]),
											1 - pf(summary(t3)$fstatistic[1], summary(t3)$fstatistic[2], summary(t3)$fstatistic[3]),
											1 - pf(summary(t4)$fstatistic[1], summary(t4)$fstatistic[2], summary(t4)$fstatistic[3]),
											1 - pf(summary(t5)$fstatistic[1], summary(t5)$fstatistic[2], summary(t5)$fstatistic[3])),


								class=c("springMax-fallMin", "springMin-fallMax", "springMax-fallMax", "springMin-fallMin", "all"),
								spring.win=j,
								fall.win=k,
								order=i)
				}
			}

		}


		### plot
			ggplot(spring.tMax, aes(x=win, y=spring.tMax, group=pop, color=pop)) + geom_line()

			setkey(spring.tMax, pop, win)
			setkey(fall.tMin, pop, win)

			tMax.tMin <- merge(spring.tMax, fall.tMin)

			ggplot(data=tMax.tMin, aes(x=spring.tMax, y=fall.tMin, group=pop, color=win)) +
			geom_line() +
			facet_wrap(~pop, ncol=4, nrow=5)

			ggplot() +
			geom_line(data=tMax.tMin, aes(x=win, y=spring.tMax), color="red") +
			geom_line(data=tMax.tMin, aes(x=win, y=fall.tMin), color="blue") +
			facet_wrap(~pop, ncol=4, nrow=5)




			ggplot(spring.tMax, aes(x=win, y=spring.tMax, group=pop, color=pop)) + geom_line()


		### make heat map plot
			ggplot(data=wwn[class!="all"], aes(x=spring.win, y=fall.win, fill=r2)) + geom_tile()	+ facet_grid(order~class) +
				scale_fill_distiller(palette = "Spectral")	+
				geom_tile()


	### find best parameter range for environmental correlation and pull out specific model

		best <- wwn[class!="all" & order=="obs"][r2>=max(wwn[class!="all" & order=="obs"]$r2)]

		spring.temp <- spring.tMax[win==best$spring.win]
		fall.temp <- fall.tMin[win==best$fall.win]

		setkey(spring.temp, pop)
		setkey(fall.temp, pop)

		sf.temp <- merge(spring.temp, fall.temp)

		setkey(sf.temp, pop)
		setkey(pbs, pop)

		temp <- merge(pbs[season=="spring"], sf.temp)[pop!="OUK_13"][dayTrue==T]

		t4 <- lm(beta~spring.tMax + fall.tMin, temp)

		temp[,pred := predict(t4)]

		pred.full <- data.table(spring.tMax = rep(c(0:350), each=length(0:350)),
													fall.tMin = rep(c(0:350), length(0:350)))
		pred.full[,pred := predict(t4, newdata=pred.full)]

		### plot
			ggplot() +
			geom_tile(data=pred.full[fall.tMin>=min(temp$fall.tMin) & fall.tMin<=max(temp$fall.tMin) &
									spring.tMax>=min(temp$spring.tMax) & spring.tMax<=max(temp$spring.tMax)]
					, aes(x=spring.tMax, y=fall.tMin, fill=pred)) +
			geom_point(data=temp, aes(x=spring.tMax, y=fall.tMin, fill=beta), colour="black",pch=21, size=5) +
			scale_colour_distiller(palette = "Spectral") +
			scale_fill_distiller(palette="Spectral")

					+
			geom_point(data=temp, aes(x=spring.tMax, y=fall.tMin, color=beta))

		ggplot() +
		geom_point(data=temp, aes(x=spring.tMax, y=fall.tMin, color=beta)) +
		scale_colour_distiller(palette = "Spectral")













## attemp #3
	### take the single spring.max & fall.min & get r^2 of 0.57. Not bad
		o <- foreach(sp.win = c(1:5), .combine="rbind")%do%{
			foreach(fall.win = c(1:5), .combine="rbind")%do%{
				print(paste(sp.win, fall.win, sep=" / "))

				spring.tMax <- foreach(j=0:90, .combine="rbind")%dopar%{
					foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{


						weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
											 date=seq(from=pbs[i]$date-(j), to=pbs[i]$date-(j-sp.win), by=1),
											 key="pop,date"))]

						weather.pop[,list(spring.tMax=mean(tmax, na.rm=T),
											spring.tMin=mean(tmin, na.rm=T),
										  date=pbs[i]$date,
										  win=j),
									list(pop)]

					}
				}

				fall.tMin <- foreach(j=0:90, .combine="rbind")%dopar%{
					foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{

						weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
											 date=seq(from=pbs[i]$date-(j), to=pbs[i]$date-(j-fall.win), by=1),
											 key="pop,date"))]

						weather.pop[,list(fall.tMin=mean(tmin, na.rm=T),
										  fall.tMax=mean(tmax, na.rm=T),
										  date=pbs[i]$date,
										  win=j),
									list(pop)]

					}
				}

				setkey(spring.tMax, pop, win)
				setkey(fall.tMin, pop, win)

				tMax.tMin <- merge(spring.tMax, fall.tMin)


				maxmin <- tMax.tMin[,list(spring.tMax.win=win[which.max(spring.tMax)],
										  fall.tMin.win=win[which.min(fall.tMin)],
										  spring.tMax=max(spring.tMax),
										  fall.tMin=min(fall.tMin),
										  spring.tMin=min(spring.tMin),
										  fall.tMax=max(fall.tMax)), pop]
				setkey(maxmin, pop)

				m <- merge(pbs[season=="spring"], maxmin)[pop!="OUK_13"][dayTrue==T]

				tm <- lm(beta~spring.tMax + fall.tMin , m)
				stm <- summary(tm)

				data.table(sp.win=sp.win, fall.win=fall.win,
							r2=stm$r.squared)
		}
	}




### attempt #4
	### calculate # days above and below threshold
	registerDoMC(20)
	o <- foreach(sp.th=seq(from=0, to=400, by=5), .combine="rbind")%dopar%{
		foreach(fall.th=seq(from=0, to=400, by=5), .combine="rbind")%do%{
			print(paste(sp.th, fall.th, sep= " / "))

			spring.tMax <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="spring"], .combine="rbind", .errorhandling="remove")%do%{


					weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
										 date=seq(from=pbs[i]$date-(90), to=pbs[i]$date, by=1),
										 key="pop,date"))]

					weather.pop[,list(spring.tMax=sum(tmax>sp.th, na.rm=T),
									  date=pbs[i]$date),
								list(pop)]

			}

			fall.tMin <- foreach(i=c(1:dim(pbs)[1])[pbs$season=="fall"], .combine="rbind", .errorhandling="remove")%do%{


					weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
										 date=seq(from=pbs[i]$date-(90), to=pbs[i]$date-(90-90), by=1),
										 key="pop,date"))]

					weather.pop[,list(fall.tMin=sum(tmin<fall.th, na.rm=T),
									  date=pbs[i]$date),
								list(pop)]

			}

			setkey(spring.tMax, pop)
			setkey(fall.tMin, pop)

			tMax.tMin <- merge(spring.tMax, fall.tMin)

			setkey(tMax.tMin, pop)

			m <- merge(pbs[season=="spring"], tMax.tMin, allow.cartesian=T)[pop!="OUK_13"][dayTrue==T]

			foreach(filter.set=unique(m$class), .errorhandling="remove", .combine="rbind")%do%{

				tm <- lm(beta~spring.tMax + fall.tMin , m[class==filter.set])
				stm <- summary(tm)

				data.table(sp.th=sp.th, fall.th=fall.th,
						r2=stm$r.squared,
						class=filter.set)
			}

		}
	}

	ggplot(data=o, aes(x=sp.th, y=fall.th, fill=r2)) +
	geom_tile() +
	scale_colour_distiller(palette = "Spectral") +
	scale_fill_distiller(palette="Spectral") +
	facet_grid(~class)







	pred.full <- data.table(spring.tMax = rep(c(0:400), each=length(0:400)),
													fall.tMin = rep(c(0:400), length(0:400)))
	pred.full[,pred := predict(tm, newdata=pred.full)]

	ggplot() +
	geom_tile(data=pred.full[fall.tMin>=min(m$fall.tMin, na.rm=T) & fall.tMin<=max(m$fall.tMin, na.rm=T) &
							spring.tMax>=min(m$spring.tMax, na.rm=T) & spring.tMax<=max(m$spring.tMax, na.rm=T)]
			, aes(x=spring.tMax, y=fall.tMin, fill=pred)) +
	geom_point(data=m, aes(x=spring.tMax, y=fall.tMin, fill=beta), colour="black",pch=21, size=5) +
	scale_colour_distiller(palette = "Spectral") +
	scale_fill_distiller(palette="Spectral")









		springCoef.plot <- ggplot(data=wwn[class!="all"],
							aes(x=spring.win, y=fall.win, fill=springCoef)) +
			geom_tile()	+
			facet_grid(~class) +
			scale_fill_distiller(palette = "Spectral")

		fallCoef.plot <- ggplot(data=wwn[class!="all"],
							aes(x=spring.win, y=fall.win, fill=fallCoef)) +
			geom_tile()	+
			facet_grid(~class) +
			scale_fill_distiller(palette = "Spectral")

















	### summarize along sp & dp q, conditional on cline
		registerDoMC(10)
		o <- foreach(i = dim(wins)[1]:100, .combine="rbind", .export="m")%dopar%{

			print(i)
			m.temp <- m[log.sp.q<= wins[i]$start][log.dp.q<= wins[i]$start]

					### no stratification on cline
						m.ag.all <- m.temp[,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="all"),
									  list(pop)]

					### SP: concordant/discordant with cline
						m.ag.clineConcord.sp <- m.temp[sign(sp.coef)==sign(clinal.coef)][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineConcord.sp"),
									  list(pop)]

						m.ag.clineDiscord.sp <- m.temp[sign(sp.coef)!=sign(clinal.coef)][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineDiscord.sp"),
									  list(pop)]

					### SP: concordant/discordant with cline + sig cline
						m.ag.clineConcord.sp.clineSig <- m.temp[sign(sp.coef)==sign(clinal.coef)][log.clinal.q<=wins[i]$start][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineConcord.sp.clineSig"),
									  list(pop)]

						m.ag.clineDiscord.sp.clineSig <- m.temp[sign(sp.coef)!=sign(clinal.coef)][log.clinal.q<=wins[i]$start][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineDiscord.sp.clineSig"),
									  list(pop)]

					### DP: concordant/discordant with cline
						m.ag.clineConcord.dp <- m.temp[sign(dp.coef)==sign(clinal.coef)][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineConcord.dp"),
									  list(pop)]

						m.ag.clineDiscord.dp <- m.temp[sign(dp.coef)!=sign(clinal.coef)][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineDiscord.dp"),
									  list(pop)]

					### DP: concordant/discordant with cline + sig cline
						m.ag.clineConcord.dp.clineSig <- m.temp[sign(dp.coef)==sign(clinal.coef)][log.clinal.q<=wins[i]$start][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineConcord.dp.clineSig"),
									  list(pop)]

						m.ag.clineDiscord.dp.clineSig <- m.temp[sign(dp.coef)!=sign(clinal.coef)][log.clinal.q<=wins[i]$start][,
										list(TT=sum(sign(dp.coef)==sign(sp.coef)),
											n=length(dp.coef),
											log.sp.q.th=wins[i]$start,
											log.dp.q.th=wins[i]$start,
											class="clineDiscord.dp.clineSig"),
									  list(pop)]
					### join
						rbindlist(list(m.ag.all,
									m.ag.clineConcord.sp, m.ag.clineDiscord.sp,
									m.ag.clineConcord.sp.clineSig, m.ag.clineDiscord.sp.clineSig,
									m.ag.clineConcord.dp, m.ag.clineDiscord.dp,
									m.ag.clineConcord.dp.clineSig, m.ag.clineDiscord.dp.clineSig))


		}
		o[,p := dbinom(TT, n, .5)]
		o[,frac := TT/n ]
		o[,q:= p.adjust(p, "fdr")]









	##############################################################################
	### merge with climate data: lookingat weather one month before collection ###
	##############################################################################

		setkey(weather, pop, date)
		setkey(pbs, pop, date)

		ww <- foreach(j=1:45, .combine="rbind")%dopar%{
			print(i)


			#j<-5

			weather.ag <- foreach(i=1:dim(pbs)[1], .combine="rbind", .errorhandling="remove")%do%{

				weather.pop <- weather[J(data.table(pop=pbs[i]$pop,
									 date=seq(from=pbs[i]$date-j, to=pbs[i]$date-(j-j), by=1),
									 key="pop,date"))]

				weather.pop[,list(nFreeze=sum(tmin<50, na.rm=T),
								  delta.temp = tmin[which.max(date)] - tmin[which.min(date)],
								  delta.temp.lm = coef(lm(I(tmin/2 + tmax/2)~date))[2],
								  temp.mu = mean(tmin/2 + tmax/2, na.rm=T),
								  tempmin = min(tmin, na.rm=T),
								  tempmax = max(tmax, na.rm=T),
								  tempmin.mu = mean(tmin, na.rm=T),
								  tempmax.mu = mean(tmax, na.rm=T),
								  date=pbs[i]$date),
							list(pop)]

			}

			setkey(weather.ag, pop, date)
			pbsw <- merge(pbs, weather.ag, all.x=T)

			temp <- pbsw[,list(beta=mean(beta),
							rd=mean(median.rd),
							latitude=mean(latitude),
							longitude=mean(longitude),
							spring.deltaTemp=delta.temp.lm[season=="spring"],
							fall.deltaTemp=delta.temp.lm[season=="fall"],
							temp.spring=mean(temp.mu[season=="spring"], na.rm=T),
							temp.fall=mean(temp.mu[season=="fall"], na.rm=T),
							temp.min.spring=tempmin[season=="spring"],
							temp.max.spring=tempmax[season=="spring"],
							temp.min.fall=tempmin[season=="fall"],
							temp.max.fall=tempmax[season=="fall"],
							temp.max.mu.spring = tempmin.mu[season=="spring"],
							temp.min.mu.fall = tempmin.mu[season=="fall"]),
						   list(pop)]

			t1m <- (lm(beta~temp.min.fall + temp.max.spring , temp))

			t1 <- summary(lm(beta~temp.min.fall + temp.max.spring , temp))

			data.table(win=j,
						temp.min.fall.beta=t1$coef[2,1],
						temp.max.spring.beta=t1$coef[3,1],
						temp.min.fall.se=t1$coef[2,2],
						temp.max.spring.se=t1$coef[3,2],
						temp.min.fall.p=t1$coef[2,4],
						temp.max.spring.p=t1$coef[3,4],
						r2=t1$r.squared)
	}



		layout(matrix(c(1,2,3), nrow=1))
		plot(-log10(temp.min.fall.p)~win, ww)
		plot(-log10(temp.max.spring.p)~win, ww)
		plot(r2~win, ww)




		layout(matrix(c(1,2), nrow=1))
		plot(beta~temp.min.fall, temp)
		plot(beta~temp.max.spring, temp)

		temp[,pred := predict(t1)]


		summary(t1 <- lm(beta~temp.spring + temp.fall , temp))


















#########################################################################################
	### merge with climate data: this is the simplest case of merging by collection date ####
	#########################################################################################

		setkey(weather, pop, date)
		pbsw <- merge(pbs, weather, all.x=T)

		### calculate diff in cumulative degree days
			pbswd <- pbsw[,list(beta=mean(beta),
								delta.cdd=cdd[season=="fall"] - cdd[season=="spring"],
								spring.cdd=cdd[season=="spring"],
								rd=mean(median.rd),
								latitude=mean(latitude),
								longitude=mean(longitude)),
							list(pop)]

		### a few simple plots
			delta.cdd <- ggplot(data=pbswd, aes(y=beta, x=delta.cdd, color=pop)) + geom_point()
			spring.cdd <- ggplot(data=pbswd, aes(y=beta, x=spring.cdd, color=pop)) + geom_point()
			plot_grid(spring.cdd, delta.cdd)

		### maybe sort of something but it makes me queasy.
			summary(lm(beta~spring.cdd + longitude + delta.cdd, pbswd))




	### long to wide
		o.wide <- dcast(o[p<.05], pop + log.sp.q.th + log.dp.q.th ~ class, value.var="frac")

		ggplot(data=o.wide, aes(y=log2(clineConcord.sp/clineDiscord.sp), x=log.sp.q.th, group=pop,
									color=factor(pop, levels=pop.order$pop))) +
		geom_line()




	### best
	o[log.sp.q.th>=-1.999381][log.sp.q.th<=-1.993381][class=="clineConcord.sp"][q<.05]


	m.best <- m[log.dp.q<=-2.5][log.sp.q<=-2.5][clinal.p<.05][sign(clinal.coef)==sign(sp.coef)][sign(sp.coef)==sign(dp.coef)]
	setkey(m.best, chr, pos)

	m.best.uniq <- m.best[!duplicated(m.best)]
	m.best.uniq[,rPos := round(pos, -4)]

	m.best.ag <- m.best.uniq[,list(n=length(pos)), list(rPos, chr)]


	setkey(m.best, chr, pos)
	m.best <- m.best[!duplicated(m.best)]

	hist(log10(expand.grid(dist(m.best.uniq[chr=="3L"]$pos))[,1]), breaks=100)


### plot
	ggplot(data=o[p<.05], aes(x=log.sp.q.th, y=log.dp.q.th, fill=TT/n)) + geom_tile() + facet_wrap(~pop, ncol=4, nrow=5)

	ggplot(data=o[p<.05][n>100], aes(TT/n)) + geom_histogram()






### defunct

o <- foreach(sp.i = 1:dim(wins)[1], .combine="rbind")%do%{
				foreach(dp.i = 1:dim(wins)[1], .combine="rbind")%dopar%{

					print(paste(sp.i, dp.i, sep=" / "))

					m.temp <- m[log.sp.q<= wins[sp.i]$start][log.dp.q<= wins[dp.i]$start]

					m.ag.all <- m.temp[,
									list(TT=sum(sign(dp.coef)==sign(sp.coef)),
										n=length(dp.coef),
										log.sp.q.th=wins[sp.i]$start,
										log.dp.q.th=wins[dp.i]$start,
										class="all"),
								  list(pop)]

					m.ag.clineConcord.dp <- m.temp[sign(dp.coef)==sign(clinal.coef)][,
									list(TT=sum(sign(dp.coef)==sign(sp.coef)),
										n=length(dp.coef),
										log.sp.q.th=wins[sp.i]$start,
										log.dp.q.th=wins[dp.i]$start,
										class="clineConcord.dp"),
								  list(pop)]

					m.ag.clineDiscord.dp <- m.temp[sign(dp.coef)!=sign(clinal.coef)][,
									list(TT=sum(sign(dp.coef)==sign(sp.coef)),
										n=length(dp.coef),
										log.sp.q.th=wins[sp.i]$start,
										log.dp.q.th=wins[dp.i]$start,
										class="clineDiscord.dp"),
								  list(pop)]

					rbindlist(list(m.ag.all, m.ag.clineConcord.dp, m.ag.clineDiscord.dp))

				}
			}

		o[,p := dbinom(TT, n, .5)]

		o.ag <- o[p<.05,list(mu=mean(TT/n), sd=sd(TT/n), n=mean(n)),
				list(mu.q=log.sp.q.th/2 + log.dp.q.th/2, pop, class)]
















	concord <- m[,list(TT=sum(sign(dp.coef[dp.q<.005 & sp.q<.005])==sign(sp.coef[dp.q<.005 & sp.q<.005])),
			n=length(dp.coef[dp.q<.005 & sp.q<.005])),
		list(pop)]


fisher.test(table(m[J("CUA_14")][sp.q<.005 & dp.q<.005]$sp.coef<.0, m[J("CUA_14")][sp.q<.005 & dp.q<.005]$dp.coef<0))


	dat <- fread("/mnt/pricey_1/dropPop/fisher_exactJ_AGA_14.coef_minp2.txt")


	ggplot(data=deltas[J(dat[minp2<1e-6])],
			aes(delta_sf)) + geom_histogram() + facet_wrap(~pop.y)


	a <- deltas[J(dat[minp2<1e-9]),list(TT=sum(abs(delta_sf)>.1, na.rm=T),
									 TF=sum(abs(delta_sf)<.1, na.rm=T)),
								list(pop.x)]

	b <- deltas[J(dat[minp2>1e-9]),list(TT=sum(abs(delta_sf)>.1, na.rm=T),
									 TF=sum(abs(delta_sf)<.1, na.rm=T)),
								list(pop.x)]

	setkey(a, pop.x)
	setkey(b, pop.x)

	ab <- merge(a, b)

	ab[,or := (TT.x/TF.x)/(TT.y/TF.y)]
