#############################
### analysis: concordance ###
#############################

#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### capture
	args = commandArgs(trailingOnly=TRUE)
	message(args[1])
	job <- as.numeric(args[1])

	#job <- 3

### libraries
	library(data.table)
	library(foreach)
	library(doMC)
	registerDoMC(10)
	library(feather)

### set wd on Rivanna
	setwd("/project/berglandlab/alan/drosRTEC")

### load data
	m <- as.data.table(read_feather(path="mnt/pricey_1/dropPop/dropPop.feather"))

### subset down to target region
	setkey(m, chr, pos)

### define inversions
	inv.dt <- rbindlist(list(
  data.table(chr="2L", inv="In2Lt", start=2225744, end=13154180),
  data.table(chr="2R", inv="In2RNS", start=11278659, end=16163839),
  data.table(chr="3L", inv="In3LP", start=3173046, end=16301941),
  data.table(chr="3R", inv="In3RP", start=12257931, end=20569732),
  data.table(chr="3R", inv="In3RMo", start=17232639, end=24857019),
  data.table(chr="3R", inv="In3RK", start=7576289, end=21966092)
))

### define jobs
	jobs <- as.data.table(expand.grid(target=c("breakpoint", "inside", "outside", "chr"),
																		win=c(50000, 100000, 500000, 1000000),
																		inv=inv.dt$inv))

### get sites for job
  getSites <- function(chr.x, start, stop, region, window, inv) {
    #chr.x<-inv.dt[1]$chr; start<-inv.dt[1]$start; stop<-inv.dt[1]$end; region="breakpoint"; window=1e6; inv<-inv.dt[1]$inv
    if(region=="breakpoint") {

      sites <- data.table(chr=chr.x, pos=c( (start-window):(start+window), (stop-window):(stop+window)), roi=T)

    } else if(region=="inside") {
      sites <- data.table(chr=chr.x, pos=(start+window):(stop-window), roi=T)

    } else if(region=="outside") {
      sites <- data.table(chr=chr.x, pos=c(1:(start-window), (stop+window):30e6), roi=T)
    } else if(region=="chr") {
      sites <- data.table(chr=chr.x, pos=c(1:30e6), roi=T)

    }

    setkey(sites, chr, pos)
		sites
	}

	sites <- getSites(chr.x = inv.dt[inv==jobs[job]$inv]$chr,
									  start = inv.dt[inv==jobs[job]$inv]$start,
									  stop = inv.dt[inv==jobs[job]$inv]$end,
									  region=jobs[job]$target,
									  window=jobs[job]$win,
									  inv = jobs[job]$inv)

### subset to target region
	setkey(m, chr, pos)

	m <- m[J(sites)]

### a few transformations
	m[,log.sp.q.all := log10(sp.q.all)]
	m[,log.dp.q.all := log10(dp.q.all)]
	m[,log.clinal.q.all := log10(clinal.q.all)]

	m[,log.sp.q.pass := log10(sp.q.pass)]
	m[,log.dp.q.pass := log10(dp.q.pass)]
	m[,log.clinal.q.pass := log10(clinal.q.pass)]


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

### summarize stats
	o.pass[,p := dbinom(TT, n, .5)]
	o.pass[,frac := TT/n ]
	o.pass[,q:= p.adjust(p, "fdr")]

### tack in run info
	o <- cbind(o.pass, jobs[job])

	### save object
		save(o, file=paste(getwd(), "/chromosome_inversion_predictability/out", job, ".Rdata", sep=""))
