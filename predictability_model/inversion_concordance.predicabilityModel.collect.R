
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R



### libraries
	library(data.table)
	library(foreach)

### set wd on Rivanna
	setwd("/project/berglandlab/alan/drosRTEC")

### load data
	fn <- list.files(paste(getwd(), "/chromosome_inversion_predictability", sep=""), "out", full.names=T)

	o <- foreach(i=fn)%do%{
		#i<-fn[1]
		message(i)
		load(i)
		return(o)
	}
	o <- rbindlist(o)

### save
	save(o, file="~/loo_inv.Rdata")
