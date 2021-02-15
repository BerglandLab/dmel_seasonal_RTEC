#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
	library(data.table)
	#library(ggplot2)
  #library(cowplot)
  #library(doMC)
  #registerDoMC(20)
  #library(speedglm)

### capture
	args = commandArgs(trailingOnly=TRUE)
	message(args[1])
	#args <-14

### load data
  #load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata")
	load(file="/scratch/aob2x/drosRTEC/dropPop/weather_and_pbs.fix.Rdata")

  m <- merge(pbs.small, weather.pop)[pop!="OUK_13"][dayTrue==T]

  m.ag <- m[,list(beta=mean(beta)), list(pop)]
  m[,pop.f:=factor(pop, levels=m.ag[order(beta)]$pop)]
  save(m, file="~/m.fix.Rdata")

  #scp aob2x@rivanna.hpc.virginia.edu:~/m.fix.Rdata ~/.
