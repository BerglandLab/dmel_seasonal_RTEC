# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
	library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### generate null distribution of R2
  tlm.bi.fn <- list.files("/scratch/aob2x/drosRTEC/dropPop/permFix", pattern="thermalLimitModels_bivariate_out_daysPrior", full.names=T)

  o.tlm.multi.ag.perm <- foreach(fn.i=tlm.bi.fn)%dopar%{
    message(fn.i)
    #fn.i <- tlm.bi.fn[1]
    load(fn.i)
    o.tlm.multi.ag[,bootSet:=as.numeric(gsub("perm", "", tstrsplit(last(tstrsplit(fn.i, "/")), "\\.")[[1]]))]
    setnames(o.tlm.multi.ag, "boot.num", "sub.boot.num")

    if(o.tlm.multi.ag$bootSet[1]>1) {
      o.tlm.multi.ag <- o.tlm.multi.ag[set!="Obs"]
    }


    o.tlm.multi.ag[,modelSet:="tlm_bivariate"]

    return(o.tlm.multi.ag)
  }
  o.tlm.multi.ag.perm <- rbindlist(o.tlm.multi.ag.perm)

  o.tlm.multi.ag.perm[sub.boot.num==1][x=="m4"]
  o.tlm.multi.ag.perm[sub.boot.num==2][x=="m3"]
  table(o.tlm.multi.ag.perm$x)


### assign empirical probabilities to observed 14-day "m4"
	fn.i <- tlm.bi.fn[1]
	load(fn.i)


	obs.14 <- o.tlm.multi[set=="Obs"][config=="m4"][daysPrior==14]
	exp.14 <- o.tlm.multi.ag.perm[x=="m4"][set!="Obs"][daysPrior==14]$r2


	thresh.out.obs.14 <- foreach(i=1:dim(obs.14)[1], .combine="rbind")%do%{
		message(i)
		obs.14[i,list(r2=r2, daysPrior=daysPrior, pr=mean(r2>exp.14)), list(sp.th, fall.th)]
	}




	obs.21 <- o.tlm.multi[set=="Obs"][config=="m4"][daysPrior==21]
	exp.21 <- o.tlm.multi.ag.perm[x=="m4"][set!="Obs"][daysPrior==21]$r2


	thresh.out.obs.21 <- foreach(i=1:dim(obs.21)[1], .combine="rbind")%do%{
		message(i)
		obs.21[i,list(r2=r2, daysPrior=daysPrior, pr=mean(r2>exp.21)), list(sp.th, fall.th)]
	}

	thresh.out.obs <- rbind(thresh.out.obs.14, thresh.out.obs.21)

	save(thresh.out.obs, file="~/paramContourData.fix.Rdata")
