# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
	library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)

### collect "Alternative Models" (geography, collection methods, thermal averages)
  alt.fn <- list.files("/scratch/aob2x/drosRTEC/dropPop/permFix", pattern="altModels", full.names=T)

  o.altModels.perm <- foreach(fn.i=alt.fn)%dopar%{
    message(fn.i)
    #fn.i <- alt.fn[1]
    load(fn.i)
    o.altModels[,bootSet:=as.numeric(gsub("perm", "", tstrsplit(last(tstrsplit(fn.i, "/")), "\\.")[[1]]))]
    setnames(o.altModels, "boot.num", "sub.boot.num")

    if(o.altModels$bootSet[1]>1) {
      o.altModels <- o.altModels[set!="Obs"]
    }
    return(o.altModels)
  }
  o.altModels.perm <- rbindlist(o.altModels.perm)

  o.altModels.perm[sub.boot.num==1][x=="latitude"]
  o.altModels.perm[sub.boot.num==2][x=="latitude"]
  table(o.altModels.perm$x, o.altModels.perm$daysPrior)


### Collect "Univariate Thermal Limit Models"
  tlm.uni.fn <- list.files("/scratch/aob2x/drosRTEC/dropPop/permFix", pattern="thermalLimitModels_univariate_out_daysPrior", full.names=T)

  o.tlm.uni.ag.perm <- foreach(fn.i=tlm.uni.fn)%dopar%{
    message(fn.i)
    #fn.i <- tlm.uni.fn[1]
    load(fn.i)
    o.tlm.uni.ag[,bootSet:=as.numeric(gsub("perm", "", tstrsplit(last(tstrsplit(fn.i, "/")), "\\.")[[1]]))]
    setnames(o.tlm.uni.ag, "boot.num", "sub.boot.num")

    if(o.tlm.uni.ag$bootSet[1]>1) {
      o.tlm.uni.ag <- o.tlm.uni.ag[set!="Obs"]
    }


    o.tlm.uni.ag[,modelSet:="tlm_univariate"]


    return(o.tlm.uni.ag)
  }
  o.tlm.uni.ag.perm <- rbindlist(o.tlm.uni.ag.perm)

  o.tlm.uni.ag.perm[sub.boot.num==1][x=="m3"]
  o.tlm.uni.ag.perm[sub.boot.num==2][x=="m3"]
  table(o.tlm.uni.ag.perm$x, o.tlm.uni.ag.perm$daysPrior)



### Collect "Bivariate Thermal Limit Models"
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
  table(o.tlm.multi.ag.perm$x, o.tlm.multi.ag.perm$daysPrior)


### save
  save(o.altModels.perm, o.tlm.uni.ag.perm, o.tlm.multi.ag.perm, file="~/mods.perm.fix.Rdata")
