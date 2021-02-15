### libraries
	library(data.table)
  library(doMC)
  registerDoMC(10)
  library(foreach)

### capture
	args = commandArgs(trailingOnly=TRUE)
	message(args[1])
	#args <- 1

### load data
  #load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata")
	load(file="/scratch/aob2x/drosRTEC/dropPop/weather_and_pbs.fix.Rdata")

### set number of perms
  nPerms <- 50
	daysPrior <- c(14,21)


##################################################
### alternative hypotheses model distributions ###
##################################################

  ### generate R2 distribution plot for average temperature

  ### average tmax & tmin
    weather.pop.ag.14 <- weather.pop[daysPrior<=14,list(tmin.spring.ave=mean(tmin[season=="spring"], na.rm=T),
                                        tmax.spring.ave=mean(tmax[season=="spring"], na.rm=T),
                                        tmin.fall.ave=mean(tmin[season=="fall"], na.rm=T),
                                        tmax.fall.ave=mean(tmax[season=="fall"], na.rm=T),
                                        sp.ave=mean(tmin[season=="spring"]/2 + tmax[season=="spring"]/2, na.rm=T),
                                        fall.ave=mean(tmin[season=="fall"]/2 + tmax[season=="fall"]/2, na.rm=T),
                                        max.cdd.spring=max(cdd[season=="spring"], na.rm=T),
                                        max.cdd.fall=max(cdd[season=="fall"], na.rm=T),
																				days.prior=14),
                                  list(pop)]

		weather.pop.ag.21 <- weather.pop[daysPrior<=21,list(tmin.spring.ave=mean(tmin[season=="spring"], na.rm=T),
																				tmax.spring.ave=mean(tmax[season=="spring"], na.rm=T),
																				tmin.fall.ave=mean(tmin[season=="fall"], na.rm=T),
																				tmax.fall.ave=mean(tmax[season=="fall"], na.rm=T),
																				sp.ave=mean(tmin[season=="spring"]/2 + tmax[season=="spring"]/2, na.rm=T),
																				fall.ave=mean(tmin[season=="fall"]/2 + tmax[season=="fall"]/2, na.rm=T),
																				max.cdd.spring=max(cdd[season=="spring"], na.rm=T),
																				max.cdd.fall=max(cdd[season=="fall"], na.rm=T),
																				days.prior=21),
																	list(pop)]

		weather.pop.ag <- rbind(weather.pop.ag.14, weather.pop.ag.21)

  ### experimental variables
    #load(file="/mnt/pricey_1/dropPop/mm.Rdata")
    load(file="/scratch/aob2x/drosRTEC/dropPop/mm.Rdata")


    mm[,PercentSimulans:=as.numeric(as.character(PercentSimulans))]
    mm[,cs:=as.numeric(as.factor(cs))]
    mm[,cm:=as.numeric(as.factor(cm))]

    mm.w <- dcast(mm[,c("pop", "season", "PercentSimulans", "cs", "cm"), with=F],
          pop~season, value.var=c("PercentSimulans", "cs", "cm"))


  ### merge & simplify
    setkey(weather.pop.ag, pop)
    setkey(pbs.small, pop)
    setkey(mm.w, pop)

    pbs.small <- merge(pbs.small, weather.pop.ag, all.x=T)
    pbs.small <- merge(pbs.small, mm.w, all.x=T)

  ### function
    permFun <- function(y, x, xName, nPerms, days.prior) {
      # y=pbs.small[dayTrue==T][pop!="OUK_13"]$beta; x=pbs.small[dayTrue==T][pop!="OUK_13"]$latitude; nPerms=10000; xName="lat"

      ### make matrix of betas for observed & permuted
        #set.seed(1234)
        beta.mat <- cbind(y, replicate(nPerms, sample(y, replace=F)))

      ### run models
        t4 <- summary(lm(beta.mat ~ x))

      ### extract out data
        foreach(i=1:length(t4), .combine="rbind")%dopar%{
          if(i%%1000==0) print(paste(xName, i, nPerms, sep=" / "))
          data.table(x=xName,
                boot.num =i,
                r2=(t4[[i]])$r.squared,
								daysPrior=days.prior)
        }

      }


  ### ecological model sets
    geography <- list("latitude", "longitude", c("latitude", "longitude"))
    cumulDegreeDays <- list("max.cdd.spring", "max.cdd.fall", c("max.cdd.spring", "max.cdd.fall"))
    aveTemp <- list("sp.ave", "fall.ave", c("sp.ave", "fall.ave"))
    aveThermalLimit.univariate <- list("tmin.spring.ave", "tmax.spring.ave","tmin.fall.ave", "tmax.fall.ave")
    aveThermalLimit.bivariate <- list(c("tmin.spring.ave", "tmin.fall.ave"),
                                      c("tmax.spring.ave", "tmax.fall.ave"),
                                      c("tmin.spring.ave", "tmax.fall.ave"),
                                      c("tmax.spring.ave", "tmin.fall.ave"))

    models.ecol <- list(geography=geography,
                   cumulDegreeDays=cumulDegreeDays,
                   aveTemp=aveTemp,
                   aveThermalLimit.univariate=aveThermalLimit.univariate,
                   aveThermalLimit.bivariate=aveThermalLimit.bivariate)

    ### run across different X variables
      o.ecol <- foreach(i=1:length(models.ecol))%do%{
          #i<-1
          o.temp <- foreach(x.i=models.ecol[[i]])%do%{

							foreach(days.i=daysPrior, .combine="rbind")%do%{

	              permFun(y=pbs.small[dayTrue==T][pop!="OUK_13"][days.prior==days.i]$beta,
	                    x=as.matrix(pbs.small[dayTrue==T][pop!="OUK_13"][days.prior==days.i][,x.i,with=F]),
	                    xName=paste(x.i, collapse="_"),
	                    nPerms=nPerms, days.prior=days.i)

							}
          }
          o.temp <- rbindlist(o.temp)
          o.temp[,modelSet:=names(models.ecol)[i]]
          return(o.temp)
      }
      o.ecol <- rbindlist(o.ecol)
      o.ecol[,set:=ifelse(boot.num==1, "Obs", "Perm")]
      o.ecol[,xnl:=gsub("_", "\n", x)]

  ### experimental model sets
    substrate <- list("cs_spring", "cs_fall", c("cs_spring", "cs_fall"))
    simulans <- list("PercentSimulans_spring", "PercentSimulans_fall", c("PercentSimulans_spring", "PercentSimulans_fall"))
    method <- list("cm_spring", "cm_fall", c("cm_spring", "cm_fall"))


    models.exp <- list(substrate=substrate,
                   simulans=simulans,
                   method=method)

    ### run across different X variables
      o.exp <- foreach(i=1:length(models.exp))%do%{
          #i<-1
          o.temp <- foreach(x.i=models.exp[[i]])%do%{
							foreach(days.i=daysPrior, .combine="rbind")%do%{
	              permFun(y=pbs.small[dayTrue==T][pop!="OUK_13"][days.prior==days.i]$beta,
	                    x=as.matrix(pbs.small[dayTrue==T][pop!="OUK_13"][days.prior==days.i][,x.i,with=F]),
	                    xName=paste(x.i, collapse="_"),
	                    nPerms=nPerms, days.prior=days.i)
							}
          }
          o.temp <- rbindlist(o.temp)
          o.temp[,modelSet:=names(models.exp)[i]]
          return(o.temp)
      }
      o.exp <- rbindlist(o.exp)
      o.exp[,set:=ifelse(boot.num==1, "Obs", "Perm")]
      o.exp[,xnl:=gsub("_", "\n", x)]



  ### merge
    o.altModels <- rbind(o.ecol, o.exp)
    save(o.altModels,
          file=paste("/scratch/aob2x/drosRTEC/dropPop/permFix/perm",
              args[1],
              ".altModels_out.Rdata",
              sep=""))
