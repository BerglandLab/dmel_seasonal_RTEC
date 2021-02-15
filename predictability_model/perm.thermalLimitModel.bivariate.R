### libraries
	library(data.table)
  library(doMC)
  registerDoMC(10)
  library(foreach)

### capture
	args = commandArgs(trailingOnly=TRUE)
	message(args[1])
	#args <-14

### load data
  #load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata")
	load(file="/scratch/aob2x/drosRTEC/dropPop/weather_and_pbs.fix.Rdata")

### set number of perms
  nPerms <- 50
	daysPrior <- c(14,21)

####################################################
### thermal limit model distributions: bivariate ###
####################################################

  o.tlm.multi <- foreach(daysPrior.i=daysPrior, .combine="rbind")%do%{
		foreach(sp.th=seq(from=-50, to=400, by=10), .combine="rbind", .errorhandling="remove")%dopar%{
    foreach(fall.th=seq(from=-50, to=400, by=10), .combine="rbind", .errorhandling="remove")%do%{

      #sp.th <- 330
      #fall.th <- 20

      print(paste(sp.th, fall.th, sep= " / "))

      weather.pop.ag.m1 <- weather.pop[daysPrior<=daysPrior.i,list(x1.sp.tmin=mean(tmin[season=="spring"] <= sp.th, na.rm=T),
                                                                x2.fall.tmin=mean(tmin[season=="fall"] <= fall.th, na.rm=T)),
                                                           list(pop)]

      weather.pop.ag.m2 <- weather.pop[daysPrior<=daysPrior.i,list(x1.sp.tmax=mean(tmax[season=="spring"] >= sp.th, na.rm=T),
                                                                  x2.fall.tmax=mean(tmax[season=="fall"] >= fall.th, na.rm=T)),
                                                            list(pop)]

      weather.pop.ag.m3 <- weather.pop[daysPrior<=daysPrior.i,list(x1.sp.tmin=mean(tmin[season=="spring"] <= sp.th, na.rm=T),
                                                                x2.fall.tmax=mean(tmax[season=="fall"] >= fall.th, na.rm=T)),
                                                           list(pop)]

      weather.pop.ag.m4 <- weather.pop[daysPrior<=daysPrior.i,list(x1.sp.tmax=mean(tmax[season=="spring"] >= sp.th, na.rm=T),
                                                                  x2.fall.tmin=mean(tmin[season=="fall"] <= fall.th, na.rm=T)),
                                                            list(pop)]


      weather.pop.list <- list(m1=weather.pop.ag.m1, m2=weather.pop.ag.m2, m3=weather.pop.ag.m3, m4=weather.pop.ag.m4)

      o.tmp <- foreach(k=1:length(weather.pop.list))%do%{
          #k<-4
          weather.pop.ag <- weather.pop.list[[k]]

          setkey(weather.pop.ag, pop)

          m <- merge(pbs.small, weather.pop.ag)[pop!="OUK_13"][dayTrue==T]
          orig.x <- names(m)[grepl("x1|x2", names(m))]
          setnames(m, names(m)[grepl("x1|x2", names(m))], c("X1", "X2"))

          ### make matrix of betas for observed & permuted
            #set.seed(1234)
            beta.mat <- cbind(m$beta, replicate(nPerms, sample(m$beta, replace=F)))

          ### run models
            #m4 <- lm(beta.mat ~ tmax + tmin, m)
            #p4 <- predict(m4)

          ### orig
            t4 <- summary(lm(beta.mat ~ X1 + X2, m))

          ### extract out data
            o.tmp <- foreach(i=1:length(t4), .combine="rbind")%do%{

              data.table(sp.th=sp.th, fall.th=fall.th,
                    boot.num =i, config=names(weather.pop.list)[k],
                    r2=(t4[[i]])$r.squared,
										daysPrior=daysPrior.i)
            }

          return(o.tmp)
        }

      o.tmp <- rbindlist(o.tmp)
      return(o.tmp)
  	}
  	}
	}

  o.tlm.multi[,set:=ifelse(boot.num==1, "Obs", "Perm")]
  o.tlm.multi.ag <- o.tlm.multi[,list(r2=max(r2), sp.th=sp.th[which.max(r2)], fall.th=fall.th[which.max(r2)], modelSet="tlm"), list(boot.num, set, x=config, daysPrior)]


  save(o.tlm.multi, o.tlm.multi.ag,
      file=paste("/scratch/aob2x/drosRTEC/dropPop/permFix/perm",
                args[1],
                ".thermalLimitModels_bivariate_out_daysPrior.Rdata",
                sep=""))

	#save(o.tlm.uni.ag, o.tlm.multi.ag,
  #  file=paste("/scratch/aob2x/drosRTEC/dropPop/perm/perm",
  #            args[1],
  #            ".thermalLimitModels_bivariate_out_daysPrior.aggregte.Rdata",
  #            sep=""))
#
