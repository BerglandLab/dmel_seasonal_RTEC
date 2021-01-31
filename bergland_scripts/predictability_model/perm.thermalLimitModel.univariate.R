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


#####################################################
### thermal limit model distributions: univariate ###
#####################################################

  o.tlm.uni <- foreach(daysPrior.i=daysPrior, .combine="rbind")%do%{
		foreach(th=seq(from=-50, to=400, by=10), .combine="rbind", .errorhandling="remove")%dopar%{
      #sp.th <- 330
      #fall.th <- 20

      message(paste(daysPrior.i, th, sep= " / "))

      weather.pop.ag.m1 <- weather.pop[daysPrior<=daysPrior.i,list(x1.sp.tmin=mean(tmin[season=="spring"] <= th, na.rm=T)),
                                                           list(pop)]

      weather.pop.ag.m2 <- weather.pop[daysPrior<=daysPrior.i,list(x1.sp.tmax=mean(tmax[season=="spring"] >= th, na.rm=T)),
                                                            list(pop)]

      weather.pop.ag.m3 <- weather.pop[daysPrior<=daysPrior.i,list(x1.fall.tmin=mean(tmin[season=="fall"] <= th, na.rm=T)),
                                                           list(pop)]

      weather.pop.ag.m4 <- weather.pop[daysPrior<=daysPrior.i,list(x1.fall.tmax=mean(tmax[season=="fall"] >= th, na.rm=T)),
                                                            list(pop)]

      weather.pop.list <- list(m1=weather.pop.ag.m1, m2=weather.pop.ag.m2, m3=weather.pop.ag.m3, m4=weather.pop.ag.m4)

      o.tmp <- foreach(k=1:length(weather.pop.list))%do%{
          #k<-4
          weather.pop.ag <- weather.pop.list[[k]]

          setkey(weather.pop.ag, pop)

          m <- merge(pbs.small, weather.pop.ag)[pop!="OUK_13"][dayTrue==T]
          orig.x <- names(m)[grepl("x1", names(m))]
          setnames(m, names(m)[grepl("x1", names(m))], c("X1"))

          ### make matrix of betas for observed & permuted
            #set.seed(1234)
            beta.mat <- cbind(m$beta, replicate(nPerms, sample(m$beta, replace=F)))

          ### run models
            #m4 <- lm(beta.mat ~ tmax + tmin, m)
            #p4 <- predict(m4)

          ### orig
            t4 <- summary(lm(beta.mat ~ X1, m))

          ### extract out data
            o.tmp <- foreach(i=1:length(t4), .combine="rbind")%do%{

              data.table(th=th,
                    boot.num =i, config=names(weather.pop.list)[k],
                    r2=(t4[[i]])$r.squared,
										daysPrior=daysPrior.i)
            }

          o.tmp
        }

      o.tmp <- rbindlist(o.tmp)
      o.tmp
  	}
	}
  o.tlm.uni[,set:=ifelse(boot.num==1, "Obs", "Perm")]
  o.tlm.uni.ag <- o.tlm.uni[r2>0,list(r2=max(r2), th=th[which.max(r2)], modelSet="tlm.uni"),
                                list(boot.num, set, x=config, daysPrior)]

  #makePlot(dat=o.tlm.uni.ag, mod="tlm.uni")
  save(o.tlm.uni, o.tlm.uni.ag,
        file=paste("/scratch/aob2x/drosRTEC/dropPop/permFix/perm",
                  args[1],
                  ".thermalLimitModels_univariate_out_daysPrior.Rdata",
                  sep=""))
