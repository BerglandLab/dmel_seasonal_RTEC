### libraries
	library(data.table)
	library(ggplot2)
  library(cowplot)
  library(doMC)
  registerDoMC(20)
  library(speedglm)

### capture
	args = commandArgs(trailingOnly=TRUE)
	message(args[1])
	#args <-14

### load data
  load(file="/mnt/pricey_1/dropPop/weather_and_pbs.Rdata")

### function
  makePlot <- function(dat, mod) {
    tmp <- dat[modelSet==mod]

    pr <- tmp[,list(pr=mean(r2[set=="Obs"]>r2[set=="Perm"]), set="Obs"), list(x)]
    setkey(tmp, set, x)
    setkey(pr, set, x)
    tmp <- merge(tmp, pr, all.x=T)

    ggplot() +
    geom_density(data=tmp[set=="Perm"], aes(r2), fill="lightblue", n=10000) +
    geom_vline(data=tmp[set=="Obs"], aes(xintercept=r2), color="red", size=1) +
    geom_text(data=tmp[set=="Obs"], aes(x=r2, y=3, label=pr), hjust=-.5) +
    facet_grid(modelSet~x) +
    ylim(0,1)
  }

### set number of perms
  nPerms <- 5000
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
            set.seed(1234)
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
        file=paste("/mnt/pricey_1/dropPop/thermalLimitModels_univariate_out_daysPrior.Rdata")



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
            set.seed(1234)
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
      file="/mnt/pricey_1/dropPop/thermalLimitModels_bivariate_out_daysPrior.Rdata")


	save(o.tlm.uni.ag, o.tlm.multi.ag,
		file="/mnt/pricey_1/dropPop/thermalLimitModels_bivariate_out_daysPrior.aggregate.Rdata")

##################################################
### alternative hypotheses model distributions ###
##################################################

  ### generate R2 distribution plot for average temperature

  ### average tmax & tmin
    weather.pop.ag <- weather.pop[,list(tmin.spring.ave=mean(tmin[season=="spring"], na.rm=T),
                                        tmax.spring.ave=mean(tmax[season=="spring"], na.rm=T),
                                        tmin.fall.ave=mean(tmin[season=="fall"], na.rm=T),
                                        tmax.fall.ave=mean(tmax[season=="fall"], na.rm=T),
                                        sp.ave=mean(tmin[season=="spring"]/2 + tmax[season=="spring"]/2, na.rm=T),
                                        fall.ave=mean(tmin[season=="fall"]/2 + tmax[season=="fall"]/2, na.rm=T),
                                        max.cdd.spring=max(cdd[season=="spring"], na.rm=T),
                                        max.cdd.fall=max(cdd[season=="fall"], na.rm=T)),
                                  list(pop)]

  ### experimental variables
    load(file="/mnt/pricey_1/dropPop/mm.Rdata")
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
    permFun <- function(y, x, xName, nPerms) {
      # y=pbs.small[dayTrue==T][pop!="OUK_13"]$beta; x=pbs.small[dayTrue==T][pop!="OUK_13"]$latitude; nPerms=10000; xName="lat"

      ### make matrix of betas for observed & permuted
        set.seed(1234)
        beta.mat <- cbind(y, replicate(nPerms, sample(y, replace=F)))

      ### run models
        t4 <- summary(lm(beta.mat ~ x))

      ### extract out data
        foreach(i=1:length(t4), .combine="rbind")%dopar%{
          if(i%%1000==0) print(paste(xName, i, nPerms, sep=" / "))
          data.table(x=xName,
                boot.num =i,
                r2=(t4[[i]])$r.squared)
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
              permFun(y=pbs.small[dayTrue==T][pop!="OUK_13"]$beta,
                    x=as.matrix(pbs.small[dayTrue==T][pop!="OUK_13"][,x.i,with=F]),
                    xName=paste(x.i, collapse="_"),
                    nPerms=nPerms)
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
              permFun(y=pbs.small[dayTrue==T][pop!="OUK_13"]$beta,
                    x=as.matrix(pbs.small[dayTrue==T][pop!="OUK_13"][,x.i,with=F]),
                    xName=paste(x.i, collapse="_"),
                    nPerms=nPerms)
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
    save(o.altModels, file="/mnt/pricey_1/dropPop/altModels_out.Rdata")
