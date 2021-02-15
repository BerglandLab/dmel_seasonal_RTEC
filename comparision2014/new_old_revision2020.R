

### old versus new sig sites
	### libraries
		library(data.table)
		library(foreach)
		library(doMC)
		registerDoMC(20)

	### functions
		matchControl <- function(targetId, param, paramId=NULL, exclude=NULL, nBoot=5, method=2) {
			if(method==1) {
				len <- length(targetId)
				if(nBoot>0) {
					out <- cbind(targetId, do.call("rbind", foreach(i=1:len, .export=c("param", "targetId", "nBoot", "len", "exclude"))%dopar% {
						print(paste(i, " / ", len, sep=""))
				#			potentialIds <- param[J(param[param$id==targetId[i], -dim(param)[2], with=FALSE])]$id  ## try to make this line vectorized
						potentialIds <- param[J(param[id==targetId[i], -dim(param)[2], with=FALSE]), nomatch=0]$id  ## try to make this line vectorized

						potentialIds <- potentialIds[!potentialIds%in%targetId]

						if(!is.null(exclude)) {
							potentialIds <- potentialIds[!exclude[potentialIds]]
						}

						ret <- NULL
						if(length(potentialIds)>0) {
							ret <- as.numeric(sample(as.character(potentialIds), nBoot, replace=TRUE))
						} else {
							ret <- rep(NA, nBoot)
						}
						ret
					}))
				} else {
					out <- matrix(targetId, ncol=1)
				}
			} else if(method==2) {
				origParamKey <- key(param)
				setkey(param, id)
				targetDt <- param[J(targetId)]

				setkeyv(param, origParamKey)
				setkeyv(targetDt, origParamKey)

				setnames(param, "id", "id.control")
				setnames(targetDt, "id", "id.target")

				ret <- param[targetDt, nomatch=NA, allow.cartesian=T]


				ret <- ret[!ret$id.control%in%targetId]
				ret$id[is.na(ret$id.control)]<-0

				sampRet <- ret[,as.numeric(sample(as.character(id.control), nBoot, replace=T)),id.target]

				ret <- matrix(sampRet$V1, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)


				out <- cbind(matrix(sampRet$id.target, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)[,1], ret)
				ret[ret==0] <- NA

			} else if(method==3) {
				targetDt <- paramId[targetId]

				setkeyv(targetDt, key(param))

				ret <- param[J(targetDt), nomatch=NA, allow.cartesian=TRUE]
				#ret <- param[J(targetDt), nomatch=NA]


				ret <- ret[!ret$id%in%targetId]

				if(!is.null(exclude)) {
					ret <- ret[!ret$id%in%exclude]
				}


				if(dim(ret)[1]>0) {
					ret$id[is.na(ret$id)]<-0
					sampRet <- ret[,as.numeric(sample(id,nBoot, replace=T)),id.1]

					ret <- matrix(sampRet$V1, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)


					out <- cbind(matrix(sampRet$id.1, nrow=dim(sampRet)[1]/nBoot, ncol=nBoot, byrow=T)[,1], ret)
					out[out==0] <- NA
				} else {
					out <- matrix(c(targetId, rep(NA, nBoot)), nrow=1)
				}
			}

			out
		}

	### Comparision to Core20, unswitched set; no PA
		### old 6D data
			load("~/6d_data.Rdata")
			p <- as.data.table(p)
			setkey(p, chr, pos)

			co <- as.data.table(co)
			setkey(co, chr, pos)

		### load SNPs to use
			snps <- fread("/mnt/pricey_1/dropPop/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
			setnames(snps, names(snps), c("chr", "pos"))
			setkey(snps, chr, pos)

			snps[,filter:="pass"]


		### do analysis with noPA data
			### load full 20 dataset
				#full20 <- fread("/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear.f_s.glm")
				#full20 <- fread("/mnt/pricey_1/dropPop/mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm")
				full20 <- fread("/mnt/pricey_1/dropPop/mel_all_paired20_2sample_noPA_caF_popyear_polymorphic.f_s.glm")

				setnames(full20, "chrom", "chr")
				setkey(full20, chr, pos)

				full20 <- merge(full20, snps, all.x=T)
				full20[,full.q.all := rank(seas.p)/(length(seas.p)+1)]
				full20[filter=="pass",full.q.pass := rank(seas.p)/(length(seas.p)+1)]
				setkey(full20, chr, pos)

			### merge
				m <- merge(p, full20)
				m <- m[!is.na(full.q.pass)][!is.na(sfsfsfX.q)]
				m <- m[!is.na(sfsfsfX.q)]

			### merge with co
				m <- merge(m, co)

			### set id
				m[,id:=c(1:dim(m)[1])]

			### generate matched control sites
				controlSet <- matchControl(targetId=m[full.q.pass<=.01]$id,
											param=data.table(chr=m$chr,
															 H=round(m$f.hat*(1-m$f.hat)*2, 2),
															 rr=round(m$rec),
															 id=m$id,
															 key="chr,H,rr"),
											nBoot=10000)

				nUniq <- apply(controlSet[,-1], 1, function(x) length(unique(x)))

				controlSet.o <- controlSet
				controlSet <- controlSet[nUniq>100,]

### define inversions

    ### are any of the overlaps in breakpoints?
      inv.dt <- rbindlist(list(
        data.table(chr="2L", inv="In2Lt", start=2225744, end=13154180),
        data.table(chr="2R", inv="In2RNS", start=11278659, end=16163839),
        data.table(chr="3L", inv="In3LP", start=3173046, end=16301941),
        data.table(chr="3R", inv="In3RP", start=12257931, end=20569732),
        data.table(chr="3R", inv="In3RMo", start=17232639, end=24857019),
        data.table(chr="3R", inv="In3RK", start=7576289, end=21966092)
      ))

      win <- 5e5


      bkpts <- foreach(i=1:dim(inv.dt)[1])%do%{
        data.table(chr=inv.dt[i]$chr, pos=(inv.dt[i]$start-win):(inv.dt[i]$end+win))
      }
      bkpts <- rbindlist(bkpts)
      setkey(bkpts, chr, pos)
      bkpts <- bkpts[!duplicated(bkpts)]

      bkpts <- merge(bkpts, m[,c("chr", "pos"), with=F])
      setkey(bkpts)

    ### save
      save(bkpts, m, controlSet, matchControl, file="~/seas2014_inv.Rdata")


### to rivanna
  scp ~/seas2014_inv.Rdata aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/.

#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

  ### libraries
    library(data.table)
    library(foreach)
    library(doMC)
    registerDoMC(10)

  ### load data
    load("/scratch/aob2x/seas2014_inv.Rdata")

		### do test
			setkey(m, id)

			focalSet <- m[J(controlSet[,1])]
			TT <- sum(focalSet$sfsfsfX.q<=.3)
			TF <- sum(focalSet$sfsfsfX.q>.3)

      setkey(focalSet, chr, pos)
      TT.inv <- dim(focalSet[J(bkpts)][sfsfsfX.q<=.3])[1]
      TF.inv <- dim(focalSet[!J(bkpts)][sfsfsfX.q<=.3])[1]

      #dim(controlSet)[2]

			o <- foreach(i=2:100, .combine="rbind")%dopar%{
				print(paste(i, dim(controlSet)[2], " / "))
				testSet <- m[J(controlSet[,i])]

				FT <- sum(testSet$sfsfsfX.q<=.3)
				FF <- sum(testSet$sfsfsfX.q>.3)

        setkey(testSet, chr, pos)
        FT.inv <- dim(testSet[J(bkpts)][sfsfsfX.q<=.3])[1]
        FF.inv <- dim(testSet[!J(bkpts)][sfsfsfX.q<=.3])[1]


				data.table(rep=i,
							or=(TT/TF)/(FT/FF),
              or.inv=(TT.inv/TF.inv)/(FT.inv/FF.inv),
              TT=TT,
              TF=TF,
              FT=FT,
              FF=FF,
              TT.inv=TT.inv,
              TF.inv=TF.inv,
              FT.inv=FT.inv,
              FF.inv=FF.inv)
      }
      o.ag <- o[,list(mu=mean(log2(or)),
                      var=var(log2(or)),
                      sd=sd(log2(or)),
                      pr=mean(or>1))]

      o.ag <- o[,list(mu=mean(log2(or.inv)),
                      var=var(log2(or.inv)),
                      sd=sd(log2(or.inv)),
                      pr=mean(or.inv>1))]
