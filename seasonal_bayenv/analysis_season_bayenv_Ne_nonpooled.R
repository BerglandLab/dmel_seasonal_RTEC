## Jan 2020

#library("nlme")

# The number of reads, adjusted for number of flies collected
load("../../data/mel_freqdp_042016_Ne_fixed.Rdata")

# Core 20 populations (N=40)
pops24 = read.table("../../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]

popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

####### Change this to target analysis ########################
focal = which( PY %in% pops24 & (S=="s" | S=="f") )
################################################################

freqfocal = freq[, focal]
dpfocal = dp[, focal] ### only looking at fall to frost
bothA=cbind(info[,1:2], freqfocal, dpfocal)
P1 = P[focal]
Y1 = Y[focal]
R1 = R[focal]
S1 = S[focal]
PY1 = PY[focal]
popID1 = popID[focal]
nsamples = length(S1)
filter=read.table("../../data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
both=merge(bothA, filter, by=c(1,2))
bothAuto = subset(both, X.CHROM != "X")
bothAuto[is.na(bothAuto)] = 0

# calculating alt and ref reads, putting into bayenv format
#SNPFILE The input file for the data is identical in format, but now the file must contain the count data for a single SNP. The allele counts must be in the same population order that they appeared the in the covariance matrix, i.e. the order that they appeared in SNPSFILE. An example input file is given in the the file rs316.
# Format:
#   Each row is 1 snp
#   First half of columns is number of alternate reads
#   Second half of columns is number of refernce reads
dpfocal = round(bothAuto[,43:82])
altfocal = round(bothAuto[,3:42]*dpfocal)
reffocal =  dpfocal-altfocal
bothaltref=cbind(bothAuto[,1:2], altfocal, reffocal)
write.table(bothaltref, file="bothaltref_Ne.txt", col.names = F, row.names = F, quote=F, sep="\t")



## Creating split files- Ne (done on lustre): 882 files (881 2000snp files + 1 not full).
bothaltref = read.table(file="bothaltref_Ne.txt", header=F, stringsAsFactors = F)
lines = seq(from=1, to=nrow(bothaltref), by=2000)
lines = c(lines, nrow(bothaltref))
for (i in 1:(length(lines)-1)){
  focaldata =  bothaltref[lines[i]:(lines[i+1]-1),]
  focalinfo = focaldata[,1:2]
  focalsnps = matrix(unlist(t(focaldata[,3:ncol(focaldata)])), ncol=40, byrow=T)
  write.table(focalsnps, file=paste("intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.", i, sep=""), col.names = F, row.names = F, quote=F, sep="\t")
  write.table(focalinfo[,1:2], file=paste("intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.", i, sep=""), col.names = F, row.names = F, quote=F, sep="\t")
}



# creating file for matrix estimation (10K snps)
non0 = apply(bothaltref[,3:(length(focal)+2)], MARGIN=1, FUN=sum)
non0ref = apply(bothaltref[,(length(focal)+3):(length(focal)*2+2)], MARGIN=1, FUN=sum)
bothaltrefNon0 = na.omit(bothaltref[non0>0 & non0ref>0, ])
bothaltref10K = bothaltrefNon0[sample(1:nrow(bothaltrefNon0), size=10000, replace=F), ]
snpsfile10K = matrix(unlist(t(bothaltref10K[,3:ncol(bothaltref10K)])), ncol=length(focal), byrow=T)
write.table(snpsfile10K, file="seasonal_paired20_snpsfile_Ne_Jan2020_sample10K.txt", col.names = F, row.names = F, quote=F, sep="\t")



#ENVIRONFILE is the file of environmental variables that you wish to estimate correlation statistics for. Each environmental variable should be standardized, i.e. subtract the mean and then divided through by the standard deviation of the variable across populations. 
# IS THAT TRUE FOR spring/fall? I wouldn't think so... (talk to Alan about this)
#Each line of the file is an environ- mental variable across populations, values should be tab separated. The variables should be listed in the same population order as they appear in the allele count files.
# We will make one with just one line for spring/fall. Spring=0, Fall=1
environ = S1
environ[S1=="f"] = 1
environ[S1=="s"] = -1
write(environ, file="seasonal_paired20_environ_Jan2020_springfall.txt", sep="\t", ncolumns = 40)




### June 2020 ###
## Making permuted environfiles
# switch within the popyear
# create table of randomly switching 1000 times. 1 is original, -1 is switched
permat = matrix(nrow=1000,ncol=20)
for (i in 1:1000){
  permat[i,] = sample(c(-1,1), size=20, replace=TRUE)  
}
sum(duplicated(permat)) # 1 is duplicated
permat2 = permat[!(duplicated(permat)),]
write.table(permat2, file="springfall_permuations_1000_June2020.txt")

PYswitch_list2 = list()
permat_switched2 = matrix(nrow=nrow(permat2),ncol=40)
for (j in 1:nrow(permat2)){ # for each of 1000 permutations
  focal_perm = permat2[j,]
  environ_new = environ
  toswitch = which(focal_perm==-1)
  PYswitch_list2[[j]] = unique(PY1)[toswitch]
  for (k in 1:length(toswitch)){
    focalpopyear = PYswitch_list2[[j]][k]
    environ_new[PY1==focalpopyear & S1=="s"] = 1
    environ_new[PY1==focalpopyear & S1=="f"] = -1
  }
  permat_switched2[j,] = environ_new
  write(environ_new, file=paste("permutations_June2020/perm", j,"_seasonal_paired20_environ_June2020_springfall.txt",sep=""), sep="\t", ncolumns = 40)
}
apply(permat_switched2, MARGIN=2, FUN=function(X) sum(as.numeric(X)==-1) )
hist(apply(permat_switched2, MARGIN=2, FUN=function(X) sum(as.numeric(X)==-1) ), breaks=100)
save(PYswitch_list2, permat_switched2, file="permutations_June2020/PYswitch_list_springfall_permuations_1000_June2020.Rdata")

for (j in 1:100){ # for each of 100 permutations
  focal_perm = permat2[j,]
  environ_new = environ
  toswitch = which(focal_perm==-1)
  PYswitch_list2[[j]] = unique(PY1)[toswitch]
  for (k in 1:length(toswitch)){
    focalpopyear = PYswitch_list2[[j]][k]
    environ_new[PY1==focalpopyear & S1=="s"] = 1
    environ_new[PY1==focalpopyear & S1=="f"] = -1
  }
  permat_switched2[j,] = environ_new
  write(environ_new, file=paste("permutations_June2020/perm", j,"_seasonal_paired20_environ_June2020_springfall.txt",sep=""), sep="\t", ncolumns = 40)
}

