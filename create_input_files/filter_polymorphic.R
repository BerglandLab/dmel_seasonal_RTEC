# Filtering for SNPs that are polymorphic in ALL 40 samples

load("data/mel_freqdp_042016_fixed.Rdata")
pops24 = read.table("data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]
popID = popinfo[,1]
S = popinfo[,5]
PY = popinfo[,6]
focal = which( PY %in% pops24 & (S=="s" | S=="f") )
freqfocal = freq[, focal]
dpfocal = dp[, focal] ### only looking at fall to frost
bothA=cbind(info[,1:2], freqfocal, dpfocal)
filter=read.table("data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
both=merge(bothA, filter, by=c(1,2))
t1 = apply(both[,3:42], MARGIN=1, FUN=function(X) sum(X<=0 | X>=1, na.rm=TRUE) )  # keep any
# ignore the NAs
#> mean(t1>0)
#[1] 0.5606287   ## only 44% don't have a 0 or 1 allele freq
both2 = both[ t1==0 , ]   # 774,841
write.table(both2[,1:2], file="data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)
