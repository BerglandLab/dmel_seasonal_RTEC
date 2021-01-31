## rdata2bayenvformat.R 
## Feb 2019

## Converting my data to format for bayenv

# Load in file NOT adjusted for Ne
load("/Users/hm8/stanford/nescent_melCA/data/mel_freqdp_042016_fixed.Rdata")

popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

# selecting the 4 samples for clinal analysis
focal = which(popID=="melPA_72011_SPT" | popID=="melSC_072010_SPT" | popID=="melGA_072008_SPT" | popID=="melFL_072010_SPT")
popID1 = popID[focal]
#"melFL_072010_SPT" "melGA_072008_SPT"  "melPA_72011_SPT" "melSC_072010_SPT" 
lat = c(25.5, 30.99, 39.88, 33.39)
freqfocal = freq[, focal]
dpfocal = dp[, focal] ### only looking at fall to frost
bothA=cbind(info[,1:2], freqfocal, dpfocal)
filter=read.table("/Users/hm8/stanford/nescent_melCA/data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
both=merge(bothA, filter, by=c(1,2))

# calculating alt and ref reads, putting into bayenv format
altfocal = round(both[,3:6]*both[,7:10])
reffocal =  both[,7:10]-altfocal
bothaltref=cbind(both[,1:2], altfocal, reffocal)
snpsfile = matrix(unlist(t(bothaltref[,3:ncol(bothaltref)])), ncol=4, byrow=T)
write.table(snpsfile, file="clinal_snpsfile_Feb2019.txt", col.names = F, row.names = F, quote=F, sep="\t")
write.table(both[,1:2], file="clinal_snpsfileinfo_Feb2019.txt", col.names = F, row.names = F, quote=F, sep="\t")

# creating test files
write.table(snpsfile[1:10,], file="clinal_snpsfile_Feb2019_test10.txt", col.names = F, row.names = F, quote=F, sep="\t")

# creating file for matrix estimation (10K snps)
non0 = apply(bothaltref[,3:6], MARGIN=1, FUN=sum)
non0ref = apply(bothaltref[,7:10], MARGIN=1, FUN=sum)
bothaltrefNon0 = na.omit(bothaltref[non0>0 & non0ref>0, ])
bothaltref1K = bothaltrefNon0[sample(1:nrow(bothaltrefNon0), size=1000, replace=F), ]
snpsfile1K = matrix(unlist(t(bothaltref1K[,3:ncol(bothaltref1K)])), ncol=4, byrow=T)
write.table(snpsfile1K, file="clinal_snpsfile_Feb2019_sample1K.txt", col.names = F, row.names = F, quote=F, sep="\t")

bothaltref10K = bothaltrefNon0[sample(1:nrow(bothaltrefNon0), size=10000, replace=F), ]
snpsfile10K = matrix(unlist(t(bothaltref10K[,3:ncol(bothaltref10K)])), ncol=4, byrow=T)
write.table(snpsfile10K, file="clinal_snpsfile_Feb2019_sample10K.txt", col.names = F, row.names = F, quote=F, sep="\t")


######## Removing snps without any alternate reads
#write.table(bothaltrefNon0[,1:2], file="clinal_snpsfileinfo_non0_Feb2019.txt", col.names = F, row.names = F, quote=F, sep="\t")
#write.table(bothaltrefNon0, file="clinal_allsnps_non0_unformatted_Feb2019.txt", col.names = F, row.names = F, quote=F, sep="\t")

bothaltrefNon0 = read.table("clinal_allsnps_non0_unformatted_Feb2019.txt", stringsAsFactors = F)

# one file per snp
for (i in 1:nrow(bothaltrefNon0)){
  focalsnp = bothaltrefNon0[i,]
  snpfile1 = matrix( focalsnp[,3:length(focalsnp)], ncol=4, byrow=T)
  write.table(snpfile1, file=paste("snpfiles/clinal_snpsfile_Feb2019_snp", i,".txt", sep=""), col.names = F, row.names = F, quote=F, sep="\t")
}
# 1939172 files

# ENVIRONFILE
write(lat, file="clinal_environfile_Feb2019.txt", sep="\t")


# SAMPLEFILE
# melPA_72011_SPT: 75
# melSC_072010_SPT: 48
# melGA_072008_SPT: 51
# melFL_072010_SPT: 39
ind = c(39, 51, 75, 48)
samplesizes = ind*2
write(samplesizes, file="clinal_samplefile_Feb2019.txt", sep="\t")


