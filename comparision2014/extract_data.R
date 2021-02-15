
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### 2014 dataset from https://datadryad.org/stash/dataset/doi:10.5061/dryad.v883p

### libraries
  library(data.table)

### load vcf file
  vcf <- fread("/scratch/aob2x/drosRTEC/6d_v7.3_output.vcf", skip=24)
  vcf[,sp:=as.numeric(gsub("SP=", "", tstrsplit(INFO, ";")[[9]]))]
  vcf[,sq:=as.numeric(gsub("SQ=", "", tstrsplit(INFO, ";")[[10]]))]
  setnames(vcf, c("#CHR", "POS"), c("chr", "pos"))

### extract out subset
  sig2014 <- vcf[sq<=.3, c("chr", "pos", "sq"), with=F]

### save
  save(sig2014, file="/scratch/aob2x/drosRTEC/sig2014.Rdata")
