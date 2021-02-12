### Aprili 2017
## merging the 24 pops

args <- commandArgs(trailingOnly = TRUE)
#args = "fisherlist.txt"

inlist = read.table(args[1], stringsAsFactors=FALSE)[,1]

#colnames(df1) = c("chrom","pos",paste(inlist[1],".G",sep=""),paste(inlist[1],".L",sep=""))
filelist = list()
for (i in 1:length(inlist) ){
    filelist[[i]] = na.omit(read.table(inlist[i], stringsAsFactors=FALSE, header=TRUE))
}

df1A = filelist[[1]][,1:2]
filter = read.table("../../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt")
auto = filter[filter[,1]!="X", ]
df1 = merge(df1A, auto, by=c(1,2))

for (i in 1:(length(inlist)-1)){
    df2 = filelist[[i+1]][,1:2]
    df3 = merge(df1, df2, by=c(1,2))
    df1 = df3
}

for (i in 1:length(inlist) ){
    outtab = merge(filelist[[i]], df1, by=c(1,2))
    write.table(outtab, file=paste(inlist[i], ".merged", sep=""), quote=FALSE, col.names=TRUE, row.names=FALSE)
}

