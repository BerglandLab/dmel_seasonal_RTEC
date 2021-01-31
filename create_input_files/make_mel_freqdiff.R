### Create file with allele freq difference sp/fall and fall/frost for each popyear
####### Run in proclus
args <- commandArgs(trailingOnly = TRUE)
filein=args[1]
fileout=args[2]

load(filein)

#load("mel_freqdp_012016.Rdata")
#fileout = paste("mel_freqdiff_012016_",chrom,".txt", sep="")

library(data.table)
P.Y = popinfo[,6]
S = popinfo[,5]
PYu = unique(P.Y)
freqdiff_list = list()
seasons = c("f","fr","s")
focalfreq = freq
focalinfo = info

s12names = paste( paste(seasons[1], seasons[2], sep="_"), PYu, sep=".")
s13names = paste( paste(seasons[1], seasons[3], sep="_"), PYu, sep=".")
s23names = paste( paste(seasons[2], seasons[3], sep="_"), PYu, sep=".")
header1 = c("chrom","pos", s12names, s13names, s23names)
writetofile = file(fileout, open="wt")
writeLines(paste(header1, collapse=" "),  con=writetofile)

for (s in 1:nrow(focalfreq)){ # For each SNP
    focalSNP = focalfreq[s,]
    s1 = matrix(ncol=length(PYu), nrow=length(seasons))
    out1 = as.character(focalinfo[s,1])
    out2 = as.numeric(focalinfo[s,2])
    
    for (i in 1:length(PYu)){  # for each pop/year
        focalPOP = PYu[i]
        focalPOPind = which(P.Y == focalPOP)
        
        for (j in 1:length(seasons)){
            focalSEASON = seasons[j]
            if(focalSEASON %in% S[focalPOPind]){
                s1[j,i] = as.numeric(focalSNP[ which(P.Y == focalPOP & S == focalSEASON) ])
            } else {
                s1[j,i] = NA
            }
        }
    }
    sdiff = c( s1[1,]-s1[2,], s1[1,]-s1[3,], s1[2,]-s1[3,])
    writeLines(paste(c(out1, out2, sdiff), collapse=" "),  con=writetofile)
}

close(writetofile)

