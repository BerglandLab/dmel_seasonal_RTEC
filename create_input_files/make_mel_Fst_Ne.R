### Create file with allele freq difference sp/fall and fall/frost for each popyear
####### Run in proclus
args <- commandArgs(trailingOnly = TRUE)
filein=args[1]
fileout=args[2]

load(filein)

#load("mel_freqdp_012016.Rdata")
#fileout = paste("mel_freqdiff_012016_",chrom,".txt", sep="")

popID = popinfo[,1]
P = popinfo[,2]
Y = popinfo[,3]
R = popinfo[,4]
S = popinfo[,5]
PY = popinfo[,6]
sfpair = popinfo[,7]
ffrpair = popinfo[,8]

focal = which( (S=="s" | S=="f") & (PY=="AGA_14" | PY=="BA_12" | PY=="BHM_14" | PY=="CUA_14" | PY=="CUA_15" | PY=="CWI_14" | PY=="LMA_14" | PY=="MA_12" | PY=="NY_12" | PY=="SCPA_14" | PY=="SON_15" | PY=="TKA_14" | PY=="VI_12"  | PY=="WI_12"  | PY=="WI_13"
    | PY=="PA_9"  | PY=="PA_10"  | PY=="PA_11"  | PY=="PA_14"  | PY=="PA_15") )
#focalB = which( (S=="s" | S=="fr") & (PY=="co_13" | PY=="rd_12" ) )
#focal = c(focalA, focalB)

popID2 = popinfo[focal,1]
P2 = popinfo[focal,2]
Y2 = popinfo[focal,3]
R2 = popinfo[focal,4]
S2 = popinfo[focal,5]
PY2 = popinfo[focal,6]
#S2[S2=="fr"] = "f"
PYu = unique(PY2)

seasons = c("f","s")

focalfreq = freq[,focal]
focalinfo = info
focaldp = dp[,focal]

options(digits=6)

#s12names = paste( paste(seasons[1], seasons[2], sep="_"), PYu, sep=".")
header1 = c("chrom","pos", PYu)
writetofile = file(fileout, open="wt")
writeLines(paste(header1, collapse=" "),  con=writetofile)


fst <- function(p1,p2,d1,d2) {
    
    f.hat <- p1/2 + p2/2
    
    H.tot <- 2*(f.hat)*(1-f.hat) * (d1 + d2)/(d1 + d2 - 1)
    
    H.with <- p1*(1-p1) * (d1)/(d1-1) + p2*(1-p2) * (d2)/(d2-1)
    
    ret <- (H.tot - H.with)/H.tot
    
    return(ret)
}


Fst = function(p1,p2){
    Fbetween = (p1 + p2)*( (2-p1-p2)/2)
    Fwithin = p1*(1-p1) + p2*(1-p2)
    (Fbetween - Fwithin)/Fwithin
}

for (s in 1:nrow(focalfreq)){ # For each SNP
    focalSNP = focalfreq[s,]
    focalDP = focaldp[s,]
    s1 = vector(length=length(PYu))
    out1 = as.character(focalinfo[s,1])
    out2 = as.numeric(focalinfo[s,2])
    
    for (i in 1:length(PYu)){  # for each pop/year
        focalPOP = PYu[i]
        focalS = focalSNP[PY2 == focalPOP & S2=="s"]
        focalF = focalSNP[PY2 == focalPOP & S2=="f"]
        focalSdp = focalDP[PY2 == focalPOP & S2=="s"]
        focalFdp = focalDP[PY2 == focalPOP & S2=="f"]
    if (is.na(focalS) | is.na(focalF) ){
            s1[i] = NA
        } else if (focalS==0 & focalF==0){
            s1[i] = 0
        } else {
            s1[i] = format(as.numeric(fst(focalS, focalF, focalSdp, focalFdp)), digits=4)
        }
    }
    writeLines(paste(c(out1, out2, s1), collapse=" "),  con=writetofile)
}

close(writetofile)
