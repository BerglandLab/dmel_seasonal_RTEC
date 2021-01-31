## January 2016
# Analysis of mel data (CA mel and nescent, mapped Dec 2015)

#setwd("/Users/heathermachado/drosAnalysisCA")
# Performed in proclus

## USAGE: Rscript make_melRDataFile.R ../mapping2015/RTECrun2/06_variant_calling/mel_X_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf mel_freqdp_X_042016.Rdata

#### or using a run script
## USAGE: qsub run_make_melRDataFile.R ../mapping2015/RTECrun2/06_variant_calling/mel_X_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf mel_freqdp_X_042016.Rdata

args <- commandArgs(trailingOnly = TRUE)
filein=args[1]
fileout=args[2]
#mel = read.table("/hsgs/projects/petrov/hmachado/mapping2015/mel/06_variant_calling/mel_mapping2015.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf",
#                 comment.char = "", header=TRUE, stringsAsFactors=FALSE)
mel = read.table(filein, comment.char = "", header=TRUE, stringsAsFactors=FALSE)

freq = mel[,10:84]
dp = mel[,85:159]
info = mel[,1:9]
info[,1] = as.character(info[,1])
# check for consistency with old data
pops = colnames(mel)[10:84]
S = c("f","fr",
      "f","f","f",
      "f",    "s",
      "fr","fr","fr",
      "fr","fr","f",
      "f","f","s",
      "s",    "f",
      "s","f","s",
      "f","s","f",
      "s","s",
      "s","n","f",
      "f","fr","s",
      "s","s","f",
          "s","f",
      "f",    "s",
      "f","s","fr","s",
    "f","s","f",
    "s","s","f",
    "s","f","f",
    "s","s","f",
    "s","s","f",
    "f","s","f",
    "s","f","s",
    "f","s","s",
    "f","s","f",
    "s","s","s",
    "f","s"
)

Y = c(12,12,
      13,13,13,
      13,   12,
      12,12,12,
      13,13,12,
      12,12,13,
      13,   12,
      12,12,12,
      09,09,12,
      12,10,
      08,03,09,
      10,11,10,
      11,10,11,
         12,12,
      12,   12,
      12,12,12,12,
    13,13,13,
    13,13,14,
    14,14,14,
    14,14,14,
    14,14,14,
    14,14,14,
    14,14,14,
    14,14,14,
    15,15,15,
    15,15,15,
    15,15
)

P=c("gb","co",
    "gb","co","fb",
    "rd",     "rd",
    "fb","gb","rd",
    "co","gb","fb",
    "rd","co","rd",
    "co",     "BA",
    "BA","MA","MA",
    "ME","PA","NY",
    "NY","FL",
    "GA","NC","PA",
    "PA","PA","PA",
    "PA","SC","PA",
         "PA","WI",
    "PA",     "VA",
    "VI","VI","WI","WI",
    "OUK","WI","WI",
    "PA","OUK","AGA",
    "AGA","BHMf1","BHMwild",
    "BHMf1","BHMwild","CUA",
    "CUA","CWI","CWI",
    "LMA","LMA","PA",
    "PA","SCPA","SCPA",
    "TKA","TKA","TWA",
    "CUA","CUA","PA",
    "PA","SCPA","SON",
    "SON","UIL"
)

R = c(rep("ca", times=17), rep("other", times=5), "pa", rep("other", times=5), rep("pa", times=5),
      "other",rep("pa", times=2), "other","pa",rep("other", times=8),"pa",rep("other", times=13),rep("pa", times=2),rep("other", times=7),rep("pa", times=2),rep("other", times=4))

popinfo = cbind(pops, P, Y, R, S)
P.Y = apply(popinfo[,2:3], MARGIN=1, FUN=paste, collapse="_")
popinfo = cbind(pops, P, Y, R, S, P.Y)
PYseas = tapply(popinfo[,5], INDEX=P.Y, FUN=unlist)
sf = cbind(names(PYseas), lapply(PYseas, FUN=function(X) any(is.element(X,"f")) ),
           lapply(PYseas, FUN=function(X) any(is.element(X,"s")) ) )
sf_paired = unlist(sf[sf[,2]==TRUE & sf[,3]==TRUE ,1])
ffr = cbind(names( PYseas), lapply(PYseas, FUN=function(X) any(is.element(X,"f")) ),
            lapply(PYseas, FUN=function(X) any(is.element(X,"fr")) ) )
ffr_paired = unlist(ffr[ffr[,2]==TRUE & ffr[,3]==TRUE ,1])
sf_pair = popinfo[,6] %in% sf_paired
ffr_pair = popinfo[,6] %in% ffr_paired
popinfo = cbind(popinfo, sf_pair, ffr_pair)

save(freq, dp , info, popinfo, file=fileout)

