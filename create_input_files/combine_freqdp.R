## June 2016

## Combining the freqdp files across chromosomes into one file
## USAGE: prefix, suffix
#args = c("data/mel_freqdp", "042016_fixed.Rdata")

args <- commandArgs(trailingOnly = TRUE)
library(data.table, lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
prefix = args[1]
suffix = args[2]
outfile = paste(prefix, suffix, sep="_")
chroms = c("2L","2R","3L","3R","X")

freqlist = list()
dplist = list()
infolist = list()

for (i in 1:length(chroms)){
    load(paste(prefix,chroms[i], suffix, sep="_"))
    freqlist[[i]] = freq
    dplist[[i]] = dp
    infolist[[i]] = info
}

freq = data.frame(rbindlist(freqlist))
dp = data.frame(rbindlist(dplist))
info = data.frame(rbindlist(infolist))

save(popinfo, freq, dp, info, file=outfile)