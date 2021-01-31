## Oct 2016
# fisher's method using fisher's exact test

#setwd("/Users/heathermachado/nescent_melCA/fishers_method_Oct2016")
#
#load("../data/mel_freqdp_042016_2L_fixed.Rdata")
#focalPopYear = "co_13"
#outprefix = "co_13_fishersmethod_testout"
#fileout = paste(outprefix, ".txt", sep="")

## Usage: spring_fall_pvalue.R data.Rdata focalPopulationYear outfile
#while read line; do
#  sbatch run_fishers_exact_pvalue_rank_Oct2016.sh $line.spring_fall_fisherexact.txt $line.spring_fall_fisherexact_pvalue_rank.txt
#done < paired_spring_fall_populations.txt
## In run_fishers_exact_pvalue_rank_Oct2016.sh:
#   Rscript fishers_exact_pvalue_rank_Oct2016.R $line.spring_fall_fisherexact.txt $line.spring_fall_fisherexact_pvalue_rank.txt
#args = c("equalbins_Oct10/co_13.spring_fall_fisherexact.txt","co_13.spring_fall_fisherexact_pvalue_rank.txt")


#install.packages("data.table", repos="http://cran.r-project.org", lib="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")
require(data.table,lib.loc="/home/hmachado/R/x86_64-unknown-linux-gnu-library/")

args = commandArgs(trailingOnly=TRUE)

focal = read.table(args[1], header=TRUE, stringsAsFactors=FALSE)
fileout = args[2]

#dpbins = 10  ## the number of dp bins 
#quant1 = quantile(focal$dp.total, probs=seq(from=0, to=1, by=1/dpbins))
#quant1[1] = quant1[1]-1 ## change the first bin so that it is LOWER than the minimum DP
#altbins = 10  ## the number of dp bins 

######### binning
outlist = list()
count = 0

dpstep=10
altstep=5
mindp = min(focal$dp.total)
maxdp = max(focal$dp.total)
dpbins = seq(from=mindp, to=maxdp, by=dpstep)
dpbins[length(dpbins)+1] = maxdp+1

for (i in 1:(length(dpbins)-1)) { # for each dp bin
  focaldp = focal[focal$dp.total>=dpbins[i] & focal$dp.total<dpbins[i+1] & !is.na(focal$greater) & !is.na(focal$less),  ]
  if (nrow(focaldp)==0){
      next
  }
  minalt = min(focaldp$alt.total)
  maxalt = max(focaldp$alt.total)
  altbins = seq(from=minalt, to=maxalt, by=altstep)
  altbins[length(altbins)+1] = maxalt+1
  
  for (j in 1:(length(altbins)-1)) {
    count = count + 1
    focalaltA = focaldp[focaldp$alt.total>=altbins[j] & focaldp$alt.total<altbins[j+1],  ]
    if (nrow(focalaltA)==0){
        next
    }
    Nsites = nrow(focalaltA)
    rand1 = runif(Nsites,0,1) # random number between 0 and 1
    rand2 = runif(Nsites,0,1) # random number between 0 and 1
    focalaltB = focalaltA[sample(1:Nsites, size=Nsites), ]
    focalaltG = focalaltB[order(focalaltB$greater),]
    focalaltG$greater.p = ((1:Nsites) - rand1)/Nsites 
    focalaltL = focalaltG[order(focalaltG$less),]
    focalaltL$lesser.p = ((1:Nsites) - rand2)/Nsites 
    #greaterRank = order(as.numeric(focalalt$greater))-rand1
    #lesserRank = order(as.numeric(focalalt$less))-rand2
    #focalalt$greater.p = greaterRank/Nsites
    #focalalt$lesser.p = lesserRank/Nsites
    focalout = data.frame(focalaltL)
    outlist[[count]] = focalout
  }
}

count = count+1
focalna = focal[is.na(focal$greater) | is.na(focal$less),  ]
if (nrow(focalna)>0){
    focalna$greater.p = "NA"
    focalna$lesser.p = "NA"
    outlist[[count]] = focalna
}
out_tableA = data.frame(rbindlist(outlist))
out_table = out_tableA[order(out_tableA[,1], out_tableA[,2]), ]
out_table2 = out_table[!duplicated(out_table[,c(1,2)]), ]

write.table(out_table2, file=fileout, row.names=FALSE, quote=FALSE, col.names=TRUE)



