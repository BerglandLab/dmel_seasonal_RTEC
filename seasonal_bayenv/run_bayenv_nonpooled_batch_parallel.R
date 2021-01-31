
## Jan 2020

args = commandArgs(trailingOnly=TRUE)

.libPaths(c("/lustre/scratch116/casm/cgp/users/hm8/R-3.6-modules","/nfs/users/nfs_h/hm8/R/x86_64-pc-linux-gnu-library/3.6","/software/R-3.6.1/lib/R/library") )


library(foreach)
library(doMC)

# Test args
#args = c("intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.1", "intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.1", "seasonal_paired20_environ_Jan2020_springfall.txt", "results/results_nonpooled_iter100K")
#args = c("intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.1", "intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.1", "permutations/perm1_seasonal_paired20_environ_Jan2020_springfall.txt", "testout_perm1_parallel_March10")

# Test files
#snpsfile = read.table("intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.1", header=F, stringsAsFactors = F)
#infofile = read.table("intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.1", header=F, stringsAsFactors = F)



snpsfile = read.table(args[1], header=F, stringsAsFactors = F)
infofile = read.table(args[2], header=F, stringsAsFactors = F)
environfile = args[3]
outprefix = args[4]
resultslist = list()
outfile = paste(outprefix, ".",unlist(strsplit(args[1], split="/", fixed=T))[2], sep="")

randomID = sample(100000:900000, size=1)
system(paste("mkdir -p input/",randomID, sep="") )
system(paste("mkdir -p output/",randomID, sep="") )

# Run bayenv for each SNP in the file
registerDoMC(4)
  
  
output = foreach(i=1:nrow(infofile), .combine = 'rbind') %dopar% {
  # select snp infp
  focaldata = snpsfile[((2*i)-1):(2*i),]
  focalinfo = infofile[i,]

  # run bayenv function
  write.table(focaldata, file=paste("input/",randomID,"/",focalinfo[1],".",focalinfo[2],".input", sep="" ), quote=FALSE, col.names = F, row.names = F, sep="\t")
  system(paste("/lustre/scratch116/casm/cgp/users/hm8/software/bayenv2 -i input/",randomID,"/",focalinfo[1],".",focalinfo[2],".input -e ",environfile," -m matrix_10Kbsub100Kiter_seasonal_paired20.nonpooled.100Kiter.mat -p 40 -n 1 -k 100000 -t -c -r 73480 -o output/",randomID, "/",focalinfo[1],".",focalinfo[2],".output", sep="") )
  # R now waits for the bayenv command to finish

  # create vector of results
  chrompos = unlist(focalinfo)
  myoutfile = paste("output/",randomID, "/",chrompos[1],".",chrompos[2],".output.bf", sep="")
  if(file.exists(myoutfile)){
      focalresult = read.table(myoutfile, stringsAsFactors = F, header=F)
      return( c(chrompos[1],chrompos[2],unlist(focalresult[2:4]) ))
  } else {
      return( c(chrompos[1],chrompos[2],"NA","NA","NA"))
  }
}
write.table(output, file=outfile, quote=F, col.names = F, row.names = F, sep = "\t")
system(paste( "chmod 775 ", outfile, sep="") )
system(paste("rm -r input/",randomID, "/", sep="" ))
system(paste("rm -r output/",randomID, "/", sep="" ))

