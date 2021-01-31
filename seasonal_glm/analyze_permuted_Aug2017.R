## Aug 2017


filter = read.table("../data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
perms = list()

for (i in 1:100){
    glm = read.table(paste("permute",i, "_mel_all_paired20_2sample_caF_popyear.glm", sep=""), header=TRUE, stringsAsFactors=FALSE)
    glm2A = na.omit(merge(glm, filter, by=c(1,2)))
    perms[[i]] = glm2A[glm2A[,4]<1 & glm2A[,4]>0  & glm2A[,1]!="X",]
}

## Getting the distribution of P-values for the permuted dataset
make_bins = function(x, size){
    x = na.omit(x)
    my_seq = seq(from=0, to=1, by=size)
    my_sum=vector()
    for (i in 1:(length(my_seq)-1) ){
        my_sum[i] = sum(x>=my_seq[i] & x<my_seq[i+1])
    }
    list(my_sum, my_seq)
}


perm_test = make_bins(perms[[1]][,4], 0.001)

perm_bins = matrix(ncol=length(perm_test[[2]])-1, nrow=100)
colnames(perm_bins) = perm_test[[2]][1:(length(perm_test[[2]])-1)]

for (i in 1:100){
    perm_bins[i,] = make_bins(perms[[i]][,4], 0.001)[[1]]
}


write.table(perm_bins, file="permutation_prop_by_bin001.txt", row.names=FALSE, quote=FALSE)

#trimmed = mybins[[1]][1:(length(mybins[[1]])-1)]
#plot(mybins[[2]][2:(length(mybins[[2]])-1)], trimmed, type="l", log="x")


