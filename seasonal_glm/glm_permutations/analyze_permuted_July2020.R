## july 2020

## for 500 permutations



perm = read.table(paste("results/permute1_mel_all_paired20_2sample_caF_popyear.glm", sep=""), header=TRUE, stringsAsFactors=FALSE)


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


perm_test = make_bins(perm[,4], 0.001)

perm_bins = matrix(ncol=length(perm_test[[2]])-1, nrow=500)
colnames(perm_bins) = perm_test[[2]][1:(length(perm_test[[2]])-1)]


for (i in 1:500){
    focalperm = read.table(paste("results/permute",i, "_mel_all_paired20_2sample_caF_popyear.glm", sep=""), header=TRUE, stringsAsFactors=FALSE)
    perm_bins[i,] = make_bins(focalperm[,4], 0.001)[[1]]
}


write.table(perm_bins, file="permutation_prop_by_bin001_perm500.txt", row.names=FALSE, quote=FALSE)
