## June 2017

# Permute final seasonal glms
# Folder: /scratch/users/hmachado/nescent_melCA/glm_permutations

# 1) create list of populations to permute spring/fall
# 2) run seasonal regression switching spring/fall for those populations
# 3) use these regressions in calculating:

## in R
pops20 =  read.table("../../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]


switching = list()
for (i in 1:500){
    toswitch = vector()
    for (j in 1:20){
        toswitch[j] = sample(0:1, size=1)
    }
    switching[[i]] = pops20[which(toswitch==1)]
}
save(switching, file="switching_pop20_permutations.Rdata")


# run permutations
for i in {1..500}; do
    Rscript permutation_seasonal_analysis_all_paired20_2sample_caF_popyear.R ../data/mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata permute$i\_mel_all_paired20_2sample_caF_popyear.glm $i
done
