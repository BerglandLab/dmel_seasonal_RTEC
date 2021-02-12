# Jan 2020
# Bayenv analysis of seasonal data

# create input files
Rscript analysis_season_bayenv_Ne_nonpooled.R

## calculate matrix for 10K sites, 100K iterations (non-pooled mode)
bayenv2 -i seasonal_paired20_snpsfile_Ne_Jan2020_sample10K.txt -p 40 -k 100000 -r 63478 > matrix_10Kbsub100Kiter_seasonal_paired20.nonpooled.out



#################   running for each SNP (just autosomes)
# submitting as an array job
# Genome has been split into 882 chunks for quicker processing
for i in {1..882}; do
    Rscript run_bayenv_nonpooled_batch.R intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.${LSB_JOBINDEX} intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.$i seasonal_paired20_environ_Jan2020_springfall.txt results/results_nonpooled_iter100K
done
cat results/results_nonpooled_iter100K.seasonal_paired20_snpsfile_Ne_Jan2020.txt.* > results/results_nonpooled_iter100K.seasonal_paired20_snpsfile_Ne_Jan2020.txt
rm results/results_nonpooled_iter100K.seasonal_paired20_snpsfile_Ne_Jan2020.txt.*


############ running 100 permutations- June 2020
for j in {1..100}; do
    for i in {1..882}; do
        Rscript run_bayenv_nonpooled_batch_parallel.R intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.$i intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.$i bayenv_permutations/perm$j\_seasonal_paired20_environ_June2020_springfall.txt results/perm$j\_nonpooled
    done
done
# Note, full permutation results too large to store in github

# For plotting, including Figure S4
#   object for plotting provided (full permutation results too large)
analysis_season_bayenv_Ne_nonpooled_results.Rmd
