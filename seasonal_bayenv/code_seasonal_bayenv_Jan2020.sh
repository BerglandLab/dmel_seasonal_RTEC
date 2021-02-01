# Jan 2020
# Bayenv analysis of seasonal data

# create input files
Rscript analysis_season_bayenv_Ne_nonpooled.R

## calculate matrix for 10K sites, 100K iterations (non-pooled mode)
bsub -J bay_mat -n2 -R "select[mem>16000]" -R "rusage[mem=16000]" -M 16000 -R "span[hosts=1]" -q basement -e bsub.error.%J -o bsub.output.%J "/lustre/scratch116/casm/cgp/users/hm8/software/bayenv2 -i seasonal_paired20_snpsfile_Ne_Jan2020_sample10K.txt -p 40 -k 100000 -r 63478 > matrix_10Kbsub100Kiter_seasonal_paired20.nonpooled.out"



#################   running for each SNP (just autosomes)
# submitting as an array job
bsub -J 'bammixbatch[1-882]' -n1 -R "select[mem>4000]" -R "rusage[mem=4000]" -M 4000 -R "span[hosts=1]" -q long -e error/bsub.error.%J -o error/bsub.output.%J 'Rscript run_bayenv_nonpooled_batch.R intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.${LSB_JOBINDEX} intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.${LSB_JOBINDEX} seasonal_paired20_environ_Jan2020_springfall.txt results/results_nonpooled_iter100K'

cat results/results_nonpooled_iter100K.seasonal_paired20_snpsfile_Ne_Jan2020.txt.* > results/results_nonpooled_iter100K.seasonal_paired20_snpsfile_Ne_Jan2020.txt
rm results/results_nonpooled_iter100K.seasonal_paired20_snpsfile_Ne_Jan2020.txt.*



############ running 100 permutations- June 2020
for j in {1..100}; do
        bsub -J 'bammixbatch[1-882]' -n8 -R "select[mem>2000]" -R "rusage[mem=2000]" -M 2000 -R "span[hosts=1]" -q normal -e error/bsub.error.%J -o error/bsub.output.%J '/software/R-3.6.1/bin/Rscript run_bayenv_nonpooled_batch_parallel.R intfiles/seasonal_paired20_snpsfile_Ne_Jan2020.txt.${LSB_JOBINDEX} intfiles/seasonal_paired20_snpsfileinfo_Ne_Jan2020.txt.${LSB_JOBINDEX} bayenv_permutations/perm'$j'_seasonal_paired20_environ_June2020_springfall.txt results/perm'$j'_nonpooled'
done



# For plotting, including Figure S4
analysis_season_bayenv_Ne_nonpooled_results.Rmd
