## GLM permutations

## Run glm permutations
# output: 500 glm permutations
./code_permute_paired20_2sample_caF_popyear.sh


## combine permutations (N sites per quantile: Figure 2A)
# output: permutation_prop_by_bin001_perm500.txt
Rscript analyze_permuted_July2020.R


## number of snps per permutation at pvalue < 0.00389
## dependency: combineperms.sh
# output: perm500_n_less_p00389.txt
./code_combine_glm_permute_paired20_2sample_caF_June2020.sh
    

## combine permutations (N sites per p-value, p-value per quantile: Figure 2B)
# output: pvalueN_mean_obs_less_perm_chrom_v2_pvalue004.txt, obs_perm500_enrich_chrom_v2_pvalue004.Rdata
bsub -J "permpart[1-500]" -q normal -R 'select[mem>=16000] rusage[mem=16000]' -M16000 -n1  -o log/bsub.%J.out -e log/bsub.%J.err "Rscript analyze_permuted_July2020_suppfigures_permsplitsv2_pvalue004.R \$LSB_JOBINDEX"
Rscript analyze_permuted_July2020_suppfigures_aggregate_permsplitsv2_pvalue004.R


## for each inversion, 500Kb around breakpoints
# output: obs_perm_filter_pvalue004N_byregion_inversionbrkpts_500Kb.txt
bsub -J "permpart[1-500]" -q normal -R 'select[mem>=16000] rusage[mem=16000]' -M16000 -n1  -o log/bsub.%J.out -e log/bsub.%J.err "Rscript analyze_permuted_July2020_suppfigures_permsplitsv2_pvalue004_inversionbrkpts_500Kb.R \$LSB_JOBINDEX"
Rscript analyze_permuted_July2020_suppfigures_aggregate_permsplitsv2_pvalue004_inversionbrkpts_500Kb.R

# output: chr_inv_enrichment_pvalues_Jan2021.txt (used in Figure 2D)
Rscript analyze_permuted_July2020_suppfigures.R

