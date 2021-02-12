#  code_samplesizeinflation_glm_vs_rfm.sh
#  Created by Heather Machado


# For supplemental figure
# script: analysis_test_samplesize_overestimates_poisson_Mar2018.R (see end)
# dependencies: test_samplesize_obs_exp_poisson_ggplot.Rdata, test_samplesize_obs_exp_ggplot.Rdata


## Artificially increasing the sample size for the glm
for i in 0 10 100; do
    Rscript seasonal_analysis_all_paired20_2sample_caF_popyear_adddp.R mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata mel_all_paired20_2sample_caF_popyear.f_s.glm mel_all_paired20_2sample_caF_popyear_dpplus$i.f_s.glm $i
    Rscript seasonal_analysis_all_paired20_2sample_caF_popyear_adddp_poisson.R mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata mel_all_paired20_2sample_caF_popyear_dpplus$i\_poisson.f_s.glm $i
done

## Artificially increasing the sample size for the rfm
while read line; do
    for i in 0 10 100; do
        Rscript fishers_exact_adddp_Mar2018.R ../data/mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata $line rank_fisher_exact_dpplus$i\_$line.txt $i
        Rscript fishers_exact_adddp_poisson_Mar2018.R ../data/mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata $line rank_fisher_exact_dpplus$i\_$line\_poisson.txt $i
        done
done < ../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt

## Analyzing the results (and plotting)
Rscript analysis_test_samplesize_overestimates_Mar2018.R
# Creates supplemental figure
Rscript analysis_test_samplesize_overestimates_poisson_Mar2018.R
