#  code_samplesizeinflation_glm_vs_rfm.sh
#  Created by Heather Machado


## Artificially increasing the sample size for the glm
Rscript seasonal_analysis_all_paired20_2sample_caF_popyear_adddp.R
Rscript seasonal_analysis_all_paired20_2sample_caF_popyear_adddp_poisson.R

## Artificially increasing the sample size for the rfm
Rscript fishers_exact_adddp_Mar2018.R
Rscript fishers_exact_adddp_poisson_Mar2018.R

## Analyzing the results (and plotting)
Rscript analysis_test_samplesize_overestimates_Mar2018.R
Rscript analysis_test_samplesize_overestimates_poisson_Mar2018.R
