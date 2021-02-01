#!/bin/sh

#  code_seasonal_glm.sh
#  
#
#  Created by Heather Machado on 01/02/2021.
#  

###### Seasonal analyses: GLM
Rscript seasonal_analysis_all_paired20_2sample_caF_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_paired20_2sample_caF_popyear.f_s.glm
Rscript seasonal_analysis_ca_f_s_popyear_paired.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_ca_popyear_paired.f_s.glm
Rscript seasonal_analysis_eur_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_eur_popyear_paired.f_s.glm

# seasonal regression without clinal populations (PA)
Rscript seasonal_analysis_all_nonclinal_paired20_2sample_caF_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm

# glm permutations
glm_permutations/code_glm_permutations.sh
