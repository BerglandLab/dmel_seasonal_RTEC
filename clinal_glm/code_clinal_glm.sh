#!/bin/sh

#  code_clinal_glm.sh
#  
#
#  Created by Heather Machado on 01/02/2021.
#  

##### Seasonal / latitudinal concordance
# Latitudinal glm regression
Rscript clinal_analysis_uniquepops_springPA_noMA.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_clinal_uniquepops_springPA_noMA.glm

# Run bootstrap matching of SNPs for seasonal/latitudinal concordance analysis
# Normal bootstrap sampling (done below), with all SNPs. Have to use version of R 3.2.2
for chrom in 2L 2R 3L 3R; do
    Rscript ../create_input_files/match_snps_dp_ch_SNPbySNP_recomb_chromposFilter.R recombinationRate/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.polymorphic.txt 100 data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt
done

# Calculating seasonal / latitudinal concordance compared to matched controls
Rscript run_analysis_seasonal_clinal_uniquepopsPAnoMA_parallelism_polymorphic_generic.R seasonal_glm/seasonal_analysis_all_paired20_2sample_caF_popyear.R clinal_glm/clinal_uniquepopsPAnoMA_parallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata

Rscript run_analysis_seasonal_clinal_uniquepopsPAnoMA_parallelism_polymorphic_generic.R seasonal_glm/seasonal_analysis_ca_f_s_popyear_paired.R clinal_glm/clinal_uniquepopsPAnoMA_CAparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata

Rscript run_analysis_seasonal_clinal_uniquepopsPAnoMA_parallelism_polymorphic_generic.R seasonal_glm/seasonal_analysis_eur_popyear.R clinal_glm/clinal_uniquepopsPAnoMA_EUparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata

# seasonal / latitudinal concordance by inversion status
# output: concordDF_Jan2021.txt (for Figure 2F)
Rscript analyze_Jan2021_inversionbrkpts_500Kb_clinalconcordance.R
