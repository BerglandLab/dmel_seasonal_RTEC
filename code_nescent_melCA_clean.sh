#### Code for Machado et al 2021
# Broad geographic sampling reveals predictable, pervasive, and strong seasonal adaptation in Drosophila

#### Master script referring to analysis scripts
#### Heather Machado



###### Data processing
## Chromosomes analyzed (loop through)
chrom=2R
chrom=2L
chrom=3L
chrom=3R

## Create RDataFile from vcf files (after Varscan SNP calling)
Rscript create_input_files/make_melRDataFile.R ../mapping2015/RTECrun2/06_variant_calling/mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf mel_freqdp_$chrom\_042016.Rdata

## Calculate Ne
Rscript create_input_files/Ne.R mel_freqdp_$chrom\_042016.Rdata mel_nescent2_samplesize.txt mel_freqdp_$chrom\_042016_Ne.Rdata auto

## Combine the wild and lines of BHM
Rscript create_input_files/combine_BHM.R mel_freqdp_$chrom\_042016_Ne.Rdata mel_freqdp_$chrom\_042016_Ne_fixed.Rdata

## Create one file for all chromosomes
Rscript create_input_files/combine_freqdp.R data/mel_freqdp 042016_fixed.Rdata  ## created mel_freqdp_042016_fixed.Rdata

## Calculate means for spring and fall for each group
Rscript create_input_files/scripts/make_mel_means_dpfreq_pops_allmean.R data/mel_freqdp_042016_Ne_fixed.Rdata data/means_dpfreq.paired20.allmean.Rdata data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt

## Filtering for SNPs that are polymorphic in ALL 40 samples (of the 20 paired sample dataset)
Rscript create_input_files/filter_polymorphic.R

## Calculate difference per season
Rscript create_input_files/make_mel_freqdiff.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_freqdiff_042016_$chrom.txt

## calculating mean dp per pop
Rscript calculate_mean_dp.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_meandp.$chrom.txt
Rscript calculate_mean_dp.R data/mel_freqdp_$chrom\_042016_fixed.Rdata data/mel_meandp_N.$chrom.txt

##  merge vcf file with SNPeff (SnpEff 4.2) (dm5.48.genome)
java -jar snpEff.jar download dm5.48
java -Xmx4G -jar snpEff.jar dm5.48.reference mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf

##  calculating Fst for samples that have a paired seasonal site (and climate data)
Rscript create_input_files/make_mel_Fst_Ne.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_Fst_Ne.seasonalpaired.txt
Rscript create_input_files/make_mel_Fst_Ne_core20.R data/mel_freqdp_042016_Ne_fixed.Rdata $chrom data/mel_Fst_Ne.core20.sp_fall.txt



###### PCA analysis
Rscript scripts_misc/PCA_clean.R



###### Seasonal analyses: GLM
Rscript seasonal_glm/seasonal_analysis_all_paired20_2sample_caF_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_paired20_2sample_caF_popyear.f_s.glm
Rscript seasonal_glm/seasonal_analysis_ca_f_s_popyear_paired.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_ca_popyear_paired.f_s.glm
Rscript seasonal_glm/seasonal_analysis_eur_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_eur_popyear_paired.f_s.glm

# seasonal regression without clinal populations (PA)
Rscript seasonal_glm/seasonal_analysis_all_nonclinal_paired20_2sample_caF_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm

# glm permutations
seasonal_glm/glm_permutations/code_glm_permutations.sh



###### Seasonal analyses: RFM (rank Fisher's method)
seasonal_fishersmethod/code_fishers_exactJ.sh

# permutations
seasonal_fishersmethod/rfm_permutations/anaysis_allseasonal_permutations_June2020.R



###### Seasonal analyses: bayenv
# includes permutations
seasonal_bayenv/code_seasonal_bayenv_Jan2020.sh



##### Seasonal / latitudinal concordance
# Latitudinal glm regression
Rscript clinal_glm/clinal_analysis_uniquepops_springPA_noMA.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_clinal_uniquepops_springPA_noMA.glm

# Run bootstrap matching of SNPs for seasonal/latitudinal concordance analysis
# Normal bootstrap sampling (done below), with all SNPs. Have to use version of R 3.2.2
for chrom in 2L 2R 3L 3R; do
    Rscript create_input_files/match_snps_dp_ch_SNPbySNP_recomb_chromposFilter.R recombinationRate/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.polymorphic.txt 100 data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt
done

# Calculating seasonal / latitudinal concordance compared to matched controls
Rscript clinal_glm/run_analysis_seasonal_clinal_uniquepopsPAnoMA_parallelism_polymorphic_generic.R seasonal_glm/seasonal_analysis_all_paired20_2sample_caF_popyear.R clinal_glm/clinal_uniquepopsPAnoMA_parallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata

Rscript clinal_glm/run_analysis_seasonal_clinal_uniquepopsPAnoMA_parallelism_polymorphic_generic.R seasonal_glm/seasonal_analysis_ca_f_s_popyear_paired.R clinal_glm/clinal_uniquepopsPAnoMA_CAparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata

Rscript clinal_glm/run_analysis_seasonal_clinal_uniquepopsPAnoMA_parallelism_polymorphic_generic.R seasonal_glm/seasonal_analysis_eur_popyear.R clinal_glm/clinal_uniquepopsPAnoMA_EUparallelism_vs_controlsDec2017_quant_polymorphic_ext.Rdata

# seasonal / latitudinal concordance by inversion status
# output: concordDF_Jan2021.txt (for Figure 2F)
Rscript clinal_glm/analyze_Jan2021_inversionbrkpts_500Kb_clinalconcordance.R

# latitudinal bayenv analysis (including Figure S5)
clinal_bayenv/code_bayenv.sh
