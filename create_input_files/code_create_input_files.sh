#!/bin/sh

#  code_create_input_files.sh
#  
#
#  Created by Heather Machado on 01/02/2021.
#

###### Data processing
## Chromosomes analyzed (loop through)
chrom=2R
chrom=2L
chrom=3L
chrom=3R

## Create RDataFile from vcf files (after Varscan SNP calling)
# vcf files available on datadryad: https://doi.org/10.5061/dryad.4r7b826
Rscript make_melRDataFile.R mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf mel_freqdp_$chrom\_042016.Rdata

## Calculate Ne
Rscript Ne.R mel_freqdp_$chrom\_042016.Rdata mel_nescent2_samplesize.txt mel_freqdp_$chrom\_042016_Ne.Rdata auto

## Combine the wild and lines of BHM
Rscript combine_BHM.R mel_freqdp_$chrom\_042016_Ne.Rdata mel_freqdp_$chrom\_042016_Ne_fixed.Rdata

## Create one file for all chromosomes
Rscript combine_freqdp.R data/mel_freqdp 042016_fixed.Rdata  ## created mel_freqdp_042016_fixed.Rdata

## Correct the Spain and Austria switch- only affects PCA analysis
Rscript correct_BA_VI.R  ## mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata

## Calculate means for spring and fall for each group
Rscript make_mel_means_dpfreq_pops_allmean.R data/mel_freqdp_042016_Ne_fixed.Rdata data/means_dpfreq.paired20.allmean.Rdata data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt

## Filtering for SNPs that are polymorphic in ALL 40 samples (of the 20 paired sample dataset)
Rscript filter_polymorphic.R

## Calculate difference per season
Rscript make_mel_freqdiff.R ../data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_freqdiff_042016_$chrom.txt

## calculating mean dp per pop
Rscript calculate_mean_dp.R ../data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_meandp.$chrom.txt
Rscript calculate_mean_dp.R ../data/mel_freqdp_$chrom\_042016_fixed.Rdata ../data/mel_meandp_N.$chrom.txt

##  merge vcf file with SNPeff (SnpEff 4.2) (dm5.48.genome)
java -jar snpEff.jar download dm5.48
java -Xmx4G -jar snpEff.jar dm5.48.reference mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf

##  calculating Fst for samples that have a paired seasonal site (and climate data)
Rscript make_mel_Fst_Ne.R ../data/mel_freqdp_042016_Ne_fixed.Rdata mel_Fst_Ne.seasonalpaired.txt
Rscript make_mel_Fst_Ne_core20.R ../data/mel_freqdp_042016_Ne_fixed.Rdata $chrom ../data/mel_Fst_Ne.core20.sp_fall.txt

## Creating sets of matched controls per SNP
for chrom in 2L 2R 3L 3R; do
    Rscript match_snps_dp_ch_SNPbySNP_recomb_chromposFilter.R ../data/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.polymorphic.txt 100 ../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt
done
