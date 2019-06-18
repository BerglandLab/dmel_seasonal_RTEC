########### April 2016

chrom=X
chrom=2R
chrom=2L
chrom=3L
chrom=3R

### After SNP calling, create RDataFile
module load R
#qsub run_make_melRDataFile.sh ../mapping2015/RTECrun2/06_variant_calling/mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf mel_freqdp_$chrom\_042016.Rdata
Rscript make_melRDataFile.R ../mapping2015/RTECrun2/06_variant_calling/mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf mel_freqdp_$chrom\_042016.Rdata

### Calculate Ne
Rscript Ne.R mel_freqdp_$chrom\_042016.Rdata mel_nescent2_samplesize.txt mel_freqdp_$chrom\_042016_Ne.Rdata auto

### Combine the wild and lines of BHM
Rscript combine_BHM.R mel_freqdp_$chrom\_042016_Ne.Rdata mel_freqdp_$chrom\_042016_Ne_fixed.Rdata

### Create one file for all chromosomes
Rscript combine_freqdp.R data/mel_freqdp 042016_fixed.Rdata  ## created mel_freqdp_042016_fixed.Rdata

### Calculate means for spring and fall for each group
#qsub run_make_mel_means_dpfreq.sh mel_freqdp_$chrom\_042016_Ne_fixed.Rdata means_dpfreq.$chrom.Rdata
Rscript scripts/make_mel_means_dpfreq.R mel_freqdp_$chrom\_042016_Ne_fixed.Rdata means_dpfreq.$chrom.Rdata
Rscript scripts/make_mel_means_dpfreq_otherFreqBin.R mel_freqdp_$chrom\_042016_Ne_fixed.Rdata means_dpfreq_otherFreqBin.$chrom.Rdata
Rscript scripts/make_mel_means_dpfreq_pops.R data/mel_freqdp_042016_Ne_fixed.Rdata data/means_dpfreq.paired20.Rdata data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt
Rscript scripts/make_mel_means_dpfreq_pops_allmean.R data/mel_freqdp_042016_Ne_fixed.Rdata data/means_dpfreq.paired20.allmean.Rdata data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt
Rscript scripts/make_mel_means_dpfreq_pops_allmean_4switch.R data/mel_freqdp_042016_Ne_fixed.Rdata data/means_dpfreq.paired20.allmean_4switch.Rdata data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt


########### Filtering for SNPs that are polymorphic in ALL 40 samples (of the 20 paired sample dataset)
# in R
load("data/mel_freqdp_auto_042016_fixed.Rdata")
pops24 = read.table("data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]
popID = popinfo[,1]
S = popinfo[,5]
PY = popinfo[,6]
focal = which( PY %in% pops24 & (S=="s" | S=="f") )
freqfocal = freq[, focal]
dpfocal = dp[, focal] ### only looking at fall to frost
bothA=cbind(info[,1:2], freqfocal, dpfocal)
filter=read.table("data/chrom_pos_medfreq01_RRgrt0.txt", stringsAsFactors=FALSE)
both=merge(bothA, filter, by=c(1,2))
t1 = apply(both[,3:42], MARGIN=1, FUN=function(X) sum(X<=0 | X>=1, na.rm=TRUE) )  # keep any
# ignore the NAs
#> mean(t1>0)
#[1] 0.5606287   ## only 44% don't have a 0 or 1 allele freq
both2 = both[ t1==0 , ]   # 774,841
write.table(both2[,1:2], file="data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt", col.names=FALSE, quote=FALSE, row.names=FALSE)


### Run seasonal regression
Rscript seasonal_analysis_all_paired20_2sample_caF_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_paired20_2sample_caF_popyear.f_s.glm
Rscript seasonal_analysis_all_paired20_2sample_caF_popyear_4switch.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_paired20_2sample_caF_popyear_4switch.f_s.glm
Rscript seasonal_analysis_ca_f_s_popyear_paired.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_ca_popyear_paired.f_s.glm
Rscript seasonal_analysis_ca_f_s_popyear_paired_switch.R data/mel_freqdp_042016_Ne_fixed.Rdata glm/mel_ca_popyear_paired_switch.f_s.glm
Rscript seasonal_analysis_eur_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_eur_popyear_paired.f_s.glm


## clinal regression
Rscript clinal_analysis_uniquepops_springPA_noMA.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_clinal_uniquepops_springPA_noMA.glm


## seasonal regression without clinal populations (PA)
Rscript seasonal_analysis_all_nonclinal_paired20_2sample_caF_popyear.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_nonclinal_paired20_2sample_caF_popyear.f_s.glm
Rscript seasonal_analysis_all_paired20_2sample_caF_popyear_4switch_noPA.R data/mel_freqdp_042016_Ne_fixed.Rdata mel_all_caF_switch_nonclinal_uniquepopsSpringPA_noMA.glm


## extra (non-primary) analyses
# seasonal_analysis_all_paired20_2sample_caF_popyear_4switch_noBerlgand2014.R


## Seasonal permutations
for j in {1..200}; do
    for chrom in 2L 2R 3L 3R; do
        #sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/permute_seasonal_analysis_all_caF.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata permute$j\_mel_all_caF.seas_pop_year.f_s.$chrom.glm
        #sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/permute_seasonal_analysis_other.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata permute$j\_mel_other.seas_pop_year.f_s.$chrom.glm
        #sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/permute_seasonal_analysis_pa_nosum.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata permute$j\_mel_pa.seas_pop_year.f_s.$chrom.glm
        #sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/permute_seasonal_analysis_ca.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata permute$j\_mel_ca.seas_pop_year.f_s.$chrom.glm
        #sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/permute_seasonal_analysis_other_noyear.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata permute$j\_mel_other.seas.f_s.$chrom.glm # actually does have pop
        sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/permute_seasonal_analysis_all_caF_noyear.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata permute$j\_mel_all_caF.seas_pop.f_s.$chrom.glm
    done
done


for p in {4..50}; do
    for i in 2L 2R 3L 3R; do
    #awk -v chrom=$i '{if ($1==chrom) print $2}' data/chrom_pos_medfreq01_RRgrt0.txt | sort -k1,1b - > $i.tmp
        awk '{print $2, $3, $4}' permute$p\_mel_all_caF.seas_pop.f_s.$i.glm | sort -k1,1b - | join $i.tmp - | awk -v chrom=$i '{print chrom, $0}' -  >> glm_permutations/permute$p\_mel_all_caF.seas_pop.f_s.glm.medfreq01_RRgrt0.txt
        #awk '{print $2, $4}' permute$p\_mel_other.seas_pop_year.f_s.$i.glm | sort -k1,1b - | join $i.tmp - | awk -v chrom=$i '{print chrom, $0}' -  >> permute$p\_mel_other.seas_pop_year.f_s.glm.medfreq01.txt
        #awk '{print $2, $4}' permute$p\_mel_pa.seas_pop_year.f_s.$i.glm | sort -k1,1b - | join $i.tmp - | awk -v chrom=$i '{print chrom, $0}' -  >> permute$p\_mel_pa.seas_pop_year.f_s.glm.medfreq01.txt
        #awk '{print $2, $4}' permute$p\_mel_ca.seas_pop_year.f_s.$i.glm | sort -k1,1b - | join $i.tmp - | awk -v chrom=$i '{print chrom, $0}' -  >> permute$p\_mel_ca.seas_pop_year.f_s.glm.medfreq01.txt
        #awk '{print $2, $3, $4}' permute$p\_mel_other.seas.f_s.$i.glm | sort -k1,1b - | join $i.tmp - | awk -v chrom=$i '{print chrom, $0}' -  >> permute$p\_mel_other.seas.f_s.glm.medfreq01.txt
    done
done

Rscript code_combine_glm_permute.R
Rscript ../scripts/analyze_permute_seasonal_na.R 200 perm200_nescentcompare_seasonal_na_concordance.Rdata




### Calculate difference per season
Rscript scripts/make_mel_freqdiff.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_freqdiff_042016_$chrom.txt
Rscript scripts/make_mel_meanfreqdiff_paired20_4switch.R data/mel_freqdiff_paired20_042016.Rdata data/mel_meanfreqdiff_paired20_4switch_042016.txt



### Run bootstrap matching of seasonal SNPs
# For block bootstrap, sample your significant SNPs by region (requiring the bootstrap to be in the region).
# for CA analysis, use SNP matches for PA population
#### OR, normal bootstrap sampling (done below), with all SNPs. Have to use version of R 3.2.2
# module load R/3.2.2
# Rscript match_snps_dp_ch_SNPbySNP_recomb.R 1 means_dpfreq.2L.Rdata bootstrap_fmean_dp.mel_2L.txt
#qsub run_seasonal_analysis.sh match_snps_dp_ch_SNPbySNP_recomb.R means_dpfreq.$chrom.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.txt 100
#for chrom in 2L 2R 3L 3R X; do
#sbatch scripts/run_seasonal_analysis.sh scripts/match_snps_dp_ch_SNPbySNP_recomb.R data/means_dpfreq_otherFreqBin.$chrom.Rdata bootstrap/bootstrap_otherfmean_dp.mel_$chrom.txt 100
#done
for chrom in 2L 2R 3L 3R; do
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/match_snps_dp_ch_SNPbySNP_recomb.R recombinationRate/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.txt 100
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/match_snps_dp_ch_SNPbySNP_recomb.R data/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.noNA_clinaluniquepopsPA_seasonalnonclinalPA.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.noNA_clinaluniquepopsPA_seasonalnonclinalPA.txt 100
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/match_snps_dp_ch_SNPbySNP_recomb.R recombinationRate/means_dpfreq_allpop_MAF10.$chrom.recRate.Rdata bootstrap/bootstrap_allmean_dp.mel_$chrom.meanfreq10.recRate.txt 100
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/match_snps_dp_ch_SNPbySNP_recomb_chromposFilter.R recombinationRate/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.fisher_exactJ.merged.20pop.txt 100 fishers_method_Oct2016/jamie_fisher_method/results/chrom_pos_L_rank_fisher_exactJ.merged.20pop.txt
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/match_snps_dp_ch_SNPbySNP_recomb_chromposFilter.R recombinationRate/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.polymorphic.txt 100 data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt
done


for chrom in 2L 2R 3L 3R; do
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/match_snps_dp_ch_SNPbySNP_recomb_chromposFilter.R recombinationRate/means_dpfreq.$chrom.medfreq01_RRgrt0.recRate.Rdata bootstrap/bootstrap_fmean_dp.mel_$chrom.medfreq01_RRgrt0.recRate.melsimSameAlleles 100 simulans/melsimSameAlleles_melchrompos.txt
done

## merge file with SNPeff (SnpEff 4.2) (dm5.48.genome) (DONE IN BIO-DAP20)
java -jar snpEff.jar download dm5.48
java -Xmx4G -jar snpEff.jar dm5.48.reference /hsgs/projects/petrov/hmachado/mapping2015/RTECrun2/06_variant_calling/mel_$chrom\_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.freq.vcf




## calculating Fst for samples that have a paired seasonal site (and climate data)
qsub scripts/run_seasonal_analysis.sh scripts/make_mel_Fst.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_Fst.seasonalpaired.$chrom.txt
qsub scripts/run_seasonal_analysis.sh scripts/make_mel_Fst_ca_sp_fall.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_Fst.ca.sp_fall.$chrom.txt
qsub scripts/run_seasonal_analysis.sh scripts/make_mel_Fst_Ne.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_Fst_Ne.seasonalpaired.$chrom.txt
qsub scripts/run_seasonal_analysis.sh scripts/make_mel_Fst_Ne.ca_sp_fall.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_Fst_Ne.ca.sp_fall.$chrom.txt

for chrom in 2L 2R 3L 3R; do
    sbatch scripts/run_seasonal_analysis.sherlock.sh scripts/make_mel_Fst_Ne_core20.R data/mel_freqdp_042016_Ne_fixed.Rdata $chrom data/mel_Fst_Ne.core20.sp_fall.$chrom.txt
done


## caculating mean dp per pop
Rscript calculate_mean_dp.R data/mel_freqdp_$chrom\_042016_Ne_fixed.Rdata mel_meandp.$chrom.txt
Rscript calculate_mean_dp.R data/mel_freqdp_$chrom\_042016_fixed.Rdata data/mel_meandp_N.$chrom.txt




############ Output new vcf
cd data
sbatch ../scripts/run_Rscript_v3.3.0.sh format_vcf.R


##### Fisher's Method seasonal analysis

## Observed data
seasonal_fishersmethod/code_fishers_exactJ.sh


## Simulations
seasonal_fishersmethod/simulations/code_neutral_from_realdata.sh


########### PCA analysis
Rscript scripts_misc/PCA_clean.R
