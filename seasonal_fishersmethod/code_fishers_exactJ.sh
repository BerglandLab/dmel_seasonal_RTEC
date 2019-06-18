### Oct 31, 2016
# HEM

### Performing seasonal rank fishers method analysis on observed data

while read line; do
    Rscript make_inputfiles.R ../data/mel_freqdp_042016_Ne_fixed.Rdata $line input_files/$line.alt_dp
done < ../data/paired_spring_fall_populations.txt
ls *.alt_dp.txt > input_files/inlist.txt


## run Jamie's fisher exact
while read line; do
    python jamie_fisher.chrompos_data.py input_files/$line fisher_exactJ/fisher_exactJ_$line
done < input_files/inlist.txt


###### Now doing for 20 populations
mkdir fisher_exactJ_pop20/
cd fisher_exactJ_pop20/
Rscript ../merge_rewrite.R ../inlist_pop20.txt
cd ../
while read line; do
    Rscript fishers_exact_pvalue_rank_Oct2016.R fisher_exactJ_pop20/$line.merged fisher_exactJ_pop20/rank_$line.merged
done < inlist_pop20.txt

while read line; do
    awk '{print $15}' fisher_exactJ_pop20/rank_$line.merged > grt.$line.merged.tmp
    awk '{print $16}' fisher_exactJ_pop20/rank_$line.merged > les.$line.merged.tmp
done < fisher_exactJ_pop20/inlist_pop20.txt
paste grt.*.merged.tmp > fisher_exactJ_pop20/greater_rank_fisher_exactJ.merged.20pop.txt
paste les.*.merged.tmp > fisher_exactJ_pop20/lesser_rank_fisher_exactJ.merged.20pop.txt
rm grt*.tmp
rm les*.tmp


## Perform fishers method
## Calculate expected vs observed
mkdir -p results
Rscript calculate_fishers_method_Jan2017.R fisher_exactJ_pop20/greater_rank_fisher_exactJ.merged.20pop.txt fisher_exactJ_pop20/lesser_rank_fisher_exactJ.merged.20pop.txt L_rank_fisher_exactJ.merged.20pop.txt


## subsampling (have to use the ranked p-values, in order to reflect the simulated datasets)
# keeping chrom_pos
for j in 1 5 10; do
    Rscript block_sampling_1perBin_grt_les_chrom_pos.R fisher_exactJ_pop20/chrom_pos_greater_rank_fisher_exactJ.merged.20pop.txt fisher_exactJ_pop20/chrom_pos_lesser_rank_fisher_exactJ.merged.20pop.txt fisher_exactJ_pop20/greater_rank_fisher_exactJ.merged.20pop.chrom_pos.sample$j\K.txt fisher_exactJ_pop20/lesser_rank_fisher_exactJ.merged.20pop.chrom_pos.sample$j\K.txt $j\000 100
    for i in {1..100}; do
        Rscript ../../scripts/calculate_fishers_method_chrom_pos_Jan2017.R fisher_exactJ_pop20/greater_rank_fisher_exactJ.merged.20pop.chrom_pos.sample$j\K.txt$i fisher_exactJ_pop20/lesser_rank_fisher_exactJ.merged.20pop.chrom_pos.sample$j\K.txt$i L_rank_fisher_exactJ.merged.20pop.chrom_pos.sample$j\K.txt$i
    done
done


## Comparing those block sampling with simulated datasets
Rscript ../../scripts/average_samples_simulated_observed_realdata_20pop_dist2_normalized_July2017.R results/enrich_L_rank_fisher_exactJ.merged.20pop.chrom_pos.sample1K.txt 100 ../../simulated_data/one_script_from_realdata inlist_rep1_pop20.txt 20pop_sample1K
Rscript ../../scripts/average_samples_simulated_observed_realdata_20pop_dist2_normalized_July2017.R results/enrich_L_rank_fisher_exactJ.merged.20pop.chrom_pos.sample5K.txt 100 ../../simulated_data/one_script_from_realdata inlist_rep1_pop20.txt 20pop_sample5K
Rscript ../../scripts/average_samples_simulated_observed_realdata_20pop_dist2_normalized_July2017.R results/enrich_L_rank_fisher_exactJ.merged.20pop.chrom_pos.sample10K.txt 100 ../../simulated_data/one_script_from_realdata inlist_rep1_pop20.txt 20pop_sample10K



###### calculating the mean spring and mean fall allele freq for the 20 pops, and the mean sp/fall diff
R
inpops = read.table("../../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt", stringsAsFactors=FALSE)[,1]
in1 = read.table("input_files/co_13.alt_dp.txt", stringsAsFactors=FALSE)[,1:2]

for (i in 1:length(inpops)){
    focal = na.omit(read.table(paste("input_files/", inpops[i], ".alt_dp.txt", sep=""), stringsAsFactors=FALSE))
    focal$freqS = focal[,3]/focal[,5]
    focal$freqF = focal[,4]/focal[,6]
    focal$SF = focal$freqS - focal$freqF
    out1 = focal[,c(1,2,7,8,9)]
    out2 = merge(in1, out1, by=c(1,2))
    in1 = out2
}

freqS = in1[,seq(from=3, to=ncol(in1), by=3)]
freqF = in1[,seq(from=4, to=ncol(in1), by=3)]
SF = in1[,seq(from=5, to=ncol(in1), by=3)]
meanfreqS = apply(freqS, MARGIN=1, FUN=mean, na.rm=TRUE)
meanfreqF = apply(freqF, MARGIN=1, FUN=mean, na.rm=TRUE)
meanSF = apply(SF, MARGIN=1, FUN=mean, na.rm=TRUE)
outtab = data.frame(chrom=in1[,1], pos=in1[,2], meanfreqS, meanfreqF, meanSF)
write.table(outtab, "meanfreq_SF_pop20.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)



