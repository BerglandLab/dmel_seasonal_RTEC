### Oct 31, 2016
# HEM

### Performing seasonal rank fishers method analysis on observed data

while read line; do
    Rscript make_inputfiles.R ../data/mel_freqdp_042016_Ne_fixed.Rdata $line input_files/$line.alt_dp
done < ../data/paired_spring_fall_populations.txt
ls *.alt_dp.txt > input_files/inlist.txt


## run fisher exact
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

