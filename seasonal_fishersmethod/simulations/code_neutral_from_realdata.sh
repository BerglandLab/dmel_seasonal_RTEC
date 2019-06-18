## March 2017
# HEM

## Simulating neutral data from the real data
# then perform fisher's exact tests, and run fisher's method

# Real paramenters
# 1) # individuals sampled: per population
# 2) # reads sampled: per SNP, per season
# 3) spring freq (use as true spring freq): per SNP

# Output:
# 1) file with simulated data
# 2) file with fisher's exact test
# 3) file with fisher's exact test rank
# 4) file with fisher's method (and enrichment over expected?)

# Simulating neutral sites:
#  Have already done this, so use set file of fisher's exact test results:
#       path: /scratch/users/hmachado/nescent_melCA/simulated_data/jamie_fisher_exact/
#       files: fisher_exactJ_pop1_neutral_1.5Msnps.txt (pops 1-25)

# Input file:
# 1 row per snp: spring alt reads, fall alt reads, spring total reads, fall total reads
# eg. 0 3 96 122


for i in {1..10}; do
    while read line; do
        Rscript simulated_datasets_from_realdata.R $i $line
    done < ../../data/paired_spring_fall_populations_noPA12.txt
done


#mkdir input_files
#mv neutralsample* input_files/


## Runnging fishers exact test
# USAGE: sbatch ../scripts/run_generic_program.sh python /scratch/users/hmachado/nescent_melCA/scripts/jamie_fisher.simulated_data.py inputfile.txt outputfile.txt
# eg. sbatch ../scripts/run_generic_program.sh python /scratch/users/hmachado/nescent_melCA/scripts/jamie_fisher.simulated_data.py datasets/pop$i\_neutral_1.5Msnps.txt jamie_fisher_exact/fisher_exactJ_pop$i\_neutral_1.5Msnps.txt


mkdir jamie_fisher_exact/

r=1
for i in neutralsample.rep$r.*; do
    python jamie_fisher.chrompos_data.py $i  jamie_fisher_exact/fisher_exactJ_$i
done

r=2
while read line; do
    python jamie_fisher.chrompos_data.py int_files/neutralsample.rep$r.$line.txt  jamie_fisher_exact/fisher_exactJ_neutralsample.rep$r.$line.txt
done < ../raw_data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt

ls jamie_fisher_exact/fisher_exactJ_neutralsample.rep$r.* > fisherlist.txt
Rscript merge_rewrite.R fisherlist.txt

cd jamie_fisher_exact/

for i in fisher_exactJ_neutralsample.rep$r.*.merged; do
   Rscript fishers_exact_pvalue_rank_Oct2016.R $i rank_$i
done

cd ../
while read i; do
    awk '{print $15}' jamie_fisher_exact/rank_fisher_exactJ_neutralsample.rep$r.$i.txt.merged > grt.$i.neutralsample.rep$r.tmp
    awk '{print $16}' jamie_fisher_exact/rank_fisher_exactJ_neutralsample.rep$r.$i.txt.merged > les.$i.neutralsample.rep$r.tmp
done < ../raw_data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt
paste grt.*.neutralsample.rep$r.tmp > jamie_fisher_exact/greater_rank_fisher_exactJ.neutralsample.rep$r.txt
paste les.*.neutralsample.rep$r.tmp > jamie_fisher_exact/lesser_rank_fisher_exactJ.neutralsample.rep$r.txt
rm grt*.neutralsample.rep$r.tmp
rm les*.neutralsample.rep$r.tmp
#rm jamie_fisher_exact/rank_fisher_exactJ_pop*.$outfile.txt
Rscript calculate_fishers_method_Jan2017.R jamie_fisher_exact/greater_rank_fisher_exactJ.neutralsample.rep$r.txt jamie_fisher_exact/lesser_rank_fisher_exactJ.neutralsample.rep$r.txt L_rank_fisher_exactJ.neutralsample.rep$r.txt


cd jamie_fisher_exact/
ls rank_fisher_exactJ_neutralsample.rep$r* > inlist.txt

module load R
R
inlist = read.table("inlist.txt", stringsAsFactors=FALSE)[,1]
df1 = na.omit(read.table(inlist[1], stringsAsFactors=FALSE, header=TRUE))[,c(1,2,13,14)]
colnames(df1) = c("chrom","pos",paste(inlist[1],".G",sep=""),paste(inlist[1],".L",sep=""))
for (i in 1:(length(inlist)-1)){
df2 = na.omit(read.table(inlist[i+1], stringsAsFactors=FALSE, header=TRUE)[,c(1,2,13,14)])
colnames(df2) = c("chrom","pos",paste(inlist[i+1],".G",sep=""),paste(inlist[i+1],".L",sep=""))
df3 = merge(df1, df2, by=c(1,2))
df1 = df3
}
write.table(df1, file="fisher_exactJ_all24pops_neutralsample.rep1.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)

df1 = na.omit(read.table(inlist[1], stringsAsFactors=FALSE, header=TRUE))[,c(1,2,15,16)]
colnames(df1) = c("chrom","pos",paste(inlist[1],".G",sep=""),paste(inlist[1],".L",sep=""))
for (i in 1:(length(inlist)-1)){
df2 = na.omit(read.table(inlist[i+1], stringsAsFactors=FALSE, header=TRUE)[,c(1,2,15,16)])
colnames(df2) = c("chrom","pos",paste(inlist[i+1],".G",sep=""),paste(inlist[i+1],".L",sep=""))
df3 = merge(df1, df2, by=c(1,2))
df1 = df3
}
write.table(df1, file="rank_fisher_exactJ_all24pops_neutralsample.rep1.txt", quote=FALSE, col.names=TRUE, row.names=FALSE)

grt = df1[,seq(from=3, to=ncol(df1), by=2)]
grtL = apply(grt, MARGIN=1, FUN=function(X) (-1)*sum(log(X)) )
les = df1[,seq(from=4, to=ncol(df1), by=2)]
lesL = apply(les, MARGIN=1, FUN=function(X) (-1)*sum(log(X)) )
L = data.frame(df1[,1:2], greater=grtL, lesser=lesL)
save(L, file="L_rank_fisher_exactJ_all24pops_neutralsample.rep1.Rdata")



############## May 5, 2017
######## testing the effect of changes in bin sizes for dp and freq
## running on new script that folds the allele freq spectrum
## run rank P-value
## also input dp and freq bin sizes
cd /scratch/users/hmachado/nescent_melCA/simulated_data/neutral_from_realdata/jamie_fisher_exact/
module load R
r=1

Ndpbin=10
Nfreqbin=10

for i in fisher_exactJ_neutralsample.rep$r.*.merged; do
    Rscript /scratch/users/hmachado/nescent_melCA/scripts/fishers_exact_pvalue_rank_folded_May2017.R $i rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_$i $Ndpbin $Nfreqbin
done


cd ../
while read i; do
    awk '{print $16}' jamie_fisher_exact/rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ_neutralsample.rep$r.$i.txt.merged > grt.dpbin$Ndpbin\_freqbin$Nfreqbin\_$i.neutralsample.rep$r.tmp
    awk '{print $17}' jamie_fisher_exact/rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ_neutralsample.rep$r.$i.txt.merged > les.dpbin$Ndpbin\_freqbin$Nfreqbin\_$i.neutralsample.rep$r.tmp
done < ../../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt
paste grt.dpbin$Ndpbin\_freqbin$Nfreqbin\_*.neutralsample.rep$r.tmp > jamie_fisher_exact/greater_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.neutralsample.rep$r.txt
paste les.dpbin$Ndpbin\_freqbin$Nfreqbin\_*.neutralsample.rep$r.tmp > jamie_fisher_exact/lesser_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.neutralsample.rep$r.txt
rm grt.dpbin$Ndpbin\_freqbin$Nfreqbin\_*.neutralsample.rep$r.tmp
rm les.dpbin$Ndpbin\_freqbin$Nfreqbin\_*.neutralsample.rep$r.tmp
#rm jamie_fisher_exact/rank_fisher_exactJ_pop*.$outfile.txt
Rscript ../../scripts/calculate_fishers_method_Jan2017.R jamie_fisher_exact/greater_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.neutralsample.rep$r.txt jamie_fisher_exact/lesser_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.neutralsample.rep$r.txt L_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.neutralsample.rep$r.txt



### new folded script
r=1
for i in fisher_exactJ_neutralsample.rep$r.*.merged; do
    Rscript /scratch/users/hmachado/nescent_melCA/fishers_method_Oct2016/fishers_exact_pvalue_rank_folded_Oct2016.R $i rank_folded_$i
done

cd ../
while read i; do
    awk '{print $16}' jamie_fisher_exact/rank_folded_fisher_exactJ_neutralsample.rep1.$i.txt.merged > grt.rank_folded_fisher_exactJ_neutralsample.rep1.$i.txt.merged.tmp
    awk '{print $17}' jamie_fisher_exact/rank_folded_fisher_exactJ_neutralsample.rep1.$i.txt.merged > les.rank_folded_fisher_exactJ_neutralsample.rep1.$i.txt.merged.tmp
done < ../../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt
paste grt.rank_folded_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp > jamie_fisher_exact/greater_rank_folded_fisher_exactJ_neutralsample.rep1.txt
paste les.rank_folded_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp > jamie_fisher_exact/lesser_rank_folded_fisher_exactJ_neutralsample.rep1.txt
rm grt.rank_folded_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp
rm les.rank_folded_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp
#rm jamie_fisher_exact/rank_fisher_exactJ_pop*.$outfile.txt
Rscript ../../scripts/calculate_fishers_method_Jan2017.R jamie_fisher_exact/greater_rank_folded_fisher_exactJ_neutralsample.rep1.txt jamie_fisher_exact/lesser_rank_folded_fisher_exactJ_neutralsample.rep1.txt L_rank_folded_freq5_fisher_exactJ_neutralsample.rep1.txt

Rscript ../../scripts/calculate_fishers_method_Jan2017.R jamie_fisher_exact/greater_rank_folded_fisher_exactJ_neutralsample.rep1.txt jamie_fisher_exact/lesser_rank_folded_fisher_exactJ_neutralsample.rep1.txt L_rank_folded_fisher_exactJ_neutralsample.rep1.txt


### no bins
r=1
for i in fisher_exactJ_neutralsample.rep$r.*.merged; do
    Rscript /scratch/users/hmachado/nescent_melCA/fishers_method_Oct2016/fishers_exact_pvalue_rank_nobin_Oct2016.R $i rank_nobin_$i
done

cd ../
while read i; do
awk '{print $15}' jamie_fisher_exact/rank_nobin_fisher_exactJ_neutralsample.rep1.$i.txt.merged > grt.rank_nobin_fisher_exactJ_neutralsample.rep1.$i.txt.merged.tmp
awk '{print $16}' jamie_fisher_exact/rank_nobin_fisher_exactJ_neutralsample.rep1.$i.txt.merged > les.rank_nobin_fisher_exactJ_neutralsample.rep1.$i.txt.merged.tmp
done < ../../data/paired_spring_fall_populations_noWI13PA9PA12PA14PA15.txt
paste grt.rank_nobin_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp > jamie_fisher_exact/greater_rank_nobin_fisher_exactJ_neutralsample.rep1.txt
paste les.rank_nobin_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp > jamie_fisher_exact/lesser_rank_nobin_fisher_exactJ_neutralsample.rep1.txt
rm grt.rank_nobin_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp
rm les.rank_nobin_fisher_exactJ_neutralsample.rep1.*.txt.merged.tmp
#rm jamie_fisher_exact/rank_fisher_exactJ_pop*.$outfile.txt
Rscript ../../scripts/calculate_fishers_method_Jan2017.R jamie_fisher_exact/greater_rank_nobin_fisher_exactJ_neutralsample.rep1.txt jamie_fisher_exact/lesser_rank_nobin_fisher_exactJ_neutralsample.rep1.txt L_rank_nobin_fisher_exactJ_neutralsample.rep1.txt



while read line; do
Rscript /scratch/users/hmachado/nescent_melCA/scripts/fishers_exact_pvalue_rank_folded_May2017.R fisher_exactJ_pop20/$line fisher_ex
actJ_pop20/rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_$line $Ndpbin $Nfreqbin
done < fisher_exactJ_pop20/inlist_pop20_merged.txt

while read line; do
awk '{print $16}' fisher_exactJ_pop20/rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_$line > grt.rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_$line.tm
p
awk '{print $17}' fisher_exactJ_pop20/rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_$line > les.rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_$line.tm
p
done < fisher_exactJ_pop20/inlist_pop20_merged.txt
paste grt.rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_*.tmp > fisher_exactJ_pop20/greater_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.me
rged.20pop.txt
paste les.rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_*.tmp > fisher_exactJ_pop20/lesser_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.mer
ged.20pop.txt
rm grt.rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_*.tmp
rm les.rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_*.tmp

mkdir -p results
Rscript ../../scripts/calculate_fishers_method_Jan2017.R fisher_exactJ_pop20/greater_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.
merged.20pop.txt fisher_exactJ_pop20/lesser_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.merged.20pop.txt L_rank_dpbin$Ndpbin\_fre
qbin$Nfreqbin\_fisher_exactJ.merged.20pop.txt
rm fisher_exactJ_pop20/greater_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.merged.20pop.txt
rm fisher_exactJ_pop20/lesser_rank_dpbin$Ndpbin\_freqbin$Nfreqbin\_fisher_exactJ.merged.20pop.txt

