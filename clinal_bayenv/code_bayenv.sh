#!/bin/sh

#  code_bayenv_Feb2019.sh
#  
#
#  Created by Heather Machado on 2/18/19.
#


## Populations used for clinal regression
# "melPA_72011_SPT" | popID=="melSC_072010_SPT" | popID=="melGA_072008_SPT" | popID=="melFL_072010_SPT"
# lat = c(25.5, 30.99, 39.88, 33.39)
#"melFL_072010_SPT": 25.5
#"melGA_072008_SPT": 30.99
#"melPA_72011_SPT": 39.88
#"melSC_072010_SPT": 33.39

## Data
# data/mel_freqdp_X_042016_Ne_fixed.Rdata

## Convert data:
rdata2bayenvformat.R

## Estimate covariance matrix
# calculate matrix for 10K sites (pooled mode)
/Users/hm8/Software/bayenv2/bayenv2 -i clinal_snpsfile_Feb2019_sample10K.txt -s clinal_samplefile_Feb2019.txt -p 4 -k 100000 -x -r 63479 > matrix_10Kbsub100Kiter.out
tail -n5 matrix_10Kbsub100Kiter.out | head -n4 - > covmatrix_10Kbsub100Kiter.txt

## Running the env test, individually per snp
mkdir -p snpresults
cd snpresults
## total of 1939173 snps (cannot run for all snps, as some have NAs or no alt or no ref reads)
for i in {1..1939173}; do
    mkdir focalsnp$i
    cd focalsnp$i
    ln -s ../snpfiles/clinal_snpsfile_Feb2019_snp$i.txt ./
    ln -s ../clinal_samplefile_Feb2019.txt ./
    ln -s ../clinal_environfile_Feb2019.txt ./
    ln -s ../covmatrix_10Kbsub100Kiter.txt ./
    bayenv2 -i clinal_snpsfile_Feb2019_snp$i.txt -s clinal_samplefile_Feb2019.txt -e clinal_environfile_Feb2019.txt -p 4 -n 1 -k 100000 -x -t -r 73480 -m covmatrix_10Kbsub100Kiter.txt
    cd ../
done
cat focalsnp*/environ_corr.clinal_environfile_Feb2019.txt > environ_corr.clinal_environfile_allnon0snps_Feb2019.txt

# Calculates seasonal and latitudinal concordance. Plotting Figure S5
# input:
#    environ_corr.clinal_environfile_allnon0snps_Feb2019.txt
#    clinal_allsnps_non0_unformatted_Feb2019.txt
Rscript analyze_bayenv_Nov2020_suppfigure.R
