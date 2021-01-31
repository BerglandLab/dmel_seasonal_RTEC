#!/bin/bash


i=$1
awk -F' ' 'NR==FNR{e[$1$2]=1;next};e[$1$2]{print $4}' ../data/chrom_pos_polymorphic_medfreq01_RRgrt0.txt results/permute$i\_mel_all_paired20_2sample_caF_popyear.glm | grep -v "NA" | sort -g > perm_paste_tmp/$i.txt
