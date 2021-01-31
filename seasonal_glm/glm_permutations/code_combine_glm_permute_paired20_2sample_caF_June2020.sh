## June 2020

mkdir perm_paste_tmp
for i in {1..500}; do
    bsub -J combperm -n1 -R "select[mem>4000]" -R "rusage[mem=4000]" -M 4000 -R "span[hosts=1]" -q normal -e error.%J -o output.%J ./combineperms.sh $i
done

## each analysis has a few NA's. Keep best 774637 p-values per permutation (since only counting low p-values)
paste perm_paste_tmp/*.txt > permute_1_500_polymorphic_mel_all_paired20_2sample_caF_popyear.glm.tmp
head -n 774637 permute_1_500_polymorphic_mel_all_paired20_2sample_caF_popyear.glm.tmp > permute_1_500_polymorphic_mel_all_paired20_2sample_caF_popyear.glm
rm permute_1_500_polymorphic_mel_all_paired20_2sample_caF_popyear.glm.tmp

for i in {1..500}; do
    cut -f $i permute_1_500_polymorphic_mel_all_paired20_2sample_caF_popyear.glm | awk '{if ($1 <0.003894319) print; else exit}' - | wc -l >> perm500_n_less_p00389.txt
done
