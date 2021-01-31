#!/bin/bash

tar -czvf /mnt/pricey_4/prediction_model_data.tar.gz \
/mnt/spicy_4/pigmentation/inputData/dat.Rdata \
/mnt/pricey_1/dropPop/chrom_pos_polymorphic_medfreq01_RRgrt0.txt \
/mnt/pricey_1/dropPop/popSS.Rdata \
/mnt/pricey_1/dropPop/popTranslate.delim \
/mnt/pricey_1/dropPop/popSS.translate.Rdata \
/mnt/pricey_1/dropPop/fisher_* \
/mnt/pricey_1/dropPop/mel_all* \
/mnt/pricey_1/dropPop/mel_clinal_uniquepops.glm \
/mnt/pricey_1/dropPop/dropPop.feather \
/mnt/pricey_1/dropPop/concordance.Rdata \
/mnt/pricey_1/dropPop/19_drop_1_concordance_data.Rdata \
/mnt/pricey_1/dropPop/concordance_o19_chr.Rdata \
/mnt/pricey_1/dropPop/19_drop_1_enrichment.Rdata \
/mnt/pricey_1/dropPop/newPopPairs.csv \
/mnt/pricey_1/dropPop/Z_*.coef_minp2.txt \
/mnt/pricey_1/dropPop/fisher_exactJ_PA_12.coef_minp2.txt \
/mnt/pricey_1/dropPop/fisher_exactJ_PA_9.coef_minp2.txt \
/mnt/pricey_1/dropPop/fisher_exactJ_PA_14.coef_minp2.txt \
/mnt/pricey_1/dropPop/fisher_exactJ_WI_13.coef_minp2.txt \
/mnt/pricey_1/dropPop/fisher_exactJ_PA_15.coef_minp2.txt \
/mnt/pricey_1/dropPop/concordance_full20_add1.orig.swap.FET.Rdata

cp /mnt/pricey_4/prediction_model_data.tar.gz /mnt/sammas_storage/bergland-lab/alan/.

#/mnt/icy_3/nescent/data/all_popinfo.csv

scp /Users/alanbergland/Documents/GitHub/drosRTEC_revisions/predictability_model/backupData.sh \
bergland@bergland-lab.bio.virginia.edu:~/.


cp /nv/vol186/bergland-lab/alan/prediction_model_data.tar.gz /project/berglandlab/alan/drosRTEC/.
