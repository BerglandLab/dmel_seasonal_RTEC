Scripts to evaluate the prediction score model. This analysis suggests that the thermal limit model

### back up data: `backupData.sh` should bundle all fo the input and output files for this analysis

### make `weather_and_pbs.fix.Rdata` object using `drosRTEC_revisions/predictability_model/dropPop.R`. That data-file has the components necessary to make Figure 3.

###

1. Generate summary statistics
### make `weather_and_pbs.fix.Rdata` object using `drosRTEC_revisions/predictability_model/dropPop.R`. That data-file has the components necessary to make Figure 3.
  `drosRTEC_revisions/predictability_model/dropPop.R`

2. Calculate Leave-one-out concordance plot stratified on Inversion
  `drosRTEC_revisions/predictability_model/inversion_concordance.predicabilityModel.R`


2. Predicability model
  ### This generates the distributions of R2 for the different potential explanatory variables, runs on whole genome.
  ### run permutations with (50 perms per job):
    ```bash
    wd=/scratch/aob2x/drosRTEC/drosRTEC_revisions
    sbatch --array=1-500 ${wd}/predictability_model/perm.runPerms.sh
    ```


  ### collect permuations
    ```bash
    wd=/scratch/aob2x/drosRTEC/drosRTEC_revisions
    sbatch ${wd}/predictability_model/perm.collect.sh
    ```



".fix" fix the Ukraine - Spain swap, drop Spain from model b/c of unknown collection date
scp aob2x@rivanna.hpc.virginia.edu:~/mods.perm.Rdata ~/.
