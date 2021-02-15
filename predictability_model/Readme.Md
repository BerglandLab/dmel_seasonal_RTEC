Scripts to evaluate the prediction score model. This analysis suggests that the thermal limit model

### back up data: `backupData.sh` should bundle all fo the input and output files for this analysis

### make `weather_and_pbs.fix.Rdata` object using `drosRTEC_revisions/predictability_model/dropPop.R`. That data-file has the components necessary to make Figure 3.


scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/drosRTEC/dropPop/weather_and_pbs.fix.Rdata ~/.
scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/drosRTEC/dropPop/mm.Rdata ~/.


1. Generate summary statistics
### make `weather_and_pbs.fix.Rdata` object using `drosRTEC_revisions/predictability_model/dropPop.R`. That data-file has the components necessary to make Figure 3.
  `drosRTEC_revisions/predictability_model/dropPop.R`

2. Predicability model
  ### This generates the distributions of R2 for the different potential explanatory variables, runs on whole genome.
  ### run permutations with (50 perms per job):
    ```bash
    wd=/scratch/aob2x/drosRTEC/drosRTEC_revisions
    sbatch --array=1-500 ${wd}/predictability_model/perm.runPerms.sh
    ```

4. Collect permutations
  ### collect permuations
    ```bash
    wd=/scratch/aob2x/drosRTEC/drosRTEC_revisions
    sbatch ${wd}/predictability_model/perm.collect.sh
    ```
