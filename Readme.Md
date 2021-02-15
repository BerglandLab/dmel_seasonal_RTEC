# Code for Machado et al 2021
## Broad geographic sampling reveals predictable, pervasive, and strong seasonal adaptation in Drosophila

#### Master script referring to analysis scripts
#### Heather Machado

# 1. Process data and generate basic output
### Creates `mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata`, needed for downstream analysis
### Available on datadryad: https://doi.org/10.5061/dryad.4r7b826
### For use with downstream scripts, place in "data" folder
`create_input_files/code_create_input_files.sh`

# 2. PCA analysis
### `Rscript scripts_misc/PCA_Oct2020_BA_VIswitch.R`

# 3. Seasonal analyses: GLM
### `seasonal_glm/code_seasonal_glm.sh`

# 4. Seasonal analyses: RFM (rank Fisher's method)
### `seasonal_rfm/code_seasonal_rfm.sh`


# 5. Seasonal analyses: bayenv
### `seasonal_bayenv/code_seasonal_bayenv_Jan2020.sh`

# 6. Seasonal / latitudinal concordance
### uses file on datadryad (https://doi.org/10.5061/dryad.4r7b826):
### `bootstrap_fmean_dp.mel.medfreq01_RRgrt0.recRate.polymorphic.txt`
### 1. latitudinal glm regression
### `Rscript clinal_glm/code_clinal_glm.sh`
### 2. latitudinal bayenv analysis (including Figure S5)
### `clinal_bayenv/code_bayenv.sh`

# 7. Inversion frequencies
### `dmel_seasonal_RTEC/InversionFreqs`

# 8. Predicability Model
### `dmel_seasonal_RTEC/predictability_model/Readme.Md`

# 9. Comparison to Bergland et al 2014
### `dmel_seasonal_RTEC/comparision2014/Readme.Md`

# 10. Evaluate decontamination with D. simulans
### `dmel_seasonal_RTEC/sim_contamination_simulation_scripts/download_map.sh`

# 11. Effect of inflated sample size estimate on glm vs rfm (Figure S8)
### `samplesizeinflation_glm_vs_rfm/code_samplesizeinflation_glm_vs_rfm.sh`
