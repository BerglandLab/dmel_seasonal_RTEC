#### Code for Machado et al 2021
# Broad geographic sampling reveals predictable, pervasive, and strong seasonal adaptation in Drosophila

#### Master script referring to analysis scripts
#### Heather Machado


###### Create files used downstream (uses vcf- available on datadryad)
create_input_files/code_create_input_files.sh
# creates mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata, needed for downstream analysis (also available on datadryad: https://doi.org/10.5061/dryad.4r7b826)
#   For use with downstream scripts, place in "data" folder

###### PCA analysis
Rscript scripts_misc/PCA_Oct2020_BA_VIswitch.R



###### Seasonal analyses: GLM
seasonal_glm/code_seasonal_glm.sh



###### Seasonal analyses: RFM (rank Fisher's method)
seasonal_rfm/code_seasonal_rfm.sh



###### Seasonal analyses: bayenv
seasonal_bayenv/code_seasonal_bayenv_Jan2020.sh



##### Seasonal / latitudinal concordance
# uses file on datadryad (https://doi.org/10.5061/dryad.4r7b826):
#          bootstrap_fmean_dp.mel.medfreq01_RRgrt0.recRate.polymorphic.txt
# latitudinal glm regression
Rscript clinal_glm/code_clinal_glm.sh
# latitudinal bayenv analysis (including Figure S5)
clinal_bayenv/code_bayenv.sh



##### Effect of inflated sample size estimate on glm vs rfm (Figure S8)
samplesizeinflation_glm_vs_rfm/code_samplesizeinflation_glm_vs_rfm.sh
