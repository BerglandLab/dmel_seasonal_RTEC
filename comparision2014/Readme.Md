1. get old data onto Rivanna

  cp \
  /Users/alanbergland/Documents/work/Projects/2005-2019/2011_6Dimensions/Paper 1 - seasonal */MS /PLoS Genet submission/  Resubmission 2/DataDryad data/6d_v7.3_output.vcf.gz \
  ~/.

  scp ~/6d_v7.3_output.vcf.gz \
  aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/.

2. Extract data from vcf file. We want to get the sites with seasonal FDR <= 0.3. Should be ~1700 of them.
  `drosRTEC_revisions/comparision_2014/extract_data.R`

3. Do enrichment test
  4
