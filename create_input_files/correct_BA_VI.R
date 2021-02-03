### Jan 2021

# This corrects the labels, specifically:
# mypop[P1=="BA"] = "AT_gr"
# mypop[P1=="VI"] = "ES_ba"
## See manuscript documenting discovery of this switch:
# Nunez et al 2019: "Updating the metadata of four misidentified samples in the DrosRTEC dataset"


## Load in data
filein = "mel_freqdp_042016_Ne_fixed.Rdata"
fileout = "mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata"
load(filein)

## identify samples that were switched (Spain and Austria)
popinfo2 = data.frame(popinfo)
subset(popinfo2, P %in% c("BA","VI"))
# pops  P  Y     R S   P.Y sf_pair ffr_pair
# X.17 melBA_102012_FAT BA 12 other f BA_12    TRUE    FALSE
# X.18 melBA_072012_SPT BA 12 other s BA_12    TRUE    FALSE
# X.39 melVI_102012_FAT VI 12 other f VI_12    TRUE    FALSE
# X.40 melVI_082012_SPT VI 12 other s VI_12    TRUE    FALSE

## make the corrected objects
dp_new = dp
freq_new = freq
popinfo_new = popinfo

dp_new[,18] = dp[,40]
freq_new[,18] = freq[,40]
popinfo_new[18,] = popinfo[40,]

dp_new[,19] = dp[,41]
freq_new[,19] = freq[,41]
popinfo_new[19,] = popinfo[41,]

dp_new[,40] = dp[,18]
freq_new[,40] = freq[,18]
popinfo_new[40,] = popinfo[18,]

dp_new[,41] = dp[,19]
freq_new[,41] = freq[,19]
popinfo_new[41,] = popinfo[19,]

popinfo[c(18,19,40,41),]
popinfo_new[c(18,19,40,41),]

dp = dp_new
freq = freq_new
popinfo = popinfo_new


save(dp, freq, info, popinfo, file=fileout)
