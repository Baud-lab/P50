#library(devtools)
#install_github("jennylsl/fast.adonis")
library(fast.adonis)
library(plyr)

#issues to install packages on the cluster so analysis ran on laptop
load('~/Downloads/collapsed_full_biomt_collapsed_clr_counts.RData')
# collapsed_clr_counts is post CLR post collapsing to higher level taxa
save_collapsed_clr_counts = t(collapsed_clr_counts)

load('~/Downloads/metadata_augmented_16S_metabo_deblur.RData')
motch = match(rownames(save_collapsed_clr_counts), metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

cc = complete.cases(metadata$study, metadata$sex)
save_collapsed_clr_counts = save_collapsed_clr_counts[cc,]
metadata = metadata[cc,]

tax_level = "f__"
collapsed_clr_counts = save_collapsed_clr_counts[,grep(tax_level, colnames(save_collapsed_clr_counts))]
D = dist(collapsed_clr_counts)
A <- -0.5*as.matrix(D)^2

asdf <- fast.adonis(A ~ sex + study, data=metadata, parallel=1, permutations=999, boot.times = 100)
print(asdf[['aov.tab']])

