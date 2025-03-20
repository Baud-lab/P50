library('beeswarm') # for dots

# Load genetic correlations
load('/users/abaud/abaud/P50_HSrats/output/VD/bivariate/all_VCs_corr_Ad1d2_zero_P50_Rn7_pruned_DGE.RData')

# Define cohort names as in paper
dict = c('NY'= 'NY',
         'MI'= 'MI',
         'TN_behavior'= 'TN1',
         'TN_breeder'= 'TN2')

# Set names for study pairs
all_VCs[,"study1"] = sapply(all_VCs[,"study1"], function(x) unname(dict[x]))
all_VCs[,"study2"] = sapply(all_VCs[,"study2"], function(x) unname(dict[x]))
all_VCs[,'study_pair'] = apply(all_VCs, 1, function(x) paste(sort(x[c("study1","study2")]),collapse="\n"))

# ordering so that decreasing mean of sample sizes
all_VCs[,'study_pair'] = factor(all_VCs$study_pair, 
                                levels = c("MI\nNY", "NY\nTN1", "MI\nTN1", "NY\nTN2", "MI\nTN2","TN1\nTN2"))


# Open pdf to save plot
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/comp_gen_corrs_across_cohorts.pdf", h= 6, w = 10.5)
par(mar=c(5.1,5.1,2.1,2.1))

# Set colours depending on p-value
cols = rep("#396C93", dim(all_VCs)[1])
cols[all_VCs$pv_chi2dof2 < 0.05] = "#E64B35FF"

# Plot boxplot
boxplot(all_VCs$corr_Ad1d2 ~ all_VCs$study_pair, 
        col="white",
        boxwex = 0.5,
        xaxt = "n", 
        xlab="", ylab="", cex.axis = 1.25,
        outline = F, 
        las = 1,  
        whisklty=3) 
# add ticks and ticks labels to x axis
axis(1, at = 1:6, labels = levels(all_VCs$study_pair), cex.axis = 1.2, padj=0.5)
# add x label
title(xlab = 'Cohort pairs', cex.lab = 1.4, line=3.8)
# add y label
title(ylab = "Same-taxon genetic correlations", cex.lab = 1.4, line=3.5)
# add dots
beeswarm(all_VCs$corr_Ad1d2 ~ all_VCs$study_pair,
         pch=16, cex=1.2, pwcol=cols, method = "swarm", add=T) 
dev.off()
