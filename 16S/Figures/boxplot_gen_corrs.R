#load genetic correlations
load('/users/abaud/abaud/P50_HSrats/output/VD/bivariate/all_VCs_corr_Ad1d2_zero_P50_Rn7_pruned_DGE.RData')

#create a variable indicating the pair of cohorts considered
studies = c('MI','NY','TN_behavior','TN_breeder')
for (i in 1:dim(all_VCs)[1]) {
  all_VCs[i,'study_pair'] = paste(sort(c(all_VCs[i,'study1'],all_VCs[i,'study2'])),collapse=' \n')
}

#calculcate median gen corr for each pair of cohorts to order the boxplot in increasing order
medians = sort(sapply(split(all_VCs$corr_Ad1d2, all_VCs$study_pair), median), decreasing = T)

#plot
all_VCs$study_pair = factor(all_VCs$study_pair, levels = names(medians))
pdf('/users/abaud/abaud/P50_HSrats/plots/comp_gen_corrs_across_cohorts.pdf', width = 20, height = 10)
par(mar = c(12,10,3,1))
boxplot(all_VCs$corr_Ad1d2 ~ all_VCs$study_pair, ylab = '', xlab = '',  boxwex=0.5, axes = FALSE)
axis(1, at = 1:6, labels = levels(all_VCs$study_pair), cex.axis = 2, padj = 1)
axis(2, padj = -1 , cex.axis = 3, las = 1)
mtext(side = 1, at = 3.5, text = 'Cohort pairs', padj = 4.2, cex = 3)
mtext(side = 2, at = 0.25, text = 'Same-taxon genetic correlations', padj = -3.5, cex = 3)
dev.off()

