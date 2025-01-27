load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData')
pca = prcomp(t(clr_counts))
vars = round(pca$sdev^2 / (sum(pca$sdev^2))*100, digits = 1)

load('/nfs/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
motch = match(colnames(clr_counts), metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]
cols = rep('white', dim(clr_counts)[2])
cols[metadata$study == 'MI'] = 'green'
cols[metadata$study == 'NY'] = 'blue'
cols[metadata$study == 'TN_behavior'] = 'red'
cols[metadata$study == 'TN_breeder'] = 'orange'

pdf('/nfs/users/abaud/abaud/P50_HSrats/plots/PCA_alltaxa.pdf')
plot(pca$x[,1], pca$x[,2], xlab = 'MDS1 (4.4%)', ylab = 'MDS2 (3.4%)', col = cols, pch = 16, las = 1, cex = 0.7)
legend('topright', fill = c('blue','green','red','orange'), legend = c('NY','MI','TN behaviour','TN breeder'), title = 'Cohort:', border = NA, box.col= NA)
dev.off()
