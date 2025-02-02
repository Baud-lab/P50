load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')
# collapsed_clr_counts is post CLR post collapsing to higher level taxa
collapsed_clr_counts = collapsed_clr_counts[grep('g__', rownames(collapsed_clr_counts)),]
pca = prcomp(t(collapsed_clr_counts))
vars = round(pca$sdev^2 / (sum(pca$sdev^2))*100, digits = 1)

load('/nfs/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
motch = match(colnames(collapsed_clr_counts), metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

cols = rep('white', dim(collapsed_clr_counts)[2])
cols[metadata$study == 'NY'] = 'blue'
cols[metadata$study == 'MI'] = 'green'
cols[metadata$study == 'TN_behavior'] = 'red'
cols[metadata$study == 'TN_breeder'] = 'orange'

pdf('/nfs/users/abaud/abaud/P50_HSrats/plots/PCA_all_genera.pdf')
plot(pca$x[,1], pca$x[,2], xlab = 'MDS1 (24%)', ylab = 'MDS2 (8%)', col = cols, pch = 16, las = 1, cex = 0.7)
legend('topright', fill = c('blue','green','red','orange'), legend = c('NY','MI','TN1','TN2'), border = NA, box.col= NA)
dev.off()

