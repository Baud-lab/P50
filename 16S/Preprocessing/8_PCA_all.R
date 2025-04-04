load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')
# collapsed_clr_counts is post CLR post collapsing to higher level taxa
save_collapsed_clr_counts = collapsed_clr_counts

pdf('/nfs/users/abaud/abaud/P50_HSrats/plots/PCA_paper.pdf')

for (tax_level in c("p__","c__","o__","f__","g__")) {

	collapsed_clr_counts = save_collapsed_clr_counts[grep(tax_level, rownames(save_collapsed_clr_counts)),]
	pca = prcomp(t(collapsed_clr_counts))
	vars = round(pca$sdev^2 / (sum(pca$sdev^2))*100, digits = 0)

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

	plot(pca$x[,1], pca$x[,2], xlab = paste('MDS1 (',vars[1],'%)', sep=''), ylab = paste('MDS2 (',vars[2],'%)', sep=''), col = cols, pch = 16, las = 1, cex = 0.7, main = tax_level)
	legend('topright', fill = c('blue','green','red','orange'), legend = c('NY','MI','TN1','TN2'), border = NA, box.col= NA)

}
dev.off()

