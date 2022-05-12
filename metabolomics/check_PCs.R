load('/homes/abaud/P50_HSrats/data/metabo/78940_palmer_amelie_metabolomics_original.RData')
lib_sizes = apply(biomt, FUN = sum, MAR = 2)

PCs = read.table('/homes/abaud/P50_HSrats/data/metabo/PCoA_108573_ordination.txt', as.is = T)
load('/homes/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo.RData')

motch = match(sub('78940.','',PCs[,1]), metadata$sample_name_metabo)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

motch = match(metadata$sample_name_metabo, names(lib_sizes))
lib_sizes = lib_sizes[motch]

w = which(metadata$phenotyping_center == 'NY' & metadata$sex == 'M')
PCs = PCs[w,]
metadata = metadata[w,]
lib_sizes = lib_sizes[w]

plot(PCs[,2], PCs[,3])

hist(PCs[,3])
#PC2 clearly bimodal, not PC1 or PC3

plot(PCs[,2], lib_sizes)
#PC1 highly negatively correlated with lib_sizes, not PC2 or PC3

#find out if a covariate explains the two groups
pdf('/homes/abaud/P50_HSrats/plots/PCs_metabo.pdf')
for (col in 1:dim(metadata)[2]) {
	trou = try(plot(PCs[,3], metadata[,col], main = colnames(metadata)[col]))
	if (inherits(trou, 'try-error')) try(boxplot(PCs[,1]~ metadata[,col], main = colnames(metadata)[col]))
}
dev.off()




