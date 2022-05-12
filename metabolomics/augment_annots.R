#this code calculates the percentile in terms of median or mean abudnance for each feature
#or something like that - double check what it does when needed
load('/homes/abaud/P50_HSrats/data/metabo/annots_115698_palmer_biom_march_2021_original.RData')
load('/homes/abaud/P50_HSrats/data/metabo/115698_palmer_biom_march_2021_original.RData')
load('/homes/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo.RData')

#blanks = grep('blank', colnames(biomt), ignore.case = T)
#biomt = biomt[,-blanks]

motch = match(colnames(biomt),metadata$sample_name_metabo)
metadata = metadata[motch,]

blanks = grep('blank', colnames(biomt), ignore.case = T)
medians = apply(biomt[,-blanks], MAR = 1, FUN = median)
percentile_medians = rank(medians, ties.method = "average") / dim(biomt)[1]*100
names(percentile_medians) = rownames(biomt)
means = apply(biomt[,-blanks], MAR = 1, FUN = mean)
percentile_means = rank(means, ties.method = "average") / dim(biomt)[1]*100
names(percentile_means) = rownames(biomt)

all(rownames(annots) == rownames(biomt))
#TRUE

save(percentile_medians, percentile_means, annots,biomt,metadata, file = '/homes/abaud/P50_HSrats/data/metabo/biomt_annots_metadata.RData')

