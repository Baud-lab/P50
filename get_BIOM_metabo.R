#in all scripts you need to set your working directly (including when you launch jobs on the cluster) to the right "P50_root" folder
#e.g. for me
#setwd("/homes/abaud/P50_HSrats")

#/pre_Apr8/no_gapFilling/65792_palmer_biom_original.biom is very first dataset put on Qiita. non gap-filled nor IIN. unnormalised. blank corrected I think
#/pre_Apr8/78940_palmer_amelie_metabolomics_original.biom on Qiita corresponds to Allegras_061419.csv. gap filled and IIN + blank corrrected I think BUT Allegras parameters which are quite different from Kellys and not recommended
#115698_palmer_biom_march_2021_original.biom gap filled but not IIN. dont know about blank corrected. Kellys parameters.
#115697_palmer_IIN_biom_march_2021_original.biom gap filled and IIN. dont know about blank corrected. Kellys parameters.


#retrieve feature annotations (come from MS2) from the biom table- which is an HDF% file
library(rhdf5)
#first dataset ever /pre_Apr8/no_gapFilling/65792_palmer_biom_original.biom
#this is the table with the large number of features
prefix = './data/metabo/65792_palmer_biom_original'
biom = h5read(paste(prefix, '.biom',sep=''),'/')
#below is another table with fewer features - features with no C13 peak have been removed
#biom = h5read('./data/metabo/115697_palmer_IIN_biom_march_2021_original.biom','/')
annots = as.data.frame(biom$observation$metadata)
rownames(annots) = paste('MZ',biom$observation$ids,sep='')
annots[annots=='nan'] = NA
my_f = function(word) {
	splot = strsplit(sub('MZ','',word),';')[[1]]
	return(c(as.numeric(splot[1]),as.numeric(splot[2])))	
}
infos = t(sapply(rownames(annots),my_f))
colnames(infos) = c('mz','RT')
annots = cbind(annots,infos)
save(annots, file = paste(prefix, '_annots.RData',sep=''))

#now quantitative measures from MS1
#biom convert -i ./data/metabo/115698_palmer_biom_march_2021_original.biom -o ./data/metabo/115698_palmer_biom_march_2021_original.txt --to-tsv
#needs RAM!
prefix = './data/metabo/65792_palmer_biom_original'
biomt = read.delim(paste(prefix, '.txt',sep=''),as.is = T, header = F)
samples = unlist(biomt[2,-1, drop = T])
biomt = biomt[-c(1,2),]
MZs = paste('MZ',biomt[,1],sep='')
biomt = biomt[,-1]
dim(biomt)
#[1] 11388  1258

biomt = lapply(biomt, FUN = as.numeric)
biomt = as.matrix(as.data.frame(biomt, stringsAsFactors = F))
colnames(biomt) = samples
rownames(biomt) = MZs

dim(biomt)

save(biomt, file = paste(prefix, '_biomt.RData',sep=''))

#equivalent to lib size in 16S but not actually from sequencing
lib_sizes = apply(biomt, FUN = sum, MAR = 2)
save(lib_sizes, file = paste(prefix,'_lib_sizes.RData',sep=''))


