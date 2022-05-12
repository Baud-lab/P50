load('/homes/abaud/P50_HSrats/data/metabo/study_spe_metabolites.RData')
#filtered_biomt_MI etc.

library(rhdf5)
fid='/homes/abaud/P50_HSrats/data/P50_rats_round8.h5'

load('/homes/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo.RData')
save_metadata = metadata

#write down to use in qsub_varianceDecomp.sh etc. 
dim(filtered_biomt_all)
#1879
dim(filtered_biomt_MI)
#1789
dim(filtered_biomt_NY)
#1977

dim(filtered_biomt_all)[1] + dim(filtered_biomt_MI)[1] + dim(filtered_biomt_NY)[1]
#5645

all_rats = unique(c(colnames(filtered_biomt_all),colnames(filtered_biomt_MI),colnames(filtered_biomt_NY)))

raw_data = matrix(ncol = dim(filtered_biomt_all)[1] + dim(filtered_biomt_MI)[1] + dim(filtered_biomt_NY)[1], nrow = length(all_rats))
colnames(raw_data)  = c(paste(rownames(filtered_biomt_all),'all',sep='_'), paste(rownames(filtered_biomt_MI),'MI',sep='_'),paste(rownames(filtered_biomt_NY),'NY',sep='_'))
rownames(raw_data) = all_rats

for (study in c('all','MI','NY')) {
	doto = get(paste('filtered_biomt_',study,sep=''))
	motch = match(all_rats, colnames(doto))
	raw_data[,paste(rownames(doto),study,sep='_')] = t(doto[,motch])

}
motch = match(rownames(raw_data), metadata$sample_name_metabo)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]
rownames(raw_data) = metadata$host_subject_id

raw_data[is.na(raw_data)] = (-999)


h5createGroup(fid,"metabo_raw_counts_new")
h5createGroup(fid,"/metabo_raw_counts_new/row_header")

max(nchar(rownames(raw_data)))
#10
h5createDataset(file=fid,dataset="metabo_raw_counts_new/row_header/sample_ID",dims=dim(raw_data)[1],storage.mode='character',size=15)
h5write(obj=rownames(raw_data),file=fid,name="metabo_raw_counts_new/row_header/sample_ID")

h5write(obj=raw_data,file=fid,name="metabo_raw_counts_new/matrix")

h5createGroup(fid,"metabo_raw_counts_new/col_header")
max(nchar(colnames(raw_data)))
#21
h5createDataset(file=fid,dataset="metabo_raw_counts_new/col_header/phenotype_ID",dims=dim(raw_data)[2],storage.mode='character',size=25)
h5write(obj=colnames(raw_data),file=fid,name="metabo_raw_counts_new/col_header/phenotype_ID")

max(nchar(metadata[,'cage']),na.rm=T)
#6
h5createDataset(file=fid,dataset="metabo_raw_counts_new/row_header/cage",dims=dim(metadata)[1],storage.mode='character',size=10)
h5write(obj=metadata[,'cage'],file=fid,name="metabo_raw_counts_new/row_header/cage")

