library(rhdf5)
fid='/users/abaud/abaud/P50_HSrats/data/P50_rats_Rn7_Helene.h5'

load('/nfs/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
save_metadata = metadata

load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/resids_qned_counts.RData')
#resids_qned_counts_all
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/resids_presence.RData')
#resids_presence_all

#write down to use in qsub_varianceDecomp.sh etc. 
dim(resids_qned_counts_all)
#296
dim(resids_qned_counts_MI)
#309
dim(resids_qned_counts_NY)
#302
dim(resids_qned_counts_TN_behavior)
#281
dim(resids_qned_counts_TN_breeder)
#320

dim(resids_qned_counts_MI)[1] + dim(resids_qned_counts_NY)[1] + dim(resids_qned_counts_TN_behavior)[1] + dim(resids_qned_counts_TN_breeder)[1]
#1212

dim(resids_presence_all)
#64 3898
dim(resids_presence_MI)
#53 1120
dim(resids_presence_NY)
#61 1179
dim(resids_presence_TN_behavior)
#59 1025
dim(resids_presence_TN_breeder)
#63 574

dim(resids_presence_all)[1] + dim(resids_presence_MI)[1] + dim(resids_presence_NY)[1] + dim(resids_presence_TN_behavior)[1] + dim(resids_presence_TN_breeder)[1]
#216

#start with colnames of bug_data are full_ids but later turned to host_subject_id
all_rats = unique(c(colnames(resids_presence_all),colnames(resids_presence_MI),colnames(resids_presence_NY), colnames(resids_presence_TN_behavior), colnames(resids_presence_TN_breeder)))

count_data = matrix(ncol = dim(resids_qned_counts_all)[1] + dim(resids_qned_counts_MI)[1] + dim(resids_qned_counts_NY)[1] + dim(resids_qned_counts_TN_behavior)[1] + dim(resids_qned_counts_TN_breeder)[1], nrow = length(all_rats))
colnames(count_data)  = c(paste(rownames(resids_qned_counts_all),'all',sep='_'), paste(rownames(resids_qned_counts_MI),'MI',sep='_'),paste(rownames(resids_qned_counts_NY),'NY',sep='_'),paste(rownames(resids_qned_counts_TN_behavior),'TN_behavior',sep='_'),paste(rownames(resids_qned_counts_TN_breeder),'TN_breeder',sep='_'))
rownames(count_data) = all_rats

presence_data = matrix(ncol = dim(resids_presence_all)[1] + dim(resids_presence_MI)[1] + dim(resids_presence_NY)[1] + dim(resids_presence_TN_behavior)[1] + dim(resids_presence_TN_breeder)[1], nrow = length(all_rats))
colnames(presence_data)  = c(paste(rownames(resids_presence_all),'all',sep='_'), paste(rownames(resids_presence_MI),'MI',sep='_'),paste(rownames(resids_presence_NY),'NY',sep='_'),paste(rownames(resids_presence_TN_behavior),'TN_behavior',sep='_'),paste(rownames(resids_presence_TN_breeder),'TN_breeder',sep='_'))
rownames(presence_data) = all_rats

for (study in c('all','MI','NY','TN_behavior','TN_breeder')) {
	resids = get(paste('resids_qned_counts_',study,sep=''))
	motch = match(all_rats, colnames(resids))
	count_data[,paste(rownames(resids),study,sep='_')] = t(resids[,motch])

	resids = get(paste('resids_presence_',study,sep=''))
	motch = match(all_rats, colnames(resids))
	presence_data[,paste(rownames(resids),study,sep='_')] = t(resids[,motch])

}
motch = match(rownames(count_data), metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]
rownames(count_data) = metadata$host_subject_id

count_data[is.na(count_data)] = (-999)

coolnames = colnames(count_data)
save(coolnames, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/coolnames_Helene.RData')

h5createGroup(fid,"phenotypes")
h5createGroup(fid,"phenotypes/deblur_counts")
h5createGroup(fid,"phenotypes/deblur_counts/row_header")

mox = max(nchar(rownames(count_data)))
#10
h5createDataset(file=fid,dataset="phenotypes/deblur_counts/row_header/sample_ID",dims=dim(count_data)[1],storage.mode='character',size=(mox+2))
h5write(obj=rownames(count_data),file=fid,name="phenotypes/deblur_counts/row_header/sample_ID")

h5write(obj=count_data,file=fid,name="phenotypes/deblur_counts/matrix")

h5createGroup(fid,"phenotypes/deblur_counts/col_header")
mox = max(nchar(colnames(count_data)))
#53
h5createDataset(file=fid,dataset="phenotypes/deblur_counts/col_header/phenotype_ID",dims=dim(count_data)[2],storage.mode='character',size=(mox+2))
h5write(obj=colnames(count_data),file=fid,name="phenotypes/deblur_counts/col_header/phenotype_ID")

#only if not added yet:

mox = max(nchar(metadata[,'cage']),na.rm=T)
#32
h5createGroup(fid,"cages")
h5createGroup(fid,"cages/cage_TNfixed2")
h5createDataset(file=fid,dataset="cages/cage_TNfixed2/array",dims=dim(count_data)[1],storage.mode='character',size=(mox + 2))
h5write(obj=metadata[,'cage'],file=fid,name="cages/cage_TNfixed2/array")
h5createDataset(file=fid,dataset="cages/cage_TNfixed2/sample_ID",dims=dim(count_data)[1],storage.mode='character',size=(mox + 2))
h5write(obj=rownames(count_data),file=fid,name="cages/cage_TNfixed2/sample_ID")


mox = max(nchar(metadata[,'dam']),na.rm=T)
#10
h5createGroup(fid,"dam")
h5createGroup(fid,"dam/dams")
h5createDataset(file=fid,dataset="dam/dams/array",dims=dim(count_data)[1],storage.mode='character',size=(mox + 2))
h5write(obj=metadata[,'dam'],file=fid,name="dam/dams/array")
h5createDataset(file=fid,dataset="dam/dams/sample_ID",dims=dim(count_data)[1],storage.mode='character',size=(mox + 2))
h5write(obj=rownames(count_data),file=fid,name="dam/dams/sample_ID")

