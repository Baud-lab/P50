code = '65792'
prefix = paste('./data/metabo/',code,'_palmer_biom_original',sep='')
load(paste(prefix, '_resids_counts.RData',sep=''))
#load(paste(prefix, '_resids_qned_counts.RData',sep=''))
#load(paste(prefix, '_resids_presence.RData',sep=''))

load('./data/metadata/metadata_augmented_16S_metabo.RData')

library(rhdf5)
fid='./data/P50_rats_round8.h5'

save_metadata = metadata #because metadata will be subsetted for counts and then for presence

#write down to use in qsub_varianceDecomp.sh etc. 
dim(resids_qned_counts_all)
#[1] 4269 1090
dim(resids_qned_counts_MI)
#[1] 4677  609
dim(resids_qned_counts_NY)
#[1] 4106  481

dim(resids_qned_counts_all)[1] + dim(resids_qned_counts_MI)[1] + dim(resids_qned_counts_NY)[1]
#[1] 789

g = grepl('MZ415', rownames(resids_qned_counts_all)) & grepl(';5.7', rownames(resids_qned_counts_all))
rownames(resids_qned_counts_all)[g]

dim(resids_presence_all)
#[1] 4473 1090
dim(resids_presence_MI)
#[1] 3587  609
dim(resids_presence_NY)
#[1] 4537  481

dim(resids_presence_all)[1] + dim(resids_presence_MI)[1] + dim(resids_presence_NY)[1]
#[1] 1857

all_rats = unique(c(colnames(resids_qned_counts_all),colnames(resids_qned_counts_MI),colnames(resids_qned_counts_NY)))
#all_rats = unique(c(colnames(resids_presence_all),colnames(resids_presence_MI),colnames(resids_presence_NY),colnames(resids_qned_counts_all),colnames(resids_qned_counts_MI),colnames(resids_qned_counts_NY)))

#rats in rows
#features in columns - actually there will be one column for a feature_MI, one for feature_NY, and one for feature_all (provided common in MI, NY and all)
#for feature_MI there will be non NA for MI rats and 0 for NY rats; feature_MI will only exist if that feature is common across MI rats
#for feature_NY there will be non NA for NY rats and 0 for MI rats; feature_NY will only exist if that feature is common across NY rats
#for feature_all there will be non NA for all rats; feature_all will only exist if that feature is common across all rats
count_data = matrix(ncol = dim(resids_qned_counts_all)[1] + dim(resids_qned_counts_MI)[1] + dim(resids_qned_counts_NY)[1], nrow = length(all_rats))
colnames(count_data)  = c(paste(rownames(resids_qned_counts_all),'all',sep='_'), paste(rownames(resids_qned_counts_MI),'MI',sep='_'),paste(rownames(resids_qned_counts_NY),'NY',sep='_'))
rownames(count_data) = all_rats

presence_data = matrix(ncol = dim(resids_presence_all)[1] + dim(resids_presence_MI)[1] + dim(resids_presence_NY)[1], nrow = length(all_rats))
colnames(presence_data)  = c(paste(rownames(resids_presence_all),'all',sep='_'), paste(rownames(resids_presence_MI),'MI',sep='_'),paste(rownames(resids_presence_NY),'NY',sep='_'))
rownames(presence_data) = all_rats

for (study in c('all','MI','NY')) {
	resids = get(paste('resids_qned_counts_',study,sep=''))
	motch = match(all_rats, colnames(resids))
	count_data[,paste(rownames(resids),study,sep='_')] = t(resids[,motch])

#	resids = get(paste('resids_presence_',study,sep=''))
#	motch = match(all_rats, colnames(resids))
#	presence_data[,paste(rownames(resids),study,sep='_')] = t(resids[,motch])

}
motch = match(rownames(count_data), metadata$sample_name_metabo)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]
#below change row label to match the labels in the phenotype table
rownames(count_data) = metadata$host_subject_id

#there were NAs in the table now they become -999, which is the code for NAs for LIMIX
count_data[is.na(count_data)] = (-999)

#not sure why we need to save those...
#coolnames = colnames(count_data)
#save(coolnames, file = './data/metabo/coolnames.RData')


h5createGroup(fid,paste("metabo_counts_notqned_",code, sep=''))
h5createGroup(fid,paste("metabo_counts_notqned_",code,"/row_header", sep=''))

max(nchar(rownames(count_data)))
#10
h5createDataset(file=fid,dataset=paste("metabo_counts_notqned_",code,"/row_header/sample_ID", sep=''),dims=dim(count_data)[1],storage.mode='character',size=15)
h5write(obj=rownames(count_data),file=fid,name=paste("metabo_counts_notqned_",code,"/row_header/sample_ID", sep=''))

h5write(obj=count_data,file=fid,name=paste("metabo_counts_notqned_",code,"/matrix",sep=''))

h5createGroup(fid,paste("metabo_counts_notqned_",code,"/col_header",sep=''))
max(nchar(colnames(count_data)))
#28
h5createDataset(file=fid,dataset=paste("metabo_counts_notqned_",code,"/col_header/phenotype_ID",sep=''),dims=dim(count_data)[2],storage.mode='character',size=30)
h5write(obj=colnames(count_data),file=fid,name=paste("metabo_counts_notqned_",code,"/col_header/phenotype_ID",sep=''))

max(nchar(metadata[,'cage']),na.rm=T)
#6
h5createDataset(file=fid,dataset=paste("metabo_counts_notqned_",code,"/row_header/cage",sep=''),dims=dim(metadata)[1],storage.mode='character',size=10)
h5write(obj=metadata[,'cage'],file=fid,name=paste("metabo_counts_notqned_",code,"/row_header/cage",sep=''))

########

metadata = save_metadata
motch = match(rownames(presence_data), metadata$sample_name_metabo)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]
rownames(presence_data) = metadata$host_subject_id

presence_data[is.na(presence_data)] = (-999)

h5createGroup(fid,paste("metabo_presence_",code, sep=''))
h5createGroup(fid,paste("metabo_presence_",code,"/row_header", sep=''))

max(nchar(rownames(count_data)))
#10
h5createDataset(file=fid,dataset=paste("metabo_presence_",code,"/row_header/sample_ID", sep=''),dims=dim(count_data)[1],storage.mode='character',size=15)
h5write(obj=rownames(count_data),file=fid,name=paste("metabo_presence_",code,"/row_header/sample_ID", sep=''))

h5write(obj=count_data,file=fid,name=paste("metabo_presence_",code,"/matrix",sep=''))

h5createGroup(fid,paste("metabo_presence_",code,"/col_header",sep=''))
max(nchar(colnames(count_data)))
#28
h5createDataset(file=fid,dataset=paste("metabo_presence_",code,"/col_header/phenotype_ID",sep=''),dims=dim(count_data)[2],storage.mode='character',size=30)
h5write(obj=colnames(count_data),file=fid,name=paste("metabo_presence_",code,"/col_header/phenotype_ID",sep=''))

max(nchar(metadata[,'cage']),na.rm=T)
#6
h5createDataset(file=fid,dataset=paste("metabo_presence_",code,"/row_header/cage",sep=''),dims=dim(metadata)[1],storage.mode='character',size=10)
h5write(obj=metadata[,'cage'],file=fid,name=paste("metabo_presence_",code,"/row_header/cage",sep=''))
