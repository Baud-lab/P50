setwd('/Users/abaud/Dropbox (Palmer Lab)/Palmer Lab/Amelie Baud/Microbiome_analysis/data/Metabolomics')
bucket = read.csv('SEED_Palmer_5E3_1E4_Filtered.csv',as.is = T)
dim(bucket)
#[1] 1515  120
#rows are features
#cols -first 3 are samples
rats  = colnames(bucket)[grep('SEED_C18p_Palmer',colnames(bucket))]
length(rats)
#96 rat samples
length(unique(rats))
#96

blanks  = colnames(bucket)[grep('Blnk',colnames(bucket))]
length(blanks)
#10

washes  = colnames(bucket)[grep('Wash',colnames(bucket))]
length(washes)
#10

submitted = read.table('submitted_metabo_64.csv',as.is = T)
#very different names for rats

#plate info
plate = read.csv('plate.csv',as.is = T, check.names = F,row.names = 1)
info = cbind(rats, NA, NA)
colnames(info) = c('mzXML','well','ID')
for (i in 1:dim(info)[1]) {
	splot = strsplit(info[i,'mzXML'],'_')[[1]]
	info[i,'well'] = sub('B','',splot[5])
	info[i,'ID'] = plate[substr(info[i,'well'],1,1),substr(info[i,'well'],2,2)]
}
info = info[grep('^0007',info[,'ID']),]

rownames(bucket) = paste(bucket[,'row.m.z'], bucket[,'row.retention.time'],sep='_') 
motch = match(colnames(bucket), info[,'mzXML'])
any(is.na(motch))
#TRUE
bucket = bucket[,!is.na(motch)]
colnames(bucket) = info[na.exclude(motch),'ID']
dim(bucket)
#[1] 1515   64
#is that reasonable number of features?
bucket = as.matrix(bucket)
save(bucket, file = 'bucket.RData')

#metadata
metadata = read.csv('metadata.csv',as.is = T,check.names = F)
metadata = metadata[1:65,]
rownames(metadata) = metadata[,'RFID']
metadata = metadata[,-1]
save(metadata,file = 'metadata.RData')




