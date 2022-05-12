#library(MASS)

prefix = './data/metabo/65792_palmer_biom_original'
load(paste(prefix,'_study_spe.RData',sep=''))

load(paste(prefix, '_lib_sizes.RData',sep=''))
load('./data/metadata/metadata_augmented_16S_metabo.RData')
#the matching will be done later

invrank= function(row) {qnorm((rank(row,na.last="keep")-0.5)/sum(!is.na(row)))}

my_regress = function(row) {
	if (study == 'all') { 
		covariates = c('sex','batch','lib_size','study_cov') 
	} else if (study %in% c('MI','NY')) {
		covariates = c('sex','batch','lib_size') 
	}

	resids =  rep(NA, length(row))	
	lm1 = lm(as.formula(paste('row ~ 1 + ',paste(covariates, collapse = '+'),sep=''))) 
	resids[!is.na(row)] = residuals(lm1)
	return(resids)
}

for (study in c('all','MI','NY')) {
	print(study)
	this_filtered_clr_counts = get(paste('filtered_clr_counts_',study,sep=''))
	this_filtered_presence =  get(paste('filtered_presence',study,sep='_'))
	if(!all(colnames(this_filtered_clr_counts) == colnames(this_filtered_presence))) stop('pb names')
	#rownames already there

	motch = match(colnames(this_filtered_clr_counts), names(lib_sizes))
	if (any(is.na(motch))) stop('pb wherzit')
	lib_size = lib_sizes[motch]

	if (any(!colnames(this_filtered_clr_counts) %in% metadata$sample_name_metabo)) stop('wherzit')
	motch = match(colnames(this_filtered_clr_counts), metadata$sample_name_metabo)
	this_metadata = metadata[motch,]

	sex = this_metadata[,'sex']
	if(any(is.na(sex))) stop('sex')
	batch = this_metadata[,'batch']
	if(any(is.na(batch))) stop('batch')
	batch = droplevels(batch)
	study_cov = rep(NA, dim(this_metadata)[1])
	study_cov[this_metadata[,'phenotyping_center'] == 'MI'] = 'MI'
	study_cov[this_metadata[,'phenotyping_center'] == 'NY'] = 'NY'
	#NAs in study_cov are possible - but want to check them by printing
	if(any(is.na(study_cov))) { w= which(is.na(study_cov)); print(table(this_metadata[w,'phenotyping_center'],this_metadata[w,'TN_track'], useNA = 'always'))}
	study_cov = factor(study_cov)
	study_cov = droplevels(study_cov)

	cc = complete.cases(sex, batch, lib_size, study_cov)
	sex = sex[cc]
	batch = batch[cc]
	lib_size = lib_size[cc]
	study_cov = study_cov[cc]
	this_filtered_clr_counts = this_filtered_clr_counts[,cc]
	this_filtered_presence = this_filtered_presence[,cc]

	qned_counts = t(apply(this_filtered_clr_counts, FUN = invrank, MAR = 1))
	resids_qned_counts = t(apply(qned_counts, FUN = my_regress, MAR = 1))
	colnames(resids_qned_counts) = colnames(this_filtered_clr_counts)
	assign(paste('resids_qned_counts',study,sep='_'),resids_qned_counts)

	resids_presence = t(apply(this_filtered_presence, FUN = my_regress, MAR = 1))
	colnames(resids_presence) = colnames(this_filtered_presence)
	assign(paste('resids_presence',study,sep='_'),resids_presence)
}

save(resids_qned_counts_all, resids_qned_counts_MI, resids_qned_counts_NY, file = paste(prefix, '_resids_qned_counts.RData',sep=''))
save(resids_presence_all, resids_presence_MI, resids_presence_NY, file = paste(prefix, '_resids_presence.RData',sep=''))



