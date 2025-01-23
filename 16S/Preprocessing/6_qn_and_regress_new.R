# does what it says: boxcox transformation (mindful of covariates) first then regress out covariates
#this script is for collapsed/taxa level data; similar script for uncollapsed/asv level data 

library(MASS)

load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData')
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/study_spe_taxa.RData')
load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')

invrank= function(row) {qnorm((rank(row,na.last="keep", ties.method = "random")-0.5)/sum(!is.na(row)))}

my_regress = function(row) {
	if (study == 'all') { 
		covariates = c('sex','prep','batch','lib_size','study_cov','weaning_suite') 
		#test covariates = c('sex','prep','batch','study_cov','weaning_suite')  #lib size fucks up even though it has normal distrib
	} else if (study %in% c('MI','NY')) {
		covariates = c('sex','prep','batch','lib_size') 
	} else if (study == 'TN_behavior') {
		covariates = c('sex','prep','batch','lib_size','weaning_suite') 
	} else if (study == 'TN_breeder') {
		covariates = c('sex','prep','batch','lib_size') 
	}

	#print(covariates)
	resids =  rep(NA, length(row))	
	lm1 = lm(as.formula(paste('row ~ 1 + ',paste(covariates, collapse = '+'),sep=''))) 
	resids[!is.na(row)] = residuals(lm1)
	return(resids)
}

for (study in c('all','MI','NY','TN_behavior','TN_breeder')) {
	print(study)
	filtered_collapsed_biomt_counts = get(paste('filtered_collapsed_clr_counts',study,sep='_'))
	filtered_collapsed_biomt_presence =  get(paste('filtered_presence',study,sep='_')) #is also collapsed even though doesnt say in name

	if (any(!colnames(filtered_collapsed_biomt_counts) %in% metadata$deblur_rooname)) stop('wherzit')
	motch = match(colnames(filtered_collapsed_biomt_counts), metadata$deblur_rooname)
	this_metadata = metadata[motch,]

	sex = this_metadata[,'sex']
	if(any(is.na(sex))) stop('sex')
	prep = this_metadata[,'prep_name_16S']
	if(any(is.na(prep))) stop('prep_name_16S')
	prep = droplevels(prep)
	batch = this_metadata[,'batch']
	if(any(is.na(batch))) stop('batch')
	batch = droplevels(batch)
	lib_size = this_metadata[,'nb_seqs_16S']
	if(any(is.na(lib_size))) stop('lib_size')
	study_cov = rep(NA, dim(this_metadata)[1])
	study_cov[this_metadata[,'phenotyping_center'] == 'MI'] = 'MI'
	study_cov[this_metadata[,'phenotyping_center'] == 'NY'] = 'NY'
	study_cov[this_metadata[,'phenotyping_center'] == 'TN' & this_metadata[,'TN_track'] == 'breeder'] = 'TN_breeder'
	study_cov[this_metadata[,'phenotyping_center'] == 'TN' & this_metadata[,'TN_track'] == 'behavior'] = 'TN_behavior'
	if(any(is.na(study_cov))) { w= which(is.na(study_cov)); print(table(this_metadata[w,'phenotyping_center'],this_metadata[w,'TN_track'], useNA = 'always'))}
	study_cov = factor(study_cov)
	study_cov = droplevels(study_cov)
	weaning_suite = this_metadata[,'weaning_suite']
	if(any(is.na(weaning_suite))) stop('weaning_suite')
	weaning_suite = droplevels(weaning_suite)

	cc = complete.cases(sex, prep, batch, lib_size, study_cov, weaning_suite)
	sex = sex[cc]
	prep = prep[cc]
	batch = batch[cc]
	lib_size = lib_size[cc]
	study_cov = study_cov[cc]
	weaning_suite = weaning_suite[cc]
	filtered_collapsed_biomt_counts = filtered_collapsed_biomt_counts[,cc]
	filtered_collapsed_biomt_presence = filtered_collapsed_biomt_presence[,cc]

	qned_counts = t(apply(filtered_collapsed_biomt_counts, FUN = invrank, MAR = 1))
	colnames(qned_counts) = colnames(filtered_collapsed_biomt_counts)
	assign(paste('qned_counts',study,sep='_'),qned_counts)

	resids_qned_counts = t(apply(qned_counts, FUN = my_regress, MAR = 1)) 
	colnames(resids_qned_counts) = colnames(filtered_collapsed_biomt_counts)
	assign(paste('resids_qned_counts',study,sep='_'),resids_qned_counts)

	resids_presence = t(apply(filtered_collapsed_biomt_presence, FUN = my_regress, MAR = 1))
	colnames(resids_presence) = colnames(filtered_collapsed_biomt_presence)
	assign(paste('resids_presence',study,sep='_'),resids_presence)

}

lost = ls()
to_save = lost[grepl('qned_counts_',lost)]
save(list = to_save, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/NONresids_qned_counts.RData')
to_save = lost[grepl('resids_qned_counts_',lost)]
save(list = to_save, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/resids_qned_counts.RData')
to_save = lost[grepl('resids_presence_',lost)]
save(list = to_save, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/resids_presence.RData')



