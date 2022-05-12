# does what it says: boxcox transformation (mindful of covariates) first then regress out covariates
#plot distribution of residuals 

library(MASS)
load('/homes/abaud/P50_HSrats/data/metabo/biomt_annots_metadata.RData')
load('/homes/abaud/P50_HSrats/data/metabo/center_spe_biomt.RData')

my_boxcox_regress = function(row) {
	m=min(row,na.rm=T)
	if (m<=0) stop('pb')

	bc=boxcox(row ~ 1 + sex + date_dissected,plotit=F,lambda=seq(-lambda_limit,lambda_limit,by=0.1))
	lambda=bc$x[which.max(bc$y)]
	resids =  rep(NA, length(row))
	if (abs(lambda) == lambda_limit) {print(paste('excluded because cannot be transformed to gaussian'))} else {
		if (lambda==0) row=log(row) else row=(row^lambda-1)/lambda
		lm1 = lm(row ~ 1 + sex + date_dissected)
		resids[!is.na(row)] = residuals(lm1) + lm1$coefficients['(Intercept)']
	}
	return(resids)
}

my_all_isna = function(row) {
	return(all(is.na(row)))
}

lambda_limit=2

for (suffix in c('MI','NY')) {
	biomt = get(paste('biomt',suffix,sep='_'))

	if (any(!colnames(biomt) %in% metadata$sample_name_metabo)) stop('wherzit')
	motch = match(colnames(biomt), metadata$sample_name_metabo)
	this_metadata = metadata[motch,]

	biomt[biomt == 0] = NA

	sex = this_metadata[,'sex']
	if(any(is.na(sex))) stop('sex')
	date_dissected = this_metadata[,'collection_timestamp']
	if(any(is.na(date_dissected))) stop('collection_timestamp')

	resids = t(apply(biomt, FUN = my_boxcox_regress, MAR = 1))
	rownames(resids) = rownames(biomt)
	colnames(resids) = colnames(biomt)
	resids = resids[!apply(resids,FUN = my_all_isna, MAR = 1),]

	assign(paste('resids',suffix,sep='_'),resids)
}

save(resids_MI, resids_NY, file = '/homes/abaud/P50_HSrats/data/metabo/boxcoxed_residuals.RData')

d = dist(resids_NY)
cluster = hclust(d)
ordre = cluster$order
pdf('/homes/abaud/P50_HSrats/plots/boxcox_distr_resids_NY_metabo.pdf')
par(mfrow = c(2,2))
for (i in seq(1,length(ordre), by = 100)) {
  hist(resids_NY[i,])
}
dev.off()



