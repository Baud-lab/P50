# calculates the P value, logP and estimate of the correlation between each bug and each phenotype, centers done separately
# uses spearman correlation at the moment
# then looks at the range of the value of interest (P, logP or cor)
# can plot a heatmap where axes (bugs and phenotypes) clustered
# can say which correlations are significant at FDR level or after Bonferroni correction

load('/homes/abaud/P50_HSrats/data/metabo/boxcox_residuals.RData')
load('/homes/abaud/P50_HSrats/data/phenotypes/all_phenos.RData')
load('/homes/abaud/P50_HSrats/data/general/sample_metadata.RData')
load('/homes/abaud/P50_HSrats/data/general/subset2include.RData')

my_cor = function(a,b, method, ret) {
	test = cor.test(a,b, use = "pairwise.complete.obs", method = method)
	return(c(estimate = test[['estimate']], pvalue = test[['p.value']], logP = -log10(test[['p.value']])))
}


for (phenotyping_center in c('MI','NY')) {
	pdf(paste('/homes/abaud/P50_HSrats/plots/corrs_oksana_',phenotyping_center,'.pdf',sep=''))
	resids = get(paste('resids',phenotyping_center,sep='_'))
	# has OTUs in rows and rats in columns

	# subset of rats
	if (phenotyping_center %in% c('MI','NY')) rats_this_center = rownames(sample_metadata)[sample_metadata[,'phenotyping_center'] == phenotyping_center] else rats_this_center = rownames(sample_metadata)[sample_metadata[,'phenotyping_center'] %in% c('MI','NY')]
	inter_rats = intersect(rats_this_center,colnames(resids))	
	inter_rats = intersect(inter_rats,rownames(all_phenos))
	inter_rats = intersect(inter_rats,good_samples)
	print(paste(length(inter_rats),'rats for',phenotyping_center))

	corr_list = list()
	for (pheno_name in colnames(all_phenos)) {
		if (all(is.na(all_phenos[inter_rats,pheno_name]))) next
		cors = apply(resids[,inter_rats],FUN = my_cor, MAR = 1, b = all_phenos[inter_rats,pheno_name],  method = "spearman", ret = 'p.value')
		corr_list [[pheno_name]] = as.data.frame(t(cors))
	}
	estimates = as.matrix(do.call('cbind',lapply(corr_list, '[', 1)))
	colnames(estimates) = names(corr_list)
	assign(paste('estimates',phenotyping_center,sep='_'), estimates)
	pvalues = as.matrix(do.call('cbind',lapply(corr_list, '[', 2)))
	colnames(pvalues) = names(corr_list)
	assign(paste('pvalues',phenotyping_center,sep='_'), pvalues)
	logPs = as.matrix(do.call('cbind',lapply(corr_list, '[', 3)))
	colnames(logPs) = names(corr_list)
	assign(paste('logPs',phenotyping_center,sep='_'), logPs)

	if (phenotyping_center == 'NY') th = 4 else th = 2
	for (k in 1:dim(which(logPs>th,arr.ind = T))[1]) {
		i = which(logPs>th,arr.ind = T)[k,1]
		j = which(logPs>th,arr.ind = T)[k,2]
		metabo = rownames(logPs)[i]
		pheno = colnames(logPs)[j]
		lm1 = lm(all_phenos[inter_rats,pheno] ~ resids[metabo,inter_rats])
		plot(resids[metabo,inter_rats],all_phenos[inter_rats,pheno], main = paste(estimates[i,j],logPs[i,j],sep=' '), las = 1, pch = 16, xlab = metabo, ylab = pheno)
		abline(lm1)
	}
	dev.off()
}

range(estimates_MI, na.rm = T)
range(logPs_MI, na.rm = T)
range(estimates_NY, na.rm = T)
range(logPs_NY, na.rm = T)
range(estimates_inter, na.rm = T)
range(logPs_inter, na.rm = T)
#estimates for cor range btw about -0.65 and +0.37
#logPs up to 14 for MI, lower for NY and 17 for inter!

unique(colnames(logPs_NY)[which(logPs_NY>5,arr.ind = T)[,2]])
#physiological_physiological_fasting_glucose
load('/homes/abaud/P50_HSrats/data/metabo/annots_unnormalized_metabo.RData')
interesting_metabos = unique(rownames(logPs_NY)[which(logPs_NY>5,arr.ind = T)[,1]])
annots[interesting_metabos,]
#Cholic acid corrs with fasting glycemia in both MI and NY

save(estimates_MI, pvalues_MI, logPs_MI,estimates_NY, pvalues_NY, logPs_NY,estimates_inter, pvalues_inter, logPs_inter, annots, file = '/homes/abaud/P50_HSrats/output/metabo/metabo_pheno_corrs.RData')

for (k in 1:dim(which(logPs_NY>4,arr.ind = T))[1]) {
	i = which(logPs_NY>4,arr.ind = T)[k,1]
	j = which(logPs_NY>4,arr.ind = T)[k,2]
	metabo = annots[rownames(logPs_NY)[i],'Compound_Name']
	pheno = colnames(logPs_NY)[j]
	if (!is.na(metabo)) print(paste(metabo,pheno, collapse = ' '))
}
for (k in 1:dim(which(logPs_MI>2,arr.ind = T))[1]) {
	i = which(logPs_MI>2,arr.ind = T)[k,1]
	j = which(logPs_MI>2,arr.ind = T)[k,2]
	metabo = annots[rownames(logPs_MI)[i],'Compound_Name']
	pheno = colnames(logPs_MI)[j]
	if (!is.na(metabo)) print(paste(metabo,pheno, collapse = ' '))
}

library(gap)
pdf('/homes/abaud/P50_HSrats/plots/QQplot_metabo_phenos_corrs.pdf',width = 10)
par(mfrow = c(1,2))
qqunif(pvalues_MI,ci=T,col='black',las=1, main = 'MI')
qqunif(pvalues_NY,ci=T,col='black',las=1, main = 'NY')
dev.off()


