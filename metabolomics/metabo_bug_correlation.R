# calculates the P value, logP and estimate of the correlation between each bug and each phenotype, centers done separately
# uses spearman correlation at the moment
# then looks at the range of the value of interest (P, logP or cor)
# can plot a heatmap where axes (bugs and phenotypes) clustered
# can say which correlations are significant at FDR level or after Bonferroni correction

load('/homes/abaud/P50_HSrats/data/metabo/boxcox_residuals.RData')
metabo_resids_inter = resids_inter
metabo_resids_MI = resids_MI
metabo_resids_NY = resids_NY
load('/homes/abaud/P50_HSrats/data/16S/boxcox_residuals.RData')
bug_resids_inter = resids_inter
bug_resids_MI = resids_MI
bug_resids_NY = resids_NY

load('/homes/abaud/P50_HSrats/data/general/sample_metadata.RData')
load('/homes/abaud/P50_HSrats/data/general/subset2include.RData')

my_cor = function(a,b, method, ret) {
	test = cor.test(a,b, use = "pairwise.complete.obs", method = method)
	return(c(estimate = test[['estimate']], pvalue = test[['p.value']], logP = -log10(test[['p.value']])))
}


for (phenotyping_center in c('MI','NY','inter')) {

	metabo_resids = get(paste('metabo_resids',phenotyping_center,sep='_'))
	bug_resids = get(paste('bug_resids',phenotyping_center,sep='_'))
	# has OTUs or metabolites in rows and rats in columns

	# subset of rats
	if (phenotyping_center %in% c('MI','NY')) rats_this_center = rownames(sample_metadata)[sample_metadata[,'phenotyping_center'] == phenotyping_center] else rats_this_center = rownames(sample_metadata)[sample_metadata[,'phenotyping_center'] %in% c('MI','NY')]
	inter_rats = intersect(rats_this_center,colnames(metabo_resids))	
	inter_rats = intersect(inter_rats,colnames(bug_resids))
	inter_rats = intersect(inter_rats,good_samples)
	print(paste(length(inter_rats),'rats for',phenotyping_center))

	corr_list = list()
	for (bug in rownames(bug_resids)) {
		cors = apply(metabo_resids[,inter_rats],FUN = my_cor, MAR = 1, b = bug_resids[bug,inter_rats],  method = "spearman", ret = 'p.value')
		corr_list [[bug]] = as.data.frame(t(cors))
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
}
range(estimates_MI, na.rm = T)
range(logPs_MI, na.rm = T)
range(estimates_NY, na.rm = T)
range(logPs_NY, na.rm = T)
range(estimates_inter, na.rm = T)
range(logPs_inter, na.rm = T)
#estimates for cor range btw about -0.65 and +0.37
#logPs up to 14 for MI, lower for NY and 17 for inter!

load('/homes/abaud/P50_HSrats/data/16S/taxonomy.RData')
interesting_bugs = unique(colnames(logPs_NY)[which(logPs_NY>10,arr.ind = T)[,2]])
taxonomy[interesting_bugs]
load('/homes/abaud/P50_HSrats/data/metabo/annots_unnormalized_metabo.RData')
interesting_metabos = unique(rownames(logPs_NY)[which(logPs_NY>10,arr.ind = T)[,1]])
annots[interesting_metabos,]

for (k in 1:dim(which(logPs_NY>10,arr.ind = T))[1]) {
	i = which(logPs_NY>10,arr.ind = T)[k,1]
	j = which(logPs_NY>10,arr.ind = T)[k,2]
	metabo = annots[rownames(logPs_NY)[i],'Compound_Name']
	bug = taxonomy[colnames(logPs_NY)[j]]
	if (!is.na(metabo)) print(paste(metabo,bug, collapse = ' '))
}
for (k in 1:dim(which(logPs_MI>10,arr.ind = T))[1]) {
	i = which(logPs_MI>10,arr.ind = T)[k,1]
	j = which(logPs_MI>10,arr.ind = T)[k,2]
	metabo = annots[rownames(logPs_MI)[i],'Compound_Name']
	bug = taxonomy[colnames(logPs_MI)[j]]
	if (!is.na(metabo)) print(paste(metabo,bug, collapse = ' '))
}


pdf('/homes/abaud/P50_HSrats/plots/QQplot_metabo_bug_corrs.pdf',width = 10)
par(mfrow = c(1,2))
qqunif(pvalues_MI,ci=T,col='black',las=1, main = 'MI')
qqunif(pvalues_NY,ci=T,col='black',las=1, main = 'NY')
dev.off()


save(estimates_MI, pvalues_MI, logPs_MI,estimates_NY, pvalues_NY, logPs_NY,estimates_inter, pvalues_inter, logPs_inter, annots, file = '/homes/abaud/P50_HSrats/output/metabo/metabo_bug_corrs.RData')

#find metabolites that are highly correlated with OTU261590 and OTU214919
dim(logPs_MI)
dim(logPs_NY)
inter = intersect(rownames(logPs_MI),rownames(logPs_NY))
length(inter)
#169 - try
sub_logPs_MI = logPs_MI[match(inter,rownames(logPs_MI)),]
sub_logPs_NY = logPs_NY[match(inter,rownames(logPs_NY)),]
all(rownames(sub_logPs_MI) == rownames(sub_logPs_NY))
#TRUE
my_check = function(row, th) {
	return(all(row >= th))
}
passes_MI = apply(sub_logPs_MI[,c('OTU261590','OTU214919')], FUN = my_check, MARGIN = 1, th = 3)
passes_NY = apply(sub_logPs_NY[,c('OTU261590','OTU214919')], FUN = my_check, MARGIN = 1, th = 3)
sum(passes_MI & passes_NY)
annots[sort(rownames(sub_logPs_MI)[passes_MI & passes_NY]),]
#MZ110.0087;0.2262  <NA>                 <NA>   <NA>
#MZ128.0191;0.2245  <NA>                 <NA>   <NA>

annots[grep('MZ1.*;0.22.*',rownames(annots)),]





