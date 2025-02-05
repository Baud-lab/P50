root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'

deblur_counts_uncollapsed_dir ='deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_IGE_IEE_cageEffect_maternalEffect/'
load(paste(root_dir,deblur_counts_uncollapsed_dir,'all_estNste.Rdata',sep=''))
asv_VCs = VCs

deblur_counts_dir ='deblur_counts/P50_Rn7_pruned_DGE_IGE_IEE_cageEffect_maternalEffect/'
load(paste(root_dir,deblur_counts_dir,'all_estNste.Rdata',sep=''))
tax_levels_VCs = VCs

community_traits_dir ='community_traits3/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(paste(root_dir,community_traits_dir,'all_estNste.Rdata',sep=''))
community_VCs = VCs

all(colnames(tax_levels_VCs) == colnames(asv_VCs))
all(colnames(tax_levels_VCs) == colnames(community_VCs))

all_VCs_full = rbind(tax_levels_VCs, asv_VCs, community_VCs)
#all_VCs_full = all_VCs_full[!grepl('_all', all_VCs_full$trait1),]
rm(list = c('tax_levels_VCs','asv_VCs','community_VCs'))

deblur_counts_uncollapsed_dir ='deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_IEE_cageEffect_maternalEffect/'
load(paste(root_dir,deblur_counts_uncollapsed_dir,'all_estNste.Rdata',sep=''))
asv_VCs = VCs

deblur_counts_dir ='deblur_counts/P50_Rn7_pruned_DGE_IEE_cageEffect_maternalEffect/'
load(paste(root_dir,deblur_counts_dir,'all_estNste.Rdata',sep=''))
tax_levels_VCs = VCs

community_traits_dir ='community_traits3/P50_Rn7_pruned_cageEffect_maternalEffect/'
load(paste(root_dir,community_traits_dir,'all_estNste.Rdata',sep=''))
community_VCs = VCs

all(colnames(tax_levels_VCs) == colnames(asv_VCs))
all(colnames(tax_levels_VCs) == colnames(community_VCs))

all_VCs_null = rbind(tax_levels_VCs, asv_VCs, community_VCs)
#all_VCs_null = all_VCs_null[!grepl('_all', all_VCs_null$trait1),]

inter=intersect(all_VCs_null$trait1,all_VCs_full$trait1)
length(inter)
#[1] 2322
all_VCs_full=all_VCs_full[match(inter,all_VCs_full$trait1),]
all_VCs_null=all_VCs_null[match(inter,all_VCs_null$trait1),]
 #both TRUE

source('/users/abaud/abaud/P50_HSrats/code/variance_decomposition/felipes_deblur/annotate_VCs_pvalues_function.R')
all_VCs_full = annotate(all_VCs_full)

## calculate T^2 using n = 2, which is conservative since in NY for example rats were 3 per cage. 
all_VCs_full$total_heritability = all_VCs_full$prop_Ad1 + 2*(2-1)*all_VCs_full$corr_Ad1s1*sqrt(all_VCs_full$prop_Ad1*all_VCs_full$prop_As1) + (2-1)^2*all_VCs_full$prop_As1 


#df = 0.7 for IGE
k = 0.7
all_VCs_full$pvalue_DGE = k*pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=1, lower.tail = FALSE) + (1-k)*pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=2, lower.tail = FALSE)
#df = 1 for DGE
#all_VCs_full$pvalue_DGE = pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=1, lower.tail = FALSE)
#all_VCs_full$sw_bonf_pvalue_DGE = all_VCs_full$pvalue_DGE*dim(all_VCs_full)[1] #sw = study wide, ie all 4 cohorts
library(qvalue)
all_VCs_full$sw_qvalue_DGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue

all_VCs_full$cw_bonf_pvalue_DGE = NA
all_VCs_full$cw_qvalue_DGE = NA
#for (study in c('all')) {
for (study in c('all','MI','NY','TN_breeder','TN_behavior')) {
	w = which(all_VCs_full$study1 == study)
	all_VCs_full[w,'cw_bonf_pvalue_DGE'] = all_VCs_full[w,'pvalue_DGE']*length(w) #cw = cohort (ie study) wide
	all_VCs_full[w,'cw_qvalue_DGE'] = qvalue(all_VCs_full[w,'pvalue_DGE'])$qvalue
}

all_VCs_full = all_VCs_full[order(all_VCs_full$pvalue_DGE),]
save(all_VCs_full, file =  paste(root_dir,'augmented_IGE_VC.RData',sep=''))

all_VCs_full = all_VCs_full[!is.na(all_VCs_full$full_taxon),]
for (study in c('all','NY','MI','TN_behavior','TN_breeder')) {
	print(study)
	w = which(all_VCs_full$study1 == study)
	print(paste('Prop FDR sig',sum(all_VCs_full[w,'cw_qvalue_DGE'] < 0.1)/ length(w)))
	subset = all_VCs_full[w,]
#	tax_level = 'g__'
#	subset = subset[grep(tax_level, subset$trait1),]
	print(paste('FDR sig taxa',sum(subset$cw_qvalue_DGE < 0.1),'out of',dim(subset)[1]))
	#print(subset[subset$cw_qvalue_DGE < 0.1,])
}

root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(paste(root_dir,'augmented_IGE_VC.RData',sep=''))

all_VCs_full = all_VCs_full[all_VCs_full$study1 == 'all',]
all_VCs_full = all_VCs_full[order(all_VCs_full$pvalue_DGE, decreasing = T),] #leave to TRUE: it's for QQ plot colours

cols = rep('black', dim(all_VCs_full)[1])
cols[all_VCs_full$cw_qvalue_DGE < 0.1] = 'orange'
cols[all_VCs_full$cw_bonf_pvalue_DGE < 0.05] = 'red'

pchs = rep(1, dim(all_VCs_full)[1])
pchs[all_VCs_full$cw_qvalue_DGE < 0.1] = 16

library(gap)
mox = max(-log10(all_VCs_full$pvalue_DGE), na.rm = T)
pdf('/users/abaud/abaud/P50_HSrats/plots/QQplot_pvalues_IGE_Helenes.pdf',bg='white')
#par(mfrow = c(2,2), mar = c(2,2,3,1))
#qqunif(all_VCs_full[grepl('_all',all_VCs_full$trait1),'pvalue_DGE'],ci=T,col = 'black',las=1, main ='All', ylim = c(0,mox))
# clear inflation there
#g_NY = grepl('_NY',all_VCs_full$trait1)
#qqunif(all_VCs_full[g_NY,'pvalue_DGE'],ci=T,las=1, main ='NY behavior: 1,167 rats', ylim = c(0,mox), col = cols[g_NY], pch = pchs[g_NY])
#g_MI = grepl('_MI',all_VCs_full$trait1)
#qqunif(all_VCs_full[g_MI,'pvalue_DGE'],ci=T,las=1, main ='MI behavior: 1,112 rats', ylim = c(0,mox), col = cols[g_MI], pch = pchs[g_MI])
#g_TN_behavior = grepl('_TN_behavior',all_VCs_full$trait1)
#qqunif(all_VCs_full[g_TN_behavior,'pvalue_DGE'],ci=T,las=1, main ='TN behavior: 950', ylim = c(0,mox), col = cols[g_TN_behavior], pch = pchs[g_TN_behavior])
#g_TN_breeder = grepl('_TN_breeder',all_VCs_full$trait1)
#qqunif(all_VCs_full[g_TN_breeder,'pvalue_DGE'],ci=T,las=1, main ='TN breeders: 555 rats', ylim = c(0,mox), col = cols[g_TN_breeder], pch = pchs[g_TN_breeder])
qqunif(all_VCs_full[,'pvalue_DGE'],ci=T,las=1, main ='All cohorts mega-analysis: 3,767 rats', ylim = c(0,mox), col = cols, pch = pchs)
dev.off()

all_VCs_full = all_VCs_full[all_VCs_full$cw_qvalue_DGE < 0.1,]
all_VCs_full = all_VCs_full[order(all_VCs_full$pvalue_DGE, decreasing = F),]
all_VCs_full

#root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
#load(paste(root_dir,'significance_DGE_VC.RData',sep=''))
#all_VCs_full


##### need to relativise per taxonomic level otherwise wrong + leads to invariant means below
relativise = function(col) { #per sample
	col = col/sum(col)
}

root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(paste(root_dir,'augmented_VC.RData',sep=''))
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/study_spe_taxa.RData')
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/study_spe_uncollapsed_taxa.RData')
pdf('/users/abaud/abaud/P50_HSrats/plots/relationship_abund_pvalues_DGE_Helenes.pdf',bg='white')
par(mfrow = c(2,2), mar = c(4,4,4,1))
for (study in c('NY','MI','TN_behavior','TN_breeder')) {
	print(study)
	#this_biomt = get(paste('filtered_full_biomt_',study,sep=''))
	this_biomt = get(paste('filtered_collapsed_full_biomt_',study,sep=''))
	this_biomt = apply(this_biomt, FUN = relativise, MAR = 2)

	grop = grep(paste('_',study, sep=''),all_VCs_full$trait1)
	motch = match(all_VCs_full[grop,'trait1'],paste(rownames(this_biomt),'_',study,sep=''))
	#7 for compositions
	this_biomt = this_biomt[motch,]
	means = apply(this_biomt, FUN = mean, MAR = 1)

	plot(means, all_VCs_full[grop,'prop_Ad1'], xlab = 'Mean relative abundance', ylab = 'Heritability', main = study, pch = 16, las = 1)
	print(cor.test(means, all_VCs_full[grop,'prop_Ad1'], method = 'pearson'))
	print(cor.test(means, all_VCs_full[grop,'prop_Ad1'], method = 'spearman'))
	try(print(t.test(means~(all_VCs_full[grop,'cw_qvalue_DGE'] < 0.1))))
}
dev.off()


load('/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_DGE_VC_wALL.RData')
taxa_sigs_all = all_VCs_full[which(all_VCs_full$study1 == 'all' & all_VCs_full$cw_qvalue_DGE < 0.1),'taxon1']
nbs = c()
for (taxon in taxa_sigs_all) {
#	w = which(all_VCs_full$taxon1 == taxon & all_VCs_full$study1 != 'all' & all_VCs_full$cw_qvalue_DGE < 0.1)
	w = which(all_VCs_full$taxon1 == taxon & all_VCs_full$study1 != 'all' & all_VCs_full$prop_Ad1 > 0.1)
	nbs = c(nbs, length(w))
}
names(nbs) = taxa_sigs_all
table(nbs)
  0   1   2 
164  49   1 




# use new mixture chi square to get final IGE P values
load("/nfs/users/abaud/abaud/outputs/P50_HSrats/VD/univariate/augmented_IGE_VC.RData")
remove = colnames(all_VCs_full)[grep('value_DGE', colnames(all_VCs_full))]
remove
all_VCs_full = all_VCs_full[,-match(remove, colnames(all_VCs_full))]
k = 0.7
all_VCs_full$pvalue_IGE = k*pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=1, lower.tail = FALSE) + (1-k)*pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=2, lower.tail = FALSE)

library(qvalue)
all_VCs_full$sw_qvalue_IGE = qvalue(all_VCs_full$pvalue_DGE)$qvalue

all_VCs_full$cw_bonf_pvalue_IGE = NA
all_VCs_full$cw_qvalue_IGE = NA
#for (study in c('all')) {
for (study in c('all','MI','NY','TN_breeder','TN_behavior')) {
	w = which(all_VCs_full$study1 == study)
	all_VCs_full[w,'cw_bonf_pvalue_IGE'] = all_VCs_full[w,'pvalue_IGE']*length(w) #cw = cohort (ie study) wide
	all_VCs_full[w,'cw_qvalue_IGE'] = qvalue(all_VCs_full[w,'pvalue_IGE'])$qvalue
}
save(all_VCs_full, file =  '/nfs/users/abaud/abaud/outputs/P50_HSrats/VD/univariate/augmented_IGE_VC.RData')







