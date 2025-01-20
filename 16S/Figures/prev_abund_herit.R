#load microbiome data
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData')
#load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')

#load metadata
load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
if (any(grepl('collapsed_clr_counts', ls()))) motch = match(colnames(collapsed_clr_counts),metadata$deblur_rooname) else motch = match(colnames(clr_counts),metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

#go from raw abundance to relative abundance
my_relativise = function(row) {
	row = row / sum(row)
	return(row)
}

#depending on whether CLRed or non CLRed microbiome data were uploaded to R, grab the right object and filter based on cohort (study)
for (study in c('MI','NY','TN_behavior','TN_breeder')) {
	print(study)
	if (!any(grepl('collapsed_clr_counts', ls()))) this_clr_counts = clr_counts[,which(metadata[,"study"] == study)]   else this_clr_counts = collapsed_clr_counts[,which(metadata[,"study"] == study)]  
	if (!any(grepl('collapsed_full_biomt', ls()))) this_full_biomt = full_biomt[,which(metadata[,"study"] == study)]   else this_full_biomt = collapsed_full_biomt[,which(metadata[,"study"] == study)]  

	#this_full_biomt has samples in columns

	#calculate prevalence
	presence = (this_full_biomt != 0)
	props_presence = apply(presence, FUN = sum, MAR = 1) / dim(presence)[2]
	if (any(grepl('collapsed_clr_counts', ls()))) assign(paste("collapsed_prevalence", study, sep='_'), props_presence) else assign(paste("prevalence", study, sep='_'), props_presence)

	#calculate median relative abundance
	rel_abundances = t(apply(this_full_biomt, MAR = 2, FUN = my_relativise)) #samples in rows now
	medians = apply(rel_abundances, FUN = median, MAR =  2)
	if (any(grepl('collapsed_clr_counts', ls()))) assign(paste("collapsed_median", study, sep='_'), medians) else assign(paste("median", study, sep='_'), medians)
}

#save prevalence and median relative abundance data to disl
lost = ls() #lost is just a name for the list retrieved by ls()
toSave = lost[grepl('prevalence', lost) | grepl('median', lost)]
save(list = toSave, file = '/users/abaud/abaud/P50_HSrats/output/prev_abund_asvs_biomt.RData')
#save(list = toSave, file = '/users/abaud/abaud/P50_HSrats/output/prev_abund_taxa_biomt.RData')

#### now will draw the figure for the paper

#load both ASV level and taxa level microbiome data, and assign common variable name for subsequent use
load('/users/abaud/abaud/P50_HSrats/output/prev_abund_asvs_biomt.RData')
load('/users/abaud/abaud/P50_HSrats/output/prev_abund_taxa_biomt.RData')
for (study in c('MI','NY','TN_behavior','TN_breeder')) {
	assign(paste("prevs", study, sep='_'), c(get(paste("prevalence", study, sep='_')), get(paste("collapsed_prevalence", study, sep='_'))))
	assign(paste("meds", study, sep='_'), c(get(paste("median", study, sep='_')), get(paste("collapsed_median", study, sep='_'))))
}

#load heritability data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(paste(root_dir,'augmented_VC.RData',sep=''))

library(corrplot)
all_VCs_full$color_Ad1 = NA
uniq_herits = sort(unique(all_VCs_full$prop_Ad1))
uniq_logPs = sort(unique(-log10(all_VCs_full$pvalue_DGE)))

#to set colour based on heritability estimate
#colours = heat.colors(length(uniq_herits))[length(uniq_herits):1]
#motch = match(all_VCs_full$prop_Ad1, uniq_herits)

#to set colour based on heritability significance
colours = heat.colors(length(uniq_logPs))[length(uniq_logPs):1]
motch = match(-log10(all_VCs_full$pvalue_DGE), uniq_logPs)

#now sets color; depends on what is uncommented above
all_VCs_full$color_Ad1 = colours[motch] 

#all_VCs_full = all_VCs_full[grepl('ASV', all_VCs_full$trait1),] #if want to plot ASVs only (not taxa)

#create legend
legend_herits = c()
legend_cols = c()
k_max = 3 #how many items we want in the legend
for (k in 1:k_max) {
	#if plot based on heritability estimate
	#w = which.min(abs(all_VCs_full$prop_Ad1 - quantile(all_VCs_full$prop_Ad1, probs = seq(0,1,length.out = k_max))[k]))
	#legend_herits = c(legend_herits, all_VCs_full[w,'prop_Ad1'])

	#if plot based on heritability significance
	w = which.min(abs(-log10(all_VCs_full$pvalue_DGE) - quantile(-log10(all_VCs_full$pvalue_DGE), probs = seq(0,1,length.out = k_max))[k]))
	legend_herits = c(legend_herits, -log10(all_VCs_full[w,'pvalue_DGE']))
	legend_cols = c(legend_cols, all_VCs_full[w,'color_Ad1'])
}
legend_herits = round(legend_herits, digits = 3)

#now plot
pdf('/users/abaud/abaud/P50_HSrats/plots/prev_abund_herit_logP_biomt.pdf')
par(mar = c(4,4,1,1))
all_prevs = c() #previously used to have all cohorts in one plot
all_meds = c() #previously used to have all cohorts in one plot
all_cols = c() #previously used to have all cohorts in one plot
for (study in c('MI','NY','TN_behavior','TN_breeder')) { #to have one plot per cohort
	prevs = get(paste("prevs", study, sep='_'))
	motch = match(paste(names(prevs), study, sep='_'), all_VCs_full$trait1)
	cols = all_VCs_full[na.exclude(motch),'color_Ad1']
	#all_cols = c(all_cols, cols) #previously used to have all cohorts in one plot
	#all_prevs = c(all_prevs, prevs[!is.na(motch)]) #previously used to have all cohorts in one plot
	meds = get(paste("meds", study, sep='_'))
	#all_meds = c(all_meds, meds[!is.na(motch)]) #previously used to have all cohorts in one plot
	#	plot(prevs, meds, col= cols) #previously used to have all cohorts in one plot

	plot(prevs, meds, col= cols, xlab = "Prevalence in the cohort", ylab = "Median relative abundance in the cohort", las = 1, cex.lab = 1.4, cex.axis = 1.2, main = paste('Cohort:', study), pch = 16)
	#legend("topleft", legend = legend_herits, fill = legend_cols, bty = 'n', cex = 1.4, title = 'Heritability:') #if colour based on heritability estimate
	legend("topleft", legend = legend_herits, fill = legend_cols, bty = 'n', cex = 1.4, title = '-logP:') #if colour based on heritability significance
}
dev.off()




