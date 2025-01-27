
#### Data were prepared in 'prev_abund_herit_dataPrep.R'

####  now will draw the figure for the paper

#load both ASV level and taxa level microbiome data, and assign common variable name for subsequent use
load('/users/abaud/abaud/P50_HSrats/output/prev_abund_asvs_biomt.RData')
load('/users/abaud/abaud/P50_HSrats/output/prev_abund_taxa_biomt.RData')
for (study in c('MI','NY','TN_behavior','TN_breeder')) {
	assign(paste("prevs", study, sep='_'), c(get(paste("prevalence", study, sep='_')), get(paste("collapsed_prevalence", study, sep='_'))))
	assign(paste("meds", study, sep='_'), c(get(paste("median", study, sep='_')), get(paste("collapsed_median", study, sep='_'))))
}

#load heritability data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))

# Preparing for plot
library(corrplot) # needed for colorlegend and COL1
# Choose 'estimate' - if based on heritability estimate; 
# Choose 'pval' if based on heritability significance
type=c("estimate") 
#type=c("pval") 

###### Helene: #######
#all_VCs_full$color_Ad1 = NA
if(type == "estimate"){
  val_oi = all_VCs_full$prop_Ad1
  lg_title = "Heritability"
}else if(type == "pval"){
  val_oi = -log10(all_VCs_full$pvalue_DGE)
  lg_title = "-logP"
}else{
  stop("choose between 'estimate' or 'pval'")
}
uniq_val = sort(unique(val_oi)) 
# here define the colors
#colours = heat.colors(length(uniq_val))[length(uniq_val):1] 
#colours = COL1('Blues',n = length(uniq_val)) # warm blue
colours = colorRampPalette(c("#ffffff","#3C5488FF"))(length(uniq_val)) # npg 
motch = match(val_oi, uniq_val)

# now sets color
all_VCs_full[,"color_Ad1"] = colours[motch] 

labels = unname(quantile(val_oi)[c(1,3,5)])



#now plot
####pdf('/users/abaud/abaud/P50_HSrats/plots/prev_abund_herit_logP_biomt.pdf')
outpdf = paste0("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/prev_abund_herit_",type,"_biomt.pdf"); cat("saving pdf to: ", outpdf, "\n")
pdf(outpdf, h=6, w = 6)
par(mar = c(5.1,5.1,2.1,2.1))
#all_prevs = c() #previously used to have all cohorts in one plot
#all_meds = c() #previously used to have all cohorts in one plot
#all_cols = c() #previously used to have all cohorts in one plot
#par(mfrow=c(2,2))
for (study in c('MI','NY','TN_behavior','TN_breeder')) { #to have one plot per cohort
	prevs = get(paste("prevs", study, sep='_'))
	motch = match(paste(names(prevs), study, sep='_'), all_VCs_full$trait1)
	cols = all_VCs_full[na.exclude(motch),'color_Ad1']
	#all_cols = c(all_cols, cols) #previously used to have all cohorts in one plot
	#all_prevs = c(all_prevs, prevs[!is.na(motch)]) #previously used to have all cohorts in one plot
	meds = get(paste("meds", study, sep='_'))
	#all_meds = c(all_meds, meds[!is.na(motch)]) #previously used to have all cohorts in one plot
	#	plot(prevs, meds, col= cols) #previously used to have all cohorts in one plot
	
	lg_xlim = c(0.01,0.1) # this should be the same for all
	ystart = round(max(meds), 2)-max(meds)/25
	lg_ylim = c(ystart-max(meds)/5, ystart)
	plot(prevs, meds, col= cols, 
	     xlab = paste0("Prevalence in ",study," cohort"), ylab = "", 
	     las = 1, cex.lab = 1.4, cex.axis = 1.25, #main = paste('Cohort:', study), 
	     pch = 16)
	title(ylab = paste0("Median relative abundance in ",study," cohort"), cex.lab = 1.4,
	      line = 3.5)

	colorlegend(colbar = rev(colours), labels = rev(round(labels, 3)),
	            ratio.colbar = 0.4, cex = 0.9,
	            xlim = lg_xlim, ylim = lg_ylim, vertical = TRUE,
	            align = "l")
	text(x=lg_xlim[2], 
	     #lg_ylim[1]-max(meds)/25, 
	     lg_ylim[2]+max(meds)/20, 
	     labels=paste0(lg_title,":"), adj=0.5)
	
}
dev.off()



########### Amelie #######
###### Amelie: ########
## all_VCs_full$color_Ad1 = NA
## uniq_herits = sort(unique(all_VCs_full$prop_Ad1))
## uniq_logPs = sort(unique(-log10(all_VCs_full$pvalue_DGE)))
## 
## #to set colour based on heritability estimate
## colours = heat.colors(length(uniq_herits))[length(uniq_herits):1]
## motch = match(all_VCs_full$prop_Ad1, uniq_herits)
## 
## #to set colour based on heritability significance
## #colours = heat.colors(length(uniq_logPs))[length(uniq_logPs):1]
## #motch = match(-log10(all_VCs_full$pvalue_DGE), uniq_logPs)

## # now sets color; depends on what is uncommented above
## #all_VCs_full$color_Ad1 = colours[motch] 

## #all_VCs_full = all_VCs_full[grepl('ASV', all_VCs_full$trait1),] #if want to plot ASVs only (not taxa)

## #create legend
## legend_herits = c()
## legend_cols = c()
## k_max = 3 #how many items we want in the legend
## for (k in 1:k_max) {
##   #if plot based on heritability estimate
##   w = which.min(abs(all_VCs_full$prop_Ad1 - quantile(all_VCs_full$prop_Ad1, probs = seq(0,1,length.out = k_max))[k]))
##   legend_herits = c(legend_herits, all_VCs_full[w,'prop_Ad1'])
##   
##   #if plot based on heritability significance
##   #w = which.min(abs(-log10(all_VCs_full$pvalue_DGE) - quantile(-log10(all_VCs_full$pvalue_DGE), probs = seq(0,1,length.out = k_max))[k]))
##   #legend_herits = c(legend_herits, -log10(all_VCs_full[w,'pvalue_DGE']))
##   
##   legend_cols = c(legend_cols, all_VCs_full[w,'color_Ad1'])
## }
## legend_herits = round(legend_herits, digits = 3)
## #legend_herits=c(0,0.009,0.179) # as in paper

## par(mar = c(4,4,1,1))
## #all_prevs = c() #previously used to have all cohorts in one plot
## #all_meds = c() #previously used to have all cohorts in one plot
## #all_cols = c() #previously used to have all cohorts in one plot
## for (study in c('MI','NY','TN_behavior','TN_breeder')) { #to have one plot per cohort
##   prevs = get(paste("prevs", study, sep='_'))
##   motch = match(paste(names(prevs), study, sep='_'), all_VCs_full$trait1)
##   cols = all_VCs_full[na.exclude(motch),'color_Ad1']
##   #all_cols = c(all_cols, cols) #previously used to have all cohorts in one plot
##   #all_prevs = c(all_prevs, prevs[!is.na(motch)]) #previously used to have all cohorts in one plot
##   meds = get(paste("meds", study, sep='_'))
##   #all_meds = c(all_meds, meds[!is.na(motch)]) #previously used to have all cohorts in one plot
##   #	plot(prevs, meds, col= cols) #previously used to have all cohorts in one plot
##   
##   plot(prevs, meds, col= cols, xlab = "Prevalence in the cohort", ylab = "Median relative abundance in the cohort", las = 1, cex.lab = 1.4, cex.axis = 1.2, main = paste('Cohort:', study), pch = 16)
##   legend("topleft", legend = legend_herits, fill = legend_cols, bty = 'n', cex = 1.4, title = 'Heritability:') #if colour based on heritability estimate
##   #legend("topleft", legend = legend_herits, fill = legend_cols, bty = 'n', cex = 1.4, title = '-logP:') #if colour based on heritability significance
## }

