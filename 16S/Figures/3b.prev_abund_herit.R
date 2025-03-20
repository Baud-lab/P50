library(corrplot) # needed for colorlegend and COL1 when plotting

#### Data were prepared in 'prev_abund_herit_dataPrep.R'

# Load both ASV level and taxa level microbiome data, and assign common variable name for subsequent use
load('/users/abaud/abaud/P50_HSrats/output/prev_abund_asvs_biomt.RData') # ASVs
load('/users/abaud/abaud/P50_HSrats/output/prev_abund_taxa_biomt.RData') # Taxa
for (study in c('MI','NY','TN_behavior','TN_breeder')) {
	assign(paste("prevs", study, sep='_'), c(get(paste("prevalence", study, sep='_')), get(paste("collapsed_prevalence", study, sep='_'))))
	assign(paste("meds", study, sep='_'), c(get(paste("median", study, sep='_')), get(paste("collapsed_median", study, sep='_'))))
}

# Load heritability data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))

# Choose 'estimate' - if based on heritability estimate; 
type=c("estimate") 
## TODO: Choose 'pval' if based on heritability significance
#type=c("pval") 

# Prepare for legend
# title depending on what plotting
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
colours = COL1('Blues',n = length(uniq_val)) # warm blue
motch = match(val_oi, uniq_val)
# now set color
all_VCs_full[,"color_Ad1"] = colours[motch] 
# now set the tick labels
labels = unname(quantile(val_oi)[c(1,3,5)])


# Open pdf to save plot
outpdf = paste0("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/prev_abund_herit_",type,"_biomt.pdf"); cat("saving pdf to: ", outpdf, "\n")
pdf(outpdf, h=6, w = 6)
par(mar = c(5.1,5.1,2.1,2.1))
#all_prevs = c() #previously used to have all cohorts in one plot
#all_meds = c() #previously used to have all cohorts in one plot
#all_cols = c() #previously used to have all cohorts in one plot
#par(mfrow=c(2,2))

# Define cohorts names as in paper
dict = c("NY" = "NY", "MI"="MI", "TN_behavior"="TN1", "TN_breeder"="TN2")

for (study in c('MI','NY','TN_behavior','TN_breeder')) { #to have one plot per cohort
  # selecting cohort name as for title
  studytitle = dict[study]
  
  # - Amelie comments - 
  prevs = get(paste("prevs", study, sep='_'))
	motch = match(paste(names(prevs), study, sep='_'), all_VCs_full$trait1)
	cols = all_VCs_full[na.exclude(motch),'color_Ad1']
	#all_cols = c(all_cols, cols) #previously used to have all cohorts in one plot
	#all_prevs = c(all_prevs, prevs[!is.na(motch)]) #previously used to have all cohorts in one plot
	meds = get(paste("meds", study, sep='_'))
	#all_meds = c(all_meds, meds[!is.na(motch)]) #previously used to have all cohorts in one plot
	#	plot(prevs, meds, col= cols) #previously used to have all cohorts in one plot
	# - end Amelie comments - 
	
	# Plot
	plot(prevs, meds, col= cols, 
	     xlab = paste0("Prevalence in ",studytitle," cohort"), ylab = "", 
	     las = 1, cex.lab = 1.4, cex.axis = 1.25,
	     pch = 16)
	# Add y lab
	title(ylab = paste0("Median relative abundance in ",studytitle," cohort"), cex.lab = 1.4,
	      line = 3.5)

	# Legend
	# define position
	lg_xlim = c(0.01,0.1) # this should be the same for all
	ystart = round(max(meds), 2)-max(meds)/25
	lg_ylim = c(ystart-max(meds)/5, ystart)
	# add legend bar
	colorlegend(colbar = rev(colours), labels = rev(round(labels, 3)),
	            ratio.colbar = 0.4, cex = 0.9,
	            xlim = lg_xlim, ylim = lg_ylim, vertical = TRUE,
	            align = "l")
	# add legend title
	text(x=lg_xlim[2], 
	     lg_ylim[2]+max(meds)/20, 
	     labels=paste0(lg_title,":"), adj=0.5)
	
}
dev.off()
