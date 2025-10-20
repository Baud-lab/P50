
# Load heritability data
load('augmented_DGE_VC_wALL.RData')
#filtering out results for "all" - focus on different centers
all_VCs_full = all_VCs_full[all_VCs_full$study1 != "all",]

#### Data were prepared in 'prev_abund_herit_dataPrep.R'

# Load both ASV level and taxa level microbiome data, and assign common variable name for subsequent use
load('prev_abund_asvs_biomt.RData') # ASVs
load('prev_abund_taxa_biomt.RData') # Taxa
#all (names(prevs) == names (meds))
#TRUE

for (study in c('MI','NY','TN_behavior','TN_breeder')) {
  print(all (names(paste("prevs", study, sep='_')) == names (paste("prevs", study, sep='_'))))
  # TRUE all
  assign(paste("prevs", study, sep='_'), c(get(paste("prevalence", study, sep='_')), get(paste("collapsed_prevalence", study, sep='_'))))
  assign(paste("meds", study, sep='_'), c(get(paste("median", study, sep='_')), get(paste("collapsed_median", study, sep='_'))))
}

# saving objects to plot
#save(list= c(paste("prevs", c('MI','NY','TN_behavior','TN_breeder'), sep='_'), 
#     paste("meds", c('MI','NY','TN_behavior','TN_breeder'), sep='_'), 
#     "all_VCs_full"), file = "source_files/fig3b_suppFig5.RData")

# loading objects to plot from source
#load("source_files/fig3b_suppFig5.RData")

# Choose 'estimate' - if based on heritability estimate; 
type=c("estimate") 
## TODO: Choose 'pval' if based on heritability significance
#type=c("pval") 

# Open pdf to save plot
library(corrplot) # needed for colorlegend and COL1 when plotting
outpdf = paste0("prev_abund_herit_",type,"_biomt.pdf"); cat("saving pdf to: ", outpdf, "\n")
pdf(outpdf, h=6, w = 6)
par(mar = c(5.1,5.1,2.1,2.1))
#all_prevs = c() #previously used to have all cohorts in one plot
#all_meds = c() #previously used to have all cohorts in one plot
#all_cols = c() #previously used to have all cohorts in one plot
#par(mfrow=c(2,2))

# Define cohorts names as in paper
dict = c("NY" = "NY", "MI"="MI", "TN_behavior"="TN1", "TN_breeder"="TN2")

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

#check = all_VCs_full[gsub("_MI|_NY|_TN_behavior|_TN_breeder","", all_VCs_full$trait1) %in% names(meds_NY[which(meds_NY > 0.10)]),c("trait1","prop_Ad1")]
#check[grep("_NY", check$trait1),]

for (study in c('MI','NY','TN_behavior','TN_breeder')) { #to have one plot per cohort
  # selecting cohort name as for title
  studytitle = dict[study]
  
  prevs = get(paste("prevs", study, sep='_'))
	motch = match(paste(names(prevs), study, sep='_'), all_VCs_full$trait1)
	cols = all_VCs_full[na.exclude(motch),'color_Ad1']
	prevs = prevs[!is.na(motch)]
	meds = get(paste("meds", study, sep='_'))
	meds = meds[!is.na(motch)]
	
	# Plot 
	plot(prevs, meds, col= cols, 
	     xlab = paste0("Prevalence in ",studytitle," cohort"), ylab = "", 
	     las = 1, cex.lab = 1.4, cex.axis = 1.25,
	     pch = 16)
	# Add y lab
	title(ylab = paste0("Median relative abundance in ",studytitle," cohort"), cex.lab = 1.4,
	      line = 3.5)

	# Legend
	# define position on the x axis
	lg_xlim = c(0.51,0.55) # this should be the same for all
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
