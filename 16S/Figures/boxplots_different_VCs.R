#change ASV etc at bottom

#Need a df with columns Vcage, Vmom, VA, VR, heritability (=VA / (VA + Vmom + Vcage + VR))
#change trait to phenotype 
#have heritability column
#Rows should be single taxa and something like "pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc", "asv_richness", "asv_shannon_h"
#single taxa should be as p_dhjsfgsh
#remove kingdom level if there!

library(dplyr) # to rearrange data to plot
library(reshape2) # to rearrange data to plot
library(ggsci) # for the colours of pal_npg
library(vioplot) # for violin plot

# Load heritability data 
### should start from augmented_DGE_VC_wALL.RData instead, keeping only ASV and taxa level phenotypes (no community phenotype)
### - Amelie comments - keep or remove the above? - end Amelie comments -
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
# For ASVs
deblur_counts_uncollapsed_dir ='deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(file.path(root_dir,deblur_counts_uncollapsed_dir,'all_estNste.Rdata'))
asv_VCs = VCs
# For taxa
deblur_counts_dir ='deblur_counts/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(file.path(root_dir,deblur_counts_dir,'all_estNste.Rdata'))
tax_levels_VCs = VCs
# Merge
all(colnames(tax_levels_VCs) == colnames(asv_VCs))
all_VCs = rbind(tax_levels_VCs, asv_VCs)

# Assign cohort
### instead should load heritability data from augmented_DGE_VC_wALL.RData and it would already be there
### - Amelie comments - keep or remove the above? - end Amelie comments -
all_VCs$study = NA
for (i in 1:dim(all_VCs)[1]) {
  # Naming cohorst as in paper
  if(grepl('_MI$', all_VCs[i,'trait1'])) all_VCs[i,'study'] = 'MI'
  if(grepl('_NY$', all_VCs[i,'trait1'])) all_VCs[i,'study'] = 'NY'
  if(grepl('_TN_behavior$', all_VCs[i,'trait1'])) all_VCs[i,'study'] ="TN1" #'TN_behavior'
  if(grepl('_TN_breeder$', all_VCs[i,'trait1'])) all_VCs[i,'study'] = "TN2" #'TN_breeder'
}

options(warn = 2)

# Change some variable names
colnames(all_VCs)[colnames(all_VCs) == 'prop_Ad1'] = 'VA'
colnames(all_VCs)[colnames(all_VCs) == 'prop_C1'] = 'Vcage'
colnames(all_VCs)[colnames(all_VCs) == 'prop_Dm1'] = 'Vmom'
colnames(all_VCs)[colnames(all_VCs) == 'prop_Ed1'] = 'VR'
colnames(all_VCs)[colnames(all_VCs) == 'trait1'] = 'phenotype'
all_VCs$heritability = all_VCs[,'VA']

# Filter out a few microbial phenotypes that are too high level
all_VCs$phenotype = sub('__','_', all_VCs$phenotype)
all_VCs = all_VCs[!grepl('^k_',all_VCs$phenotype),] #remove the kingdom level phenotypes
all_VCs = all_VCs[!grepl('^d_', all_VCs$phenotype),] #remove the domain level phenotypes

# Code reused from code used for Suppl fig (originally from Grieneisen et al)
mods <- all_VCs
mods3temp <- mods %>%
  mutate(Cage = 100*(Vcage/(Vcage+Vmom+VA+VR)),
         Maternal = 100*(Vmom/(Vcage+Vmom+VA+VR)),
         `Additive genetic` = 100*(VA/(Vcage+Vmom+VA+VR)),
  )
mods3 = mods3temp[,c('phenotype','study', 'Cage', 'Maternal', 'Additive genetic')]

# Melt the data to get the Components as a column
mods3 <- reshape2::melt(mods3)
colnames(mods3) <- c("phenotype",'study', "component", "est")

plotdat <- mods3


# Reorder levels so they plot in the wanted order
comp_order = c('Additive genetic','Maternal','Cage')
study_order = c('NY','MI','TN1','TN2') # order depending on cohort sample size

plotdat$component <- factor(plotdat$component, levels = comp_order)
plotdat$study <- factor(plotdat$study, levels =study_order)

# Preparing to plot
# define colours for different variance components - will correspond to order in comp_order
pal =  c("#3C5488FF","#00A087FF","#91D1C2FF") # npg darkblue and greens  
# define violin plot positions on x axis
ats = c(1:3, 5:7, 9:11, 13:15) # where to put boxplots
# define formula
formula = as.formula(plotdat$est ~ plotdat$component + plotdat$study)

### Option violin plot + boxplot
# Open pdf to save plot
pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/VCs_merged_viopl_col.pdf', width = 7, h=6)
vio_col = rep(pal, length(study_order))
box_col = "white"

par(mar=c(5.1,5.1,2.1,2.1))
# Plot violinplot - coloured
# then boxplot - white
vioplot(formula, col = vio_col,
        at = ats, names = rep('', length(ats)), 
        las = 1, xlab = '', xaxt = "n", ylab='',
        cex.axis = 1.25,  
        border = "black", 
        drawRect=F,
        wex=1,
        areaEqual = F
); boxplot(formula, col = box_col,
        boxwex = 0.3,
        xaxt = "n", yaxt="n", box="n", 
        outline = F, 
        at = ats, 
        names = rep('', length(ats)), 
        las = 1,  
        add=T,
        whisklty=3,
        frame=F)
# Add ticks and ticks labels
axis(side = 1, at = c(2,6,10,14), labels = study_order,tick = FALSE, cex.axis=1.25)
# Add y label
title(ylab = "Proportion variance explained", cex.lab = 1.4,
      line = 3.5)
# Add legend
legend(x = 'topleft',  legend = comp_order, fill = pal, cex = 1.2, bg="white", bty="n", inset = c(-0.01,-0.01))

dev.off()


## ### Option boxplot only
## # Open pdf to save plot
## pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/VCs_merged_bplot.pdf', width = 7, h=6)
## par(mar=c(5.1,5.1,2.1,2.1))
## 
## box_col  = c(rep(pal, length(study_order)))
## # boxplot
## boxplot(formula, col = box_col, 
##         xaxt = "n" , 
##         outline = F, 
##         range = 0.00000000001, 
##         at = ats, ylim = c(0,25), 
##         names = rep('', length(ats)), 
##         las = 1, ylab = 'Proportion variance explained', xlab = '', cex.lab = 1.4,
##         cex.axis=1.25) 
## # Add ticks and ticks labels
## axis(side = 1, at = c(2,6,10,14), labels = study_order,tick = FALSE, cex.axis = 1.25)
## # Add legend
## legend(x = 'topleft',  legend = comp_order, fill = pal, border = NULL, cex = 1.25)
## 
## dev.off()
