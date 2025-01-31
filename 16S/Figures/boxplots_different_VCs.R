#change ASV etc at bottom

#Need a df with columns Vcage, Vmom, VA, VR, heritability (=VA / (VA + Vmom + Vcage + VR))
#change trait to phenotype 
#have heritability column
#Rows should be single taxa and something like "pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc", "asv_richness", "asv_shannon_h"
#single taxa should be as p_dhjsfgsh
#remove kingdom level if there!

library(dplyr)
library(reshape2)
library(ggsci) # for the colours of pal_npg
#library(ggplot2)
#library(cowplot)

#load heritability data - should start from augmented_DGE_VC_wALL.RData instead, keeping only ASV and taxa level phenotypes (no community phenotype)
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
#for ASVs
deblur_counts_uncollapsed_dir ='deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(file.path(root_dir,deblur_counts_uncollapsed_dir,'all_estNste.Rdata'))
asv_VCs = VCs
#for taxa
deblur_counts_dir ='deblur_counts/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(file.path(root_dir,deblur_counts_dir,'all_estNste.Rdata'))
tax_levels_VCs = VCs
#merge
all(colnames(tax_levels_VCs) == colnames(asv_VCs))
all_VCs = rbind(tax_levels_VCs, asv_VCs)

#assign cohort. instead should load heritability data from augmented_DGE_VC_wALL.RData and it would already be there
all_VCs$study = NA
for (i in 1:dim(all_VCs)[1]) {
  if(grepl('_MI$', all_VCs[i,'trait1'])) all_VCs[i,'study'] = 'MI'
  if(grepl('_NY$', all_VCs[i,'trait1'])) all_VCs[i,'study'] = 'NY'
  if(grepl('_TN_behavior$', all_VCs[i,'trait1'])) all_VCs[i,'study'] ="TN1" #'TN_behavior'
  if(grepl('_TN_breeder$', all_VCs[i,'trait1'])) all_VCs[i,'study'] = "TN2" #'TN_breeder'
}

options(warn = 2)

#change some variable names
colnames(all_VCs)[colnames(all_VCs) == 'prop_Ad1'] = 'VA'
colnames(all_VCs)[colnames(all_VCs) == 'prop_C1'] = 'Vcage'
colnames(all_VCs)[colnames(all_VCs) == 'prop_Dm1'] = 'Vmom'
colnames(all_VCs)[colnames(all_VCs) == 'prop_Ed1'] = 'VR'
colnames(all_VCs)[colnames(all_VCs) == 'trait1'] = 'phenotype'
all_VCs$heritability = all_VCs[,'VA']

#filter out a few microbial phenotypes that are too high level
all_VCs$phenotype = sub('__','_', all_VCs$phenotype)
all_VCs = all_VCs[!grepl('^k_',all_VCs$phenotype),] #remove the kingdom level phenotypes
all_VCs = all_VCs[!grepl('^d_', all_VCs$phenotype),] #remove the domain level phenotypes

#code reused from code used for Suppl fig (originally from Grieneisen et al)
mods <- all_VCs
mods3temp <- mods %>%
  mutate(Cage = 100*(Vcage/(Vcage+Vmom+VA+VR)),
         Maternal = 100*(Vmom/(Vcage+Vmom+VA+VR)),
         `Additive genetic` = 100*(VA/(Vcage+Vmom+VA+VR)),
  )
mods3 = mods3temp[,c('phenotype','study', 'Cage', 'Maternal', 'Additive genetic')]

#melt the data to get the Components as a column
mods3 <- reshape2::melt(mods3)
colnames(mods3) <- c("phenotype",'study', "component", "est")

plotdat <- mods3

# reorder levels so they plot in the wanted order
comp_order = c('Additive genetic','Maternal','Cage')
#study_order = c('TN_behavior','MI','TN_breeder','NY')
#study_order = c('TN_breeder','TN_behavior','MI','NY') # order as in panel A of the same figure
study_order = c('NY','MI','TN1','TN2') # order in order of cohort sample size

plotdat$component <- factor(plotdat$component,levels = comp_order)
plotdat$study <- factor(plotdat$study,levels =study_order)

#plot
pal = #c("#396C93", "#E64B35FF","#F39B7FFF") # darkblue[2] from panel A and npg reds
  #c("#4DBBD5FF","#E64B35FF","#F39B7FFF") # npg lightblue and reds  
  c("#3C5488FF","#00A087FF","#91D1C2FF") # npg darkblue and greens  
  #c("#3C5488FF","#E64B35FF","#F39B7FFF") # npg darkblue and reds  # my favorite
#c("#8491B4FF", "#E64B35FF","#F39B7FFF") # npg blue[2] and reds  


ats = c(1:3, 5:7, 9:11, 13:15) #where to put boxplots

#pdf(paste('/nfs/users/abaud/abaud/P50_HSrats/plots/VCs_merged.pdf',sep=''), width = 10)
#pdf(paste('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/VCs_merged.pdf',sep=''), width = 10, h=6)

library(vioplot)
# Option boxplot with violinplot - which has a meaning in showing all values
# - can be adapted removing whiskers etc...
formula = as.formula(plotdat$est ~ plotdat$component + plotdat$study)

#pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/VCs_merged_viopl_white.pdf', width = 7, h=6)
#vio_col = "white"
#box_col  = c(rep(pal, length(study_order)))

pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/VCs_merged_viopl_col.pdf', width = 7, h=6)
vio_col = rep(pal, length(study_order))
box_col = "white"

par(mar=c(5.1,5.1,2.1,2.1))
vioplot(formula, col = vio_col,#col = rep(pal, 3),
        at = ats, names = rep('', length(ats)), 
        las = 1, xlab = '', xaxt = "n", ylab='',
        cex.axis = 1.25,  
        border = "black", 
        drawRect=F,
        wex=1,
        areaEqual = F
        #ylim = c(0,25)
); boxplot(formula, col = box_col, #col="white",
        boxwex = 0.3,
        xaxt = "n", yaxt="n", box="n", 
        outline = F, 
        #border = NULL,
        #range = 0.00000000001, 
        at = ats, #ylim = c(0,25), 
        names = rep('', length(ats)), 
        las = 1,  
        add=T,
        whisklty=3,
        frame=F)
axis(side = 1, at = c(2,6,10,14), labels = study_order,tick = FALSE, cex.axis=1.25)
title(ylab = "Proportion variance explained", cex.lab = 1.4,
      line = 3.5)
legend(x = 'topleft',  legend = comp_order, fill = pal, cex = 1.2, bg="white", bty="n", inset = c(-0.01,-0.01))
dev.off()


# Option boxplot only
pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/VCs_merged_bplot.pdf', width = 7, h=6)
par(mar=c(5.1,5.1,2.1,2.1))
box_col  = c(rep(pal, length(study_order)))
boxplot(formula, col = box_col, #col="white",
        #boxwex = 0.3,
        xaxt = "n" , 
        outline = F, 
        #border = NULL, # don't see any difference
        range = 0.00000000001, 
        at = ats, ylim = c(0,25), 
        names = rep('', length(ats)), 
        las = 1, ylab = 'Proportion variance explained', xlab = '', cex.lab = 1.4,
        cex.axis=1.25) #, 
        #add=T,
        #whisklty=3)
axis(side = 1, at = c(2,6,10,14), labels = study_order,tick = FALSE, cex.axis = 1.25)
legend(x = 'topleft',  legend = comp_order, fill = pal, border = NULL, cex = 1.25)

dev.off()





####### Amelie #######
## pal <- c('#225ea8','#016450','#a1dab4') #use same colors throughout the manuscript; here same as Grieneisen et al
## ats = c(1:3,5:7,9:11, 13:15) #where to put boxplots
## #pdf(paste('/nfs/users/abaud/abaud/P50_HSrats/plots/VCs_merged.pdf',sep=''), width = 10)
## boxplot(plotdat$est ~ plotdat$component + plotdat$study, col = rep(pal, 3), 
##         xaxt = "n" , outline = FALSE, border = NULL,range = 0.00000000001, 
##         at = ats, ylim = c(0,25), names = rep('', length(ats)), las = 1, ylab = 'Proportion variance explained', xlab = '')
## axis(side = 1, at = c(2,6,10,14), labels = study_order,tick = FALSE)
## legend(x = 'top',  legend = comp_order, fill = pal, border = NULL)
## #dev.off()


