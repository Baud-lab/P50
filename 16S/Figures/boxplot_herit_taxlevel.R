library(dplyr)
library(reshape2)
library(scales)

# Loading and arranging data 
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))
save_all_VCs = all_VCs_full

data_prep = function(study){
  # Data operations:
  all_VCs = save_all_VCs
  all_VCs = all_VCs[grep(study,all_VCs$trait1),]
  
  colnames(all_VCs)[colnames(all_VCs) == 'prop_Ad1'] = 'VA'
  colnames(all_VCs)[colnames(all_VCs) == 'prop_C1'] = 'Vcage'
  colnames(all_VCs)[colnames(all_VCs) == 'prop_Dm1'] = 'Vmom'
  colnames(all_VCs)[colnames(all_VCs) == 'prop_Ed1'] = 'VR'
  colnames(all_VCs)[colnames(all_VCs) == 'trait1'] = 'phenotype'
  all_VCs$heritability = all_VCs[,'VA']
  
  all_VCs$phenotype = sub('__','_', all_VCs$phenotype)
  all_VCs$phenotype = sub(paste('_',study,sep=''),'', all_VCs$phenotype)
  all_VCs = all_VCs[!grepl('^k_',all_VCs$phenotype),]
  all_VCs = all_VCs[!grepl('^d_', all_VCs$phenotype),]
  all_VCs[all_VCs$phenotype == 'shannon_entropy','phenotype'] = 'shannon'
  all_VCs[all_VCs$phenotype == 'observed_features','phenotype'] = 'richness'
  
  # 283 single-taxa phenotypes + 7 community phenotypes
  mods <- all_VCs
  mods3temp <- mods %>%
    mutate(Cage = 100*(Vcage/(Vcage+Vmom+VA+VR)),
           Maternal = 100*(Vmom/(Vcage+Vmom+VA+VR)),
           `Additive genetic` = 100*(VA/(Vcage+Vmom+VA+VR)),
    )
  mods3 = mods3temp[,c('phenotype', 'Cage', 'Maternal', 'Additive genetic')]
  
  # Melt the data to get the Components as a column
  mods3 <- reshape2::melt(mods3) #the reshaping leads to the loss of any column that isn't one of the 3 above
  colnames(mods3) <- c("phenotype", "component", "est")
  
  #add on the p value column
  #mods2 <- subset(mods, select = c("phenotype", "p_adjust"))
  #mods3 <- merge(mods2, mods3, by.x = "phenotype")
  #colnames(mods3) <- c("Model", "pedigree_p", "component", "est")
  
  plotdat <- mods3
  
  # Need to reorder levels so they plot in the right order left to right: Cage, maternal, baboon id, additive genetic
  plotdat$component <- as.factor(plotdat$component)
  plotdat$component <- factor(plotdat$component,levels(plotdat$component))
  
  # Add a column for if it's taa or not
  not_tax <- c("PC1", "PC2", "PC3", "PC4", "PC5", "richness", "shannon")
  
  # Reorder phenotypes so it's by highest heritability, and group by taxa vs community phenotype
  h_order <- mods
  h_order$status <- ifelse(h_order$phenotype %in% not_tax, "Non-taxa phenotype", "Taxa")
  h_order$tax_level <- substr(h_order$phenotype, 0,2)
  # Add in level, order levels
  h_order <- h_order %>%
    mutate(tax_level2 = case_when(status == "Non-taxa phenotype" ~ "composition",
                                  tax_level == "p_" ~ "phylum",
                                  tax_level == "c_" ~ "class",
                                  tax_level == "o_" ~ "order",
                                  tax_level == "f_" ~ "family",
                                  tax_level == "g_" ~ "genus",
                                  tax_level == "s_" ~ "species",
                                  tax_level == "AS" ~ "ASV")) %>%
    mutate(tax_level_numeric = case_when(status == "Non-taxa phenotype" ~ 7,
                                         tax_level == "p_" ~ 7,
                                         tax_level == "c_" ~ 6,
                                         tax_level == "o_" ~ 5,
                                         tax_level == "f_" ~ 4,
                                         tax_level == "g_" ~ 3,
                                         tax_level == "s_" ~ 2,
                                         tax_level == "AS" ~ 1)) %>%
    arrange(tax_level_numeric, heritability) #heritability here is var_Ad so good
  
  
  h_order = group_by(h_order, status) 
  
  h_order = h_order %>% dplyr::mutate(variable_order = 1:n()) %>% dplyr::select(phenotype, variable_order, status, tax_level2, heritability, Vmom, Vcage,cw_qvalue_DGE) 
  
  return(list(h_order = h_order, plotdat = plotdat))
}

taxon_box <- function(h_order){
  # Box plot for supp fig 3
  # Define cohort name as in paper
  dict = c("NY" = "NY", "MI" = "MI", "TN_behavior" ="TN1", "TN_breeder" = "TN2")
  # define order - as in plot
  h_order$tax_level2 = factor(h_order$tax_level2, levels = c('phylum','class','order','family','genus','species','ASV'))
  
  # Plot 
  # define plot limits by creating a blank boxplot
  bp = boxplot(h_order$heritability ~ h_order$tax_level2, 
               col = alpha("white", 0), border=alpha("white", 0), 
               ylab = "",
               xlab = '', 
               main = dict[study], 
               varwidth = TRUE,
               las = 1, cex.axis = 1.2, cex.main = 1.5, type="n", )
  # define x coordinates for points
  at.dict = seq_along(levels(h_order$tax_level2)); names(at.dict) = levels(h_order$tax_level2)
  # add points with jitter
  set.seed(20); points(x = jitter(unname(at.dict[h_order$tax_level2]), factor=0.5),
                       y = h_order$heritability, #jitter(counts_oi[,1], factor=1), 
                       pch=16, cex=1, 
                       col = "#4DBBD5FF")
  # add boxplots
  bp = boxplot(h_order$heritability ~ h_order$tax_level2, 
               col = alpha("white", 0.5),
               outline=F, 
               ylab = "",
               xlab = '', 
               main = study, 
               varwidth = TRUE,
               las = 1, cex.axis = 1.2, cex.main = 1.5, add=T)
  # add y lab
  title(ylab ="Heritability", cex.lab=1.4, line=3.8)
}

# Open pdf to save plots
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/herit_per_tax_level.pdf", h = 12,w = 16)
par(mar=c(5.1,5.1,2.1,2.1), mfrow=c(2,2))

for (study in c("NY","MI","TN_behavior","TN_breeder")){
  res = data_prep(study)
  taxon_box(res$h_order)
}
dev.off()