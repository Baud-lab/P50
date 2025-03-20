#change ASV etc at bottom

#### Figure 2C ####

#Need a df with columns Vcage, Vmom, VA, VR, heritability (=VA / (VA + Vmom + Vcage + VR))
#change trait to phenotype 
#have heritability column
#Rows should be single taxa and something like "pc1_bc", "pc2_bc", "pc3_bc", "pc4_bc", "pc5_bc", "asv_richness", "asv_shannon_h"
#single taxa should be as p_dhjsfgsh
#remove kingdom level if there!

options(warn = 2)

library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)

root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))
save_all_VCs = all_VCs_full

my_f = function(signif_taxon) {
   g = grep(signif_taxon,plotdat$phenotype)[1]
  annotate("text", x = plotdat[g,'variable_order']-1, y = sum(plotdat[g,c('heritability','Vmom','Vcage')]*100)+1, size = 5, label = "*")
}

#pdf(paste('/nfs/users/abaud/abaud/P50_HSrats/plots/herit_per_tax_level.pdf',sep=''), width = 15, height = 10)
par(mfrow = c(2,2))
for (study in c('MI','NY','TN_behavior','TN_breeder')) {
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
  
  #melt the data to get the Components as a column
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
  
  #add a column for if it's taa or not
  not_tax <- c("PC1", "PC2", "PC3", "PC4", "PC5", "richness", "shannon")
  
  #reorder phenotypes so it's by highest heritability, and group by taxa vs community phenotype
  h_order <- mods
  h_order$status <- ifelse(h_order$phenotype %in% not_tax, "Non-taxa phenotype", "Taxa")
  h_order$tax_level <- substr(h_order$phenotype, 0,2)
  #add in level, order levels
  
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
  
  
  ### now merge plotdat and h_order; last 3 cols of h_order will not be used; what will be used is order, status etc in h_order
  plotdat <- merge(plotdat, h_order, by.x = "phenotype", by.y = "phenotype")
  
  #plot(1:4,1:4,pch = 16, col = c('#a1dab4','#016450','#41b6c4','#225ea8'))
  #pal <- c('#a1dab4','#016450','#225ea8')
  pal <- c("#3C5488FF","#00A087FF","#91D1C2FF") 
  
  ## no legend
  th2 <- theme(
    axis.line.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.x = element_text(color="black", size=15),
    axis.text.x = element_text(color="black", size=15),
    axis.ticks.y = element_blank(),
    strip.background = element_blank() #facet laabel background white, not grey
    , strip.text.y =  element_blank() #remove teh facet labels
    , strip.text.x = element_text()
  )
  
  #community phenotype plot
  y_max1 = max(apply(subset(plotdat, status == "Non-taxa phenotype")[,c('heritability','Vmom','Vcage')], FUN = sum, MAR = 1))*100
  y_max1 = round_any(y_max1, 20, f = ceiling)
  y_max2 = max(apply(subset(plotdat, status != "Non-taxa phenotype")[,c('heritability','Vmom','Vcage')], FUN = sum, MAR = 1))*100
  y_max2 = round_any(y_max2, 20, f = ceiling)
  y_max = max(y_max1, y_max2)
  
  
  pC_1 <- ggplot(subset(plotdat, status == "Non-taxa phenotype"), aes(x = variable_order, y = est, fill=component)) +
    geom_bar(stat = "identity",colour="white", linewidth=0.25) +
    scale_fill_manual(values = pal) +
    guides(fill = guide_legend(reverse = TRUE, title = "Variance component")) +  #make additive genetic top
    coord_flip(clip = "off") + #wont' clip text labels
    labs(y=" ", x=" ") +
    annotate("text", x = 7+5, y = 25,label = "Community phenotype (n=7)" , color="black", size=5) +
    scale_x_discrete(expand=c(0,0)) +
    ##add in dummy lines so everything lines up
    annotate("segment", x = 2, xend = 4, y = -1, yend = -1, colour = "white") + 
    annotate("text", x = 3, y = -3, size = 3, angle = 90, label = "b", color = "white") + 
    th2 + 
    ylim(-5,y_max) + 
    lapply(unique(plotdat[plotdat$cw_qvalue_DGE <= 0.1 & plotdat$status == "Non-taxa phenotype",'phenotype']), my_f)
  
  
  line_info <- h_order %>%
    group_by(tax_level2) %>%
    summarise(line_min = min(variable_order), line_max = max(variable_order), line_mid = mean(variable_order))
  
  #taxa plot + legend
  th3 <- theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
               #axis.line.x = element_blank(), 
               axis.line.y = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               legend.title = element_text(size=14),
               legend.text = element_text(size=12),
               legend.position = c(0.7, 0.075),
  #             legend.position = c(0.6, 0.04), 
               axis.text.y = element_blank(),
               axis.title.y = element_text(color="black", size=20),
               axis.title.x = element_text(color="black", size=15),
               axis.text.x = element_text(color="black", size=15),
               axis.ticks.y = element_blank(),
               strip.background = element_blank() #facet laabel background white, not grey
               , strip.text.y =  element_blank() #remove teh facet labels
               , strip.text.x = element_text()
  )
  
  range_p = range(plotdat[grep('^p_',plotdat$phenotype),'variable_order'])
  range_c = range(plotdat[grep('^c_',plotdat$phenotype),'variable_order'])
  range_o = range(plotdat[grep('^o_',plotdat$phenotype),'variable_order'])
  range_f = range(plotdat[grep('^f_',plotdat$phenotype),'variable_order'])
  range_g = range(plotdat[grep('^g_',plotdat$phenotype),'variable_order'])
  range_s = range(plotdat[grep('^s_',plotdat$phenotype),'variable_order'])
  range_asv = range(plotdat[grep('^ASV_',plotdat$phenotype),'variable_order'])
  
  pC_2 <- ggplot(subset(plotdat, status != "Non-taxa phenotype"), aes(x = variable_order, y = est, fill=component)) +
    geom_bar(stat = "identity",colour="white", linewidth=0.25) +
    scale_fill_manual(values = pal) +
    guides(fill = guide_legend(reverse = TRUE, title = "Variance component")) +  #make additive genetic top
    coord_flip(clip = "off") + #wont' clip text labels
    labs(y="Percent variance explained", x="Microbial phenotype") +
    annotate("text", x = range_p[2]+5, y = 30,label = paste("Single-taxon phenotype (n=",range_p[2],")",sep=''), color="black", size=5) +
    scale_x_discrete(expand=c(0,0)) +
  
    annotate("segment", x = range_p[1], xend = range_p[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = range_p[2], y = -3, size = 5, angle = 90, label = "Phylum") +
    
    annotate("segment", x = range_c[1], xend = range_c[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = mean(range_c), y = -3, size = 5, angle = 90, label = "Classy") +
    
    annotate("segment", x = range_o[1], xend = range_o[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = mean(range_o), y = -3, size = 5, angle = 90, label = "Order") +
    
    annotate("segment", x = range_f[1], xend = range_f[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = mean(range_f), y = -3, size = 5, angle = 90, label = "Family") +
    
    annotate("segment", x = range_g[1], xend = range_g[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = mean(range_g), y = -3, size = 5, angle = 90, label = "Genus") +
  
    annotate("segment", x = range_s[1], xend = range_s[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = mean(range_s), y = -3, size = 5, angle = 90, label = "Species") +
    
    annotate("segment", x = range_asv[1], xend = range_asv[2], y = -1, yend = -1, colour = "darkgrey") + 
    annotate("text", x = mean(range_asv), y = -3, size = 5, angle = 90, label = "ASV") +
    
    th3 + 
    ylim(-5,y_max) + 
    
    lapply(unique(plotdat[plotdat$cw_qvalue_DGE <= 0.1 & plotdat$status != "Non-taxa phenotype",'phenotype']), my_f)
  
  
  
  #pdf(paste('/nfs/users/abaud/abaud/P50_HSrats/plots/baboon_fig_',study,'.pdf',sep=''), width = 5, height = 20)
  #print(plot_grid(pC_1, pC_2, align = "hv", ncol = 1, rel_heights = c(1,18)))
  #dev.off()
  
  boxplot(h_order$heritability ~ factor(h_order$tax_level2, levels = c('phylum','class','order','family','genus','species','ASV')), ylab = 'Heritability', xlab = '', varwidth = TRUE, las = 1, main = study)
  # boxplot(h_order$Vmom ~ factor(h_order$tax_level2, levels = c('phylum','class','order','family','genus','species','ASV')), ylab = 'Maternal effects', xlab = '',, varwidth = TRUE, las = 1)
  # boxplot(h_order$Vcage~ factor(h_order$tax_level2, levels = c('phylum','class','order','family','genus','species','ASV')), ylab = 'Cage effects', xlab = '',, varwidth = TRUE, las = 1)
  # boxplot(h_order$Vcage~ factor(h_order$tax_level2, levels = c('composition','phylum','class','order','family','genus','species','ASV')), ylab = 'Cage effects', xlab = '',, varwidth = TRUE, las = 1)

}
#dev.off()



###### Trying new thingss
study = "MI"

h_order$tax_level2 = factor(h_order$tax_level2, levels = c('phylum','class','order','family','genus','species','ASV'))
bp = boxplot(h_order$heritability ~ h_order$tax_level2, 
             col = alpha("white", 0), border=alpha("white", 0), #alpha("#4DBBD5", 0.2),
             #outline=F, #pch=16,outline=T,
             ylab = "",#paste0(ylob, ' (raw counts)'), 
             xlab = '', 
             main = study, 
             varwidth = TRUE,
             #main = "St6galnac1", 
             las = 1, cex.axis = 1.2, cex.main = 1.5, type="n", )

at.dict = seq_along(levels(h_order$tax_level2)); names(at.dict) = levels(h_order$tax_level2)
set.seed(20); points(x = jitter(unname(at.dict[h_order$tax_level2]), factor=0.5),
                     y = h_order$heritability, #jitter(counts_oi[,1], factor=1), 
                     pch=16, cex=1, 
                     col = "#4DBBD5FF")
bp = boxplot(h_order$heritability ~ h_order$tax_level2, 
             col = alpha("white", 0.5),
             outline=F, #pch=16,outline=T,
             ylab = "",#paste0(ylob, ' (raw counts)'), 
             xlab = '', 
             main = study, 
             varwidth = TRUE,
             #main = "St6galnac1", 
             las = 1, cex.axis = 1.2, cex.main = 1.5, add=T)

title(ylab ="Heritability", cex.lab=1.4, line=3.8)

