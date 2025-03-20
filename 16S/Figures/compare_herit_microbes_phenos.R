#load phenotype data
phenos_dir ='/users/abaud/abaud/P50_HSrats/output/VD/univariate/residuals_02Jan2020/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(file.path(phenos_dir,'all_estNste.Rdata'))
all_VCs_phenos = VCs

#load microbiome data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))
all_VCs_16S = all_VCs_full

#determine if phenotype is behavioural or not
non_behav = c('adams','baculum','bicknell','insulin','physiological')
behav = c('ccc','ccp','crf','delay_discounting','elevated_plus_maze',
          'light_reinforcement','locomotor_testing','nicsa','novel_object',
          'novelty_seeking','open_field','pavca','reaction_time','social_')

# for each row (MARGIN=1)
# return true if any of the behav is found in "trait1"; 
# paste(behav, collapse="|") - is the string grepped, looking for any of the words in between |
all_VCs_phenos$behavioural_pheno = apply(all_VCs_phenos, 
                                         MARGIN=1, 
                                         function(x) { grepl(paste(behav, collapse="|"), x["trait1"]) } )

#visual check
sort(all_VCs_phenos[which(all_VCs_phenos[,'behavioural_pheno']),'trait1'])
sort(all_VCs_phenos[-which(all_VCs_phenos[,'behavioural_pheno']),'trait1'])


# Define function to obtain counts for each set of phenotypes to plot
#     subdf: sub-dataframe with set of phenotypes of interest
#     var: proportional variance to plot - ideally could do this for other 'prop_'
#     br: breaks in histogram

prepare_hist = function(subdf, var = "prop_Ad1", br){
  hist_values = hist(subdf[, var], breaks = br, plot = F)
  counts = hist_values$counts/sum(nrow(subdf))
  return(counts)
}

# defining breaks; if want to compare different subset in the same plot, need breaks to be all the same
br = seq(0, 0.65, by = 0.05)

# creating dataframe to plot
toplot = rbind("NY" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'NY',], br=br),
               "MI" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'MI',], br=br),
               "TN1" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_behavior',], br=br), 
               "TN2" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_breeder',], br=br),
               "Behaviour" = prepare_hist(all_VCs_phenos[all_VCs_phenos$behavioural_pheno,], br=br),
               "Physiology" = prepare_hist(all_VCs_phenos[ ! all_VCs_phenos$behavioural_pheno,], br=br))
colnames(toplot) = br[-1]

n = 4 # numbers of microbiome studies
cols = c(colorRampPalette(c("#003A6B","#ACD0E6"))(n), "grey40", "grey75") 
#cols = c("#003A6B", "#ACD0E6", "#925E9F","#D1C4E9", "grey40", "grey75")  # another option of colours
#scales::show_col(cols) # to see colors


# Starting with plot
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/barplots_herits_studies_pheno.pdf", width = 8, h = 6)
# Set plot margin
par(mar = c(5.1, 5.1, 2.1, 2.1)) # default: c(5.1, 4.1, 4.1, 2.1)

# Set ylim 
ylimi = c(0,1)
#ylimi = range(toplot) # use this if want to be flexible depending on the data set

# plot
barplot(toplot[,1:11],
        beside = T, names = colnames(toplot[,1:11]), 
        col = cols, las = 1, 
        ylim = ylimi, 
        cex.names = 1.25, cex.axis = 1.25, cex.lab=1.4,
        ylab = 'Proportion of heritable traits', xlab = 'Heritability')
# legend
legend (x = 'topright', 
        legend = rownames(toplot), 
        fill = cols, border = NA, cex =1.2, bty="o", inset=c(0.04,0)) 
# inset so that box of legend aligns to last bin - with margins par(mar = c(5.1, 5.1, 2.1, 2.1))  and width = 8
dev.off()
