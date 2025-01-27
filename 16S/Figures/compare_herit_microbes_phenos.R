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

### Helene: apply instead of for loop:
#     for each row (MARGIN=1)
#     return true if any of the behav is found in "trait1"; 
#     paste(behav, collapse="|") - is the string grepped, looking for any of the words in between |
all_VCs_phenos$behavioural_pheno = apply(all_VCs_phenos, 
                                         MARGIN=1, 
                                         function(x) { grepl(paste(behav, collapse="|"), x["trait1"]) } )

#visual check
sort(all_VCs_phenos[which(all_VCs_phenos[,'behavioural_pheno']),'trait1'])
sort(all_VCs_phenos[-which(all_VCs_phenos[,'behavioural_pheno']),'trait1'])


########### Helene plot ######## 
# define function to obtain counts for each set of phenotypes to plot
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
toplot = rbind("Microbiome TN_breeder" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_breeder',], br=br),
               "Microbiome TN_behavior" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_behavior',], br=br), 
               "Microbiome MI" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'MI',], br=br),
               "Microbiome NY" = prepare_hist(all_VCs_16S[all_VCs_16S$study1 == 'NY',], br=br),
               "Behavioural phenotypes" = prepare_hist(all_VCs_phenos[all_VCs_phenos$behavioural_pheno,], br=br),
               "Non behavioural phenotypes" = prepare_hist(all_VCs_phenos[ ! all_VCs_phenos$behavioural_pheno,], br=br))
colnames(toplot) = br[-1]

n = 4 # numbers of microbiome studies
#show_col(colorRampPalette(c("#003A6B","#ACD0E6"))(n)) # To see colours
cols = c(colorRampPalette(c("#003A6B","#ACD0E6"))(n), "grey40", "grey75") 
#show_col(cols) # need library(scales) for this

pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/barplots_herits_studies_pheno_HT.pdf", width = 8, h = 6)
# Set plot margin
par(mar = c(5.1, 5.1, 2.1, 2.1)) # default: c(5.1, 4.1, 4.1, 2.1)

# Set ylim 
ylimi = c(0,1)
#ylimi = range(toplot) # to be flexible depending on the data set

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
        fill = cols, border = NA, cex=1.3)
dev.off()



### Amelie: for loop:
# Checks: - TRUE
#H_behav = apply(all_VCs_phenos, MARGIN=1, function(x) {grepl(paste(behav, collapse="|"), x["trait1"])} )
#all(H_behav == all_VCs_phenos$behavioural_pheno) # TRUE

## all_VCs_phenos$behavioural_pheno = NA
## for (i in 1:dim(all_VCs_phenos)[1]) {
##     behavioural = FALSE
##     for (cate in behav) {
##         if (grepl(cate, all_VCs_phenos[i,'trait1'])) {behavioural = TRUE ; next} 
##     }
##     all_VCs_phenos[i,'behavioural_pheno'] = behavioural
## }

###### Amelie plot ######
### Checks:  - ALL TRUE
#all(toplot["Microbiome MI",] == hist_MI$counts/sum(all_VCs_16S$study1 == 'MI')) # TRUE
#all(toplot["Microbiome NY",] == hist_NY$counts/sum(all_VCs_16S$study1 == 'NY')) # TRUE
#all(toplot["Microbiome TN_behavior",] == hist_TN_behavior$counts/sum(all_VCs_16S$study1 == 'TN_behavior')) # TRUE
#all(toplot["Microbiome TN_breeder",] == hist_TN_breeder$counts/sum(all_VCs_16S$study1 == 'TN_breeder')) # TRUE
#all(toplot["Behavioural phenotypes",] == hist_behavioural_phenos$counts/sum(all_VCs_phenos$behavioural_pheno)) # TRUE
#all(toplot["Non behavioural phenotypes",] == hist_non_behavioural_phenos$counts/sum(!all_VCs_phenos$behavioural_pheno)) # TRUE

## #gather values for plotting
## range(c(all_VCs_phenos$prop_Ad1, all_VCs_16S$prop_Ad1))
## hist_MI = hist(all_VCs_16S[all_VCs_16S$study1 == 'MI','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
## hist_NY = hist(all_VCs_16S[all_VCs_16S$study1 == 'NY','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
## hist_TN_behavior = hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_behavior','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
## hist_TN_breeder = hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_breeder','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
## hist_behavioural_phenos = hist(all_VCs_phenos[all_VCs_phenos$behavioural_pheno,'prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
## hist_non_behavioural_phenos = hist(all_VCs_phenos[!all_VCs_phenos$behavioural_pheno,'prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
## 
## #plot
## #pdf('~/P50_HSrats/plots/barplots_herits_studies_pheno_new.pdf', width = 10)
## cols = gray.colors(n=4, start = 0.9, end = 0.3) # AB
## barplot(rbind(hist_TN_breeder$counts/sum(all_VCs_16S$study1 == 'TN_breeder'), 
##               hist_TN_behavior$counts/sum(all_VCs_16S$study1 == 'TN_behavior'), 
##               hist_MI$counts/sum(all_VCs_16S$study1 == 'MI'),
##               hist_NY$counts/sum(all_VCs_16S$study1 == 'NY'), 
##               hist_behavioural_phenos$counts/sum(all_VCs_phenos$behavioural_pheno), 
##               hist_non_behavioural_phenos$counts/sum(!all_VCs_phenos$behavioural_pheno)), 
##         beside = T, names = seq(0.05,0.65,by = 0.05), col = c(cols, "orange","red"), las = 1, 
##         ylim = c(0,0.8), 
##         ylab = 'Proportion of microbiome traits in the given cohort', xlab = 'Heritability')
## legend (x = 'topright', 
##         legend = c('Microbiome TN_breeders','Microbiome TN_behavior','Microbiome MI','Microbiome NY','Behavioural phenotypes','Non behavioural phenotypes'), 
##         fill = c(cols, "orange","red"), border = NA)
## #dev.off()
