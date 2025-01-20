#load phenotype data
phenos_dir ='/users/abaud/abaud/P50_HSrats/output/VD/univariate/residuals_02Jan2020/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
load(paste(phenos_dir,'all_estNste.Rdata',sep=''))
all_VCs_phenos = VCs

#load microbiome data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(paste(root_dir,'augmented_VC.RData',sep=''))
all_VCs_16S = all_VCs_full

#determine if phenotype is behavioural or not
non_behav = c('adams','baculum','bicknell','insulin','physiological')
behav = c('ccc','ccp','crf','delay_discounting','elevated_plus_maze','light_reinforcement','locomotor_testing','nicsa','novel_object','novelty_seeking','open_field','pavca','reaction_time','social_')
all_VCs_phenos$behavioural_pheno = NA
for (i in 1:dim(all_VCs_phenos)[1]) {
    behavioural = FALSE
    for (cate in behav) {
        if (grepl(cate, all_VCs_phenos[i,'trait1'])) {behavioural = TRUE ; next} 
    }
    all_VCs_phenos[i,'behavioural_pheno'] = behavioural
}
#visual check
sort(all_VCs_phenos[which(all_VCs_phenos[,'behavioural_pheno']),'trait1'])
sort(all_VCs_phenos[-which(all_VCs_phenos[,'behavioural_pheno']),'trait1'])

#gather values for plotting
range(c(all_VCs_phenos$prop_Ad1, all_VCs_16S$prop_Ad1))
hist_MI = hist(all_VCs_16S[all_VCs_16S$study1 == 'MI','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
hist_NY = hist(all_VCs_16S[all_VCs_16S$study1 == 'NY','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
hist_TN_behavior = hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_behavior','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
hist_TN_breeder = hist(all_VCs_16S[all_VCs_16S$study1 == 'TN_breeder','prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
hist_behavioural_phenos = hist(all_VCs_phenos[all_VCs_phenos$behavioural_pheno,'prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)
hist_non_behavioural_phenos = hist(all_VCs_phenos[!all_VCs_phenos$behavioural_pheno,'prop_Ad1'], breaks = seq(0,0.65,by = 0.05), plot = F)

#plot
pdf('~/P50_HSrats/plots/barplots_herits_studies_pheno_new.pdf', width = 10)
cols = gray.colors(n=4, start = 0.9, end = 0.3)
barplot(rbind(hist_TN_breeder$counts/sum(all_VCs_16S$study1 == 'TN_breeder'), hist_TN_behavior$counts/sum(all_VCs_16S$study1 == 'TN_behavior'), hist_MI$counts/sum(all_VCs_16S$study1 == 'MI'),hist_NY$counts/sum(all_VCs_16S$study1 == 'NY'), hist_behavioural_phenos$counts/sum(all_VCs_phenos$behavioural_pheno), hist_non_behavioural_phenos$counts/sum(!all_VCs_phenos$behavioural_pheno)), beside = T, names = seq(0.05,0.65,by = 0.05), col = c(cols, 'orange','red'), las = 1, ylim = c(0,0.8), ylab = 'Proportion of microbiome traits in the given cohort', xlab = 'Heritability')
legend (x = 'topright', legend = c('Microbiome TN_breeders','Microbiome TN_behavior','Microbiome MI','Microbiome NY','Behavioural phenotypes','Non behavioural phenotypes'), fill = c(cols,'orange','red'), border = NA)
dev.off()




