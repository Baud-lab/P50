library(rhdf5)
library(parallel) # required for mclapply

# Loading 'unpruned_bug_QTLs'
load('/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/QTLs_alpha1e-04_unpruned.RData')
unpruned_bug_QTLs = unpruned_bug_QTLs[unpruned_bug_QTLs$tax_level != 'community_trait',]
DGE_QTLs = unpruned_bug_QTLs
DGE_QTLs = DGE_QTLs[DGE_QTLs$logP > 5.8,]

# - Amelie comments -
#leave in as there are still some _all
#       VCs$taxon1 = sub('_all','',VCs$taxon1)

#taxon = 'ASV_33692'
#alpha = 0.01
#load(paste(pvalues_dir_DGE,'/QTLs_alpha', alpha, '_unpruned_porcupine_', taxon,'.RData',sep=''))
#files2plot = files_DGE[grep(taxon, files_DGE)]
#pheno_names = sub('.h5','',files2plot)

#alternative1
DGE_QTLs = DGE_QTLs[grep('g__',DGE_QTLs$full_taxon),] #will not show sig QTLs that are only at the family level for example
pheno_names = unique(DGE_QTLs$measure) # has study extensions, includes both taxa and ASVs (but not community levels phenotypes)
#### bug?? see chr10 
my_extract_genus = function(full_taxon) {
  splot = strsplit(full_taxon,';')[[1]]
  splot = splot[grep('g_',splot)]
  genus = strsplit(splot, '__')[[1]][2]
}
DGE_QTLs$genus = sapply(DGE_QTLs$full_taxon, my_extract_genus) #no NA due to prior exclusion of phenotypes without genus taxonomy
taxa = DGE_QTLs[match(pheno_names, DGE_QTLs$measure),'genus'] #genera corresponding to unique measures (after excluding those mapping to higher level than taxon)
# - end Amelie comments -

# Setting mock colours to set them later - col1:coln
n <- length(unique(taxa))
colours = paste0("col", 1:n)
names(colours) = unique(taxa)


# Loading cumpos
load('/users/abaud/abaud/P50_HSrats/data/cumpos_P50_rats_Rn7.RData')

my_f = function(k) {
  # - Amelie comments -
  measure = pheno_names[k]
  taxon = taxa[k]
  if (grepl('ASV', measure)) pvalues_dir_DGE = '/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/' else
    pvalues_dir_DGE = '/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
  
  cat("starting with measure", measure, "\n")
  DGE_h5 = h5read(paste(pvalues_dir_DGE,'/',measure,'.h5',sep=''),'/')
  
  my_pvalues = c()
  my_cumposs = c()
  my_chrs = c()
  my_poss = c()
  my_betas = c()
  for (chr in 1:20) {
    if (!paste('pvalues_chr',chr,sep='') %in% names(DGE_h5)) stop()
    my_pvalues = c(my_pvalues, DGE_h5[[paste('pvalues_chr',chr,sep='')]])
    my_chrs = c(my_chrs, DGE_h5[[paste('chrs_chr',chr,sep='')]])
    my_poss = c(my_poss, DGE_h5[[paste('poss_chr',chr,sep='')]])
    my_betas = c(my_betas, DGE_h5[[paste('betas_chr',chr,sep='')]])
    my_cumposs = c(my_cumposs, cumpos[[paste('cumpos_chr',chr,sep='')]])
    if (length(my_pvalues) != length(my_cumposs)) stop('pb')
  }
  DGE_h5$chr = my_chrs
  DGE_h5$pvalues = my_pvalues
  DGE_h5$cumpos = my_cumposs
  DGE_h5$pos = my_poss
  DGE_h5$betas = my_betas
  DGE_h5$trait1 = measure
  
  DGE_h5 = data.frame(DGE_h5[['chr']],DGE_h5[['pos']],DGE_h5[['cumpos']],DGE_h5[['pvalues']],DGE_h5[['betas']],DGE_h5[['trait1']])
  colnames(DGE_h5) = c('chr','pos','cumpos','pvalues','betas', 'trait1')
  DGE_h5$logP=-log10(DGE_h5[,'pvalues'])
  
  DGE_h5$col = 'darkgrey'
  w = which(DGE_QTLs[,'measure'] == measure)
  for (i in w) {
    quels = which(DGE_h5[,'chr'] == DGE_QTLs[i,'chr'] & DGE_h5[,'pos'] == DGE_QTLs[i,'pos'])
    #		if (length(quels) != 1) print(k)
    DGE_h5[quels,'col'] = colours[taxon] 
  }
  #	for (chr in 1:20) {
  #		this_min_logP = min(DGE_h5[DGE_h5$col != 'darkgrey','logP'])
  #		if (any)
  #	}
  
  #	if (any(DGE_h5[,'logP']> 7 & DGE_h5[,'chr'] == 1) ) print(k)
  ret = DGE_h5[DGE_h5[,'logP']>=3,]
  # - end Amelie comments -
  
  cat("done with measure ", measure, "\n")
  return(ret)
}
res = mclapply(1:length(pheno_names), my_f, mc.cores = 14)

res = do.call('rbind',res)
#unique(res$col)

outfile = "/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/QTLs_alpha1e-04_unpruned_DGE_CE_MaE_toPlot.RData"
cat("saving file to", outfile, "\n")
save(res, file = outfile)
