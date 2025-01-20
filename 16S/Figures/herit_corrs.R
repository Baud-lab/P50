#load herit data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(paste(root_dir,'augmented_VC.RData',sep=''))

#build center_spe_herits table with one row per microbiome phenotype and 4 columns corresponding to 4 cohorts
union = unique(all_VCs_full$taxon)
center_spe_herits = matrix(nrow = length(union), ncol = 4, NA)
rownames(center_spe_herits) = union
colnames(center_spe_herits) = c('NY', 'MI', 'TN_behavior','TN_breeder')
for (i in 1:dim(all_VCs_full)[1]) {
    center_spe_herits[all_VCs_full[i,'taxon1'],all_VCs_full[i,'study1']] = all_VCs_full[i,'prop_Ad1']
} 

#plot using R's pairs plot
range(all_VCs_full$prop_Ad1)
pdf('~/P50_HSrats/plots/compare_herits_Helenes_diff_centers.pdf')
cols = rep('black')
my_cor <- function(x, y,...) {
    cor = cor.test(x, y, use = 'pairwise.complete.obs')
    txt <- paste('cor = ',format(cor$estimate, digits = 2)[1],sep='')
    if (cor[['p.value']] < (0.05/6)) col = 'red' else col = 'black' #red if passes Bonferroni correction
         text(0.1, 0.1, txt, cex = 2, col = col)
}
my_points = function(x,y,...) {
    points(x,y, ylim = c(0,0.2), pch = 19, las = 1, cex.lab = 2)
}
my_main = function(x,y,...) {
    points(x,y)
}
my_diag = function(x,y,...) {
    txt = colnames(center_spe_herits)
    text(0.14, 0.14, txt, cex = 3, col = 'blue')
}
pairs(center_spe_herits, panel = points,xlim = c(0,0.2),ylim = c(0,0.2), las = 1, upper.panel = my_points, lower.panel = my_cor, xaxt='n', yaxt='n')
dev.off()
