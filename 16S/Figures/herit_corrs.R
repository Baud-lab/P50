#load herit data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))

#build center_spe_herits table with one row per microbiome phenotype and 4 columns corresponding to 4 cohorts
union = unique(all_VCs_full$taxon)
center_spe_herits = matrix(nrow = length(union), ncol = 4, NA)
rownames(center_spe_herits) = union
colnames(center_spe_herits) = c('NY', 'MI', 'TN_behavior','TN_breeder')

# this could be done with apply
for (i in 1:dim(all_VCs_full)[1]) {
    center_spe_herits[all_VCs_full[i,'taxon1'],all_VCs_full[i,'study1']] = all_VCs_full[i,'prop_Ad1']
} 

# ordering to align as the rest of the plots in figure
center_spe_herits = center_spe_herits[,c('TN_breeder', 'TN_behavior',  'MI', 'NY')]

#NY cohort  (N = 1,167 rats), 
#MI cohort (N = 1,112 rats), 
#TN behaviour cohort (N = 950), 
#TN breeder cohort (N = 555).
Ns = c("NY\n(N = 1,167)", "MI\n(N = 1,112)", "TN_behavior\n(N = 950)", "TN_breeder\n(N = 555)")

#plot using R's pairs plot
range(all_VCs_full$prop_Ad1)
#pdf('~/P50_HSrats/plots/compare_herits_Helenes_diff_centers.pdf')
pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/compare_herits_Helenes_diff_centers.pdf', w=6,h=6)
#names(par())
par(mar=c(5.1,5.1,2.1,2.1))
# Colnames for plotting
#colnames(center_spe_herits) = c('NY\n(N=1,167)', 'MI\n(N=1,112)', 'TN_behavior\n(N=950)','TN_breeder\n(N=555)')
colnames(center_spe_herits) = Ns[match(colnames(center_spe_herits), unlist(lapply(strsplit(Ns, "\n"),"[[", 1)))]

cols = rep('black')
#colr = "#00A087FF" # "red"
colr="#E64B35FF"
my_cor <- function(x, y,...) {
    cor = cor.test(x, y, use = 'pairwise.complete.obs')
    txt <- paste('cor = ',format(cor$estimate, digits = 2)[1],sep='')
    if (cor[['p.value']] < (0.05/6)) col = colr else col = 'black' # color = colr and font = bold if passes Bonferroni correction
    text(0.1, 0.1, txt, cex = 1.4, col = col)
}
my_points = function(x,y,...) {
    points(x,y, ylim = c(0,0.2), las = 1, cex.lab = 2, 
           pch = 16) #, #19 
           #cex=0.8)
}
my_main = function(x,y,...) {
    points(x,y)
}

#my_diag = function(x,y,...) {
#    txt = colnames(center_spe_herits)
#    text(0.14,0.14, txt, cex = 3, col = 'blue', font=2)
#}

pairs(center_spe_herits, panel = points, xlim = c(0,0.2), ylim = c(0,0.2), #cex.main=2,
      las = 1, upper.panel = my_points, lower.panel = my_cor, xaxt='n', yaxt='n', font.labels=2, 
      gap=0.5) # gap between panels

dev.off()


