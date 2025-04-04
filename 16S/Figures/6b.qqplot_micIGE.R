library(gap) # qqunif function - to plot
# Loading VC data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_IGE_VC.RData'))

# Selecting only the ones corresponding to all cohorts together
all_VCs_full = all_VCs_full[all_VCs_full$study1 == 'all',]
all_VCs_full = all_VCs_full[order(all_VCs_full$pvalue_DGE, decreasing = T),] #leave to TRUE: it's for QQ plot colours

# Setting colours depending on the p-value
cols = rep('black', dim(all_VCs_full)[1])
cols[all_VCs_full$cw_qvalue_DGE < 0.1] = 'orange' # FDR < 10%
cols[all_VCs_full$cw_bonf_pvalue_DGE < 0.05] = 'red' # bonferroni p < 0.05
# Setting pch depending on p-value
pchs = rep(1, dim(all_VCs_full)[1])
pchs[all_VCs_full$cw_qvalue_DGE < 0.1] = 16 # FDR < 10%


# Open pdf to save plot
pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/QQplot_pvalues_IGE_Helenes.pdf', h=6, w=7)
par(mar=c(5.1,5.1,2.1,2.1))

# Define y limit max
mox = max(-log10(all_VCs_full$pvalue_DGE), na.rm = T)
# Plot
qqunif(all_VCs_full[,'pvalue_DGE'],ci=T,
       las=1, main ='All cohorts mega-analysis: 3,767 rats', 
       ylim = c(0,mox), col = cols, pch = pchs, 
       cex.axis = 1.25, cex.lab=1.4, cex = 1.2)
# legend
legend("topleft", pch=16, col=c("red","orange"), 
       legend=c(expression(Bonferroni~italic(p)~"<"~0.05), "FDR < 10%"), 
       cex = 1.2, title = "Significance of Mic-IGE:", title.font = 2)
dev.off()


###
## #doing plots per cohort 
## par(mfrow = c(2,2), mar = c(2,2,3,1))
## qqunif(all_VCs_full[grepl('_all',all_VCs_full$trait1),'pvalue_DGE'],ci=T,col = 'black',las=1, main ='All', ylim = c(0,mox))
##  clear inflation there
## g_NY = grepl('_NY',all_VCs_full$trait1)
## qqunif(all_VCs_full[g_NY,'pvalue_DGE'],ci=T,las=1, main ='NY behavior: 1,167 rats', ylim = c(0,mox), col = cols[g_NY], pch = pchs[g_NY])
## g_MI = grepl('_MI',all_VCs_full$trait1)
## qqunif(all_VCs_full[g_MI,'pvalue_DGE'],ci=T,las=1, main ='MI behavior: 1,112 rats', ylim = c(0,mox), col = cols[g_MI], pch = pchs[g_MI])
## g_TN_behavior = grepl('_TN_behavior',all_VCs_full$trait1)
## qqunif(all_VCs_full[g_TN_behavior,'pvalue_DGE'],ci=T,las=1, main ='TN behavior: 950', ylim = c(0,mox), col = cols[g_TN_behavior], pch = pchs[g_TN_behavior])
## g_TN_breeder = grepl('_TN_breeder',all_VCs_full$trait1)
## qqunif(all_VCs_full[g_TN_breeder,'pvalue_DGE'],ci=T,las=1, main ='TN breeders: 555 rats', ylim = c(0,mox), col = cols[g_TN_breeder], pch = pchs[g_TN_breeder])
