load('P50_Rn7_pruned_DGE_IGE_IEE_cageEffect_maternalEffect_all_estNste.Rdata')
all_VCs_full = VCs
load('P50_Rn7_pruned_DGE_IEE_cageEffect_maternalEffect_all_estNste.Rdata')
all_VCs_null = VCs

inter=intersect(all_VCs_null$trait1,all_VCs_full$trait1)
length(inter)
all_VCs_full=all_VCs_full[match(inter,all_VCs_full$trait1),]
all_VCs_null=all_VCs_null[match(inter,all_VCs_null$trait1),]

library(gap)

pdf('QQplot_IGE_pvalues_bootstrap_dfs.pdf',bg='white', width = 20, height = 15)
par(mfrow = c(3,4))

	#https://stat.ethz.ch/pipermail/r-help/2013-June/354527.html
for (k in seq(0,1,by = 0.1)) {
	all_VCs_full$pvalue_IGE = k*pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=1, lower.tail = FALSE) + (1-k)*pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=2, lower.tail = FALSE)
	qqunif(all_VCs_full[,'pvalue_IGE'],ci=T,las=1, main ='P values null simulations', sub = paste('Mixture parameter:',k), pch = 16)
}
dev.off()
