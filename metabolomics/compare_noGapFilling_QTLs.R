my_strsplit = function(mot) {
	splot = strsplit(mot, ';', fixed = T)[[1]]
	mz = as.numeric(sub('MZ','',splot[1]))
	rt = as.numeric(strsplit(splot[2], '_', fixed = T)[[1]][1])
	return(c(mz,rt))
}


load('/nfs/leia/research/stegle/abaud/P50_HSrats/output/old/pvalues/metabo7/round8_DGE_cageEffect/QTLs_alpha0.001.RData')
old_QTLs = all_QTLs
old_QTLs = cbind(old_QTLs, do.call('rbind', lapply(old_QTLs$measure, FUN= my_strsplit)))
colnames(old_QTLs) = c(colnames(all_QTLs), c('mz','rt'))

load('/nfs/leia/research/stegle/abaud/P50_HSrats/output/pvalues_LOCO/metabo_counts_new3/pruned_dosagesDGE/QTLs_alpha1e-04.RData')
new_QTLs = all_QTLs
new_QTLs = cbind(new_QTLs, do.call('rbind', lapply(new_QTLs$measure, FUN= my_strsplit)))
colnames(new_QTLs) = c(colnames(all_QTLs), c('mz','rt'))

#head(old_QTLs[,c('measure','chr','pos','logP')])
#head(new_QTLs[,c('measure','chr','pos','logP')])

par(mfrow = c(1,2))
hist(old_QTLs[old_QTLs$logP > 5, 'logP'], breaks = 10, ylim = c(0, 100))
hist(new_QTLs[new_QTLs$logP > 5, 'logP'], breaks = 10, ylim = c(0, 100))

logP_th = 6
sum(old_QTLs$logP > logP_th & grepl('_all',old_QTLs$measure))
sum(new_QTLs$logP > logP_th & grepl('_all',new_QTLs$measure))

sum(old_QTLs$logP > logP_th & grepl('_MI',old_QTLs$measure))
sum(new_QTLs$logP > logP_th & grepl('_MI',new_QTLs$measure))

sum(old_QTLs$logP > logP_th & grepl('_NY',old_QTLs$measure))
sum(new_QTLs$logP > logP_th & grepl('_NY',new_QTLs$measure))


top_old_QTLs = old_QTLs[old_QTLs$logP > 10,]
top_new_QTLs = new_QTLs[new_QTLs$logP > 7,]
top_new_QTLs[,c('measure','chr','pos','logP')]

for (i in 1:dim(top_old_QTLs)[1]) {
	w = which(new_QTLs$chr == top_old_QTLs[i,'chr'] & abs(new_QTLs$pos - top_old_QTLs[i,'pos']) < 500000 & new_QTLs$logP > 5 & abs(new_QTLs$mz - top_old_QTLs[i,'mz']) < 1 & abs(new_QTLs$rt - top_old_QTLs[i,'rt']) <= 0.1)
#	if (length(w) == 0 | !grepl(';5.7', top_old_QTLs[i,'measure'])) next
	if (length(w) == 0 | !grepl(';5.7', top_old_QTLs[i,'measure'])) next
	print(top_old_QTLs[i,c('measure','chr','pos','logP')])
	print(new_QTLs[w,c('measure','chr','pos','logP')])
}

for (i in 1:dim(top_new_QTLs)[1]) {
	w = which(old_QTLs$chr == top_new_QTLs[i,'chr'] & abs(old_QTLs$pos - top_new_QTLs[i,'pos']) < 500000 & old_QTLs$logP > 5 & abs(old_QTLs$mz - top_new_QTLs[i,'mz']) < 1 & abs(old_QTLs$rt - top_new_QTLs[i,'rt']) <= 0.1)
#	if (length(w) == 0 | !grepl(';5.7', top_old_QTLs[i,'measure'])) next
	if (length(w) == 0) next
	print(top_new_QTLs[i,c('measure','chr','pos','logP')])
	print(old_QTLs[w,c('measure','chr','pos','logP')])
}

old_QTLs[which(abs(old_QTLs$mz - 433) < 2 & abs(old_QTLs$rt - 6.8) <= 0.2),c('measure','chr','pos','logP')]
new_QTLs[which(abs(new_QTLs$mz - 433) < 2 & abs(new_QTLs$rt - 6.8) <= 0.2),c('measure','chr','pos','logP')]

all_QTLs[all_QTLs$chr == 14 & abs(all_QTLs$pos - 92672944) < 1000000,]
check = all_QTLs[all_QTLs$chr == 16 & abs(all_QTLs$pos - 90610504 ) < 1000000 & grepl(';5.7',all_QTLs$measure),]
check = check[order(check$measure),]

load('/homes/abaud/P50_HSrats/data/metabo/annots_115698_palmer_biom_march_2021_original.RData')
