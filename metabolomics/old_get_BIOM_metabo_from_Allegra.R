#quantitative measures
abundances = read.csv('/homes/abaud/P50_HSrats/data/metabo/adducts/Amelie_BlankFiltered_061419.csv', as.is = T, check.names = F)
dim(abundances)
#[1] 2119 1406
rownames(abundances) = paste('MZ',round(abundances[,'row m/z'],digits = 4),';',round(abundances[,'row retention time'],digits = 4),';',abundances[,'row ID'],sep='')
g = regexpr('0007',colnames(abundances))
unique(g)
#-1 1
#meaning that either not present or present at beginning of colname
all(grepl('Blank',colnames(abundances)[g==(-1)]) | grepl('Blk',colnames(abundances)[g==(-1)]) | grepl('Std_MIx',colnames(abundances)[g==(-1)]) | colnames(abundances)[g==(-1)] %in% c('row ID','row m/z','row retention time'))
# TRUE
abundances = abundances[,g==1]
colnames(abundances) = sub('.mzXML filtered Peak area','',colnames(abundances),fixed = T)
dim(abundances)
#[1]  2119 1152
load('/homes/abaud/P50_HSrats/data/metabo/unnormalized_metabo.RData')
biomt = biomt[,!grepl('Blank',colnames(biomt))]
dim(biomt)
#[1] 7412 1150
inter = intersect(colnames(abundances), colnames(biomt))
length(inter)
#1150
colnames(abundances)[!colnames(abundances) %in% inter]
#[1] "00077E756F_180412183056" "00077E8185_180413221941"
match(c('00077E756F','00077E756F_180412183056'),colnames(abundances))
#[1] 15 16
match(c('00077E8185','00077E8185_180413221941'),colnames(abundances))
#[1] 284 285
# what are 00077E8185_180413221941 and 00077E756F_180412183056? seem like replicates
colnames(biomt)[!colnames(biomt) %in% colnames(abundances)]
#none
abundances = as.matrix(abundances)
#plot(abundances[,'00077E756F'],abundances[,'00077E756F_180412183056'])
#plot(abundances[,'00077E8185'],abundances[,'00077E8185_180413221941'])
abundances = abundances[,-match(c('00077E756F_180412183056','00077E8185_180413221941'),colnames(abundances))]
save(abundances, file = '/homes/abaud/P50_HSrats/data/metabo/Allegra_metabo.RData')

####### compare with previous BIOM

load('/homes/abaud/P50_HSrats/data/metabo/unnormalized_metabo.RData')

my_parse = function(feature_name) {
	splot = strsplit(feature_name,';')[[1]]
	mz = round(as.numeric(sub('MZ','',splot[1])),digits = 2)
	rt = round(as.numeric(splot[2]),digits = 2)
	row_ID = as.numeric(splot[3])
	return(c(mz,rt,row_ID))
}
biomt_features = t(sapply(rownames(biomt),my_parse))
abundances_features = t(sapply(rownames(abundances),my_parse))
motch = list()
for (i in 1:dim(abundances_features)[1]) {
	w = which(biomt_features[,1] == abundances_features[i,1] & biomt_features[,2] == abundances_features[i,2])
	motch[[i]] = unname(w)
}
table(unlist(lapply(motch,length)))
#   0    1    2 
#1606  618    9 
remove = which(unlist(lapply(motch,length)) != 1)
abundances = abundances[-remove,colnames(abundances) %in% inter]
motch = motch[-remove]
biomt = biomt[unlist(motch),match(colnames(abundances),colnames(biomt))]
s = sample(1:dim(abundances)[1],32)
pdf('/homes/abaud/P50_HSrats/plots/compare_old_new_feature_tables.pdf',width = 15, height = 15)
par(mfrow = c(4,4))
for (i in s) {
	plot(biomt[i,], abundances[i,], xlab = 'old, unnormalised feature table from Qiita', ylab = 'new, (un)normalised? feature table (061419.csv)', main = rownames(biomt)[i])
	abline(0,1,col = 'grey')
}
dev.off()

###### collapse adducts
load('/homes/abaud/P50_HSrats/data/metabo/Allegra_metabo.RData')
my_parse = function(feature_name) {
	splot = strsplit(feature_name,';')[[1]]
	mz = round(as.numeric(sub('MZ','',splot[1])),digits = 2)
	rt = round(as.numeric(splot[2]),digits = 2)
	row_ID = as.numeric(splot[3])
	return(c(mz,rt,row_ID))
}
infos_abundances = t(sapply(rownames(abundances),my_parse))
colnames(infos_abundances) = c('mz','rt','row_ID')
infos_abundances = as.data.frame(infos_abundances, stringsAsFactors = F)

edges = read.csv('/homes/abaud/P50_HSrats/data/metabo/adducts/ms1_corr_edges_msannotation.csv', as.is = T, check.names = F)

in_adducts = unique(c(edges$ID1,edges$ID2))
new_abundances = NULL
infos_abundances$map = NA
count = 0
done = c()
for (i in in_adducts) {
	w = which(edges$ID1 == i)
	if (length(w) == 0) next
 	collapse = c(i,edges[w,'ID2'])
	#if (length(collapse) != length(unique(collapse))) stop()
	#collapse will be value in infos_abundances$row_ID that are collapsed
	collapse = unique(collapse)
	if (any(collapse %in% done)) {if (!all(collapse %in% done)) stop('pb code') ; next}
	done = unique(c(done,collapse))
	count = count + 1
	#infos_abundances$map has row in new_abundances
	infos_abundances[infos_abundances$row_ID %in% collapse,'map'] = count
	new_abundance = apply(abundances[infos_abundances$row_ID %in% collapse,],MAR = 2, FUN = sum)
	new_abundance = matrix(new_abundance, nrow = 1)
	rownames(new_abundance) = paste('collapsed_abundances_row',count,sep='')
	new_abundances = rbind(new_abundances,new_abundance)
}
k = dim(new_abundances)[1]
infos_abundances[!infos_abundances$row_ID %in% in_adducts,'map'] = ((k+1):(k+sum(!infos_abundances$row_ID %in% in_adducts)))
new_abundances = rbind(new_abundances,abundances[!infos_abundances$row_ID %in% in_adducts,])
save(infos_abundances,new_abundances, file = '/homes/abaud/P50_HSrats/data/metabo/adducts/collapsed_Allegra_metabo.RData')



