abundances = read.csv('/homes/abaud/P50_HSrats/data/metabo/adducts/061419.csv', as.is = T, check.names = F)
all(is.na(abundances[,dim(abundances)[2]]))
#TRUE
abundances = abundances[,-dim(abundances)[2]]
edges = read.csv('/homes/abaud/P50_HSrats/data/metabo/adducts/ms1_corr_edges_msannotation.csv', as.is = T, check.names = F)
unique(edges$EdgeType)
#"MS1 annotation"
table(edges$Score)
# 3  4  5  7  9 13 
#39 43 10 21 36 78 
my_parse = function(annot) {
	outs = strsplit(annot,' ')[[1]]
	from = outs[1]
	to = outs[2]
	dmz = outs[3]
	dmz=sub('dm/z=','',dmz)
	return(c(from,to,dmz))
}
add = t(sapply(edges$Annotation, my_parse))
add = as.data.frame(add,stringsAsFactors = F)
row.names(add) = NULL
edges = cbind(edges,add)
colnames(edges)[(dim(edges)[2]-2):dim(edges)[2]] = c('from_mol','to_mol','dmz')
toble_from = table(edges[,'from_mol'])
toble_to = table(edges[,'to_mol'])
uni = unique(c(names(toble_from),names(toble_to)))
uni_toble_from = toble_from[match(uni,names(toble_from))]
uni_toble_from[is.na(uni_toble_from)] = 0
names(uni_toble_from) = uni
uni_toble_to = toble_to[match(uni,names(toble_to))]
uni_toble_to[is.na(uni_toble_to)] = 0
names(uni_toble_to) = uni

#plot(as.vector(uni_toble_from),as.vector(uni_toble_to))
#abline(0,1,col = 'red')
#somewhat correlated

in_adducts = c(edges$ID1,edges$ID2)
table(table(in_adducts))
# 2  3  4  6  8 12 
#39 26  7  7  9 13 
# what are those?

edges = cbind(edges, abundances[match(edges$ID1,abundances[,'row ID']),c('row m/z','row retention time')],abundances[match(edges$ID2,abundances[,'row ID']),c('row m/z','row retention time')])
colnames(edges)[(dim(edges)[2]-3):dim(edges)[2]] = c('MZ1','RT1','MZ2','RT2')
output = edges[,c('Score','Annotation','MZ1','RT1','MZ2','RT2')]
#write.table(output, file = '/homes/abaud/P50_HSrats/plots/metabo_edges.txt',sep='\t',col.names = T, row.names = F, quote = F)

in_adducts = unique(in_adducts)
new_abundances = NULL
all_infos = NULL
all(edges$ID2 > edges$ID1)
#TRUE
for (i in unique(edges$ID1)) {
	w = which(edges$ID1 == i)
	collapse = c(i,edges[w,'ID2'])
	new_abundance = apply(abundances[abundances[,'row ID'] %in% collapse,4:dim(abundances)[2]],MAR = 2, FUN = sum)
	infos = c(paste(abundances[abundances[,'row ID'] %in% collapse,1],collapse = ','),paste(abundances[abundances[,'row ID'] %in% collapse,2],collapse = ','),paste(abundances[abundances[,'row ID'] %in% collapse,3],collapse = ','))	
	new_abundances = rbind(new_abundances,new_abundance)
	all_infos = rbind(all_infos,infos)
}
colnames(all_infos) = colnames(abundances)[1:3]
add = data.frame(all_infos,new_abundances,stringsAsFactors = F, row.names = NULL, check.names = F)
new_abundances = rbind(abundances[!abundances[,'row ID'] %in% in_adducts,],add)
save(new_abundances, file = '/homes/abaud/P50_HSrats/data/metabo/adducts/collapsed_abundances_061419.RData')


