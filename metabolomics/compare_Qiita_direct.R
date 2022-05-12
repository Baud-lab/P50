load('/homes/abaud/P50_HSrats/data/metabo/78940_palmer_amelie_metabolomics_original.RData')
#biomt
load('/homes/abaud/P50_HSrats/data/metabo/Amelie_BlankFiltered_09092019_correct.RData')
#abundances

colnames_biomt = colnames(biomt)
corresp = cbind(colnames(abundances),NA)
rownames(corresp) = colnames(abundances)
for (colname_abundances in colnames(abundances)) {
	g = grep(colname_abundances, colnames_biomt)
	if (length(g)!=1) print(colname_abundances) else corresp[colname_abundances,2] = colnames_biomt[g]
}

biomt = biomt[,gs_in_biomt]

cors = c()
for (i in 1:dim(biomt)[2]) {
	cors = c(cors, cor(biomt[,i], abundances[,i]))
}

#which(cors<0.1)
#[1]   15  310  312 1079
plot(biomt[,15], abundances[,15])
abline(0,1,col = 'red')
#lots of 0s, but that are also found in original MS1 table from Allegra: Allegras_061419.csv
#00077E756F.mzXML filtered Peak area

par(mfrow = c(1,2))
plot(abundances[,'00077E756F_180412183056'],biomt[,"11479.RATCECUM.MI.00077E756F"])
abline(0,1,col = 'red')
plot(abundances[,'00077E756F'],biomt[,"11479.RATCECUM.MI.00077E756F"])
abline(0,1,col = 'red')

par(mfrow = c(1,2))
plot(abundances[,'00077E8185_180413221941'],biomt[,"11479.RATCECUM.NY.00077E8185"])
abline(0,1,col = 'red')
plot(abundances[,'00077E8185'],biomt[,"11479.RATCECUM.NY.00077E8185"])
abline(0,1,col = 'red')

#ok so one of two duplicates kept in BIOM