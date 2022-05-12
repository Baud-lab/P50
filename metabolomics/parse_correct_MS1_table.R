#quantitative measures - imputed by Allegra and blank filtered by Kelly
abundances = read.csv('/homes/abaud/P50_HSrats/data/metabo/Amelie_BlankFiltered_09092019_correct.csv', as.is = T, check.names = F)
dim(abundances)
#[1] 2119 1406
#rownames(abundances) = paste('MZ',round(abundances[,'row m/z'],digits = 4),';',round(abundances[,'row retention time'],digits = 4),';',abundances[,'row ID'],sep='')
rownames(abundances) = paste('MZ',round(abundances[,'row m/z'],digits = 4),';',round(abundances[,'row retention time'],digits = 4),sep='')
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
abundances = as.matrix(abundances)
save(abundances, file = '/homes/abaud/P50_HSrats/data/metabo/Amelie_BlankFiltered_09092019_correct.RData')
