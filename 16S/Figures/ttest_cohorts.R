
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')
abundances = t(collapsed_clr_counts)
#taxa in cols

#need to choose tax level (eg. genus) otherwise not correct to get relative abundances (as same reads count towards genus, order...)
tax_level = 'f__'
abundances = abundances[,grep(tax_level, colnames(abundances))]

load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo.RData')
my_strsplit = function(mot, code) {
    splot = strsplit(mot, '.', fixed = T)[[1]]
    if (length(splot) == 3) return(splot[code])
    if (grepl('blank',mot, ignore.case =T)) {
        if (code %in% c(1,2)) return(splot[code]) else return('BLANK')
    } 
    if (length(splot) == 4) return(splot[code])
    if (length(splot) == 5) {
        if (code %in% c(1,2)) return(splot[code]) else return(splot[code+2])
    }
    stop('pb')
}
rats = unname(unlist(lapply(rownames(abundances),my_strsplit,code =3)))
motch = match(rats, metadata$clean_ids )
metadata = metadata[motch,]
abundances = abundances[!is.na(metadata$study),]
metadata = metadata[!is.na(metadata$study),]


studies = c('MI','TN_breeder','NY','TN_behavior')

high_families = c("f__Bacteroidaceae","f__Muribaculaceae","f__Oscillospiraceae_88309","f__Lachnospiraceae")

sigs = c()
for(fam in colnames(abundances)) {
    ttest = anova(lm(abundances[,fam] ~ metadata$study))
    sigs = c(sigs, ttest[,'Pr(>F)'][1] < (0.05/dim(abundances)[2]))
}
names(sigs) = colnames(abundances)

table(sigs)
#130 TRUE
#264 FALSE

sigs[high_families]
#all TRUE






