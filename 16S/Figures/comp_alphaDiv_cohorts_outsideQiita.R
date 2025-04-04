
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')
#use collapsed_full_biomt; not CLR transformed; collapsed at taxonomic levels
abundances = t(collapsed_full_biomt)
#taxa in cols

#need to choose tax level (eg. genus) otherwise not correct to get relative abundances (as same reads count towards genus, order...)
tax_level = 'g__'
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
metadata$new_study = metadata$study
metadata[which(metadata$study == 'TN_behavior'),'new_study'] = 'TN1'
metadata[which(metadata$study == 'TN_breeder'),'new_study'] = 'TN2'

my_count = function(row) {
    sum(row != 0)
}
nb_obs_taxa = apply(abundances, MAR = 1, FUN = my_count)

metadata$new_study = factor(metadata$new_study, levels = c('NY','MI','TN1','TN2'))
anova = anova(lm(nb_obs_taxa ~ metadata$study))

pdf('/users/abaud/abaud/P50_HSrats/plots/alpha_cohorts_g.pdf')
boxplot(nb_obs_taxa ~ metadata$new_study, las = 1, xlab = 'Cohort', ylab = 'Number of different bacterial genera observed', cex.lab = 1.5, cex.axis = 1.5)
dev.off()
