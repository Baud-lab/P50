
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')
#use collapsed_full_biomt; not CLR transformed; collapsed at taxonomic levels
abundances = t(collapsed_full_biomt)
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

my_relativise = function(row) {
	row = row / sum(row)
	return(row)
}
rel_abundances = t(apply(abundances, MAR = 1, FUN = my_relativise)) #per sample

studies = c('MI','TN_breeder','NY','TN_behavior')
means = matrix(NA, ncol = length(studies), nrow = dim(rel_abundances)[2])
colnames(means) = studies
rownames(means) = colnames(rel_abundances)
for (study in studies) {
	w = which(metadata$study == study)
	means[,study] = apply(rel_abundances[w,], FUN = mean, MAR = 2) #### means is actually sometimes medians now
}

threshold = 0.1

#for condition 1
my_any = function(row) {
	any(row >= threshold)
}
any_above_threshold = apply(means, FUN = my_any, MAR = 1)

#for condition 2
#overall_means = apply(means, MAR = 1, FUN = mean)

condition = any_above_threshold
#condition = (overall_means >= threshold)

sums = apply(means[!condition,], MAR = 2, FUN = sum)

means = means[condition,]
overall_means = apply(means, MAR = 1, FUN = mean)
means = means[order(overall_means),]
means = rbind(means, sums)
rownames(means)[dim(means)[1]] = 'Sum of lower abundance families'

library(RColorBrewer)
n <- dim(means)[1]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
if (n > length(col_vector)) rep = TRUE else rep = FALSE
cols = sample(col_vector, n, replace = rep)

pdf('/users/abaud/abaud/P50_HSrats/plots/average_genera_barplots__f.pdf', width = 8)
barplot(means, col= cols, las = 1, border = NA, ylab = 'Mean relative abundance across all samples in the cohort')
plot(1,1, col = 'white')
legend(x = 'topleft', legend = rownames(means)[seq(length(cols),1, by = -1)], fill = cols[seq(length(cols),1, by = -1)], border = 'white', bty = 'n', cex = 0.6)
dev.off()

