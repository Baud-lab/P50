# Loading counts (not CLR transformed)
load('collapsed_full_biomt_collapsed_clr_counts.RData')
# use collapsed_full_biomt; not CLR transformed; collapsed at taxonomic levels
abundances = t(collapsed_full_biomt)
# taxa in cols

# Loading metadata 
load('metadata_16Spaper.RData')

# TODO: need to choose tax level (eg. genus) otherwise not correct to get relative abundances (as same reads count towards genus, order...)
tax_level = 'f__'
abundances = abundances[,grep(tax_level, colnames(abundances))]

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


# To collapse taxa in "Sum of lower abundance families"
# TODO: change the threshold to see more or less fragmentation in the bars
#       lower the threshold, higher the fragmentation of bars
threshold = 0.1 # 5 sections in the bars with our data

## for condition 1
my_any = function(row) {
	any(row >= threshold)
}
any_above_threshold = apply(means, FUN = my_any, MAR = 1)
condition = any_above_threshold

## for condition 2
#overall_means = apply(means, MAR = 1, FUN = mean)
#condition = (overall_means >= threshold)

# For both conditions now
sums = apply(means[!condition,], MAR = 2, FUN = sum)

means = means[condition,]
overall_means = apply(means, MAR = 1, FUN = mean)
means = means[order(overall_means),]
means = rbind(means, sums)
rownames(means)[dim(means)[1]] = 'Sum of lower abundance families'

# Setting things up to plot
# defining cohort names as in paper
dict = c("NY"="NY", "MI"="MI", "TN_behavior"="TN1", "TN_breeder"="TN2")
colnames(means) = unname(dict[colnames(means)])
colord = c("NY", "MI", "TN1", "TN2")
means = means[,colord]
n <- dim(means)[1]

# Setting colors specifically 
coolors = c("#4A6990FF", "#7AA6DCFF", "#00A087FF", "#91D1C2FF", "#EFC000FF") # colors in paper, exact number of colours as in `means` 
cols = coolors[1:n]

## if want to have a look at colors
#scales::show_col(cols)

# Saving objects to plot
#save(means, cols, tax_level, file=paste0("source_files/fig2a__",gsub("__","",tax_level),".RData"))

# Loading objects to plot
#load("source_files/fig2a__f.RData")

# Open pdf to save plot
pdf(paste0("average_genera_barplots__",gsub("__","",tax_level),".pdf"), h= 7, w = 10)
par(mar=c(5.1,5.1,2.1,7.1))

# barplot
barplot(means, col= cols, las = 1, border = NA, 
        ylab = '', # define y lab separately
        cex.axis=1.25)
# add y lab
title( ylab = 'Mean relative abundance across all samples in the cohort', 
       cex.lab = 1.4, line=3.5)
# legend
legend(x = 'topright', legend = rownames(means)[seq(length(cols),1, by = -1)], 
       fill = cols[seq(length(cols),1, by = -1)],
       border = NA, 
       bty = 'n', cex = 0.6, 
       xpd=T, inset=c(-0.2,0))

dev.off()

