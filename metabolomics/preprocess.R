load('/homes/abaud/P50_HSrats/data/metabo/biomt_annots_metadata.RData')

# filter out OTUs based on prop_rats_with and median_c
# will use > median_c_th
median_c_th = 0
commons = list()
pdf('/homes/abaud/P50_HSrats/plots/filter_metabo.pdf')
par(mfrow = c(2,2))
#no metabo samples from TN
for (phenotyping_center in c('NY','MI')) {
  print(phenotyping_center)
  this_counts = biomt[,metadata[,"phenotyping_center"] == phenotyping_center]
  pce_counts = (this_counts!=0)
  prop_rats_with = apply(pce_counts, FUN = mean, MAR = 1)
  median_c = apply(this_counts, FUN = median, MAR = 1)
  hist(prop_rats_with, xlab = 'Proportion of rats with the bug (pce/absence)', main = phenotyping_center)
  hist(median_c, xlab = 'Median bug abundance', main = phenotyping_center)
  plot(prop_rats_with,median_c,xlab = 'Proportion of rats with the bug (pce/absence)', ylab = 'Median bug abundance', main = phenotyping_center)
  plot(prop_rats_with,median_c,xlab = 'Proportion of rats with the bug (pce/absence)', ylab = 'Median bug abundance',ylim = c(0,10), main = phenotyping_center)
  tests = (median_c>median_c_th)
  print(sum(tests))
  commons[[phenotyping_center]] = rownames(this_counts)[tests]
  somple = sample(commons[[phenotyping_center]], size = 8)
  for (bug in somple) {
    boxplot(this_counts[bug,], outline = F, xlab = bug, main = paste('Median c is:',median_c[bug]))
  }
}
dev.off()
#[1] "NY"
#[1] 1977
#[1] "MI"
#[1] 1789
inter1 = intersect(commons[['MI']],commons[['NY']])
length(inter1)
#1763

# plot expression of these bugs in two centers
cols = rep('grey', dim(biomt)[1])
cols[which(rownames(biomt) %in% commons[['MI']])] = 'red'
cols[which(rownames(biomt) %in% commons[['NY']])] = 'blue'
cols[which(rownames(biomt) %in% inter1)] = 'purple'
pdf('/homes/abaud/P50_HSrats/plots/medians_inter_centers_metabo.pdf', width = 10)
par(mfrow = c(1,2))
median_MI = apply(biomt[,metadata[,"phenotyping_center"] == 'MI'], FUN = median, MAR = 1)
median_NY = apply(biomt[,metadata[,"phenotyping_center"] == 'NY'], FUN = median, MAR = 1)

plot(median_MI,median_NY, col = cols, xlab = 'Median in MI', ylab = 'Median in NY')
abline(0,1,col = 'grey')
# !! all but purple won't show in plot below as log10 not defined for 0
plot(log10(median_MI),log10(median_NY), col = cols, xlab = 'log10(median in MI)', ylab = 'log10(median in NY)')
abline(0,1,col = 'grey')

dev.off()

#create center specific biomt: filter out rare OTUs and few samples with all 0s for those common bugs
biomt_MI = biomt[commons[['MI']],metadata[,"phenotyping_center"] == 'MI']
exclude_MI = colnames(biomt_MI)[apply(biomt_MI, FUN = sum, MAR = 2)==0]
length(exclude_MI)
#[1] 2 rat whose bugs are not in common for that phenotyping center
if (length(exclude_MI)!=0) biomt_MI = biomt_MI[,-match(exclude_MI,colnames(biomt_MI))]

biomt_NY = biomt[commons[['NY']],metadata[,"phenotyping_center"] == 'NY']
exclude_NY = colnames(biomt_NY)[apply(biomt_NY, FUN = sum, MAR = 2)==0]
length(exclude_NY)
#[1] 2 rat whose bugs are not in common for that phenotyping center
if (length(exclude_NY)!=0) biomt_NY = biomt_NY[,-match(exclude_NY,colnames(biomt_NY))]

save(biomt_MI,biomt_NY, metadata, file = '/homes/abaud/P50_HSrats/data/metabo/center_spe_biomt.RData')
