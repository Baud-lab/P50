#filter collapsed taxa level data 
#similar script for uncollapsed/asv level data 

#for the abundance phenotypes reported in the paper, only taxa present in at least 50% of the rats in the considered cohort will be kept
#during development, presence/absence phenotype also considered and in that case only taxa presence in more than 20% and less than 80% of the rats kept

load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')

load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
motch = match(colnames(collapsed_clr_counts),metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

all(colnames(collapsed_clr_counts) == colnames(collapsed_full_biomt))
#TRUE

my_all_is_NA = function(x) {
  return(all(is.na(x)))
}

non_varying = function(x) {
  return(length(unique(x)) == 1)
}

calculate_max_prop_ties = function(x) {
  return(max(table(x)) / length(x))
}

# filter out OTUs FOR COUNTS PHENOTYPE
commons = list()
pdf('/users/abaud/abaud/P50_HSrats/plots/filter_bugs_new_deblur.pdf')
#par(mfrow = c(2,2))
for (study in c('all','MI','NY','TN_behavior','TN_breeder')) {
  print(study)
  if (study=='all') {
    this_collapsed_clr_counts = collapsed_clr_counts
    this_collapsed_full_biomt = collapsed_full_biomt
  } else if (study == 'TN_behavior') {
    this_collapsed_clr_counts = collapsed_clr_counts[,which(metadata[,"phenotyping_center"] == 'TN' & metadata[,'TN_track'] == 'behavior')]  
    this_collapsed_full_biomt = collapsed_full_biomt[,which(metadata[,"phenotyping_center"] == 'TN' & metadata[,'TN_track'] == 'behavior')]  
  } else if (study == 'TN_breeder') {
    this_collapsed_clr_counts = collapsed_clr_counts[,which(metadata[,"phenotyping_center"] == 'TN' & metadata[,'TN_track'] == 'breeder')]       
    this_collapsed_full_biomt = collapsed_full_biomt[,which(metadata[,"phenotyping_center"] == 'TN' & metadata[,'TN_track'] == 'breeder')]  
  } else {
    this_collapsed_clr_counts = collapsed_clr_counts[,which(metadata[,"phenotyping_center"] == study)]  
    this_collapsed_full_biomt = collapsed_full_biomt[,which(metadata[,"phenotyping_center"] == study )]  
  }

  presence = (this_collapsed_full_biomt != 0)
  props_presence = apply(presence, FUN = sum, MAR = 1) / dim(presence)[2]
  #hist(props_presence, xlab = 'Proportion of rats with the bug', main = study)
  #almost all around 0 and 1!!

  medians = apply(this_collapsed_full_biomt, FUN = median, MAR =  1)
  sums = apply(this_collapsed_full_biomt, FUN = sum, MAR =  1)
  mox = max(sums)
  w = which(sums > quantile(sums, probs = seq(0,1,by=0.1))['90%'])
  pdf('/users/abaud/abaud/P50_HSrats/plots/filter_bugs_new.pdf')
  plot(props_presence,medians)
  plot(props_presence,medians, ylim = c(0,10))
  plot(props_presence,sums)
  plot(props_presence,sums, ylim = c(0,100000))
  plot(props_presence[w],sums[w], ylim = c(80,mox))
  dev.off()

  w = which(props_presence > 0.2 & props_presence < 0.8)
  print(paste('# taxa present/absent in at least 20%: ', length(w), paste = ''))  
  if (length(w) == 0) stop('no pass presence test FOR PRESENCE PHENOTYPE')
  filtered_presence = presence[w,]
 
  ###### main filtering step here
  tests = (props_presence > 0.5)
  print(paste('# taxa present in at least 50%: ', sum(tests), paste = ''))  
  commons[[study]] = rownames(this_collapsed_clr_counts)[tests]

  #somple = sample(commons[[study]], size = 16)
  #for (bug in somple) {
  #  hist(this_collapsed_clr_counts[bug,], xlab = 'CLR count', ylab = 'Number of rats', freq = TRUE, main = paste('Median relative abundance (%) is:',median_clr_abund[bug]), sub = paste('% of 0s =',round(sum(is.na(this_collapsed_full_biomt[bug,]))/length(this_collapsed_clr_counts[bug,]*100))))
  #}

  filtered_collapsed_clr_counts = this_collapsed_clr_counts[commons[[study]],]
  filtered_collapsed_full_biomt = this_collapsed_full_biomt[commons[[study]],]

  exclude_rats = colnames(filtered_collapsed_clr_counts)[apply(filtered_collapsed_clr_counts, FUN = my_all_is_NA, MAR = 2)]
  if (length(exclude_rats)!=0) {
    print('some rats whose bugs are not in common for that phenotyping center were excluded')
    filtered_collapsed_clr_counts = filtered_collapsed_clr_counts[,-match(exclude_rats,colnames(filtered_collapsed_clr_counts))]
    filtered_collapsed_full_biomt = filtered_collapsed_full_biomt[,-match(exclude_rats,colnames(filtered_collapsed_full_biomt))]
  }

  exclude_taxa = rownames(filtered_collapsed_full_biomt)[apply(filtered_collapsed_full_biomt, FUN = non_varying, MAR = 1)]
  if (length(exclude_taxa)!=0) {
    print('some non varying taxa excluded')
    filtered_collapsed_clr_counts = filtered_collapsed_clr_counts[-match(exclude_taxa,rownames(filtered_collapsed_clr_counts)),]
    filtered_collapsed_full_biomt = filtered_collapsed_full_biomt[-match(exclude_taxa,rownames(filtered_collapsed_full_biomt)),]
  }

  assign(paste('filtered_collapsed_full_biomt',study,sep='_'), filtered_collapsed_full_biomt)
  assign(paste('filtered_collapsed_clr_counts',study,sep='_'), filtered_collapsed_clr_counts)
  assign(paste('filtered_presence',study,sep='_'), filtered_presence)

}
dev.off()

inter1 = intersect(commons[['MI']],commons[['NY']])
inter2 = intersect(inter1,commons[['TN_behavior']])
length(inter2)
inter3 = intersect(inter2,commons[['TN_breeder']])
length(inter3)

# plot expression of these bugs in two centers
cols = rep('grey', dim(collapsed_clr_counts)[1])
cols[which(rownames(collapsed_clr_counts) %in% commons[['MI']])] = 'red'
cols[which(rownames(collapsed_clr_counts) %in% commons[['NY']])] = 'blue'
cols[which(rownames(collapsed_clr_counts) %in% commons[['TN_behavior']])] = 'yellow'
cols[which(rownames(collapsed_clr_counts) %in% intersect(commons[['MI']],commons[['NY']]))] = 'purple'
cols[which(rownames(collapsed_clr_counts) %in% intersect(commons[['MI']],commons[['TN_behavior']]))] = 'orange'
cols[which(rownames(collapsed_clr_counts) %in% intersect(commons[['NY']],commons[['TN_behavior']]))] = 'green'
cols[which(rownames(collapsed_clr_counts) %in% inter2)] = 'brown'
pdf('/users/abaud/abaud/P50_HSrats/plots/medians_inter_centers_taxa_deblur.pdf', width = 15)
par(mfrow = c(1,3))
median_MI = apply(collapsed_clr_counts[,which(metadata[,"phenotyping_center"] == 'MI')], FUN = median, MAR = 1)
median_NY = apply(collapsed_clr_counts[,which(metadata[,"phenotyping_center"] == 'NY')], FUN = median, MAR = 1)
median_TN = apply(collapsed_clr_counts[,which(metadata[,"phenotyping_center"] == 'TN')], FUN = median, MAR = 1)

plot(median_MI,median_NY, col = cols, xlab = 'Median in MI', ylab = 'Median in NY', main = 'MI red, NY blue, TN_behavior yellow')
abline(0,1,col = 'grey')
# !! all but purple won't show in plot below as log10 not defined for 0
#plot(log10(median_MI),log10(median_NY), col = cols, xlab = 'log10(median in MI)', ylab = 'log10(median in NY)')
#abline(0,1,col = 'grey')

plot(median_MI,median_TN, col = cols, xlab = 'Median in MI', ylab = 'Median in TN')
abline(0,1,col = 'grey')
# !! all but purple won't show in plot below as log10 not defined for 0
#plot(log10(median_MI),log10(median_TN), col = cols, xlab = 'log10(median in MI)', ylab = 'log10(median in TN)')
#abline(0,1,col = 'grey')

plot(median_NY,median_TN, col = cols, xlab = 'Median in NY', ylab = 'Median in TN')
abline(0,1,col = 'grey')
# !! all but purple won't show in plot below as log10 not defined for 0
#plot(log10(median_NY),log10(median_TN), col = cols, xlab = 'log10(median in NY)', ylab = 'log10(median in TN)')
#abline(0,1,col = 'grey')

dev.off()
#shows that median levels of expression correlate between two centers - whenever both have median != 0

save(filtered_presence_all, filtered_presence_MI, filtered_presence_NY, filtered_presence_TN_behavior, filtered_presence_TN_breeder, filtered_collapsed_clr_counts_MI, filtered_collapsed_clr_counts_NY, filtered_collapsed_clr_counts_TN_behavior, filtered_collapsed_clr_counts_TN_breeder, filtered_collapsed_clr_counts_all, filtered_collapsed_full_biomt_MI, filtered_collapsed_full_biomt_NY, filtered_collapsed_full_biomt_TN_breeder, filtered_collapsed_full_biomt_TN_behavior, filtered_collapsed_full_biomt_all, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/study_spe_taxa.RData')

