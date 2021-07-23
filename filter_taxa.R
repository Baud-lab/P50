prefix = './data/metabo/65792_palmer_biom_original'
load(paste(prefix,'_sample_filtered_clred_counts.RData',sep=''))

#need to match on rats
load(paste(prefix,'_sample_filtered_biomt.RData',sep=''))
motch = match(colnames(clr_counts), colnames(biomt))
any(is.na(motch))
#FALSE
biomt = biomt[,motch]

#need to match on rats
load('./data/metadata/metadata_augmented_16S_metabo.RData')
motch = match(colnames(clr_counts), metadata$sample_name_metabo)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

my_all_is_NA = function(x) {
  return(all(is.na(x)))
}

non_varying = function(x) {
  return(length(unique(x)) == 1)
}

#commons below tracks the name of the features that are common, as defined above 
commons = list()
#pdf('./plots/filter_bugs.pdf')
#par(mfrow = c(2,2))
for (study in c('all','MI','NY')) {
  print(study)
  if (study=='all') {
    this_clr_counts = clr_counts[,which(metadata[,"phenotyping_center"] %in% c('MI','NY'))]  
    this_biomt = biomt[,which(metadata[,"phenotyping_center"] %in% c('MI','NY'))]  
  } else {
    this_clr_counts = clr_counts[,which(metadata[,"phenotyping_center"] == study)]  
    this_biomt = biomt[,which(metadata[,"phenotyping_center"] == study )]  
  }

  #create a new table that indicates whether the feature is present or absent in the rat - this is where it is important to have removed the offset from the biom table in clr_counts.R
  this_presence = (this_biomt != 0)
  #calculate the fraction of rats in which the feature is present (remember that rats are in columns):
  props_presence = apply(this_presence, FUN = sum, MAR = 1) / dim(this_presence)[2]
  #hist(props_presence, xlab = 'Fraction of rats with the feature', main = study)
  #almost all around 0 and 1!!
  #now FOR THE PRESENCE PHENOTYPE keep only features that are present in at least 20% of the rats AND absent in at least 20% of the rats
  #### !!! HARDCODED
  w = which(props_presence > 0.2 & props_presence < 0.8)
  print(paste('# metabolites present/absent in at least 20%: ', length(w), paste = ''))  
  if (length(w) == 0) stop('no pass presence test FOR PRESENCE PHENOTYPE')
  #filter down the presence table to keep only those features satisfying the requirement
  this_filtered_presence = this_presence[w,]
   print(dim(this_filtered_presence))

  ###### now for the COUNTS PHENOTPE keep only features that are present in at least 50% of the rats
  # filter out features FOR COUNTS PHENOTYPE: requirement is that feature is present in at least presence_th of the samples (here 50%, called "common")
  #### !!! HARDCODED
  tests = (props_presence > 0.5)
  print(paste('# metabolites present in at least 50%: ', sum(tests), paste = ''))  
  this_filtered_clr_counts = this_clr_counts[tests,]
  print(dim(this_filtered_clr_counts))
  
  #old below needs checking
  #somple = sample(commons[[study]], size = 16)
  #for (bug in somple) {
  #  hist(this_clr_counts[bug,], xlab = 'CLR count', ylab = 'Number of rats', freq = TRUE, main = paste('Median relative abundance (%) is:',median_clr_abund[bug]), sub = paste('% of 0s =',round(sum(is.na(this_biomt[bug,]))/length(this_clr_counts[bug,]*100))))
  #}

  #can we do without this_filtered_biomt?
  assign(paste('filtered_clr_counts',study,sep='_'), this_filtered_clr_counts)
  assign(paste('filtered_presence',study,sep='_'), this_filtered_presence)

}
#dev.off()

inter1 = intersect(rownames(filtered_clr_counts_MI),rownames(filtered_clr_counts_NY))
length(inter1)

# plot expression of these metabolites in two centers
cols = rep('grey', dim(clr_counts)[1])
cols[which(rownames(clr_counts) %in% rownames(filtered_clr_counts_MI))] = 'red'
cols[which(rownames(clr_counts) %in% rownames(filtered_clr_counts_NY))] = 'blue'
cols[which(rownames(clr_counts) %in% inter1)] = 'purple'

pdf('./plots/medians_inter_centers_metabolites.pdf', width = 15)
par(mfrow = c(1,2))
median_MI = apply(clr_counts[,which(metadata[,"phenotyping_center"] == 'MI')], FUN = median, MAR = 1)
median_NY = apply(clr_counts[,which(metadata[,"phenotyping_center"] == 'NY')], FUN = median, MAR = 1)
plot(median_MI,median_NY, col = cols, xlab = 'Median in MI', ylab = 'Median in NY', main = 'MI red, NY blue')
plot(log10(median_MI),log10(median_NY), col = cols, xlab = 'log10 median in MI', ylab = 'log10 median in NY', main = 'MI red, NY blue')
# !! in plot above all but purple won't show in plot below as log10 not defined for 0
abline(0,1,col = 'grey')
dev.off()
#shows that median levels of expression correlate between two centers

save(filtered_presence_all, filtered_presence_MI, filtered_presence_NY, filtered_clr_counts_MI, filtered_clr_counts_NY, filtered_clr_counts_all, file = paste(prefix,'_study_spe.RData',sep=''))
