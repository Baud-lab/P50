#this script filters the data collapsed at the different taxonomic level 
#similar script used for uncollapsed/asv level data 

#for the abundance phenotypes reported in the paper, only taxa present in at least 50% of the rats in the considered cohort will be kept

load('collapsed_full_biomt_collapsed_clr_counts.RData')

load('metadata_augmented_16S_metabo_deblur.RData')
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

commons = list()
for (study in c('all','MI','NY','TN_behavior','TN_breeder')) { #TN_behavior is TN1 in paper; TN_breeder is TN2
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
 
  ###### main filtering step here
  tests = (props_presence > 0.5)
  print(paste('# taxa present in at least 50%: ', sum(tests), paste = ''))  
  commons[[study]] = rownames(this_collapsed_clr_counts)[tests]

  filtered_collapsed_clr_counts = this_collapsed_clr_counts[commons[[study]],]
  filtered_collapsed_full_biomt = this_collapsed_full_biomt[commons[[study]],]

  assign(paste('filtered_collapsed_full_biomt',study,sep='_'), filtered_collapsed_full_biomt)
  assign(paste('filtered_collapsed_clr_counts',study,sep='_'), filtered_collapsed_clr_counts)

}

save(filtered_collapsed_clr_counts_MI, filtered_collapsed_clr_counts_NY, filtered_collapsed_clr_counts_TN_behavior, filtered_collapsed_clr_counts_TN_breeder, filtered_collapsed_clr_counts_all, filtered_collapsed_full_biomt_MI, filtered_collapsed_full_biomt_NY, filtered_collapsed_full_biomt_TN_breeder, filtered_collapsed_full_biomt_TN_behavior, filtered_collapsed_full_biomt_all, file = 'study_spe_taxa.RData')

