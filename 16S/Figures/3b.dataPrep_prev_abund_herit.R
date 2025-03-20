# Load microbiome data
# for ASVs
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData') # loads 'clr_counts' and 'full_biomt' - 2.7GB each 
outfile = '/users/abaud/abaud/P50_HSrats/output/prev_abund_asvs_biomt.RData'

## for taxa 
# TODO: uncomment following lines to get file for taxa
#load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData') # loads 'collapsed_clr_counts' 'collapsed_full_biomt' 'full_taxa' - 200mb all
#outfile = '/users/abaud/abaud/P50_HSrats/output/prev_abund_taxa_biomt.RData'

# Load metadata
load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
if (any(grepl('collapsed_clr_counts', ls()))) motch = match(colnames(collapsed_clr_counts),metadata$deblur_rooname) else motch = match(colnames(clr_counts),metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

# Function to go from raw abundance to relative abundance
my_relativise = function(row) {
  row = row / sum(row)
  return(row)
}

# Depending on whether CLRed or not CLRed data were uploaded to R, grab the right object and filter based on cohort (study)
for (study in c('MI','NY','TN_behavior','TN_breeder')) {
  print(study)
  if (!any(grepl('collapsed_clr_counts', ls()))) this_clr_counts = clr_counts[,which(metadata[,"study"] == study)]   else this_clr_counts = collapsed_clr_counts[,which(metadata[,"study"] == study)]  
  if (!any(grepl('collapsed_full_biomt', ls()))) this_full_biomt = full_biomt[,which(metadata[,"study"] == study)]   else this_full_biomt = collapsed_full_biomt[,which(metadata[,"study"] == study)]  
  
  #this_full_biomt has samples in columns
  
  # calculate prevalence
  presence = (this_full_biomt != 0)
  props_presence = apply(presence, FUN = sum, MAR = 1) / dim(presence)[2]
  if (any(grepl('collapsed_clr_counts', ls()))) assign(paste("collapsed_prevalence", study, sep='_'), props_presence) else assign(paste("prevalence", study, sep='_'), props_presence)
  
  # calculate median relative abundance
  rel_abundances = t(apply(this_full_biomt, MAR = 2, FUN = my_relativise)) #samples in rows now
  medians = apply(rel_abundances, FUN = median, MAR =  2)
  if (any(grepl('collapsed_clr_counts', ls()))) assign(paste("collapsed_median", study, sep='_'), medians) else assign(paste("median", study, sep='_'), medians)
}

# save prevalence and median relative abundance data to files
# selecting the object to save from the list of objects
lost = ls() #lost is just a name for the list retrieved by ls()
toSave = lost[grepl('prevalence', lost) | grepl('median', lost)]
save(list = toSave, file = outfile)
