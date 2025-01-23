load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData')
all(rownames(clr_counts) == names(full_biomt))
#TRUE
load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/taxonomy.RData')

# now get new full_biomt with collapses at various taxonomic levels
all(rownames(clr_counts) == taxonomy[,'asv_id'])
#TRUE
parse_taxonomy = function(tax_word) {
  splot = strsplit(tax_word,';')[[1]]
  #treat deblur taxonomic calls
  missing = 7 - length(splot)
  if (missing>0) splot = c(splot,rep(NA,missing))
  #treat OTU taxonomic calls (old)
  splot[splot %in% c("k__","p__","c__","o__","f__","g__","_s__")] = NA
  return(splot)
}
parsed_taxonomy = t(sapply(taxonomy[,'full_taxon'],FUN = parse_taxonomy))
rownames(parsed_taxonomy) = taxonomy[,'asv_id']

#parsed_taxonomy = cbind(taxonomy, parsed_taxonomy)
#save(parsed_taxonomy, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/parsed_taxonomy.RData')

my_f = function(taxon) {
  w = which(parsed_taxonomy[,k] == taxon)
  full_taxon = paste(parsed_taxonomy[w[1],1:k], collapse = ';')
  ASVs = rownames(parsed_taxonomy)[which(parsed_taxonomy[,k] == taxon)]

  this_clr_count = apply(clr_counts[ASVs,,drop = F], FUN = sum, MARGIN = 2, na.rm = T)
  this_raw_count = apply(full_biomt[ASVs,,drop = F], FUN = sum, MARGIN = 2)

  return(list(this_clr_count, this_raw_count, taxon, full_taxon))
}

collapsed_clr_counts = NULL
collapsed_full_biomt = NULL
full_taxa = c()
#dim(parsed_taxonomy)[2] = 7 different levels
for (k in 1:dim(parsed_taxonomy)[2]) {
  #all choices at given level

  res = lapply(na.exclude(unique(parsed_taxonomy[,k])), my_f)
  
  this_collapsed_clr_counts = do.call('rbind',lapply(res,'[[',1))
  rownames(this_collapsed_clr_counts) =  unlist(lapply(res,'[[',3))
  collapsed_clr_counts = rbind(collapsed_clr_counts, this_collapsed_clr_counts)
 
  this_collapsed_full_biomt = do.call('rbind',lapply(res, '[[',2))
  rownames(this_collapsed_full_biomt) =  unlist(lapply(res,'[[',3))
  collapsed_full_biomt = rbind(collapsed_full_biomt, this_collapsed_full_biomt)

  this_full_taxa = unlist(lapply(res,'[[',4))
  names(this_full_taxa) = rownames(this_collapsed_full_biomt)
  full_taxa = c(full_taxa, this_full_taxa)
}

colnames(collapsed_clr_counts) = colnames(clr_counts)
colnames(collapsed_full_biomt) = colnames(clr_counts)

save(collapsed_clr_counts, collapsed_full_biomt, full_taxa, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')

