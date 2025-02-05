#get full taxon from taxonomy. previously run
#load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/taxonomy.RData')
#parse_taxonomy = function(tax_word) {
#  splot = strsplit(tax_word,';')[[1]]
#  #treat deblur taxonomic calls
#  missing = 7 - length(splot)
#  if (missing>0) splot = c(splot,rep(NA,missing))
#  #treat OTU taxonomic calls
#  splot[splot %in% c("d__","p__","c__","o__","f__","g__","s__")] = NA
#  return(splot)
#}
#parsed_taxonomy = t(sapply(taxonomy[,'full_taxon'],FUN = parse_taxonomy))
#parsed_taxonomy = cbind(taxonomy, parsed_taxonomy)
#save(parsed_taxonomy, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/parsed_taxonomy.RData')


#call with source('/users/abaud/abaud/P50_HSrats/code/variance_decomposition/felipes_deblur/annotate_VCs_pvalues_function.R')
#will retrieve annotate

annotate = function(VCs) {
#expected VCs dataframe with column trait1 

	#create taxon from trait
	if ('trait1' %in% colnames(VCs)) VCs$taxon1 = VCs$trait1 else if ('measure' %in% colnames(VCs)) VCs$taxon1 = VCs$measure
	VCs$taxon1 = sub('_MI','',VCs$taxon1)
	VCs$taxon1 = sub('_NY','',VCs$taxon1)
	VCs$taxon1 = sub('_TN_breeder','',VCs$taxon1)
	VCs$taxon1 = sub('_TN_behavior','',VCs$taxon1)
#leave in as there are still some _all
	if (any(grepl('_all',VCs$taxon1))) VCs$taxon1 = sub('_all','',VCs$taxon1)

	#create study from trait and taxon
	VCs$study1 = NA
	for (i in 1:dim(VCs)[1]) {
		if ('trait1' %in% colnames(VCs)) VCs[i,'study1'] = sub(paste(VCs[i,'taxon1'],'_',sep=''),'',VCs[i,'trait1'], fixed = T) else #fixed = TRUE added later, check ok
			if ('measure' %in% colnames(VCs)) VCs[i,'study1'] = sub(paste(VCs[i,'taxon1'],'_',sep=''),'',VCs[i,'measure'], fixed = T)
	}
	if ('trait2' %in% colnames(VCs)) {
		VCs$study2 = NA
		for (i in 1:dim(VCs)[1]) {
			VCs[i,'study2'] = sub(paste(VCs[i,'taxon1'],'_',sep=''),'',VCs[i,'trait2'], fixed = T) #fixed = TRUE added later, check ok
		}	
	}

	load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/parsed_taxonomy.RData')

	VCs$full_taxon = NA
	for (i in 1:dim(VCs)[1]) {
	    w1 = which(parsed_taxonomy == VCs[i,'taxon1'], arr.ind = TRUE)
	    col = unique(w1[,'col'])
	    if (length(col) == 0) next
	    if (col == 3) col = 11 #col == 3 for ASVs; return 11th col which is s__
	    example_row = w1[1,'row']
	    VCs[i,'full_taxon'] = paste(parsed_taxonomy[example_row,5:col],collapse = ';')
	}

	VCs$tax_level = NA
	VCs[grepl('ASV', VCs$taxon1),'tax_level'] = 'ASV' 
	my_strsplit = function (mot) {
		tax_level = strsplit(mot,'__')[[1]][1]
	}
	VCs[grepl('__', VCs$taxon1),'tax_level'] = sapply(VCs[grepl('__', VCs$taxon1),'taxon1'], my_strsplit)
	VCs[!grepl('ASV', VCs$taxon1) & !grepl('__', VCs$taxon1),'tax_level'] = 'community_trait' 

	return(VCs)
}

# [1] "measure"    "marker"     "chr"        "pos"        "pvalue"    
# [6] "logP"       "ci_starts"  "ci_stops"   "merge"      "tax_level" 
#[11] "taxon1"     "study1"     "full_taxon"

