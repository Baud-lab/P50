#change all full_biomt into collapsed_full_biomt and vice versa
load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
load(paste('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/deblur_genotyped_2sd.RData',sep=''))

motch = match(rownames(deblur),metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

full_biomt = t(deblur)
#[1] 93090  3886

#cols must be rats
#transformation done by column/rat
calculate_gm_mean = function(x, na.rm=TRUE){
	exp(sum(log(x), na.rm=na.rm) / length(x)) #double checked equivalence with equation (1) in https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full?report=reader; remember log (a^b) = blog(a) 
}

do_clr_default = function(data){
	transform_col = function(x){
    	gMean = calculate_gm_mean(x)
    	transformed_col = log(x / gMean)
    	return(transformed_col)
	}
	transformed_data = apply(data, MAR = 2, FUN = transform_col)
	return(transformed_data)
}

#absolutely needs to be done on the ASV rather than collapsed taxon table as otherwise reads are counted twice!! (or one per taxonomic level)
offset = 0.00001
full_biomt = full_biomt + offset

geom_means = apply(full_biomt, FUN = calculate_gm_mean, MAR = 2)
plot(geom_means, metadata[,'nb_seqs_16S'])
cor(geom_means, metadata[,'nb_seqs_16S'])
#cor = 0.604672

clr_counts = do_clr_default(full_biomt)
#bugs in rows, indivs in cols

#### double check implementation of CLR
#library(rgr)
#counts2 = t(clr(t(full_biomt + min(full_biomt[full_biomt>0])/2)))
#any(abs(counts - counts2) > 0.00001)
#FALSE so the same

#### For info
#means_raw = apply(full_biomt, FUN = mean, MAR = 2)
#sds_raw = apply(full_biomt, FUN = sd, MAR = 2)
#means_tr = apply(counts, FUN = mean, MAR = 2)
#sds_tr = apply(counts, FUN = sd, MAR = 2)
#pdf('/homes/abaud/P50_HSrats/plots/mean_sd_relationship_befAf_CLR.pdf')
#plot(means_raw, sds_raw)
#plot(means_tr, sds_tr)
#dev.off()

#take back offset as use full_biomt (not in clr table) to track presence/absence from full_biomt table later
full_biomt = full_biomt - offset

save(full_biomt, clr_counts, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData')

pdf('/users/abaud/abaud/P50_HSrats/plots/behaviour_CLR_deblur.pdf')
somple = sample(1:dim(full_biomt)[2], size = 50)
for(i in somple) {
	plot(full_biomt[,i], clr_counts[,i])
}
dev.off()


