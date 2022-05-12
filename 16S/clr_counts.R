prefix = './data/metabo/65792_palmer_biom_original'
load(paste(prefix,'_sample_filtered_biomt.RData',sep=''))

#cols must be rats
#transformation done by column/rat
calculate_gm_mean = function(x, na.rm=TRUE){
	exp(sum(log(x), na.rm=na.rm) / length(x)) #double checked equivalence with equation (1) in https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full?report=reader; remember log (a^b) = blog(a) 
}

do_clr_default = function(data){
	#used to do this below but now add a little something to the entire table to begin with
	#if(any(data==0)) data = data + min(data[data>0])/2 #no negative values in data (counts); adds a little something that's half smallest non 0 count

	transform_col = function(x){
    	gMean = calculate_gm_mean(x)
    	transformed_col = log(x / gMean)
    	return(transformed_col)
	}
	transformed_data = apply(data, MAR = 2, FUN = transform_col)
	return(transformed_data)
}

offset = 0.00001
biomt = biomt + offset

clr_counts = do_clr_default(biomt)
#metabolits in rows, indivs in cols

#take back offset from biomt table BUT NOT CLR TABLE as use biomt to track presence later so need 0 for absence
# unnecessary as don't save biomt here: biomt = biomt - offset

save(clr_counts, file = paste(prefix,'_sample_filtered_clred_counts.RData',sep=''))

