###### Saving - this will go in another file
processing_dir = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/'
load(file.path(processing_dir,'full_biomt_clr_counts.RData')) # 10 GB no sufficient # with 15 are enough; loading "clr_counts", "full_biomt"
# Saving the two objects separately so that need less computing mem
save(full_biomt, file = "/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_full_biomt.RData")
save(clr_counts, file = "/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_clr_counts.RData")

# print(object.size(c(clr_counts, full_biomt)), units="Gb")
# dim(clr_counts)


#### Small file to test
### data_type = "raw_counts"
### asv = "ASV_5095"
### 
### if(data_type == "raw_counts"){
###   load("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_full_biomt.RData") # loading "full_biomt"
###   stopifnot("full_biomt" %in% ls())
###   stopifnot(!"clr_counts" %in% ls()) # to avoid any possible confusion
###   # will give error if loaded the wrong file
### }else if( data_type == "clr_counts" ){
###   load("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_clr_counts.RData") # loading "clr_counts"
###   stopifnot("clr_counts" %in% ls())
###   stopifnot(!"full_biomt" %in% ls()) # to avoid any possible confusion
###   # will give error if loaded the wrong file
### }else{
###   stop("data_type need to be 'raw_counts' or 'clr_counts'")
### }
### 
### 
### if(data_type == "raw_counts"){
###  counts = t(full_biomt)[,asv, drop = F] #all rats together, original data
###  # will give error of not found if loaded the wrong file
### }else if( data_type == "clr_counts" ){
###  counts = t(clr_counts)[,asv, drop = F] #all rats together
### }else{
###  stop("data_type need to be 'raw_counts' or 'clr_counts'")
### }
### 
### save(data_type,asv,counts, file = paste0("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_",data_type,"_",asv,".Rdata"))  ### 