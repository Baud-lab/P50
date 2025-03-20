## # Files to load split in smaller ones - less memory req - in dataPrep script
## processing_dir = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/'
## load(file.path(processing_dir,'full_biomt_clr_counts.RData')) # need to ask for 15GB of mem on cluster to do this; loading "clr_counts", "full_biomt"
## # Saving the two objects separately so that need less computing mem
## save(full_biomt, file = "/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_full_biomt.RData")
## save(clr_counts, file = "/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_clr_counts.RData")

suppressMessages(library("scales")) # needed for alpha() - box transparency

# Loading counts - raw or clr - files saved from dataPrep script
# TODO: change here depending if want to plot using raw counts or clr counts
data_type = "raw_counts"
cat("Loading",data_type,"\n")

## # Loading small file to test - cam be created using code in dataPrep script
## #asv = "ASV_5095"
## #load(paste0("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_",data_type,"_",asv,".Rdata")) # small file to test

if(data_type == "raw_counts"){
  load("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_full_biomt.RData") # loading "full_biomt"
  stopifnot("full_biomt" %in% ls())
  stopifnot(!"clr_counts" %in% ls()) # to avoid any possible confusion
  # will give error if loaded the wrong file
}else if( data_type == "clr_counts" ){
  load("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_clr_counts.RData") # loading "clr_counts"
  stopifnot("clr_counts" %in% ls())
  stopifnot(!"full_biomt" %in% ls()) # to avoid any possible confusion
  # will give error if loaded the wrong file
}else{
  stop("data_type need to be 'raw_counts' or 'clr_counts'")
}
ls()

# - Amelie comments -
# Reading and preparing genos 
cat("working on genos\n")
local_genos = read.table('/users/abaud/abaud/P50_HSrats/data/dosages/P50_Rn7_chr10qtl_allSNPS.raw', as.is = T, header = T, check.names = F)
# print(object.size(local_genos), units="Gb")
sample_names = local_genos[,2]
local_genos = local_genos[,-c(1:6)]
#lead SNP is in local_genos
dim(local_genos)
#[1] 17812  5103
length(sample_names)
#[1] 17812
my_strsplit = function(mot, div) {
	strsplit(mot, div, fixed = T)[[1]]
}
res1 = t(sapply(colnames(local_genos), my_strsplit, div = '_'))
res2 = t(sapply(res1[,1], my_strsplit, div = ':'))
snp_infos = data.frame(paste(res2[,1],res2[,2],sep='_'),as.numeric(sub('chr','',res2[,1])),as.numeric(res2[,2]), stringsAsFactors = F)
colnames(snp_infos) = c('id','chr','pos')
#w = which(snp_infos$pos == 101974959)
w = which(snp_infos$pos == 102006150)

genos = local_genos[,w]
names(genos) = sample_names
#some are NA

#genos[is.na(genos)] = "Missing"
chr_genos = as.character(genos)
chr_genos[genos == 0] = 'Single copy'
chr_genos[genos == 2] = 'Partial dupl./tripl.' #chr_genos[genos == 2] = 'Partial dupl./tripl.'
chr_genos[genos == 1] = 'Het'
chr_genos = factor(chr_genos, levels = c('Single copy','Het','Partial dupl./tripl.'))
names(chr_genos) = names(genos)

save_genos = chr_genos
# - end Amelie comments -


# Loading metadata
cat("Loading metadata\n")
load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData') # loads 'metadata'
# DEfine cohorts names as in paper
dict = c("NY"="NY", "MI"="MI", "TN_behavior"="TN1", "TN_breeder"="TN2")
metadata$study = unname(dict[metadata[,"study"]])

# Function to plot boxplot per each asv
plot_asv = function(asv){
  cat("doing asv ", asv,"\n")
  if(data_type == "raw_counts"){
    counts = t(full_biomt)[,asv, drop = F] #all rats together, original data
    # will give error of not found if loaded the wrong file
  }else if( data_type == "clr_counts" ){
    counts = t(clr_counts)[,asv, drop = F] #all rats together
  }else{
    stop("data_type need to be 'raw_counts' or 'clr_counts'")
  }
  dim(counts)
  
  for (study in c('NY','TN1','TN2','MI')) {
    if (asv == 'ASV_5163' & study == 'MI') next
    if (asv %in% c('ASV_5095','ASV_2821','ASV_17008') & study != 'MI') next

    # Selecting metadata related to study
    meta_oi = metadata[which(metadata$study == study),] #to select rats when full_biomt or clr_counts are used
    # Selecting rows of counts related to study and aligning meta_oi to it
    motch = match(rownames(counts), meta_oi$deblur_rooname)
    length(motch)
    counts_oi = counts[!is.na(motch),, drop = F] # selecting only the ones in the study
    meta_oi = meta_oi[na.exclude(motch),] 
    # Selecting genos related to study and aligning meta_oi and counts_oi to it
    motch = match(meta_oi$host_subject_id, names(save_genos))
    genos_oi = save_genos[na.exclude(motch)]
    meta_oi = meta_oi[!is.na(motch),]
    counts_oi = counts_oi[!is.na(motch),, drop = F]
    # Look at results
    table(genos_oi, meta_oi$sex, useNA = 'always')
    
    # x-axis label with sample size
    sizes = table(genos_oi)
    #xlabs = paste0(names(sizes), "\n(N=", sizes,")")
    xlabs = names(sizes)
    
    #boxplot(counts[,1] ~ meta_oi$sex + genos, varwidth = T, notch = T, outline = F, ylab = 'Raw counts Paraprevotella ASV_5163', xlab = 'Sex:Genotype')
    #males have fewer Paraprevotella - opposite of positive correlation between testosterone and Para in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7612624/
    #but consistent with slightly higher Para in WT females in https://www.frontiersin.org/files/Articles/330013/fmicb-09-01008-HTML/image_m/fmicb-09-01008-g005.jpg (Fig 5)
    
    #before QN, F0 M0 F1 high and M1 F2 M2 0
    if (asv == 'ASV_5095' | asv == 'ASV_5163') ylob = paste('Paraprevotella',asv,study,sep=' ')
    if (asv == 'ASV_2821') ylob = paste('A. muciniphila',asv,study,sep=' ')
    if (asv == 'ASV_17008') ylob = paste('Muribaculaceae',asv,study,sep=' ')
    #boxplot(counts[,1] ~ genos_oi, varwidth = T, notch = F, outline = F, ylab = 'CLR transformed counts', xlab = 'Genotype', main = moin, las = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
    #boxplot(counts_oi[,1] ~ meta_oi$sex + genos_oi, varwidth = F, notch = T, outline = T, ylab = 'Raw counts', xlab = 'Sex:Genotype', main = paste(asv,study,sep='_'))
    #par(font.main=4)
    
    formula = as.formula(counts_oi[,1] ~ genos_oi) 
    # formula = as.formula(counts_oi[,1] ~ meta_oi$sex + genos_oi) # other possibility 
    
    # Plot 
    # prepare base plot (no boxplot yet)
    bp = boxplot(formula, 
                 col = alpha("white", 0), border=alpha("white", 0), outline=F, ylab = "",xlab = '', xaxt="n",
                  main = "St6galnac1", drawRect=F, varwidth = TRUE, las = 1, cex.axis = 1.2, cex.main = 1.5)
    # add ticks and ticks labels to x axis
    axis(1, at=1:length(xlabs), labels=xlabs, cex.axis=1.25)
    # add y label
    title(ylab = paste0(ylob, ' (raw counts)'), cex.lab=1.4, line=3.8)
    
    # Define x positions for points
    at.dict = seq_along(sizes); names(at.dict) = names(sizes)
    nooutl = !names(counts_oi[,1]) %in% names(bp$out)
    Nout = table(at.dict[unname(genos_oi[!nooutl])]); names(Nout) = names(sizes)
    set.seed(20); 
    # add points
    points(x = jitter(unname(at.dict[unname(genos_oi[nooutl])]), factor=1), #c(1,0.9,0.5)),
           y = counts_oi[nooutl, 1], 
           pch=16, cex=1, 
           col = "#4DBBD5")
    # add boxes
    bp = boxplot(formula, 
                 col = alpha("white", 0.5),
                 outline=F, 
                 ylab = "",
                 xlab = '', xaxt="n",
                 drawRect=F,
                 varwidth = TRUE,
                 las = 1, cex.axis = 1.2, cex.main = 1.5,
                 add=T)  
    # add number of outliers at the bottom
    axis(1, at=1:length(xlabs), labels=paste0("(",Nout, " outliers)"), tck=F, lwd=0, line=1.5, cex.axis=1, font = 3)
  }
  
}

cat("Starting plot\n")
# Open pdf to save plot
pdf(paste0("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/all_chr10_boxplots_",data_type,".pdf"), h = 6, w = 7)
par(font.main=4, mar=c(5.1,6.1,2.1,1.1))

#plot_asv("ASV_5163")
plot_asv(asv="ASV_5095")
plot_asv("ASV_2821")
plot_asv("ASV_17008")

dev.off()

