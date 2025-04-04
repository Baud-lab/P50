load('metadata_augmented_16S_metabo_deblur.RData')
load('deblur_biom.RData')
#remove from deblur blanks (which are not in metadata),  duplicates (whose secondary ids are not in metadata), bad rats and unknown rats - genotyped or not
motch = match(rownames(deblur),metadata[which(metadata$bad_rat == FALSE & metadata$unknown_rat == FALSE),'deblur_rooname'])
deblur = deblur[!is.na(motch),]
dim(deblur)
#[1]  4129 93090

motch = match(rownames(deblur),metadata$deblur_rooname)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]
#at this point no more blanks, duplicates, bad or unknown rats

#filter out non genotyped first, then filter out on nb of reads; keep in two steps to have better visibility of number of samples excluded by each step

#exclude non genotyped first
keep = which(metadata$genotyped)  
deblur = deblur[keep,]
dim(deblur)
#[1]  4107 93090
metadata = metadata[keep,]

m = mean(metadata$nb_seqs_16S)
s = sd(metadata$nb_seqs_16S)

#choose nb reads threshold
th_low = m - 2*s #6613.483 with deblur, very close to 6166 used on Qiita (from OTUs) so left
th_high = m + 2*s #34663.27, very close to 34856
#remove samples with fewer reads than threshold
keep = which(metadata$nb_seqs_16S >= th_low & metadata$nb_seqs_16S <= th_high)  
deblur = deblur[keep,]
metadata = metadata[keep,]

dim(deblur)
#[1]  3886 93090 

save(deblur, file = paste('deblur_genotyped_2sd.RData',sep=''))

