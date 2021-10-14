prefix = './data/metabo/65792_palmer_biom_original'

load(paste(prefix, '_lib_sizes.RData',sep=''))
#distribution of lib sizes to determine the threshold for excluding samples with too low or too high lib sizes
pdf('./plots/filter_samples_metabo.pdf')
hist(lib_sizes[!grepl('blank',names(lib_sizes),ignore.case = T)], xlab = 'Sum of MS1 table per sample', breaks = 20, main = '', xlim = c(0, max(lib_sizes)))
hist(lib_sizes[grepl('blank',names(lib_sizes),ignore.case = T)],border = 'red', add = T, breaks = 20)
m = mean(lib_sizes[!grepl('blank',names(lib_sizes),ignore.case = T)])
s = sd(lib_sizes[!grepl('blank',names(lib_sizes),ignore.case = T)])
abline(v = m - 2*s, col = 'grey')
abline(v = m + 2*s, col = 'grey')
legend (x = 'topright', legend = c('Rats','Blanks'), col = c('black','red'), bty = 'n', lty = 1)
dev.off()

#start filtering samples
load(paste(prefix, '_biomt.RData',sep=''))
load('./data/metadata/metadata_augmented_16S_metabo.RData')

# first remove from biomt blanks (which are not in metadata), 
motch = match(colnames(biomt),metadata[,'sample_name_metabo'])
#check col in MS1 table that don't have a metadata: blank or needs checking
check = colnames(biomt)[is.na(motch)]
all(grepl('blank',check, ignore.case = T))
#FALSE
check[!grepl('blank',check, ignore.case = T)]
# "11479.00077E756F" "11479.00077E9E8F.BAD" 
#lib_sizes[grep('00077E756F', colnames(biomt))]
#            11479.00077E756F 11479.RATCECUM.MI.00077E756F 
#                  2,647,690,046                  33,466,639,939 
#so exclude 11479.00077E756F and the one that says BAD
#in other words the two cols that are not in the metadata and that are not blanks are fine to exclude
#below exclude those two + all blanks
biomt = biomt[,!is.na(motch)]
dim(biomt)
#7412 1150
#now remove bad rats (as determined by Apurva, ie swapped samples, or inconsistent...) and unknown rats (those that used to have a row in the metabo table but no real metadata) - genotyped or not
motch = match(colnames(biomt),metadata[which(metadata$bad_rat == FALSE & metadata$unknown_rat == FALSE),'sample_name_metabo'])
biomt = biomt[,!is.na(motch)]
dim(biomt)
#7412 1146
#put cols of biomt and rows or metadata in same order
motch = match(colnames(biomt),metadata$sample_name_metabo)
any(is.na(motch))
#FALSE
metadata = metadata[motch,]

#now filter out non genotyped first, then filter out on nb of reads; keep in two steps to have better visibility of number of samples excluded by each step

#step 1: exclude non genotyped first
keep = which(metadata$genotyped)  
biomt = biomt[,keep]
metadata = metadata[keep,]
dim(biomt)

#step 2: choose nb reads threshold; first update m and s as previously calculated including blanks, bad rats etc.
motch = match(colnames(biomt), names(lib_sizes))
any(is.na(motch))
#should be FALSE
lib_sizes = lib_sizes[motch]
m = mean(lib_sizes[!grepl('blank',names(lib_sizes),ignore.case = T)])
s = sd(lib_sizes[!grepl('blank',names(lib_sizes),ignore.case = T)])
th_low = m - 2*s
th_high = m + 2*s
#remove samples outside of lib siz range (!! that's range across filtered samples!)
keep = which(lib_sizes >= th_low & lib_sizes <= th_high)  
biomt = biomt[,keep]
dim(biomt)

save(biomt, file = paste(prefix,'_sample_filtered_biomt.RData',sep=''))

