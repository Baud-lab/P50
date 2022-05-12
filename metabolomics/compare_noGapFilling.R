#needs RAM!
#biom convert -i /homes/abaud/P50_HSrats/data/metabo/no_gapFilling/65792_palmer_biom_original.biom -o /homes/abaud/P50_HSrats/data/metabo/no_gapFilling/65792_palmer_biom_original.txt --to-tsv
#biom convert -i /homes/abaud/P50_HSrats/data/metabo/no_gapFilling/65793_palmer_biom_normalized.biom -o /homes/abaud/P50_HSrats/data/metabo/no_gapFilling/65793_palmer_biom_normalized.txt --to-tsv


biomt = read.delim('/homes/abaud/P50_HSrats/data/metabo/pre_Apr8/no_gapFilling/65792_palmer_biom_original.txt',as.is = T, header = F)
samples = unlist(biomt[2,-1, drop = T])
biomt = biomt[-c(1,2),]
MZs = paste('MZ',biomt[,1],sep='')
biomt = biomt[,-1]
biomt = lapply(biomt, FUN = as.numeric)
biomt = as.matrix(as.data.frame(biomt, stringsAsFactors = F))
colnames(biomt) = samples
rownames(biomt) = MZs
unnormed = biomt
save(unnormed, file = '/homes/abaud/P50_HSrats/data/metabo/pre_Apr8/no_gapFilling/65792_palmer_biom_original.RData')

#### comparing now

load('/homes/abaud/P50_HSrats/data/metabo/pre_Apr8/no_gapFilling/65792_palmer_biom_original.RData')
original = unnormed
load('/homes/abaud/P50_HSrats/data/metabo/115698_palmer_biom_march_2021_original.RData')
gapFilled = biomt

dim(original)
#7,412
dim(gapFilled)
#11,388

my_sum_zeros = function(row) {
	return(sum(row == 0))
}
par(mfrow = c(1,2))
hist(apply(original, FUN = my_sum_zeros, MAR = 1))
#most have mostly zeros
hist(apply(gapFilled, FUN = my_sum_zeros, MAR = 1))
#has been gapFilled

###### get back here!

inter = intersect(colnames(original),colnames(gapFilled))
colnames(unnormed)[!colnames(unnormed) %in% inter]
#None
colnames(gapFilled)[!colnames(gapFilled) %in% inter]
#                   V16                  V1082                  V1245 
#    "11479.00077E756F" "11479.00077E9E8F.BAD"       "11479.Blank.93" 
#                 V1246 
#      "11479.Blank.94" 
unnormed = unnormed[,match(inter, colnames(unnormed))]
gapFilled = gapFilled[,match(inter, colnames(gapFilled))]

all(rownames(gapFilled) %in% rownames(unnormed))
all(rownames(unnormed) %in% rownames(gapFilled))

setA = sort(rownames(unnormed)[grepl('MZ105.', rownames(unnormed), fixed = T)])
setB = sort(rownames(gapFilled)[grep('MZ105.', rownames(gapFilled), fixed = T)])
pdf('/homes/abaud/P50_HSrats/plots/check_gapFilling_new.pdf', width = 20, height = 20)
par(mfrow = c(4,3))
for (j in setB) {
	for (i in setA) {
		plot(unnormed[i,], gapFilled[j,], xlab = paste('unnormed',i), ylab = paste('gapFilled',j), main = paste(i,j))
		abline(0,1,col = 'red')
	}
}
dev.off()
#either perfectly correlated or 0 in unnormed when MZ is same +/- 0.0001 (more flexible on RT)
#rest of the time not very highly correlated
#so it looks like features were dropped on the one hand, and the features that were retained were filled


setA = sort(rownames(unnormed)[grepl('MZ433.', rownames(unnormed), fixed = T) & grepl(';6.8', rownames(unnormed), fixed = T)])
setB = sort(rownames(gapFilled)[grepl('MZ433.', rownames(gapFilled), fixed = T) & grepl(';6.8', rownames(gapFilled), fixed = T)])

pdf('/homes/abaud/P50_HSrats/plots/check_gapFilling.pdf', width = 20, height = 20)
par(mfrow = c(4,3))
for (i in setA) {
	for (j in setB) {
		plot(unnormed[i,], gapFilled[j,], xlab = paste('unnormed',i), ylab = paste('gapFilled',j))
		abline(0,1,col = 'red')
	}
}
dev.off()

# here there is a case of having both 0s in unnormed, perfectly correlated values *and non perfectly correlated values*. 
# so some values borrowed from feature in original table (setA) and assigned to (slightly different) feature in gap filled table (setB)
# some values were 0 in original table and were "saved" to non 0 by gap filling 
# some values were perhaps borrowed from another feature in the original table (setA) and assigned to that feature in setB
# so those new setB features would be compound features in addition to being imputed. would need to be within MZ and RT range

#blue plot above shows that one feature has tons of 0s when the other one does not and for pairs of non 0 values,
#correlation is very high but not along diagonal (greater values for feature that has no 0s)
#looks like same true peak detected as two peaks but one falls below threshold for some samples
#"good" feature is MZ433.3313;6.8099, which used to give logP 11 QTLon chr14 but is further away in MZ from feature in gapFilling
#seems like feature in gapFilling table is centered on "bad" feature, which explains why I can't find the high logP peak.

#one feature (not the one that gave logP 11) took precedence when it came to gap filling - too bad
#black dots got imputed too - from where??
