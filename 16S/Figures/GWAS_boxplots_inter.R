local_genos = read.table('/users/abaud/abaud/P50_HSrats/data/dosages/P50_Rn7_chr10qtl_allSNPS.raw', as.is = T, header = T, check.names = F)
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
w = which(snp_infos$pos == 102006150)
genos = local_genos[,w]
names(genos) = sample_names
#some are NA
genos[is.na(genos)] = "Missing"
save_genos = genos


local_genos = read.table('/users/abaud/abaud/P50_HSrats/data/dosages/P50_Rn7_chr1qtl_allSNPS.raw', as.is = T, header = T, check.names = F)
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
w2 = which(snp_infos$pos == 196217481)
genos2 = local_genos[,w2]
names(genos2) = sample_names
genos2[is.na(genos2)] = "Missing"
save_genos2 = genos2

processing_dir = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/'
load(paste(processing_dir,'full_biomt_clr_counts.RData',sep=''))



pdf('/users/abaud/abaud/P50_HSrats/plots/all_chr10_boxplots_clr_inter.pdf', width = 10, height = 10)
par(mfrow = c(2,2))

for (asv in c('ASV_5163','ASV_5095','ASV_2821','ASV_17008')) {
    #counts = t(full_biomt)[,asv, drop = F] #all rats together, original data
    counts = t(clr_counts)[,asv, drop = F] #all rats together
    #load(paste(processing_dir,'study_spe_uncollapsed_taxa.RData',sep='')) #before QN
    #counts = t(filtered_clr_counts_NY)[,'ASV_5163', drop = F]
    #load(paste(processing_dir,'NONresids_qned_counts_uncollapsed.RData',sep='')) #before accounting for covs
    #counts = t(qned_counts_NY)[,'ASV_5163', drop = F]
    #load(paste(processing_dir,'resids_qned_counts_uncollapsed.RData',sep=''))
    #counts = t(resids_qned_counts_NY)[,'ASV_5163', drop = F]
    save_counts = counts

    for (study in c('NY','TN_behavior','TN_breeder','MI')) {

        genos = save_genos
        genos2 = save_genos2
        counts = save_counts

        load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')
        metadata = metadata[which(metadata$study == study),] #to select rats when full_biomt or clr_counts are used
        motch = match(rownames(counts), metadata$deblur_rooname)
        metadata = metadata[na.exclude(motch),]
        counts = counts[!is.na(motch),, drop = F]
        motch = match(metadata$host_subject_id,	names(genos))
        metadata = metadata[!is.na(motch),]
        counts = counts[!is.na(motch),, drop = F]
        genos = genos[na.exclude(motch)]
        genos2 = genos2[na.exclude(motch)]

        if (asv == 'ASV_5095' | asv == 'ASV_5163') moin = paste('Paraprevotella',asv,study,sep=' ')
        if (asv == 'ASV_2821') moin = paste('A. muciniphila',asv,study,sep=' ')
        if (asv == 'ASV_17008') moin = paste('Muribaculaceae',asv,study,sep=' ')
        boxplot(counts[,1] ~ genos*genos2, varwidth = T, notch = F, outline = F, ylab = 'CLR transformed counts', xlab = 'Genotype', main = moin, las = 1, cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
#        boxplot(counts[,1] ~ metadata$sex + genos, varwidth = F, notch = T, outline = T, ylab = 'Raw counts', xlab = 'Sex:Genotype', main = paste(asv,study,sep='_'))
    }
}
dev.off()

anova(lm(counts[,1] ~ metadata$sex * genos))
#with full_biomt
Analysis of Variance Table

Response: counts[, 1]
                     Df  Sum Sq Mean Sq F value    Pr(>F)    
metadata$sex          1   17472   17472  3.1806    0.0748 .  
genos                 1  232901  232901 42.3978 1.144e-10 ***
metadata$sex:genos    1     323     323  0.0588    0.8084    
Residuals          1069 5872262    5493                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#with clr_counts: both sex and genos effects more pronounced
Analysis of Variance Table

Response: counts[, 1]
                     Df Sum Sq Mean Sq F value    Pr(>F)    
metadata$sex          1   1261  1261.3 25.4315 5.381e-07 ***
genos                 1   4814  4813.9 97.0629 < 2.2e-16 ***
metadata$sex:genos    1     67    66.8  1.3466    0.2461    
Residuals          1069  53017    49.6                      
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1






