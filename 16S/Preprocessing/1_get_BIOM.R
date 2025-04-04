#biomformat only works locally
#library(biomformat)

#download artefact 175568 (left-most in the analysis pipeline) from Qiita's analysis 57950
#fid = '175568_57950_analysis_16S_FilterfeaturesagainstreferencefilterfeaturesPhylogenetictreedatabasesgg202210202210taxonomyasvnwkqza.biom' 
#biom = read_biom(fid)
#save(biom, file = '175568_57950_analysis_16S_FilterfeaturesagainstreferencefilterfeaturesPhylogenetictreedatabasesgg202210202210taxonomyasvnwkqza.RData')
#rsync 

load('175568_57950_analysis_16S_FilterfeaturesagainstreferencefilterfeaturesPhylogenetictreedatabasesgg202210202210taxonomyasvnwkqza.RData')
deblur = as.matrix(as.data.frame(biom$data))
#not intuitive but correct...
rownames(deblur) = unlist(biom$columns)
taxonomy = read.table('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/taxonomy_Greengenes2.txt', as.is = T)
motch = match(unlist(biom$rows), taxonomy[,1])
sum(is.na(motch))
#0
table(table(motch))
#    1 
#93090 
taxonomy = taxonomy[motch,]
all(unlist(biom$rows) == taxonomy[,1])
#TRUE
taxonomy[,3] = paste('ASV_',1:dim(taxonomy)[1],sep='')
colnames(taxonomy) = c('seq','full_taxon','asv_id','Geez')
colnames(deblur) = taxonomy[,'asv_id']
dim(deblur)
#93,090 ASVs for 4,422 samples 
save(deblur, file = 'deblur_biom.RData')
save(taxonomy, file = 'taxonomy.RData')
