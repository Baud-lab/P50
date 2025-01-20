#only works locally
#library(biomformat)
#fid = 
#read_biom(fid)
#save(biom, file = 'biom.RData')
#rsync -avz biom.RData abaud@ant-login.linux.crg.es:/users/abaud/abaud/P50_HSrats/data/felipes_deblur/

load('/users/abaud/abaud/P50_HSrats/data/felipes_deblur/175568_57950_analysis_16S_FilterfeaturesagainstreferencefilterfeaturesPhylogenetictreedatabasesgg202210202210taxonomyasvnwkqza.RData')
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
#93,090 ASVs for 4,422 samples - reported in paper
save(deblur, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/deblur_biom.RData')
save(taxonomy, file = '/users/abaud/data/secondary/P50_HSrats/felipes_deblur/taxonomy.RData')

my_strsplit = function(mot, code) {
    splot = strsplit(mot, '.', fixed = T)[[1]]
    if (length(splot) == 3) return(splot[code])
    if (grepl('blank',mot, ignore.case =T)) {
        if (code %in% c(1,2)) return(splot[code]) else return('BLANK')
    } 
    if (length(splot) == 4) return(splot[code])
    if (length(splot) == 5) {
        if (code %in% c(1,2)) return(splot[code]) else return(splot[code+2])
    }
    stop('pb')
}
corresp = rbind(c('4075','41941','175547'),c('4078','41937','175555'),c('6462','70124','175553'),c('6826','72629','175552'),c('6827','72564','175549'),c('6935','73651','175546'),c('7783','79794','175551'),c('7808','79988','175550'),c('7831','80252','175548'))
colnames(corresp) = c('prep','otu_biom_id','deblur_biom_id')

#get this
deblur_biom_ids = unlist(lapply(rownames(deblur),my_strsplit,1))
all(sort(corresp[,'deblur_biom_id']) == sort(unique(deblur_biom_ids)))
#TRUE
#I also checked prep IDs were correctly selected by F.
unique(unlist(lapply(rownames(deblur),my_strsplit,code =2)))
#11479
#get this too
deblur_preps = corresp[match(deblur_biom_ids, corresp[,'deblur_biom_id']),'prep']
#get this too
rats = unname(unlist(lapply(rownames(deblur),my_strsplit,code =3)))

#below characterisation
load('/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo.RData')
metadata$deblur_biom_id = NA
metadata$deblur_rooname = NA
for (i in 1:dim(deblur)[1]) {
    if (grepl('blank', rownames(deblur)[i], ignore.case = T)) next
    w = which(metadata$clean_ids == rats[i] & grepl(deblur_preps[i], metadata$prep_name_16S))
    if (length(w) == 0 & nchar(rats[i]) == 7) {
        rats[i] = paste('000',rats[i],sep='')
        w = which(metadata$clean_ids == rats[i] & grepl(deblur_preps[i], metadata$prep_name_16S))
        if (length(w) == 0) stop('pb')
    }
    if (length(w)>1) {
        print(metadata[w,c('full_ids')])
        w = w[1]
    }
    metadata[w,'deblur_biom_id'] = deblur_biom_ids[i]
    metadata[w,'deblur_rooname'] = rownames(deblur)[i]
}
check = metadata[is.na(metadata$deblur_rooname),]
table(check$in_16s_prep)
#  NO  YES 
#2814   30 - reported in paper
check[check$in_16s_prep == 'YES','nb_seqs_16S']
#all very low so dropped by deblur - ok - reported in paper
sum(grepl('blank', rownames(deblur), ignore.case = T))
#268 - reported in paper
save(metadata, file = '/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')


