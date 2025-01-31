library(RColorBrewer)
library(rhdf5)
library(parallel)

load('/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/QTLs_alpha1e-04_unpruned.RData')
unpruned_bug_QTLs = unpruned_bug_QTLs[unpruned_bug_QTLs$tax_level != 'community_trait',]
DGE_QTLs = unpruned_bug_QTLs
DGE_QTLs = DGE_QTLs[DGE_QTLs$logP > 5.8,]

#leave in as there are still some _all
#       VCs$taxon1 = sub('_all','',VCs$taxon1)

#taxon = 'ASV_33692'
#alpha = 0.01
#load(paste(pvalues_dir_DGE,'/QTLs_alpha', alpha, '_unpruned_porcupine_', taxon,'.RData',sep=''))
#files2plot = files_DGE[grep(taxon, files_DGE)]
#pheno_names = sub('.h5','',files2plot)

load('/users/abaud/abaud/P50_HSrats/data/cumpos_P50_rats_Rn7.RData')

draw_manhattan_plot=function(dat) {
	cols <- c( "lightgrey", "lightblue")
	#cols <-"lightgrey"
	dat$colour=cols[1]
	dat[which(dat[,'chr'] %in% c(2,4,6,8,10,12,14,16,18,20)),'colour']=cols[2]
	dat$cex = 0.4
	w = which(dat[,'col'] != 'darkgrey')
	if (length(w)!=0) {
		dat[w,'cex'] = 2
		dat[w,'colour'] = dat[w,'col']
	}
	
	#need both order first by decreasing logP then by colour otherwise QTL point gets overpainted so is not visible
	dat = dat[order(dat$colour != cols[1]),]
	r = range(c(0,range(dat$logP)))
	plot(dat[,c('cumpos','logP')],ylab="",col=dat[,'colour'],xlab="",main='',cex=dat[,'cex'],pch=16,axes=FALSE,xaxs="i")
	
	axis(2,cex.axis=2,las = 1)

	mtext("-logP",2,padj=-3,cex=2)

	mtext(c(1:20),1,at=tapply(dat$cumpos,dat$chr,median)+40000,adj=1,las=1,cex=1.5)
	mtext('Chromosome',1,at=median(dat$cumpos),padj=2,las=1,cex=1.5)
}

draw_legend = function(dat) {
	#otherwise QTL point gets overpainted so is not visible
	plot(0,1,col='white')
	dat = dat[dat$col != 'darkgrey',] #needed otherwise match below goes to first occurence, which is darkgrey
	w = which(dat$logP > 5.8 & paste(dat$chr, dat$pos, sep=':') %in% c(all_ld$SNP_A,all_ld$SNP_B))
	sigs = dat[w,'full_taxon']
	cols = dat[match(sigs, dat$full_taxon),'col']

	ids = paste(dat[w,'chr'], dat[w,'pos'],sep=':')
	tops = c()
	for (id in ids) {
		w = which(all_ld$SNP_A == id | all_ld$SNP_B == id)
		top = unique(all_ld[w,'top'])
		if (length(top) != 1) stop('pb') 
		tops = c(tops, top)
	}	
	sigs = sigs[order(tops)]
	cols = cols[order(tops)]

	sigs = unique(sigs) #unique because of same taxon across cohorts
	cols = unique(cols)
	legend(x = 'top', legend = sigs, col = cols, border = 'white', pch = 16, cex = 1, y.intersp = 1.2, box.col = 'white')
}


map2color<-function(x,pal,limits=NULL){
    if(is.null(limits)) limits=range(x)
    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
#colours = sample(map2color(1:length(pheno_names),rainbow(length(pheno_names)),limits=NULL))

##colours = sample(map2color(1:length(unique(taxa)),rainbow(length(unique(taxa))),limits=NULL))
library(RColorBrewer)

#alternative1
DGE_QTLs = DGE_QTLs[grep('g__',DGE_QTLs$full_taxon),] #will not show sig QTLs that are only at the family level for example
pheno_names = unique(DGE_QTLs$measure) # has study extensions, includes both taxa and ASVs (but not community levels phenotypes)
#### bug?? see chr10 
my_extract_genus = function(full_taxon) {
	splot = strsplit(full_taxon,';')[[1]]
	splot = splot[grep('g_',splot)]
	genus = strsplit(splot, '__')[[1]][2]
}
DGE_QTLs$genus = sapply(DGE_QTLs$full_taxon, my_extract_genus) #no NA due to prior exclusion of phenotypes without genus taxonomy
taxa = DGE_QTLs[match(pheno_names, DGE_QTLs$measure),'genus'] #genera corresponding to unique measures (after excluding those mapping to higher level than taxon)

#alternative2
#pheno_names = unique(DGE_QTLs$measure) # has study extensions
#taxa = sub('_MI','',pheno_names)
#taxa = sub('_NY','',taxa)
#taxa = sub('_TN_breeder','',taxa)
#taxa = sub('_TN_behavior','',taxa)


n <- length(unique(taxa))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

col_vector = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colours = sample(col_vector, n)
#plot(1:length(pheno_names), col = colours)
names(colours) = sample(unique(taxa))
#motch = match(taxa, names(colours))
#colours = colours[motch]


#my_strsplit = function(mot) {
#	splot = strsplit(mot, '_')[[1]]
#	shorty = paste(splot[1], splot[2], sep='_')
#}
#short_pheno_names = unlist(lapply(pheno_names, my_strsplit))

my_f = function(k) {

	measure = pheno_names[k]
	taxon = taxa[k]
	if (grepl('ASV', measure)) pvalues_dir_DGE = '/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/' else
		pvalues_dir_DGE = '/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
	
	print(measure)
	DGE_h5 = h5read(paste(pvalues_dir_DGE,'/',measure,'.h5',sep=''),'/')
	
	my_pvalues = c()
	my_cumposs = c()
	my_chrs = c()
	my_poss = c()
	my_betas = c()
	for (chr in 1:20) {
		if (!paste('pvalues_chr',chr,sep='') %in% names(DGE_h5)) stop()
		my_pvalues = c(my_pvalues, DGE_h5[[paste('pvalues_chr',chr,sep='')]])
		my_chrs = c(my_chrs, DGE_h5[[paste('chrs_chr',chr,sep='')]])
		my_poss = c(my_poss, DGE_h5[[paste('poss_chr',chr,sep='')]])
		my_betas = c(my_betas, DGE_h5[[paste('betas_chr',chr,sep='')]])
		my_cumposs = c(my_cumposs, cumpos[[paste('cumpos_chr',chr,sep='')]])
		if (length(my_pvalues) != length(my_cumposs)) stop('pb')
	}
	DGE_h5$chr = my_chrs
	DGE_h5$pvalues = my_pvalues
	DGE_h5$cumpos = my_cumposs
	DGE_h5$pos = my_poss
	DGE_h5$betas = my_betas
	DGE_h5$trait1 = measure
	
	DGE_h5 = data.frame(DGE_h5[['chr']],DGE_h5[['pos']],DGE_h5[['cumpos']],DGE_h5[['pvalues']],DGE_h5[['betas']],DGE_h5[['trait1']])
	colnames(DGE_h5) = c('chr','pos','cumpos','pvalues','betas', 'trait1')
	DGE_h5$logP=-log10(DGE_h5[,'pvalues'])
	
	DGE_h5$col = 'darkgrey'
	w = which(DGE_QTLs[,'measure'] == measure)
	for (i in w) {
		quels = which(DGE_h5[,'chr'] == DGE_QTLs[i,'chr'] & DGE_h5[,'pos'] == DGE_QTLs[i,'pos'])
#		if (length(quels) != 1) print(k)
		DGE_h5[quels,'col'] = colours[taxon] 
	}
#	for (chr in 1:20) {
#		this_min_logP = min(DGE_h5[DGE_h5$col != 'darkgrey','logP'])
#		if (any)
#	}

#	if (any(DGE_h5[,'logP']> 7 & DGE_h5[,'chr'] == 1) ) print(k)
	ret = DGE_h5[DGE_h5[,'logP']>=3,]
	return(ret)
}
res = mclapply(1:length(pheno_names), my_f, mc.cores = 14)

res = do.call('rbind',res)
#all_DGE_h5[all_DGE_h5[,'logP']>10,'logP'] = 10
#w = which(all_DGE_h5[,'col'] == 'black')
#if (length(w)!=0) all_DGE_h5[w,'col'] = 'darkgrey'

source('/users/abaud/abaud/P50_HSrats/code/variance_decomposition/felipes_deblur/annotate_VCs_pvalues_function.R')
res_sigs = res[res$logP > 5.8,]
uniqs = unique(res_sigs$trait1)
motch = match(uniqs, res_sigs$trait1)
res_sigs = res_sigs[motch,]
res_sigs = annotate(res_sigs)
asvs = grep('ASV', res_sigs$trait1)
res_sigs[asvs,'full_taxon'] = paste(res_sigs[asvs,'full_taxon'], res_sigs[asvs,'taxon1'], sep=';')
motch = match(res$trait1, res_sigs$trait1)
res$full_taxon = res_sigs[motch,'full_taxon']

res <- res[order(res$logP, decreasing = T),]

#jpeg(paste('/users/abaud/abaud/P50_HSrats/plots/porcupine_',taxon,'.jpg',sep=''), width=1300,height = 500)
jpeg(paste('/users/abaud/abaud/P50_HSrats/plots/porcupine_uncollapsed_genus2.jpg',sep=''), width=1300,height = 600)
par(mar = c(3,6,2,1))
draw_manhattan_plot(res)
graphics.off()



#### WRONG! no matching on position, just names
library(VennDiagram)
 
# Generate 3 sets of 200 words
set1 <- DGE_QTLs[DGE_QTLs$study1 == 'MI','genus']
set2 <- DGE_QTLs[DGE_QTLs$study1 == 'NY','genus']
set3 <- DGE_QTLs[DGE_QTLs$study1 == 'TN_behavior','genus']
set4 <- DGE_QTLs[DGE_QTLs$study1 == 'TN_breeding','genus']

# Chart
venn.diagram(
  x = list(set1, set2, set3,set4),
  category.names = c("MI" , "NY" , "TN_behavior",'TN_breeding'),
  filename = '/users/abaud/abaud/P50_HSrats/plots/venn_genus.png',
  output=TRUE
)





pvalues_dir='/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/'
target_loci = c('1:196217481','4:70834123','10:101974959')
all_ld = NULL
for (target_locus in target_loci) {
	splot = strsplit(target_locus,':')[[1]]
	target_chr = splot[1]
	target_pos = splot[2]
	load(paste0(pvalues_dir,'LD_',target_chr, '_', target_pos,'.RData'))
	ld$top = target_locus
	all_ld = rbind(all_ld, ld)
}
all_ld = all_ld[all_ld$R2 >= 0.80,]
#all_ld = all_ld[all_ld$R2 < 0.95 & all_ld$CHR_A == 1,]

# !!!! needs to be run in same iteration as real plot otherwise colours don't match
pdf(paste('/users/abaud/abaud/P50_HSrats/plots/porcupine_uncollapsed_legend.pdf',sep=''), width = 15)
draw_legend(res)
dev.off()


cor.test(res[[1]][['betas']],res[[2]][['betas']])
cor.test(res[[1]][['logP']],res[[2]][['logP']])











