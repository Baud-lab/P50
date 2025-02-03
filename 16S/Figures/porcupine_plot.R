#library(RColorBrewer)
library(rhdf5)
library(parallel)
library(paletteer)

##### Functions to plot manhattan
draw_manhattan_plot=function(dat, def.cex = 0.4, sig.cex = 2) {
	cols <- c("grey70", "grey55") #c( "grey70", "grey40") # c("#999999", "#808080")
	#cols <-"lightgrey"
	#dat = res[1:1000,]
	#dat = res[res$chr %in% c(1,2),]
	# default color for chromosomes (odd numbers)
	dat$colour=cols[1]
	# chromosome with even number are of the other color
	dat[which(dat[,'chr'] %in% c(2,4,6,8,10,12,14,16,18,20)),'colour']=cols[2]
	# default dimension and type of the dot, to change with dimensions of output image
	dat$cex = def.cex
	dat$pch = 16
	
	# selecting rows that have a different color than the default - significant dots
	w = which(dat[,'col'] != 'darkgrey')
	sig.pch = c("NY" = 15, "MI" = 16, "TN_behavior" = 17, "TN_breeder" = 18)
	
	if (length(w)!=0) {
		dat[w,'cex'] = sig.cex
		dat[w,'colour'] = dat[w,'col']
		dat[w,'pch'] = sig.pch[dat[w,"study1"]]
		dat[dat[,'pch'] == 18, "cex"] = sig.cex + 0.3
	}
	
	#need both order first by decreasing logP then by colour otherwise QTL point gets overpainted so is not visible
	dat = dat[order(dat$colour != cols[1]),]
	r = range(c(0,range(dat$logP)))
	par(mgp=c(0,0,0))
	plot(dat[,c('cumpos','logP')], ylab="", col=dat[,'colour'], 
	     xlab="", main='', cex=dat[,'cex'],
	     pch=dat[,"pch"], axes=FALSE, #xaxs="i", 
	     )
	box()
	axis(2,cex.axis=2,las = 1)

	mtext("-logP",2,padj=-3,cex=2)

	mtext(c(1:20), 1, 
	      at=tapply(dat$cumpos,dat$chr,mean),#median),#+40000,
	      adj=1,las=1,cex=1.5)
	mtext('Chromosome',1,at=median(dat$cumpos),padj=2,las=1,cex=1.5)
}

###### Function to plot legend 
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

######### Function to get a palette of colours?
#map2color<-function(x,pal,limits=NULL){
#    if(is.null(limits)) limits=range(x)
#    pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
#}
#colours = sample(map2color(1:length(pheno_names),rainbow(length(pheno_names)),limits=NULL))
##colours = sample(map2color(1:length(unique(taxa)),rainbow(length(unique(taxa))),limits=NULL))

######### Load res to plot
#load("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/QTLs_alpha1e-04_unpruned_DGE_CE_MaE_toPlot_10pheno.RData")
load("/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/QTLs_alpha1e-04_unpruned_DGE_CE_MaE_toPlot.RData")

######## Annotate peaks?
#source('/users/abaud/abaud/P50_HSrats/code/variance_decomposition/felipes_deblur/annotate_VCs_pvalues_function.R') # annotate() function 
source('/users/abaud/htonnele/git/lab/P50/16S/Figures/annotate_VCs_pvalues_function.R') # annotate() function 
# selecting significant res
res_sigs = res[res$logP > 5.8,] 
# defining which traits are significant 
uniqs = unique(res_sigs$trait1)
# selecting significant traits
motch = match(uniqs, res_sigs$trait1)
res_sigs = res_sigs[motch,] # dim(res_sigs) 111 8
# annotating significant traits with full taxon name
res_sigs = annotate(res_sigs)
# adding ASV n to taxon name
asvs = grep('ASV', res_sigs$trait1)
res_sigs[asvs,'full_taxon'] = paste(res_sigs[asvs,'full_taxon'], res_sigs[asvs,'taxon1'], sep=';')
motch = match(res$trait1, res_sigs$trait1)
res$full_taxon = res_sigs[motch,'full_taxon']
res$study1 = res_sigs[motch,'study1']
#length(unique(res$full_taxon)); # this is 97
#length(unique(res_sigs$full_taxon)); # this is 97 too
#length(uniqs) # this is 111
res <- res[order(res$logP, decreasing = T),]


###### Setting colors
colres = unique(res$col)

#library(paletteer) # moved above
#show_col(c(paletteer_d("colorBlindness::SteppedSequential5Steps"), paletteer_d("dichromat::BluetoDarkOrange_18")[14:18])) # 30 colours
coolors = as.vector(c(paletteer_d("colorBlindness::SteppedSequential5Steps"), 
            paletteer_d("dichromat::BluetoDarkOrange_18")[14:18], 
            paletteer_d("RColorBrewer::PuRd")[4:9],
            paletteer_c("grDevices::ag_GrnYl", 15))) # 42
#coolors = paletteer_d("colorBlindness::SteppedSequential5Steps")[1:length(colres)]

#coolors = as.vector(paletteer_c("grDevices::Viridis", length(colres[-1]))) # -1 because of "darkgrey"
#scales::show_col(coolors)
set.seed(1); names(coolors) = sample(colres[-which(colres == "darkgrey")])
coolors = c(coolors, "darkgrey" = "darkgrey")
scales::show_col(coolors)
#coolors[which(!names(coolors) %in% res_sigs$col)] = "black" # To test something # all the 51 colors are supposedly used

res$col = unname(coolors[res$col])
scales::show_col(unique(res$col))
#res_sigs$col = unname(coolors[res_sigs$col])
#which( ! unique(res$col) %in% unique(res_sigs$col)) 

#### TODO: Check if use all colors or only the sig ones - if only sig ones can choose better colors 
#          Can check by putting "black" to all the colors that do not appear in sig ones

# TODO: Have to add - in manhattan function - the pch different depending on the cohort - study - found in res_sig
subres = res[res$logP > 4,]

#jpeg(paste('/users/abaud/abaud/P50_HSrats/plots/porcupine_',taxon,'.jpg',sep=''), width=1300,height = 500)
#jpeg(paste('/users/abaud/abaud/P50_HSrats/plots/porcupine_uncollapsed_genus2.jpg',sep=''), width=1300,height = 600)
#jpeg(paste("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/porcupine_uncollapsed_genus2_HT.jpg",sep=''), width=1300,height = 600)
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/porcupine_uncollapsed_genus2_HT.pdf", h=8, w=20)
par(mar = c(3,6,2,1))
draw_manhattan_plot(subres, def.cex=0.6)
#graphics.off()
dev.off()


#### WRONG! no matching on position, just names
## library(VennDiagram)
##  
## # Generate 3 sets of 200 words
## set1 <- DGE_QTLs[DGE_QTLs$study1 == 'MI','genus']
## set2 <- DGE_QTLs[DGE_QTLs$study1 == 'NY','genus']
## set3 <- DGE_QTLs[DGE_QTLs$study1 == 'TN_behavior','genus']
## set4 <- DGE_QTLs[DGE_QTLs$study1 == 'TN_breeding','genus']
## 
## # Chart
## venn.diagram(
##   x = list(set1, set2, set3,set4),
##   category.names = c("MI" , "NY" , "TN_behavior",'TN_breeding'),
##   filename = '/users/abaud/abaud/P50_HSrats/plots/venn_genus.png',
##   output=TRUE
## )




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
#pdf(paste('/users/abaud/abaud/P50_HSrats/plots/porcupine_uncollapsed_legend.pdf',sep=''), width = 15)
draw_legend(res)
#dev.off()


cor.test(res[[1]][['betas']],res[[2]][['betas']])
cor.test(res[[1]][['logP']],res[[2]][['logP']])











