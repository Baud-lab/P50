load('/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData')
# collapsed_clr_counts is post CLR post collapsing to higher level taxa
save_collapsed_clr_counts = collapsed_clr_counts

load('/nfs/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData')

dict = c("NY"="NY", "MI"="MI", "TN_behavior"="TN1", "TN_breeder"="TN2")
metadata$study = unname(dict[metadata[,"study"]])
#coolors = c("#E64B35FF", "#4DBBD5FF", "#F39B7FFF", "#00A087FF") #"#F39B7FFF")
#scales::show_col(ggsci::pal_npg("nrc")(10))
#"#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"

#scales::show_col(ggsci::pal_jco()(10))
#"#0073C2FF" "#EFC000FF" "#868686FF" "#CD534CFF" "#7AA6DCFF" "#003C67FF" "#8F7700FF" "#3B3B3BFF" "#A73030FF" "#4A6990FF"
coolors= c(NY="#0073C2FF",MI="#1EA896FF", TN1="#CD534CFF", TN2="#EFC000FF")


#pdf('/nfs/users/abaud/abaud/P50_HSrats/plots/PCA_paper.pdf')
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/PCA_paper.pdf", h=7, w =8)
par(mar = c(5.1,5.1,2.1,2.1)) #5.1))

for (tax_level in c("f__")){ #c("p__","c__","o__","f__","g__")) {  # c("f__")){

	collapsed_clr_counts = save_collapsed_clr_counts[grep(tax_level, rownames(save_collapsed_clr_counts)),]
	pca = prcomp(t(collapsed_clr_counts))
	vars = round(pca$sdev^2 / (sum(pca$sdev^2))*100, digits = 0)

	motch = match(colnames(collapsed_clr_counts), metadata$deblur_rooname)
	any(is.na(motch))
	#FALSE
	meta_oi = metadata[motch,]

	cols = rep('white', dim(collapsed_clr_counts)[2])
	for (s in names(coolors)){
	  cols[meta_oi$study == s] = coolors[s] #'blue'
	}
	#cols[meta_oi$study == 'NY'] = coolors["NY"] #'blue'
	#cols[meta_oi$study == 'MI'] = coolors["MI"] #'green'
	#cols[meta_oi$study == 'TN1'] = coolors["TN1"] #'red'
	#cols[meta_oi$study == 'TN2'] = coolors["TN2"] #'orange'

	plot(pca$x[,1], pca$x[,2], 
	     xlab = paste('MDS1 (',vars[1],'%)', sep=''), 
	     ylab = "",#paste('MDS2 (',vars[2],'%)', sep=''), 
	     col = cols, pch = 16, las = 1, 
	     cex = 0.85, 
	     main = tax_level,
	     cex.axis = 1.2, 
	     cex.lab=1.4) 
	     #(mar = c(5.1,5.1,2.1,2.1))
	title(ylab = paste('MDS2 (',vars[2],'%)', sep=''), cex.lab = 1.4, line=3.5)
	
	lgpos = "topleft"
	if(tax_level == "f__"){ # because have points plotting there for "f__"
	  lgpos = "bottomright"
	}
	legend(lgpos, 
	       fill = coolors, 
	       #grconvertX(1,"npc")
	       #inset = c(-0.2,0),xpd=T,
	       legend = names(coolors), 
	       border = NA, cex=1.1, 
	       bty = "n", #bty ="o", box.col = "black",
	       x.intersp = 0.6 
	)

}
dev.off()

