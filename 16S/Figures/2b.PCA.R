load('collapsed_full_biomt_collapsed_clr_counts.RData')
# collapsed_clr_counts is post CLR post collapsing to higher level taxa
save_collapsed_clr_counts = collapsed_clr_counts

# Loading metadata
load('metadata_16Spaper.RData')

dict = c("NY"="NY", "MI"="MI", "TN_behavior"="TN1", "TN_breeder"="TN2")
metadata$study = unname(dict[metadata[,"study"]])
coolors= c(NY="#0073C2FF",MI="#1EA896FF", TN1="#CD534CFF", TN2="#EFC000FF")

tax_level = "f__" # or any of "p__","c__","o__","f__","g__"

# or can use a for loop to plot all:
#for (tax_level in c("f__", "p__","c__","o__","f__","g__")) {
  
collapsed_clr_counts = save_collapsed_clr_counts[grep(tax_level, rownames(save_collapsed_clr_counts)),]
pca = prcomp(t(collapsed_clr_counts))
vars = round(pca$sdev^2 / (sum(pca$sdev^2))*100, digits = 0)

motch = match(colnames(collapsed_clr_counts), metadata$deblur_rooname)
any(is.na(motch))
#FALSE
meta_oi = metadata[motch,]

cols = rep('white', dim(collapsed_clr_counts)[2])
for (s in names(coolors)){
  cols[meta_oi$study == s] = coolors[s] 
}

# Saving objects used to plot
#save(pca, vars, tax_level, cols, coolors, file = paste0("source_files/fig2b__",gsub("__","",tax_level),".RData"))

# Loading objects used to plot
#load("source_files/fig2b__f.RData")

pdf(paste0("PCA__",gsub("__","",tax_level),"_paper.pdf"), h=7, w =8)

par(mar = c(5.1,5.1,2.1,2.1)) #5.1))

plot(pca$x[,1], pca$x[,2], 
     xlab = paste('MDS1 (',vars[1],'%)', sep=''), 
     ylab = "", #paste('MDS2 (',vars[2],'%)', sep=''), 
     col = cols, pch = 16, las = 1, 
     cex = 0.85, 
     main = tax_level,
     cex.axis = 1.2, 
     cex.lab=1.4) 
title(ylab = paste('MDS2 (',vars[2],'%)', sep=''), cex.lab = 1.4, line=3.5)

lgpos = "topleft"
if(tax_level == "f__"){ # because have points plotting there for "f__"
  lgpos = "bottomright"
}
legend(lgpos, 
       fill = coolors, 
       legend = names(coolors), 
       border = NA, cex=1.1, 
       bty = "n", 
       x.intersp = 0.6 
)

dev.off()

#} # uncomment if using for loop
