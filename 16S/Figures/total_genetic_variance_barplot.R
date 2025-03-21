# Loading VD data 
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_IGE_VC.RData'))

colnames(all_VCs_full)
all_VCs_full$total_heritability

all_VCs_full = all_VCs_full[all_VCs_full$study1 == 'all',] # not sure I have to, but will do as it is done in the qqplot

# one bar with three different colors for:
## 1. all_VCs_full$prop_Ad1 
## 2. 2*(2-1)*all_VCs_full$corr_Ad1s1*sqrt(all_VCs_full$prop_Ad1*all_VCs_full$prop_As1) 
## 3. (2-1)^2*all_VCs_full$prop_As1

# next to this bar another bar with only one color corresponding to 
## 1. all_VCs_full$prop_Ad1 from DGE only model
load(file.path("/users/abaud/abaud/P50_HSrats/output/VD/univariate/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect", "all_estNste.Rdata"))
asv_VCs = VCs

# for the **three phenotypes most significantly affected by Mic-IGE** (as said in the text), so 6 bars total
sel = all_VCs_full[order(all_VCs_full$pvalue_DGE, decreasing = F),][1:3,]
selDGE= asv_VCs[asv_VCs$trait1 %in% sel$trait1,]
  
# total herit is 4.4, 7.35, 5 times greater than classical heritability across these 3 phenotypes
sel$total_heritability / selDGE$prop_Ad1 # 


# Creating matrix to plot te 3 most significant 
toplot = matrix(NA, ncol = 2*nrow(sel), nrow = 4)

colnames(toplot) = c(sel$trait1, paste0(sel$trait1, "_classicHerit"))
rownames(toplot) = c("Mic-DGE", "cov(Mic-DGE,Mic-IGE)", "Mic-IGE", "Mic-DGE alone")

toplot["Mic-DGE",1:nrow(sel)] = sel$prop_Ad1
toplot["cov(Mic-DGE,Mic-IGE)",1:nrow(sel)] = 2*(2-1)*sel$corr_Ad1s1*sqrt(sel$prop_Ad1*sel$prop_As1)
toplot["Mic-IGE",1:nrow(sel)]= (2-1)^2*sel$prop_As1
toplot[c("Mic-DGE alone"),1:nrow(sel)] = 0

toplot["Mic-DGE alone", (nrow(sel)+1) : (2*nrow(sel))] = selDGE$prop_Ad1
toplot[c("Mic-DGE", "cov(Mic-DGE,Mic-IGE)", "Mic-IGE"),(nrow(sel)+1) : (2*nrow(sel))] = 0

# Order as trait1, trait1DGE-alone; trait2, trait2DGE-alone  ...
ord = c(sapply(sel$trait1, function(t) grep(t, colnames(toplot), value = T)))
toplot = toplot[,ord]

coolors = c("#8491B4FF","#91D1C2FF","#F39B7FFF","#3C5488FF")

# Function to plot bar plot at different point on x axis
bars <- function(trait1, space=0, add=F, ...){
  traitplot = toplot[,grep(trait1, colnames(toplot))] 
  # add barplot
  bp = barplot(traitplot, col=coolors, 
               space=c(space, 0.1), xlim=c(0, ncol(toplot)+2), ylim=c(0,max(apply(toplot, 2,sum)+0.04 )),
               xaxt="n",
               border = NA, cex.axis = 1.25,
               las=1, add=add, ...)
  # add segment of significance on top
  segments(bp[1], max(apply(traitplot, 2,sum))+0.01 , bp[2], max(apply(traitplot, 2,sum))+0.01)
  # add stars of significance on top
  text(mean(c(bp[1], bp[2])), max(apply(traitplot, 2,sum)) + 0.015, "**", cex=1.25, font=2)
  return(bp)
}

# Open pdf to save plot
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_barplot.pdf", h = 6, w = 7)
par(mar=c(5.1,5.1,2.5,3.5))

# Bar plot
bp1 = bars(sel$trait1[1])
bp2 = bars(sel$trait1[2], space = 2+1, add=T, axes = F)
bp3 = bars(sel$trait1[3], space = 3+2+1, add=T, axes = F)

# add ticks and ticks labels for x axis
labx = gsub("_all", "",sel$trait1)
axis(1, at=c(mean(bp1), mean(bp2), mean(bp3)), labels = labx, tck=F, lwd=0, cex.axis = 1.25)
# add y label
title(ylab="total genetic variance", cex.lab=1.4, line=4)
# add legend
legend("topright", fill = coolors, 
       border = NA , bty="o", legend = c(rownames(toplot)), cex=1.1, xpd=T, #,
       inset = c(-0.1,-0.1))

dev.off()

