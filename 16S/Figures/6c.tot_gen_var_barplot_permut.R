# Loading VD data from model with DGE only
load('augmented_DGE_VC_wALL.RData')
DGEonly_VCs = all_VCs_full
rm(all_VCs_full)

# Now loading VD data from model with DGE and IGE
load('augmented_IGE_VC.RData')
# one bar with three different colors for:
## 1. all_VCs_full$prop_Ad1 
## 2. 2*(2-1)*all_VCs_full$corr_Ad1s1*sqrt(all_VCs_full$prop_Ad1*all_VCs_full$prop_As1) 
## 3. (2-1)^2*all_VCs_full$prop_As1

# Now loading VD data from scrambled 
load("scrambled_VCs_full_model4Helene.Rdata") # loading scrambled_VCs_full_model

# for the **three phenotypes most significantly affected by Mic-IGE** (as said in the text), so 6 bars total
sel = all_VCs_full[order(all_VCs_full$pvalue_DGE, decreasing = F),][1:3,]
selDGE= DGEonly_VCs[DGEonly_VCs$trait1 %in% sel$trait1,]
#selDGE= DGEonly_VCs[DGEonly_VCs$trait1 %in% c("ASV_13916_MI", "ASV_18948_MI", "ASV_17551_MI"),]

# total herit is 4.4, 7.35, 5 times greater than classical heritability across these 3 phenotypes
sel$total_heritability / selDGE$prop_Ad1 

# Creating matrix to plot the 3 most significant 
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
# Checking tot_heritability corresponds to sum of single variance components
#sum(toplot[,"ASV_13916_all"]) == sel[sel$trait1 == "ASV_13916_all", "total_heritability"] # TRUE
#sum(toplot[,"ASV_18948_all"]) == sel[sel$trait1 == "ASV_18948_all", "total_heritability"] # TRUE
#sum(toplot[,"ASV_17551_all"]) == sel[sel$trait1 == "ASV_17551_all", "total_heritability"] # TRUE

# Saving objects to plot - toplot (model wt IGE), sel (model with DGE), scrambled_VCs_full_model (permut)
save(toplot, sel, scrambled_VCs_full_model, file = "source_files/fig6c.RData")

# Loading objects to plot
#load("source_files/fig6c.RData")
coolors = c("#8491B4FF","#91D1C2FF","#F39B7FFF","#3C5488FF")

# Function to plot bar plot at different point on x axis
# need toplot and scrambled_VCs_full_model
bars <- function(trait1, space=0, add=F, ...){
  traitplot = toplot[,grep(trait1, colnames(toplot))] 
  # add barplot
  bp = barplot(traitplot, col=coolors, 
               space=c(space, 0.1), 
               xlim=c(0, ncol(toplot)+6), ylim=c(0,max(apply(toplot, 2,sum)+0.04 )),
               xaxt="n",
               border = NA, cex.axis = 1.25,
               las=1, add=add, ...)
  # add segment of significance on top
  segments(bp[1], max(apply(traitplot, 2,sum))+0.01 , bp[2], max(apply(traitplot, 2,sum))+0.01)
  # add stars of significance on top
  text(mean(c(bp[1], bp[2])), max(apply(traitplot, 2,sum)) + 0.015, "**", cex=1.25, font=2)
  
  # add points with results from permutations
  scram_herit = scrambled_VCs_full_model[grep(trait1, scrambled_VCs_full_model$trait1), "total_heritability"]
  factor = (50*0.24/(space+3))
  #print(paste0("space= ", space, "; factor= ", factor))
  set.seed(21)
  points(x = jitter(rep(space+3, length(scram_herit)), factor=factor), #c(1,0.9,0.5)),
         y = scram_herit, 
         pch=16, cex=0.5, 
         col = adjustcolor("grey50",0.5))
  # add boxplot with results from permutations
  bx <- boxplot(scram_herit, 
          at = space+3, add = TRUE, 
          col = adjustcolor("white", 0.5), border = "black", width = 0.8,
          axes = FALSE, outline = F
          )
  
  return(bp)
}


# Open pdf to save plot
pdf("tot_herit_barplot_perm.pdf", h = 6, w = 7)
par(mar=c(5.1,5.1,2.5,3.5))

# Bar plot: full model with IGE, model with DGE only and boxplot with tot heritability from permutations
bp1 = bars(sel$trait1[1])
bp2 = bars(sel$trait1[2], space = 3+1.5, add=T, axes = F)
bp3 = bars(sel$trait1[3], space = 4+3+2, add=T, axes = F)

# add ticks and ticks labels for x axis
labx = gsub("_all", "",sel$trait1)
axis(1, at=c(mean(bp1), mean(bp2), mean(bp3)) + 0.65, labels = labx, tck=F, lwd=0, cex.axis = 1.25)
# add y label
title(ylab="total genetic variance", cex.lab=1.4, line=4)
# add legend
lgd = c(coolors)
names(lgd) = c(rownames(toplot))
lgd = lgd[c(3:1,4)]

l = legend("topright", fill = lgd, 
       border = NA , bty="o", legend = names(lgd), cex=0.8, xpd=T,
       inset = c(-0.1,-0.1)) #,
       #title = "Real data", title.font = 2, title.adj = 0.1)

lgd2 = c("permuted cage mates   " = adjustcolor("grey50", 0.5))
legend(l$rect$left, l$text$y[length(coolors)] - (l$text$y[1] - l$text$y[2]),  
       fill = lgd2, legend = names(lgd2), cex=0.8, xpd=T) #, 
       #title = "Permutations", title.font = 2, title.adj = 0.1)
dev.off()


# Check boxplot
#boxplot(scrambled_VCs_full_model$total_heritability ~ scrambled_VCs_full_model$trait1)
