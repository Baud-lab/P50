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
# The three phenotypes most significantly affected by Mic-IGE were ASVs in the Muribaculaceae family
sel = all_VCs_full[order(all_VCs_full$pvalue_DGE, decreasing = F),][1:3,]
selDGE= asv_VCs[asv_VCs$trait1 %in% sel$trait1,]
  
# total herit is ... times greater than classical heritability across these 3 phenotypes
sel$total_heritability / selDGE$prop_Ad1 # 

#sel$tot_her_IGE = sel$prop_Ad1 + 2*(2-1)*sel$corr_Ad1s1*sqrt(sel$prop_Ad1*sel$prop_As1) + (2-1)^2*sel$prop_As1

### need traits on cols
### component on rows
## toplot = matrix(NA, ncol = nrow(sel), nrow = 3)
## colnames(toplot) = sel$trait1
## rownames(toplot) = c("DGE", "DGE_IGE", "IGE")
## 
## toplot["DGE",] = sel$prop_Ad1
## toplot["DGE_IGE",] = 2*(2-1)*sel$corr_Ad1s1*sqrt(sel$prop_Ad1*sel$prop_As1)
## toplot["IGE",]= (2-1)^2*sel$prop_As1
   
toplot = matrix(NA, ncol = 2*nrow(sel), nrow = 4)
#toplot = matrix(NA, ncol = 2*nrow(sel), nrow = 6)
colnames(toplot) = c(sel$trait1, paste0(sel$trait1, "_classicHerit"))
rownames(toplot) = c("Mic-DGE", "cov(Mic-DGE,Mic-IGE)", "Mic-IGE", "Mic-DGE alone")
#rownames(toplot) = c("Mic-DGE","cov(Mic-DGE,Mic-IGE)", "Mic-IGE", "Mic-DGE alone", "STE-DGE", "STE-DGE alone")
toplot["Mic-DGE",1:nrow(sel)] = sel$prop_Ad1
toplot["cov(Mic-DGE,Mic-IGE)",1:nrow(sel)] = 2*(2-1)*sel$corr_Ad1s1*sqrt(sel$prop_Ad1*sel$prop_As1)
toplot["Mic-IGE",1:nrow(sel)]= (2-1)^2*sel$prop_As1
toplot[c("Mic-DGE alone"),1:nrow(sel)] = 0

toplot["Mic-DGE alone", (nrow(sel)+1) : (2*nrow(sel))] = selDGE$prop_Ad1
toplot[c("Mic-DGE", "cov(Mic-DGE,Mic-IGE)", "Mic-IGE"),(nrow(sel)+1) : (2*nrow(sel))] = 0

# Saving STE
ste = matrix(NA, ncol = 2*nrow(sel), nrow = 2)
colnames(ste) = c(sel$trait1, paste0(sel$trait1, "_classicHerit"))
rownames(ste) = c("STE-DGE", "STE-DGE alone")
ste["STE-DGE",1:nrow(sel)] = sel$STE_Ad1
ste["STE-DGE alone",1:nrow(sel)] = 0
ste["STE-DGE alone", (nrow(sel)+1) : (2*nrow(sel))] = selDGE$STE_Ad1
ste["STE-DGE", (nrow(sel)+1) : (2*nrow(sel))] = 0

#ord = order(colnames(toplot)) # order alfabetically, lose info about most sign
ord = c(sapply(sel$trait1, function(t) grep(t, colnames(toplot), value = T)))
toplot = toplot[,ord]
ste = ste[,ord] 
#barplot(toplot[,grep(sel$trait1[1], colnames(toplot))], col=colors()[c(23,12, 25)], )
#barplot(toplot, col=colors()[c(23,12, 25)], )
#"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
#coolors = c("#8491B4FF", "#A4036FFF", "#F39B7FFF", "#3C5488FF")#,"") #"#610345" "#C59FC9"
coolors = c("#8491B4FF","#91D1C2FF","#F39B7FFF","#3C5488FF")

#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, length=0.05,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_barplot.pdf", h = 6, w = 8)
par(mar=c(5.1,5.1,2.5,3.5))
#t = sel$trait1[1]
bars <- function(trait1, space=0, add=F, ...){
  traitplot = toplot[,grep(trait1, colnames(toplot))] 
  bp = barplot(traitplot, col=coolors, 
          space=c(space, 0.1), xlim=c(0, ncol(toplot)+2), ylim=c(0,max(apply(toplot, 2,sum)+0.04 )),
          xaxt="n",
          border = NA, cex.axis = 1.25,
          las=1, add=add, ...)
  segments(bp[1], max(apply(traitplot, 2,sum))+0.01 , bp[2], max(apply(traitplot, 2,sum))+0.01)
          #angle=90, code=3, length=0.05) # this cna be used in case want to use arrows instead of segments
  text(mean(c(bp[1], bp[2])), max(apply(traitplot, 2,sum)) + 0.015, "**", cex=1.25, font=2)
  #error.bar(bp[1], toplot["Mic-DGE",grep(paste0(t,"$"), colnames(toplot))], 
  #          ste["STE-DGE",grep(paste0(t,"$"), colnames(toplot))])
  #error.bar(bp[2], toplot["Mic-DGE alone",grep(paste0(t,"_classicHerit"), colnames(toplot))], 
  #          ste["STE-DGE alone",grep(paste0(t,"_classicHerit"), colnames(toplot))])
  #bp1 = bp; rm(bp); rm(t)
  return(bp)
}
bp1 = bars(sel$trait1[1])
bp2 = bars(sel$trait1[2], space = 2+1, add=T, axes = F)
bp3 = bars(sel$trait1[3], space = 3+2+1, add=T, axes = F)

#labx = unlist(lapply(strsplit(sel$trait1, "_"), function(x) paste(x[1:2], collapse = " ")))
labx = gsub("_all", "",sel$trait1)
axis(1, at=c(mean(bp1), mean(bp2), mean(bp3)), labels = labx, tck=F, lwd=0, cex.axis = 1.25)
title(ylab="total genetic variance", cex.lab=1.8, line=3.5)
legend("topright", fill = coolors, 
       border = NA , bty="o", legend = c(rownames(toplot)), cex=1.2, xpd=T, #,
       inset = c(-0.1,-0.1))

dev.off()



###########
# checking if all DGE alone < DGE when IGE
rownames(all_VCs_full) = all_VCs_full$trait1
rownames(asv_VCs) = asv_VCs$trait1
#sum(rownames(asv_VCs) %in% rownames(all_VCs_full))
#sum(rownames(all_VCs_full) %in% rownames(asv_VCs))

all_VCs_full = all_VCs_full[rownames(all_VCs_full) %in% rownames(asv_VCs),] # 1st element in asv_VCs corresponds to 283 el in all_VCs_full
asv_VCs = asv_VCs[rownames(all_VCs_full),]
all(rownames(all_VCs_full) == rownames(asv_VCs))

sum(asv_VCs$prop_Ad1 < all_VCs_full$prop_Ad1)

plotdots <- function(x,y,xylim=T,line="xy",...){
  pch=16
  if(xylim){ axlim = range(c(x,y))}else{axlim = NULL}
  #print(axlim)}
  plot(x, y,
       xlim = axlim, ylim = axlim,type="n",
       ...)
  if(line=="xy") abline(b=1,a=0, col="grey", lty=2)
  if(line=="cor") abline(lm(y~x), col="grey", lty=2)
  points(x,y, col=dot.col,pch=pch)
  
  legend("topleft", pch=pch, col=coolors, 
         legend = lgd)
  if(line=="cor") mtext(paste0("cor(x,y) = ",round(cor(x,y), 2)), line=-1)
}

coolors=c("red","orange")
thrsh = c(0.01,0.05)

dot.col = rep("grey40", nrow(all_VCs_full))
dot.col[all_VCs_full$prop_As1 > thrsh[1]] = coolors[1]
dot.col[all_VCs_full$prop_As1 > thrsh[1] & all_VCs_full$pvalue_DGE < thrsh[2]] = coolors[2]

lgd = c(paste0("IGE > ",thrsh[1]), paste0("& p-val < ", thrsh[2]))

pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_vs_DGE.pdf")
plotdots(all_VCs_full$prop_Ad1, asv_VCs$prop_Ad1, 
         xlab="DGE(with IGE)", ylab = "DGE(alone)")
plotdots(all_VCs_full$total_heritability, asv_VCs$prop_Ad1, 
         xlab="tot_herit(with IGE)", ylab = "DGE(alone)")
dev.off()


pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_vs_DGE_corcol.pdf")
coolors=c("red","orange","#4DBBD5FF")
thrsh = c(0.01,0.05,0.9)
dot.col[all_VCs_full$prop_As1 > thrsh[1] & all_VCs_full$pvalue_DGE < thrsh[2] & all_VCs_full$corr_Ad1s1 > thrsh[3]] = coolors[3]
lgd = c(paste0("IGE > ",thrsh[1]), paste0("& p-val < ", thrsh[2]) , paste0("& cor(DGE,IGE) > ", thrsh[3]) )

plotdots(all_VCs_full$prop_Ad1, asv_VCs$prop_Ad1, 
         xlab="DGE(with IGE)", ylab = "DGE(alone)")
plotdots(all_VCs_full$total_heritability, asv_VCs$prop_Ad1, 
         xlab="tot_herit(with IGE)", ylab = "DGE(alone)")
dev.off()

##-> could you check whether the largest differences are for the most positive corr_Ad1s1?
## plot diff (asv_VCs$prop_Ad1 - all_VCs_full$prop_Ad1) or (asv_VCs$prop_Ad1 - all_VCs_full$total_heritability) VS all_VCs_full$corr_Ad1s1
pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_delta_DGE_vs_corr_corcol.pdf")

plotdots(all_VCs_full$corr_Ad1s1, (all_VCs_full$prop_Ad1 - asv_VCs$prop_Ad1 ),
         xlab = "cor(IGE,DGE)", ylab="DGE(with IGE) - DGE(alone)", 
         xylim=F, line="cor")
#cor.test(abs(all_VCs_full$corr_Ad1s1), (all_VCs_full$prop_Ad1 - asv_VCs$prop_Ad1 ))
#lm((all_VCs_full$prop_Ad1 - asv_VCs$prop_Ad1) ~ all_VCs_full$corr_Ad1s1)

plotdots(all_VCs_full$corr_Ad1s1, (all_VCs_full$total_heritability - asv_VCs$prop_Ad1 ),
         xlab = "cor(IGE,DGE)", ylab="tot_herit(with IGE) - DGE(alone)", 
         xylim=F, line="cor")
dev.off()


######### 
boxdots = function(y,...){
  set.seed(12)
  pch=16
  boxplot(y, 
          outline = F, 
          ylim=range(y),
          boxwex=0.8,
          col=alpha("white",0), 
          border=alpha("white",0),...)
          #xlab = "DGE(with IGE) - DGE(alone)")
  points(jitter(rep(1, nrow(all_VCs_full)), factor=7), 
         y, 
         ylab='', col = dot.col,
         pch=pch, cex=1.2)
  boxplot(y, 
          outline = F, 
          boxwex=0.8,
          col=alpha("white",0.5), 
          #border=alpha("white",0),
          yaxt="n",
          add=T)
}

coolors=c("red","orange","#4DBBD5FF")
thrsh = c(0.01,0.05,0.9)
dot.col = rep("grey40", nrow(all_VCs_full))
dot.col[all_VCs_full$prop_As1 > thrsh[1]] = coolors[1]
dot.col[all_VCs_full$prop_As1 > thrsh[1] & all_VCs_full$pvalue_DGE < thrsh[2]] = coolors[2]
dot.col[all_VCs_full$prop_As1 > thrsh[1] & all_VCs_full$pvalue_DGE < thrsh[2] & all_VCs_full$corr_Ad1s1 > thrsh[3]] = coolors[3]

lgd = c(paste0("IGE > ",thrsh[1]), paste0("& p-val < ", thrsh[2]) , paste0("& cor(DGE,IGE) > ", thrsh[3]))
pch=16

pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_delta_DGE_vs_corr_bplot.pdf")
par(mfrow=c(1,2))
boxdots(all_VCs_full$prop_Ad1 - asv_VCs$prop_Ad1, ylab="DGE(with IGE) - DGE(alone)")
legend("topleft", pch=pch, col=coolors, 
       legend = lgd, xpd=T, inset=c(0,-0.12), cex = 0.8, bg="white")
boxdots(all_VCs_full$total_heritability - asv_VCs$prop_Ad1, ylab="tot_herit(with IGE) - DGE(alone)")
dev.off()
