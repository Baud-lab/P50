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
  
# total herit is 2.6-3.2 times greater than classical heritability across these 3 phenotypes
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
colnames(toplot) = c(sel$trait1, paste0(sel$trait1, "_classicHerit"))
rownames(toplot) = c("Mic-DGE", "cov(Mic-DGE,Mic-IGE)", "Mic-IGE", "Mic-DGE alone")

toplot["Mic-DGE",1:nrow(sel)] = sel$prop_Ad1
toplot["cov(Mic-DGE,Mic-IGE)",1:nrow(sel)] = 2*(2-1)*sel$corr_Ad1s1*sqrt(sel$prop_Ad1*sel$prop_As1)
toplot["Mic-IGE",1:nrow(sel)]= (2-1)^2*sel$prop_As1
toplot["Mic-DGE alone",1:nrow(sel)] = 0

toplot["Mic-DGE alone",(nrow(sel)+1) : (2*nrow(sel))] = selDGE$prop_Ad1
toplot[c("Mic-DGE", "cov(Mic-DGE,Mic-IGE)", "Mic-IGE"),(nrow(sel)+1) : (2*nrow(sel))] = 0

#ord = order(colnames(toplot)) # order alfabetically, lose info about most sign
ord = c(sapply(sel$trait1, function(t) grep(t, colnames(toplot), value = T)))
toplot = toplot[,ord]
# 
#barplot(toplot[,grep(sel$trait1[1], colnames(toplot))], col=colors()[c(23,12, 25)], )
#barplot(toplot, col=colors()[c(23,12, 25)], )
#"#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF" "#B09C85FF"
#coolors = c("#8491B4FF", "#A4036FFF", "#F39B7FFF", "#3C5488FF")#,"") #"#610345" "#C59FC9"
coolors = c("#8491B4FF","#91D1C2FF","#F39B7FFF","#3C5488FF")

pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/tot_herit_barplot.pdf", h = 6, w = 8)
par(mar=c(5.1,5.1,2.5,3.5))
barplot(toplot[,grep(sel$trait1[1], colnames(toplot))], col=coolors, 
        space=0.1, xlim=c(0, ncol(toplot)+2), ylim=c(0,max(apply(toplot[,sel$trait1], 2,sum) )),
        xaxt="n",border = NA,cex.axis = 1.25,
        las=1, )

barplot(toplot[,grep(sel$trait1[2], colnames(toplot))], col=coolors, 
        space=c(2+1, 0.1), xaxt="n", axes = F, border = NA,
        add=T)

barplot(toplot[,grep(sel$trait1[3], colnames(toplot))], col=coolors, 
        space=c(3+2+1, 0.1), xaxt="n", axes = F,border = NA,
        add=T)

#labx = unlist(lapply(strsplit(sel$trait1, "_"), function(x) paste(x[1:2], collapse = " ")))
labx = gsub("_all", "",sel$trait1)
axis(1, at=c(1,4,7), labels = labx, tck=F, lwd=0, cex.axis = 1.25)
title(ylab="total genetic variance", cex.lab=1.8, line=3.5)
legend("topright", fill = coolors, 
       border = NA , bty="o", legend = c(rownames(toplot)), cex=1.2, xpd=T, #,
       inset = c(-0.1,-0.1))

dev.off()




