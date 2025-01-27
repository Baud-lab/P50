library("scales") # for alpha
library('beeswarm') # for dots

#load genetic correlations
load('/users/abaud/abaud/P50_HSrats/output/VD/bivariate/all_VCs_corr_Ad1d2_zero_P50_Rn7_pruned_DGE.RData')

#create a variable indicating the pair of cohorts considered
#studies = c('MI','NY','TN_behavior','TN_breeder')
#studies = c('TN_breeder','TN_behavior','MI','NY') #c('MI','NY','TN_behavior','TN_breeder')

# this could be done with apply
for (i in 1:dim(all_VCs)[1]) {
  all_VCs[i,'study_pair'] = paste(sort(c(all_VCs[i,'study1'],all_VCs[i,'study2'])),collapse='\n')
}

# ordering so that align to the other plots in the figure
all_VCs[,'study_pair'] = factor(all_VCs$study_pair, 
                                levels = c("MI\nTN_breeder","MI\nTN_behavior","MI\nNY","NY\nTN_breeder","NY\nTN_behavior","TN_behavior\nTN_breeder"))


#calculcate median gen corr for each pair of cohorts to order the boxplot in increasing order
#medians = sort(sapply(split(all_VCs$corr_Ad1d2, all_VCs$study_pair), median), decreasing = T)

#plot
#all_VCs$study_pair = factor(all_VCs$study_pair, levels = names(medians))

pdf("/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/comp_gen_corrs_across_cohorts.pdf", h= 6, w = 10.5)
par(mar=c(5.1,5.1,2.1,2.1))
#vio_col = alpha("#91D1C2FF", 0.3) #c("lightblue")
#vio_col = alpha("#4DBBD5FF", 0.3) #c("lightblue")
#vio_col=alpha("#F39B7FFF", 0.3)
#show_col(vio_col)
box_col = "white"
#bee_col = "#4DBBD5FF" # "#396C93"
bee_col = "#396C93" # Second colour in bar plot

formula=as.formula(all_VCs$corr_Ad1d2 ~ all_VCs$study_pair)
#vioplot(formula, col = vio_col,#col = rep(pal, 3),
#        #at = ats, 
#        #names = rep('', length(ats)), 
#        las = 1, xlab = '', xaxt = "n", ylab='',
#        cex.axis = 1.25,
#        border = "black", 
#        drawRect=F,
#        wex=1,
#        areaEqual = F
#        #ylim = c(0,25)
#);
boxplot(formula, col = box_col, #col="white",
           boxwex = 0.5,
           xaxt = "n", #yaxt="n", box="n",
           xlab="", ylab="", cex.axis = 1.25,
           outline = F, 
           #border = NULL,
           #range = 0.00000000001, 
           #at = ats, #ylim = c(0,25), 
           #names = rep('', length(ats)), 
           las = 1,  
           #add=T,
           whisklty=3) #,
           #frame=F)
beeswarm(formula, add=T, col=bee_col, pch=16, cex=1.2)

axis(1, at = 1:6, labels = levels(all_VCs$study_pair), cex.axis = 1.2, padj=0.5)
title(xlab = 'Cohort pairs', cex.lab = 1.4, line=3.8)
title(ylab = "Same-taxon genetic correlations", cex.lab = 1.4, line=3.5)
dev.off()


#beeswarm(formula, add=T)



########## Amelie ######

## #pdf('/users/abaud/abaud/P50_HSrats/plots/comp_gen_corrs_across_cohorts.pdf', width = 20, height = 10)
## #par(mar = c(12,10,3,1))
## par(mar=c(5.1,5.1,2.1,2.1))
## 
## boxplot(all_VCs$corr_Ad1d2 ~ all_VCs$study_pair, ylab = '', xlab = '',  boxwex=0.5, axes = FALSE)
## axis(1, at = 1:6, labels = levels(all_VCs$study_pair), cex.axis = 2, padj = 1)
## axis(2, padj = -1 , cex.axis = 3, las = 1)
## mtext(side = 1, at = 3.5, text = 'Cohort pairs', padj = 4.2, cex = 3)
## mtext(side = 2, at = 0.25, text = 'Same-taxon genetic correlations', padj = -3.5, cex = 3)
## 
## #dev.off()

