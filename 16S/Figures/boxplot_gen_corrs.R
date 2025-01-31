library("scales") # for alpha
library('beeswarm') # for dots

#load genetic correlations
load('/users/abaud/abaud/P50_HSrats/output/VD/bivariate/all_VCs_corr_Ad1d2_zero_P50_Rn7_pruned_DGE.RData')

#create a variable indicating the pair of cohorts considered
#studies = c('MI','NY','TN_behavior','TN_breeder')
#studies = c('TN_breeder','TN_behavior','MI','NY') #c('MI','NY','TN_behavior','TN_breeder')

dict = c('NY'= 'NY',
         'MI'= 'MI',
         'TN_behavior'= 'TN1',
         'TN_breeder'= 'TN2')

all_VCs[,"study1"] = sapply(all_VCs[,"study1"], function(x) unname(dict[x]))
all_VCs[,"study2"] = sapply(all_VCs[,"study2"], function(x) unname(dict[x]))

#for (i in 1:dim(all_VCs)[1]) {
#  all_VCs[i,'study_pair_for'] = paste(sort(c(all_VCs[i,'study1'],all_VCs[i,'study2'])),collapse='\n')
#}

# Done with apply
all_VCs[,'study_pair'] = apply(all_VCs, 1, function(x) paste(sort(x[c("study1","study2")]),collapse="\n"))
#all(all_VCs[,'study_pair'] == all_VCs[,'study_pair_for'])


# ordering so that decreasing mean of sample sizes
all_VCs[,'study_pair'] = factor(all_VCs$study_pair, 
                                #levels = c("MI\nTN_breeder","MI\nTN_behavior","MI\nNY","NY\nTN_breeder","NY\nTN_behavior","TN_behavior\nTN_breeder"))
                                levels = c("MI\nNY", "NY\nTN1", "MI\nTN1", "NY\nTN2", "MI\nTN2","TN1\nTN2"))

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




###### Plotting with colour depending on p-values
## Creating fake p-values in while waiting for true ones
all_VCs$logP = runif(nrow(all_VCs), 0, 5)
a = -log10(0.05)
sig = which(all_VCs$logP < a)
#nonsig = which(all_VCs$logP > a | all_VCs$logP == a)

## NB: it works as long as correlations with similar magnitude are similarly sig/non-sig;
#      otherwise there's a bit of overlap of dots, which might be annoying
#### 
boxplot(all_VCs$corr_Ad1d2 ~ all_VCs$study_pair, 
        col = box_col, #col="white",
        boxwex = 0.5,
        xaxt = "n", #yaxt="n", box="n",
        xlab="", ylab="", cex.axis = 1.25,
        outline = F, 
        las = 1,  
        whisklty=3) #,
axis(1, at = 1:6, labels = levels(all_VCs$study_pair), cex.axis = 1.2, padj=0.5)
title(xlab = 'Cohort pairs', cex.lab = 1.4, line=3.8)
title(ylab = "Same-taxon genetic correlations", cex.lab = 1.4, line=3.5)

bee_col = "#396C93"
beeswarm(all_VCs$corr_Ad1d2[-sig] ~ all_VCs$study_pair[-sig],
         add=T, 
         col=bee_col, bg=alpha(bee_col, 0.5), 
         pch= 21, 
         cex=1.2,
         method = "swarm") # "center", "hex", "square"

bee_col_sig = "#E64B35FF"
beeswarm(all_VCs$corr_Ad1d2[sig] ~ all_VCs$study_pair[sig], 
         add=T, 
         col=bee_col_sig, bg=alpha(bee_col_sig, 0.5), 
         pch=21, 
         cex=1.2,
         method="swarm")



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

