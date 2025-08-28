#setwd("/users/abaud/htonnele/PRJs/P50_HSrats/16S/data/")
library("plotrix")

load('augmented_DGE_VC_wALL.RData')
DGEonly_VCs = all_VCs_full
rm(all_VCs_full)

# Now loading VD data from model with DGE and IGE
load('augmented_IGE_VC.RData')

# Now loading VD data from scrambled 
load("scrambled_VCs_full_model4Helene.Rdata")

######## plots from here ########
pdf("suppN.VCs_permut_boxplot.pdf", w = 8, h = 5)
# setting parameters 
par(las = 1, 
    pch = 16, 
    cex.axis = 1.25, 
    cex.lab=1.4, 
    cex.main = 1.5)
foctor = 10
bw = 1.3
lmar = 6
rmar = 1
tmar = 2
bmar = 4
pch.cex = 1.2

### For DGE
# Create a custom layout with more control
layout_matrix <- matrix(c(2, 4, 6, 7, 1, 3, 5, 8), nrow = 2, ncol = 4, byrow = TRUE)
layout(layout_matrix, heights = c(1, 5), widths = c(1, 1, 1, 0.7))
#layout.show(6)
traitplot = c("ASV_13916", "ASV_18948", "ASV_17551")

for(t in traitplot){
  set.seed(101)# for jitter
  # select trait to plot
  toplot = scrambled_VCs_full_model[scrambled_VCs_full_model$trait1 == paste0(t, "_all"), ]
  # take DGE values from real data
  DGEonly = DGEonly_VCs[DGEonly_VCs$trait1 == paste0(t, "_all"), "prop_Ad1"]
  DGEnIGE = all_VCs_full[all_VCs_full$trait1 == paste0(t, "_all"), "prop_Ad1"]
  
  # in all three ASV there is one outlier in permutations, select it to break the plot
  outn = which.max(toplot$prop_Ad1)
  noout = toplot$prop_Ad1[-outn]
  outl = toplot$prop_Ad1[outn]
  
  ### boxplot - bottom part 
  par(mar = c(bmar, lmar, 0, rmar))
  # prepare base plot (no box yet)
  bp = boxplot(noout, 
               col = adjustcolor("white", 0), border=adjustcolor("white", 0), 
               ylim = round(range(c(range(noout), DGEonly, DGEnIGE)), 3), 
               outline=T, ylab = "", xlab = "", xaxt="n",
               drawRect=F, varwidth = TRUE, 
               frame=F)
  # add ticks and ticks labels to x axis
  axis(1, at=1, labels=t, tick = F, font = 2)
  # add y label
  title(ylab = paste0("               Mic-DGE"), line=4.5, font.lab=2)
  # add points
  points(x = jitter(rep(1, length(noout)), factor=foctor), 
         y = noout, cex = pch.cex,
         col = adjustcolor("grey50", 0.5))
  # add box
  boxplot(noout, outline = F, add=T, xaxt="n", yaxt="n", frame=F, 
          boxwex = bw, 
          col = adjustcolor("white", 0.7), border = "black")
  box(bty="c")
  # line of prop_Ad1 when IGE or when DGE only
  abline(h = DGEnIGE, col="red", lwd = 3, lty = 2)
  abline(h = DGEonly, col="#3C5488FF", lwd = 3, lty = 2)
  
  ### outlier - top part
  par(mar = c(0.5, lmar, tmar, rmar))
  plot(1, outl, 
       xaxt="n", yaxt="n",
       ylim = c(outl-0.002, outl+0.002), 
       ylab = "", cex = pch.cex,
       col =adjustcolor("grey50", 0.5), bty = "l")
  axis(2, at = round(outl, 3))
}
# add legend (as a new plot on the right)
par(mar = c(0, 0.5, tmar, 0))
plot.new()
legend("topleft", lty = 2, lwd = 2, 
       col = c("red", "#3C5488FF"), legend = c("DGE+IGE", "DGE only"), 
       cex = 1.2,
       bty = "o")
plot.new() # to fill last plot

### For IGE
layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4, byrow = TRUE)
layout(layout_matrix, widths = c(1, 1, 1, 0.7))

for(t in traitplot){
  set.seed(101) # for jitter
  # select trait to plot
  toplot = scrambled_VCs_full_model[scrambled_VCs_full_model$trait1 == paste0(t, "_all"), ]
  # take IGE values from real data
  DGEnIGE = all_VCs_full[all_VCs_full$trait1 == paste0(t, "_all"), "prop_As1"]
  
  # no need to select for outlier as it is in the range of the value from real data
  noout = toplot$prop_As1
  ylimi = round(range(c(range(noout), DGEnIGE)), 3)
  ### boxplot 
  par(mar = c(bmar, lmar, tmar, rmar))
  # prepare base plot (no box yet)
  bp = boxplot(noout, 
               col = adjustcolor("white", 0), border=adjustcolor("white", 0), 
               ylim = ylimi, 
               outline=T, ylab = "", xlab = "", xaxt="n",
               drawRect=F, varwidth = TRUE, 
               frame=F)
  # add ticks and ticks labels to x axis
  axis(1, at=1, labels=t, tick = F, font = 2)
  # add y label
  title(ylab = paste0("Mic-IGE"), line=4.5, font.lab=2)
  # add points
  points(x = jitter(rep(1, length(noout)), factor=foctor), 
         y = noout, cex = pch.cex,
         col = adjustcolor("grey50", 0.5))
  # add box
  boxplot(noout, outline = F, add=T, xaxt="n", yaxt="n", boxwex = bw,
          frame=F, col = adjustcolor("white", 0.7), border = "black")
  box(bty="l")
  # line of prop_As1 when IGE 
  abline(h = DGEnIGE, col="red", lwd = 3, lty = 2)
}
# add legend (as a new plot on the right)
par(mar = c(0, 0.5, tmar, 0))
plot.new()
legend("topleft", lty = 2, lwd = 2, 
       col = c("red"), legend = c("DGE+IGE"), 
       cex = 1.2,
       bty = "o")

dev.off()
