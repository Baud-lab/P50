# Load heritability data
load('augmented_DGE_VC_wALL.RData')

# Build center_spe_herits table with one row per microbiome phenotype and 4 columns corresponding to 4 cohorts
#filtering out results for "all" - focus on different centers
all_VCs_full = all_VCs_full[all_VCs_full$study1 != "all",]
union = unique(all_VCs_full$taxon)
center_spe_herits = matrix(nrow = length(union), ncol = 4, NA)
rownames(center_spe_herits) = union
colnames(center_spe_herits) = c('NY', 'MI', 'TN_behavior','TN_breeder')

for (i in 1:dim(all_VCs_full)[1]) {
    center_spe_herits[all_VCs_full[i,'taxon1'],all_VCs_full[i,'study1']] = all_VCs_full[i,'prop_Ad1']
} 

# Define cohorts names as in paper and order as in rest of the figures
dict = c('NY'= 'NY',
         'MI'= 'MI',
         'TN_behavior'= 'TN1',
         'TN_breeder'= 'TN2')
colnames(center_spe_herits) = unname(dict[colnames(center_spe_herits)])

#NY cohort  (N = 1,167 rats), 
#MI cohort (N = 1,112 rats), 
#TN behaviour cohort (N = 950), aka TN1 
#TN breeder cohort (N = 555). aka TN2
Ns = c("NY\n(N = 1,167)", "MI\n(N = 1,112)", "TN1\n(N = 950)", "TN2\n(N = 555)")

# Colnames for plotting
colnames(center_spe_herits) = Ns[match(colnames(center_spe_herits), unlist(lapply(strsplit(Ns, "\n"),"[[", 1)))]

# Setting dot colours for all
#dotcol = rep("grey20", nrow(center_spe_herits))
inscol = "grey50" # colour of non-sign
dotcol = rep(inscol, nrow(center_spe_herits))
names(dotcol) = rownames(center_spe_herits)

# Load results from porcupine for dot colour
load("porcupine_colors.RData")
# Selecting only significant ones
#sign = tosave[which(tosave$col != "darkgrey"),] # all significant ones
top_taxa = c("ASV_3613_NY","ASV_3613_MI",
             "ASV_3613_TN2","ASV_18566_NY",
             "ASV_18566_MI","ASV_18566_TN1",
             "ASV_5163_NY","ASV_5163_TN1") # ref to supplementary tables 2-4
sign = tosave[tosave[,"trait1"] %in% top_taxa & tosave$col != "darkgrey",] # only 3 top peaks
rm(tosave) # no need ot keep and quite big
sign[,"slim_trait1"] = gsub("_MI|_NY|_TN_breeder|_TN_behavior", "", sign$trait1)
row_sig = names(dotcol)[names(dotcol) %in% sign[,"slim_trait1"]]
motch = match(row_sig, sign$slim_trait1)
#sign[motch,"slim_trait1"] == row_sig # rownames(center_spe_herits[row_sig,]) == sign[motch,"slim_trait1"]
dotcol[row_sig] = sign[motch, "col"] 

# Function for plot on lower triangle
my_cor <- function(x, y, ...) {
  cor = cor.test(x, y, use = 'pairwise.complete.obs')
  txt <- paste('cor = ',format(cor$estimate, digits = 2)[1],sep='')
  if (cor[['p.value']] < (0.05/6)) col = colr else col = 'black' # color = colr and font = bold if passes Bonferroni correction
  text(0.1, 0.1, 
       txt, cex = 1.4, col = col)
}

# Helper function to create nice axis breaks - used in 'my_points'
create_axis_breaks <- function(x,y, lim=NULL) {
  # Get the range for the current panel
  if(is.null(lim)){
    range_vals <- range(c(x, y), na.rm = TRUE)
  }else{
    # Specify range_vals 
    range_vals <- lim
  }
  
  # Create regular breaks for ticks (more frequent)
  tick_count <- 5  # Number of intervals desired
  tick_breaks <- pretty(range_vals, n = tick_count)
  
  # Create breaks for labels (less frequent)
  label_count <- 3  # Number of intervals desired for labels
  label_breaks <- pretty(range_vals, n = label_count)  # -1 to account for endpoints
  
  # Create labels vector (empty strings for ticks without labels)
  labels <- ifelse(tick_breaks %in% label_breaks,
                   as.character(tick_breaks),
                   "")
  
  return(list(ticks = tick_breaks, labels = labels))
}

# Modified points function for upper triangle
my_points <- function(x, y, ...) {
  cexaxis = 1.25
  lim = c(0,0.2)
  # Get current plot coordinates
  mfg <- par('mfg')
  current_row <- mfg[1]
  current_col <- mfg[2]
  total_rows <- mfg[3]
  
  # Plot black points first
  black_indices <- which(pch.col == pch.bg )
  if(length(black_indices) > 0) {
    points(x[black_indices], y[black_indices], pch = 16, col =pch.bg)
  }
  
  # Then plot colored points on top
  colored_indices <- which(pch.col != pch.bg)
  if(length(colored_indices) > 0) {
    points(x[colored_indices], y[colored_indices], pch = 16, col = pch.col[colored_indices])
  }
  
  ## # Plot points (all together, no colour priority)
  ## points(x, y, pch = 16, col= pch.col)
  
  # Get axis breaks for this panel
  breaks <- create_axis_breaks(x,y, lim) # NB: change here depending on x-ylim
  
  # Add top axis only for top row
  if (current_row == 1) {
    axis(side = 3, at = breaks$ticks, labels = breaks$labels, gap.axis = 2,las = 1, cex.axis = cexaxis)
  }
  
  # Add right axis only for rightmost column
  if (current_col == total_rows) {
    axis(side = 4, at = breaks$ticks, labels = breaks$labels,las = 1, cex.axis = cexaxis)
  }
}


# Define colours for significant correlation - used in 'my_cor'
colr = "#E64B35FF"
pch.col = adjustcolor(dotcol, alpha.f = 1)
pch.bg = adjustcolor(inscol, alpha.f = 1)

# Open pdf to save plot
# NB: pdf was edited in Adobe Illustrator to do the following:
#     keep dot colour for only those dots which taxon showed a significant QTL in both cohorts 
#     increase the size of the dots which taxon showed a significant QTL in both cohorts 
pdf('compare_herits_diff_centers.pdf', w=6,h=6)

# Plot using R's pairs plot 
# calling 'my_points' - upper tri - and 'my_cor' - lower tri
pairs(center_spe_herits,
      panel = points,
      xlim = c(0,0.2), ylim = c(0,0.2), #cex.main=2,
      upper.panel = my_points, # NB: in function my_points() change the lim according to ylim and xlim and cexaxis
      lower.panel = my_cor,
      xaxt='n',yaxt='n', 
      font.labels = 2,
      gap = 0.6,
      oma = c(2.1, 2.1, 5.1, 5.1))

dev.off()

