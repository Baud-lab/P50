#load herit data
root_dir = '/users/abaud/abaud/P50_HSrats/output/VD/univariate/'
load(file.path(root_dir,'augmented_VC.RData'))

#build center_spe_herits table with one row per microbiome phenotype and 4 columns corresponding to 4 cohorts
union = unique(all_VCs_full$taxon)
center_spe_herits = matrix(nrow = length(union), ncol = 4, NA)
rownames(center_spe_herits) = union
colnames(center_spe_herits) = c('NY', 'MI', 'TN_behavior','TN_breeder')

# this could be done with apply
for (i in 1:dim(all_VCs_full)[1]) {
    center_spe_herits[all_VCs_full[i,'taxon1'],all_VCs_full[i,'study1']] = all_VCs_full[i,'prop_Ad1']
} 

# ordering to align as the rest of the plots in figure
#center_spe_herits = center_spe_herits[,c('TN_breeder', 'TN_behavior',  'MI', 'NY')]
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


#plot using R's pairs plot
range(all_VCs_full$prop_Ad1)
#pdf('~/P50_HSrats/plots/compare_herits_Helenes_diff_centers.pdf')
#names(par())
#par(mar=c(5.1,5.1,2.1,2.1))

my_cor <- function(x, y,...) {
  cor = cor.test(x, y, use = 'pairwise.complete.obs')
  txt <- paste('cor = ',format(cor$estimate, digits = 2)[1],sep='')
  if (cor[['p.value']] < (0.05/6)) col = colr else col = 'black' # color = colr and font = bold if passes Bonferroni correction
  text(0.1, 0.1, 
       txt, cex = 1.4, col = col)
}

# Helper function to create nice axis breaks
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
  
  # Plot points
  points(x, y, pch = 16)
  
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

cols = rep('black')
colr="#E64B35FF"

pdf('/users/abaud/htonnele/PRJs/P50_HSrats/16S/plot/compare_herits_Helenes_diff_centers.pdf', w=6,h=6)

# Main plotting call
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


####### Trying with ggally
#ggpairs(center_spe_herits[!apply(center_spe_herits, 1, function(x) any(is.na(x))),], 
#        upper = list(continuous = "points", na="NA"), 
#        lower = list(continuous = "cor", na="NA"), 
#        diag = colnames(center_spe_herits))

# length(unique(which(is.na(center_spe_herits), arr.ind = T)[,"row"])) # checking NAs
# table(apply(center_spe_herits, 1, function(x) any(is.na(x))))
# dim(center_spe_herits[!apply(center_spe_herits, 1, function(x) any(is.na(x))),])
