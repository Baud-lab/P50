load('/homes/abaud/P50_HSrats/data/metabo/78940_palmer_amelie_metabolomics_original.RData')

max(biomt)
#[1] 6.59e+09

pdf('/homes/abaud/P50_HSrats/plots/distributions_metabolites.pdf')
s = sample(1:dim(biomt)[1], size = 100)
nonnulls = c()
for (i in sort(s)) {
	hist(biomt[i,])
	nonnulls = c(nonnulls, sum(biomt[i,]!=0))
}
dev.off()

#lots of 0s always, very large numbers on the x axis
#sometimes nice exponential decay
#sometimes peak at 0 and then kind of normal distribution
#most of the time non 0s not visible

