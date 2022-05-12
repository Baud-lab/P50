#all_poss_annots = read.delim('/homes/abaud/P50_HSrats/data/metabo/all_poss_annots.txt',as.is = T, header = F)
#uniq = unique(c(all_poss_annots[,1],all_poss_annots[,2]))
#sort(uniq)
# [1] "[2M-2H2O+K]+" "[2M-H2O+H]+"  "[2M-H2O+K]+"  "[2M+H]+"      "[2M+K]+"     
# [6] "[2M+Na]+"     "[M-2H2O+?]+"  "[M-2H2O+?]2+" "[M-2H2O+H]+"  "[M-H2O+?]+"  
#[11] "[M-H2O+?]2+"  "[M-H2O+H]+"   "[M-H2O+K]+"   "[M-H2O+Na]+"  "[M+?]+"      
#[16] "[M+?]2+"      "[M+H]+"       "[M+K]+"       "[M+Na]+"     

exact_masses_adducts = c('H' = 1.00728, 'H2O' = 18.01056, 'Na' = 22.98922, 'Ca' = 39.96204, 'K' = 38.96316)
exact_masses_BAs = c('(cheno)deoxycholic/ursocholic acid' = 392.293,
						'cholic acid' = 408.288,
						'glycocholic acid' = 465.309,
						'taurocholic acid' = 515.292,
						'glycochenodeoxycholic acid' = 449.314,
						'taurochenodeoxycholic acid' = 499.297,
						'lithocholic acid' = 376.298,
						'12ketodeoxycholic acid' = 390.277,
						'coprocholic acid' = 450.335,
						'3bhydroxy5cholenoic acid' = 374.282,
						'9HpOTrE' = 310.214,
						'pinolenic acid' = 278.225, 
						'formic acid' = 46.005, 
						'acetic acid' = 59.013,
						'propionic acid' = 73.029,
						'butyric acid' = 87.045,
						'lactic acid' = 90.032,
						'linoleic acid/conjugated linoleic acid' = 280.24,
						'EpOME' = 296.235,
						'13-Keto-9Z,11E-octadecadienoic acid/9-Oxo-10E,12Z-octadecadienoic acid/9S-HOTrE' = 294.219,
						'3-Hydroxydodecanoic acid' = 216.173,
						'(12Z)-9,10-Dihydroxyoctadec-12-enoic acid' = 314.246,
						'Octadecatrienoic Acid ethyl ester' = 306.256,
						'Eicosatetraenoic acid' = 304.24,
						'Docosapentaenoic acid' = 330.256,
						'Glutamic acid' = 147.053,
						'Myristoleic acid' = 226.193,
						'Oleanolic acid' = 456.36,
						'9-Hexadecenoic acid' = 254.225,
						'Stearidonic acid' = 276.209,
						'Sumaresinolic acid' = 472.355,
						'Enterolactone' = 298.121,
						'glucoronic acid' = 194.043,
						'Ferulic acid' = 194.058,
						'Fertaric acid' = 326.064,
						'Caffeic acid' = 180.042,
						'Vanillin' = 152.047,
						'Secoisolariciresinol diglucoside' = 686.279,
						'Linolenic acid' = 278.225,
						'Palmitic acid' = 256.24,
						'Oleic acid' = 282.256,
						'Cortisol' = 362.209,
						'Lysophosphatidic acid' = 436.259,
						'Dopamine' = 153.079,
						'Levodopa' = 197.069,
						'Amantadine' = 151.136,
						'2-Butanone, 4-(2,6,6-trimethyl-2-cyclohexen-1-yl)' = 208.183,
						'Adipic acid' = 146.058,
						'Suberic acid' = 318.168,
						'Sebacic acid' = 318.183,
						'Mevalonic acid' = 148.074,
						'Geranyl pyrophosphate' = 314.068,
						'Phenobarbital' = 232.085,
						'Borneol/Eucalyptol/Geraniol/Linalool' = 154.136,
						'Menthol/Citronellol' = 156.151,
						'Camphor' = 152.12,
						'Carene/Sabinene/Camphene/Thujene/Limonene' = 136.125,
						'AcetylCoA' = 809.126)

calculate_derivative_mzs = 	function(target) {
	derivative_mzs = c('[M+H]+' = exact_masses_BAs[target] + exact_masses_adducts['H'],
					 '[M+Na]+' = exact_masses_BAs[target] + exact_masses_adducts['Na'], 
					 '[M+K]+' = exact_masses_BAs[target] + exact_masses_adducts['K'], 
					 '[M+Ca]+' = exact_masses_BAs[target] + exact_masses_adducts['Ca'], 
					 '[M-H2O+H]+' = exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['H'], 
					 '[M-2H2O+H]+' = exact_masses_BAs[target] - 2*exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
					 '[M-3H2O+H]+' = exact_masses_BAs[target] - 3*exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
#					 '[M-4H2O+H]+' = exact_masses_BAs[target] - 4*exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
					 '[M-H2O+Na]+' = exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['Na'], 
					 '[M-2H2O+Na]+' = exact_masses_BAs[target] - 2*exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
					 '[M-3H2O+Na]+' = exact_masses_BAs[target] - 3*exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
#					 '[M-4H2O+Na]+' = exact_masses_BAs[target] - 4*exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
					 '[M-H2O+K]+' = exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['K'], 
					 '[M-2H2O+K]+' = exact_masses_BAs[target] - 2*exact_masses_adducts['H2O'] + exact_masses_adducts['K'],
					 '[M-3H2O+K]+' = exact_masses_BAs[target] - 3*exact_masses_adducts['H2O'] + exact_masses_adducts['K'],
#					 '[M-4H2O+K]+' = exact_masses_BAs[target] - 4*exact_masses_adducts['H2O'] + exact_masses_adducts['K'],
					 '[M-H2O+Ca]+' = exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['Ca'], 
					 '[M-2H2O+Ca]+' = exact_masses_BAs[target] - 2*exact_masses_adducts['H2O'] + exact_masses_adducts['Ca'],
					 '[M-3H2O+Ca]+' = exact_masses_BAs[target] - 3*exact_masses_adducts['H2O'] + exact_masses_adducts['Ca'],

					 '[M+2H]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'])/2, 
					 '[M-H+2Na]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'])/2, 
					 '[M-H+2K]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'])/2, 
					 '[M-H+2Ca]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Ca'])/2, 
					 '[M-H+NaK]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'])/2, 
					 '[M-H+NaCa]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['Ca'])/2, 
					 '[M-H+CaK]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Ca'] + exact_masses_adducts['K'])/2, 
					 '[M-H+2Na-H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - exact_masses_adducts['H2O'])/2, 
					 '[M-H+2K-H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[M-H+2Ca-H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Ca'] - exact_masses_adducts['H2O'])/2, 
					 '[M-H+NaK-H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[M-H+NaCa-H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['Ca'] - exact_masses_adducts['H2O'])/2, 
					 '[M-H+CaK-H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Ca'] + exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[M-H+2Na-2H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-H+2K-2H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-H+2Ca-2H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Ca'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-H+NaK-2H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-H+NaCa-2H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['Ca'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-H+CaK-2H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Ca'] + exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-H+2Na-3H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-H+2K-3H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-H+2Ca-3H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Ca'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-H+NaK-3H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-H+NaCa-3H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['Ca'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-H+CaK-3H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Ca'] + exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 
#					 '[M-H+2Na-4H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - 4*exact_masses_adducts['H2O'])/2, 
#					 '[M-H+2K-4H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - 4*exact_masses_adducts['H2O'])/2, 
#					 '[M-H+2Na-4H20]2+' = (exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 4*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+3Na]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['Na'])/2, 
					 '[M-2H+3K]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['K'])/2, 
					 '[M-2H+2NaK]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] + exact_masses_adducts['K'])/2, 
					 '[M-2H+Na2K]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + exact_masses_adducts['Na'] + 2*exact_masses_adducts['K'])/2, 
					 '[M-2H+3Na-H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['Na'] - exact_masses_adducts['H2O'])/2, 
					 '[M-2H+3K-H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[M-2H+2NaK-H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] + exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[M-2H+Na2K-H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + exact_masses_adducts['Na'] + 2*exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[M-2H+3Na-2H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['Na'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+3K-2H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+2NaK-2H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+Na2K-2H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + exact_masses_adducts['Na'] + 2*exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+3Na-3H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['Na'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+3K-3H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 3*exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+2NaK-3H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[M-2H+Na2K-3H2O]2+' = (exact_masses_BAs[target] - 2*exact_masses_adducts['H'] + exact_masses_adducts['Na'] + 2*exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 

					 '[2M+H]+' = 2 * exact_masses_BAs[target] + exact_masses_adducts['H'],
					 '[2M+Na]+' = 2 * exact_masses_BAs[target] + exact_masses_adducts['Na'],
					 '[2M+K]+' = 2 * exact_masses_BAs[target] + exact_masses_adducts['K'],
					 '[2M-H2O+H]+' = 2 * exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
					 '[2M-2H2O+H]+' = 2 * exact_masses_BAs[target] - 2 *exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
					 '[2M-3H2O+H]+' = 2 * exact_masses_BAs[target] - 3 *exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
#					 '[2M-4H2O+H]+' = 2 * exact_masses_BAs[target] - 4 *exact_masses_adducts['H2O'] + exact_masses_adducts['H'],
					 '[2M-H2O+Na]+' = 2 * exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
					 '[2M-2H2O+Na]+' = 2 * exact_masses_BAs[target] - 2 *exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
					 '[2M-3H2O+Na]+' = 2 * exact_masses_BAs[target] - 3 *exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
#					 '[2M-4H2O+Na]+' = 2 * exact_masses_BAs[target] - 4 *exact_masses_adducts['H2O'] + exact_masses_adducts['Na'],
					 '[2M-H2O+K]+' = 2 * exact_masses_BAs[target] - exact_masses_adducts['H2O'] + exact_masses_adducts['K'],
					 '[2M-2H2O+K]+' = 2 * exact_masses_BAs[target] - 2 *exact_masses_adducts['H2O'] + exact_masses_adducts['K'],
					 '[2M-3H2O+K]+' = 2 * exact_masses_BAs[target] - 3 *exact_masses_adducts['H2O'] + exact_masses_adducts['K'],
#					 '[2M-4H2O+K]+' = 2 * exact_masses_BAs[target] - 4 *exact_masses_adducts['H2O'] + exact_masses_adducts['K'],

					 '[2M+2H]2+' = (2 * exact_masses_BAs[target] + 2*exact_masses_adducts['H'])/2,
					 '[2M-H+2Na]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'])/2, 
					 '[2M-H+2K]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'])/2, 
					 '[2M-H+2Na]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'])/2, 
					 '[2M-H+2Na-H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2K-H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2Na-H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2Na-2H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2K-2H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2Na-2H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 2*exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2Na-3H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2K-3H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2, 
					 '[2M-H+2Na-3H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 3*exact_masses_adducts['H2O'])/2
#					 '[2M-H+2Na-4H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['Na'] - 4*exact_masses_adducts['H2O'])/2, 
#					 '[2M-H+2K-4H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + 2*exact_masses_adducts['K'] - 4*exact_masses_adducts['H2O'])/2, 
#					 '[2M-H+2Na-4H20]2+' = (2*exact_masses_BAs[target] - exact_masses_adducts['H'] + exact_masses_adducts['Na'] + exact_masses_adducts['K'] - 4*exact_masses_adducts['H2O'])/2
					 )
	names(derivative_mzs) = sub('.','_',names(derivative_mzs),fixed = T)
	return(derivative_mzs)
}

load('/homes/abaud/P50_HSrats/data/metabo/annots_unnormalized_metabo.RData')
load("/nfs/leia/research/stegle/abaud/P50_HSrats/output/pvalues/metabo/round2_include_DGE/QTLs_alpha0.001.RData")
motch = match(sub('_NY','',sub('_MI','',all_QTLs$measure)),rownames(annots))
all_QTLs = cbind(all_QTLs,annots[motch,])
#w = which(all_QTLs$chr == 16 & abs(all_QTLs$pos - 90610504) < 100000)
w = which(all_QTLs$chr == 14 & abs(all_QTLs$pos - 92672944) < 100000)
#w = which(all_QTLs$chr == 16 & abs(all_QTLs$pos - 90610504) < 100000 & abs(all_QTLs$RT - 5.74)<0.1)
#w = which(all_QTLs$chr == 14 & abs(all_QTLs$pos - 92672944) < 100000 & abs(all_QTLs$RT - 5.11)<0.1)
check = all_QTLs[w,c('measure','pos','logP','Adduct','Compound_Name','Compound_Source','LibMZ','RT','mz')]
check = check[order(check$mz),]

all_derivative_mzs = c()
for (target in names(exact_masses_BAs)) {
	these_derivative_mzs = calculate_derivative_mzs(target)
	all_derivative_mzs = c(all_derivative_mzs,these_derivative_mzs)
}
check$derivative = NA
for (i in 1:dim(check)[1]) {
	for (j in 1:length(all_derivative_mzs)) {
		if (abs(check[i,'mz'] - all_derivative_mzs[j]) < 0.01) {
			if (is.na(check[i,'derivative'])) check[i,'derivative'] = names(all_derivative_mzs)[j] else
				check[i,'derivative'] = paste(check[i,'derivative'],names(all_derivative_mzs)[j], sep = ';')
		}
	}
}

#these_annots = annots[grepl('5.74',annots$RT) | grepl('5.11',annots$RT),]
#for (mass in these_annots[,'mz']) {
#for (mass in 280.24) {
#	exact_masses_BAs = c(exact_masses_BAs, mass)
#	names(exact_masses_BAs)[length(exact_masses_BAs)] = as.character(mass)
#	all_derivative_mzs = c(all_derivative_mzs,calculate_derivative_mzs(as.character(mass)))
#}

#masses_annots = NULL
#for (i in 1:dim(check)[1]) {
#	for (j in 1:length(all_derivative_mzs)) {
#		if (abs(check[i,'mz'] - all_derivative_mzs[j]) < 0.01) {
#			if (is.na(check[i,'derivative'])) {
#					check[i,'derivative'] = names(all_derivative_mzs)[j] 
#					splot = strsplit(names(all_derivative_mzs)[j],'_')[[1]][2]
#					ma = these_annots[these_annots$mz == as.numeric(splot),]
#					masses_annots = rbind(masses_annots,ma)
#				} else {
#				if (!grepl(names(all_derivative_mzs)[j],check[i,'derivative'], fixed =T)) {
#					check[i,'derivative'] = paste(check[i,'derivative'],names(all_derivative_mzs)[j], sep = ';')
#					splot = strsplit(names(all_derivative_mzs)[j],'_')[[1]][2]
#					ma = these_annots[these_annots$mz == as.numeric(splot),]
#					masses_annots = rbind(masses_annots,ma)
#				}
#			}
#		}
#	}
#}
#masses_annots = masses_annots[!is.na(masses_annots$Compound_Name),]
#save(check,masses_annots, file = '/nfs/leia/research/stegle/abaud/P50_HSrats/output/pvalues/metabo/round2_include_DGE/annotated_chr16_QTL.RData')
#save(check, file = '/nfs/leia/research/stegle/abaud/P50_HSrats/output/pvalues/metabo/round2_include_DGE/annotated_chr14_QTL.RData')
#write.table(check, file = '/nfs/leia/research/stegle/abaud/P50_HSrats/plots/metabo_annotated_chr16_QTL.txt',sep='\t',quote = F, col.names = T, row.names = F)
write.table(check[,c('measure','RT','mz','Compound_Name','Adduct','derivative','logP')], file = '/homes/abaud/P50_HSrats/plots/metabo_annotated_chr16_QTL.txt',sep='\t',quote = F, col.names = T, row.names = F)

#I'm missing a group of features. check what they are correlated to
load('/homes/abaud/P50_HSrats/data/metabo/unnormalized_metabo.RData')
load('/homes/abaud/P50_HSrats/data/general/sample_metadata.RData')
ny_rats = rownames(sample_metadata)[sample_metadata$phenotyping_center == 'NY']
biomt = biomt[,colnames(biomt) %in% ny_rats]
cors = c()
for (row in 1:dim(biomt)[1]) {
	cors = c(cors, cor(biomt['MZ243.2105;5.7437',],biomt[row,]))
}
names(cors) = rownames(biomt)
cors = cors[order(cors, decreasing = T)]
correlates = head(cors,15)
motch = match(names(correlates),sub('_NY','',check$measure))
cbind(check[motch,],correlates)

sort(all_derivative_mzs[abs(all_derivative_mzs-215.1794) < 1])
sort(all_derivative_mzs[abs(all_derivative_mzs-221.1534) < 1])
sort(all_derivative_mzs[abs(all_derivative_mzs-229.1947) < 1])
sort(all_derivative_mzs[abs(all_derivative_mzs-243.2105) < 1])
sort(all_derivative_mzs[abs(all_derivative_mzs-247.1687) < 1])
sort(all_derivative_mzs[abs(all_derivative_mzs-261.1846) < 1])

min(abs(all_derivative_mzs-215.1794))
min(abs(all_derivative_mzs-221.1534))
min(abs(all_derivative_mzs-229.1947))
min(abs(all_derivative_mzs-243.2105))
min(abs(all_derivative_mzs-247.1687))
min(abs(all_derivative_mzs-261.1846))

#so only mz 215 upwards map to transporter on chr16, from NY and RT 5.7
#mz 107 - 207 (which incl. [M+2H]2+_lithocholic acid and [M+2H]2+_(cheno)deoxycholic/ursocholic acid, 9(S)-HpOTrE (fatty acid) from NIST14 and Amantadine from NIST14) also map to decarboxylase
#these are mostly MI and RT 5.1

#also check what other loci features above map to
w1 = which(all_QTLs$chr == 16 & abs(all_QTLs$pos - 90610504) < 100000)
w2 = which(all_QTLs$chr == 14 & abs(all_QTLs$pos - 92672944) < 100000)
check = all_QTLs[unique(c(w1,w2)),c('measure','pos','logP','Adduct','Compound_Name','Compound_Source','LibMZ','RT','mz')]
check = check[order(check$mz),]
check_QTLs = NULL
check_QTLs = all_QTLs[all_QTLs$measure %in%  check$measure,]
check_QTLs = check_QTLs[check_QTLs$logP > 5,]
check_QTLs = check_QTLs[order(check_QTLs$chr,check_QTLs$pos),]
check_QTLs[,c('measure','chr','pos','logP','Adduct','Compound_Name')]

1576  MZ207.1378;5.1163_NY   1  55715202  5.385152    <NA>
1607  MZ229.1946;5.1163_NY   1  55715202  5.381201   229.2
2374  MZ389.2297;3.3363_NY   1 245797899  5.413284    <NA>
11640 MZ785.5899;5.7436_MI   2  25813803  5.381198    2M+H
11022  MZ357.279;5.1161_NY   2  96045478  5.474995 M+H-H2O
2157  MZ207.1378;5.1163_NY   3  99427784  5.350594    <NA>
119   MZ107.0853;5.1165_MI   4 144903030  5.021753    <NA>
1991  MZ339.2681;5.1154_MI   9  41327277  5.058056    <NA>
11062 MZ389.2297;3.3363_NY  10  64682939  5.629973    <NA>
11488 MZ433.3313;6.8099_MI  14  15616655  5.003639    <NA>
11099 MZ425.3411;6.3594_MI  14  22566597 13.968408    <NA>
11487 MZ433.3313;6.8099_MI  14  22829268 20.588923    <NA>
11512 MZ443.3506;6.3797_MI  14  22829268  9.364865    <NA>
11524 MZ455.3145;6.8139_MI  14  22829268 17.596515    <NA>
11544  MZ465.3318;6.361_MI  14  22829268  7.297769    <NA>
1998  MZ341.2835;6.7258_MI  14  24024915  5.628584    <NA>
11391 MZ359.2882;6.7258_MI  14  24024915  5.224225 M+H-H2O
1506  MZ203.1793;5.1163_NY  14  92660368  5.350089    <NA>
3380   MZ465.3318;6.361_MI  14  92660368  5.923592    <NA>
3272  MZ341.2835;6.7258_MI  14  92672944  5.169008    <NA>
3345  MZ425.3411;6.3594_MI  14  92672944  6.991395    <NA>
3376  MZ455.3145;6.8139_MI  14  92672944  9.849975    <NA>
3362  MZ433.3313;6.8099_MI  14  92724095 11.790916    <NA>
1587  MZ215.1794;5.7433_NY  16  90610504  6.108728    <NA>
1596  MZ221.1534;5.7437_NY  16  90610504  6.067211    <NA>
1609  MZ229.1947;5.7439_NY  16  90610504  5.410595   229.2
1689  MZ243.2105;5.7437_NY  16  90610504  5.940315    <NA>
1699  MZ247.1687;5.7444_NY  16  90610504  6.315550    <NA>
1959  MZ321.2568;5.7431_NY  16  90610504  6.406137    <NA>
1994  MZ339.2681;5.7435_NY  16  90610504  6.695119    <NA>
11089 MZ415.2806;5.7433_NY  16  90610504  6.751622    <NA>
11471 MZ431.2469;5.7424_NY  16  90610504  5.061301    <NA>
11646 MZ807.5734;5.7436_NY  16  90610504  5.526682    <NA>
2309  MZ341.2835;6.7258_MI  18  61609181  5.178403    <NA>
11025 MZ359.2882;6.7258_MI  18  61609181  5.364053 M+H-H2O

# [2] "Massbank:EA282304 Benzoylecgonine| 3-benzoyloxy-8-methyl-8-azabicyclo[3.2.1]octane-4-carboxylic acid"                                                                                                                                                                                                                                                 
# from cocaine
# check QTLs
annots[which(annots$Compound_Name == "Massbank:EA282304 Benzoylecgonine| 3-benzoyloxy-8-methyl-8-azabicyclo[3.2.1]octane-4-carboxylic acid"),]
#MZ290.1374;2.5238 [M+H]+  #adduct checked
all_QTLs[all_QTLs$measure %in% c('MZ290.1374;2.5238_NY','MZ290.1374;2.5238_MI'),]
#no QTLs






