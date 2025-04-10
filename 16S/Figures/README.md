Data files and intermediate results files available from https://figshare.com/account/home#/collections/7761632

## Figure 2. Variation in the HS rat gut microbiome

#### Panel A.
Code: `2a.average_taxonomy_barplots.R` 
<details>
<summary>Input:</summary>

+ Full biomatrix -> _collapsed\_full\_biomt_
```
collapsed_full_biomt_collapsed_clr_counts.RData # created by P50/16S/Preprocessing/4_merge_taxonomic_level.R
```

+ Metadata
```
metadata_16Spaper.RData # figshare
```
</details>

Output: **average\_genera\_barplots\_\_f.pdf**
<br/><br/>

#### Panel B.
Code: `2b.PCA.R`
<details>
<summary>Input:</summary>

+ CLR counts -> _collapsed\_clr\_counts_
```
collapsed_full_biomt_collapsed_clr_counts.RData # created by P50/16S/Preprocessing/4_merge_taxonomic_level.R
```

+ Metadata 2
```
metadata_16Spaper.RData # figshare
```
</details>
 
Output: **PCA_paper.pdf**
<br/><br/>

## Figure 3. Characteristics of polygenic host genetic effects

#### Panel A.
Code: `3a.compare_herit_microbes_phenos.R`
<details>
<summary>Input:</summary>
  
+ VD of phenotypes
```
phenos_all_estNste.Rdata # figshare
```

+ VD of microbiome
```
augmented_VC.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```
</details>
  
  Output: **barplots\_herits\_studies\_pheno.pdf**
<br/><br/>

#### Panel B (and Supp. Fig. 5).
Code: `3b.dataPrep_prev_abund_herit.R` + `3b.prev_abund_herit.R`
<details>
<summary>Input - dataPrep:</summary>

+ For ASVs
```
full_biomt_clr_counts.RData # created by P50/16S/Preprocessing/3_clr_counts.R
```

+ For taxa
```
collapsed_full_biomt_collapsed_clr_counts.RData # created by P50/16S/Preprocessing/4_merge_taxonomic_level.R
```

+ Heritability data
```
augmented_VC.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```
</details>

<details>
<summary>Input - plot:</summary>

+ For ASVs - intermediate output from 3b.dataPrep_prev_abund_herit.R
```
prev_abund_asvs_biomt.RData 
```

+ For taxa - intermediate output from 3b.dataPrep_prev_abund_herit.R
```
prev_abund_taxa_biomt.RData
```
</details>

Output: **prev\_abund\_herit\_estimate\_biomt.pdf**
<br/><br/>

#### Panel C.
Code: `3c.boxplots_different_VCs.R`
<details>
<summary>Input:</summary>

+ Heritability data
```
augmented_VC.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```
</details>
 
Output: **VCs\_merged\_viopl\_col.pdf**
<br/><br/>

#### Panel D.
Code: `3d.herit_corrs.R`
<details>
<summary>Input:</summary>

+ Heritability data 
```
augmented_VC.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```
</details>
  
Output: **compare\_herits\_diff\_centers.pdf**
<br/><br/>

#### Panel E.
Code: `3e.boxplot_gen_corrs.R`
<details>
<summary>Input:</summary>

+ Genetic correlations 
```
all_VCs_corr_Ad1d2_zero_P50_Rn7_pruned_DGE.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```
</details>

Output: **comp\_gen\_corrs\_across\_cohorts.pdf**
<br/><br/>

## Figure 4. Microbiome-associated loci

#### Panel.
Code: `4.dataPrep_porcupine_plot.R` + `4.porcupine_plot.R` <br/>
Source: `annotate_VCs_pvalues.R` - annotate() function

<details>
<summary>Input - dataPrep:</summary>

+ Cumulative position annotation 
```
cumpos_P50_rats_Rn7.RData # figshare
```

+ Unpruned QTLs
```
QTLs_alpha1e-04_unpruned.RData # figshare
```

</details> 

<details>
<summary>Input - plot:</summary>

+ summarised QTLs - intermediate output from dataPrep
```
QTLs_alpha1e-04_unpruned_DGE_CE_MaE_toPlot.RData
```

</details>
  
Output: **porcupine\_uncollapsed\_genus2.pdf**
<br/><br/>

## Figure 5. Association between Paraprevotella and the *St6galnac1* locus on chromosome 10

#### Panel D.
Code: `5d.GWAS_boxplots.R`
<details>
<summary>Input - plot:</summary>

+ CLR and raw counts 
```
full_biomt_clr_counts.RData # created by P50/16S/Preprocessing/3_clr_counts.R
```

+ Geno positions 
```
/users/abaud/abaud/P50_HSrats/data/dosages/P50_Rn7_chr10qtl_allSNPS.raw # figshare
```

+ Metadata 
```
metadata_16Spaper.RData # figshare
```
</details>

Output: **all\_chr10\_boxplots\_raw\_counts.pdf**
<br/><br/>

## Figure 6. Indirect (social) genetic effects on microbiome phenotypes

#### Panel B. 
Code: `6b.qqplot_micIGE.R`
<details>
<summary>Input:</summary>

+ VC data 
```
augmented_IGE_VC_allOnly.RData # figshare
```
</details>

Output: **QQplot\_pvalues\_IGE\_Helenes.pdf**
<br/><br/>

#### Panel C. 
Code: `6c.total_genetic_variance_barplot.R`
<details>
<summary>Input:</summary>

+ VD with DGE and IGE
```
augmented_IGE_VC_allOnly.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```

+ VD with DGE only (without IGE)
```
augmented_VC.RData # created by P50/16S/Preprocessing/annotate_VCs_pvalues.R from output of CoreQuantGen; also available from figshare
```
</details>

Output: **tot\_herit\_barplot.pdf**
<br/><br/>


#### Panel D (and Supp. Fig. 14). 
Code: `6d.simulations.R` <br/>
<details>
<summary>Input:</summary>
  
+ Results from simulations (MI and NY)
```
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/simulations/MI/P50_Rn7_pruned_DGE_cageEffect_None_all_estNste*.Rdata #(17) (DGEonly-MI)
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/simulations/MI/P50_Rn7_pruned_DGE_IGE_cageEffect_None_all_estNste*.Rdata #(18) (wt IGE-MI)

/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/simulations/NY/P50_Rn7_pruned_DGE_cageEffect_None_all_estNste*.Rdata #(19) (DGEonly-NY)
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/simulations/NY/P50_Rn7_pruned_DGE_IGE_cageEffect_None_all_estNste*.Rdata #(20) (wt IGE-NY)
```

+ Simulated values (MI and NY)
```
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/simulations/MI/params_uni*.txt #(21 -MI)
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/simulations/NY/params_uni*.txt #(22 -NY)
```
</details>

Output: **{MI,NY}\_DG1\_IG1\_VCs\_from\_sim\_0.9.0.0.neg0.9.pdf** 
<br/><br/>


## Supp. Figure 3. Comparison of heritability at different taxonomic levels

#### Panels. 
Code: `supp3.boxplot_herit_taxlevel.R` 
<details>
<summary>Input:</summary>

+ VD data 
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_VC.RData #(5)
```
</details>

Output: **herit\_per\_tax\_level.pdf**
<br/><br/>

## Supp. Figure 4. Decomposition of the variance of microbiome phenotypes

#### Panels. 
Code: `supp4.baboon_VD_figure.R` 
<details>
<summary>Input:</summary>

+ VD data 
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_VC.RData #(5)
```
</details>

Output: **baboon\_fig\_{study}.pdf**
