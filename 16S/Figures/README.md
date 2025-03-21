## Figure 2. Variation in the HS rat gut microbiome

#### Panel A.
Code: `2a.average_taxonomy_barplots.R` 
<details>
<summary>Input:</summary>

+ Full biomatrix -> _collapsed\_full\_biomt_
```
/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData #(1)
```

+ Metadata
```
/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo.RData #(2)
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
/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData #(1)
```

+ Metadata 2
```
/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData #(3)
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
/users/abaud/abaud/P50_HSrats/output/VD/univariate/residuals_02Jan2020/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/all_estNste.Rdata #(4)
```

+ VD of microbiome
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_VC.RData #(5)
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
/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData #(6)
```

+ For taxa
```
/users/abaud/data/secondary/P50_HSrats/felipes_deblur/collapsed_full_biomt_collapsed_clr_counts.RData #(1)
```

+ Heritability data
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_VC.RData #(5)
```
</details>

<details>
<summary>Input - plot:</summary>

+ For ASVs - intermediate output from dataPrep
```
/users/abaud/abaud/P50_HSrats/output/prev_abund_asvs_biomt.RData 
```

+ For taxa - intermediate output from dataPrep
```
/users/abaud/abaud/P50_HSrats/output/prev_abund_taxa_biomt.RData
```
</details>

Output: **prev\_abund\_herit\_estimate\_biomt.pdf**
<br/><br/>

#### Panel C.
Code: `3c.boxplots_different_VCs.R`
<details>
<summary>Input:</summary>

+ VD for ASVs
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/all_estNste.Rdata #(7)
```

+ VD for taxa
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/deblur_counts/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/all_estNste.Rdata #(8)
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
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_VC.RData #(5)
```
</details>
  
Output: **compare\_herits\_Helenes\_diff\_centers.pdf**
<br/><br/>

#### Panel E.
Code: `3e.boxplot_gen_corrs.R`
<details>
<summary>Input:</summary>

+ Genetic correlations 
```
/users/abaud/abaud/P50_HSrats/output/VD/bivariate/all_VCs_corr_Ad1d2_zero_P50_Rn7_pruned_DGE.RData #(9)
```
</details>

Output: **comp\_gen\_corrs\_across\_cohorts.pdf**
<br/><br/>

## Figure 4. Microbiome-associated loci

#### Panel.
Code: `4.dataPrep_porcupine_plot.R` + `4.porcupine_plot.R` <br/>
Source: `fun_annotate_VCs_pvalues.R` - annotate() function

<details>
<summary>Input - dataPrep:</summary>

+ Cumulative position annotation 
```
/users/abaud/abaud/P50_HSrats/data/cumpos_P50_rats_Rn7.RData #(10)
```

+ Unpruned QTLs
```
/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/QTLs_alpha1e-04_unpruned.RData #(11)
```

+ Dir with GWAS for ASVs 
```
/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/ #(12)
```

+ Dir with GWAS for taxa 
```
/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/ #(13)
```
</details> 

<details>
<summary>Input - plot:</summary>

+ summarised QTLs - intermediate output from dataPrep
```
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/QTLs_alpha1e-04_unpruned_DGE_CE_MaE_toPlot.RData
```

+ Dir to load snps in LD
```
/users/abaud/abaud/P50_HSrats/output/pvalues_LOCO_unpruned/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/ #(12)
```

+ For “annotate” function 
```
/users/abaud/data/secondary/P50_HSrats/felipes_deblur/parsed_taxonomy.RData #(14)
```
</details>
  
Output: **porcupine\_uncollapsed\_genus2.pdf**
<br/><br/>

## Figure 5. Association between Paraprevotella and the *St6galnac1* locus on chromosome 10

#### Panel D.
Code: `5d.dataPrep_GWAS_boxplots.R` + `5d.GWAS_boxplots.R`
<details>
<summary>Input - dataPrep:</summary>

+ CLR and raw counts 
```
/users/abaud/data/secondary/P50_HSrats/felipes_deblur/full_biomt_clr_counts.RData #(6)
```
</details>

<details>
<summary>Input - plot:</summary>

+ full biomatrix -> *full_biomt* - intermediate output from dataPrep
```
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_full_biomt.RData
```

+ CLR -> *clr_counts* - intermediate output from dataPrep
```
/users/abaud/htonnele/PRJs/P50_HSrats/16S/output/felipes_deblur_clr_counts.RData
```

+ Geno positions 
```
/users/abaud/abaud/P50_HSrats/data/dosages/P50_Rn7_chr10qtl_allSNPS.raw #(15)
```

+ Metadata 
```
/users/abaud/abaud/P50_HSrats/data/metadata/metadata_augmented_16S_metabo_deblur.RData #(3)
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
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_IGE_VC.RData #(16)
```
</details>

Output: **QQplot\_pvalues\_IGE\_Helenes.pdf**
<br/><br/>

#### Panel C. 
Code: `6c.total_genetic_variance_barplot.R`
<details>
<summary>Input:</summary>

+ VD data
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/augmented_IGE_VC.RData #(16)
```

+ VD without IGE 
```
/users/abaud/abaud/P50_HSrats/output/VD/univariate/deblur_counts_uncollapsed/P50_Rn7_pruned_DGE_cageEffect_maternalEffect/all_estNste.Rdata #(7)
```
</details>

Output: **tot\_herit\_barplot.pdf**
<br/><br/>


#### Panel D (and Supp. Fig. 14). 
Code: `6d.simulations.R` <br/>
Source: `fun_prepareVD_res.R` - prepare_res() function
<details>
<summary>Input:</summary>
  
+ Results from simulations - in VD folder (MI and NY)
```
/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/VD/cage_422/setUhDGhIG/DG1_IG1/ #(17 -MI)
/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/VD/cage_654/setUhDGhIG/DG1_IG1/ #(18 -NY)

+ P50_Rn7_pruned_DGE_cageEffect_None_all_estNste.Rdata #(DGE only) 
+ P50_Rn7_pruned_DGE_IGE_cageEffect_None_all_estNste.Rdata #(wt IGE) 
```

+ Simulated values - in simulation folder (MI and NY)
```
/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/mockphenos/cage_422/setUhDGhIG/DG1_IG1/ #(19 -MI)
/users/abaud/htonnele/nf_PRJs/nf-CoreQuantGen/simulations/output/simP50/mockphenos/cage_654/setUhDGhIG/DG1_IG1/ #(20 -NY)
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
