## HNSCC Manuscript Workflow

Workflow code supporting the manuscript 'Assessing Individual Head and Neck Squamous Cell Carcinoma Patient Response to Therapy Through Integration of Functional and Genomic Data'

### Required Files

Files required for the `data` folder:

Files can be downloaded from https://github.com/biodev/hnscc_data

* hnscc_clinical.txt
* hnscc_inhibitor_data.txt
* hnscc_targetome_subset.txt
* hnscc_wes_somatic_calls.txt
* hnscc_gene_cna_calls.txt
* hnscc_rnaseq_counts.txt
* hnscc_rna_annots.txt
* hnscc_rppa.txt
* hnscc_rppa_abs.txt 
* wgcna_mm.txt
* tcga_plus_hnscc_exprs.txt.gz
* hnscc_netprop_results.txt

Distributed as part of this repo:

* Derived from the [TCGA-HNSC manuscript](https://www.nature.com/articles/nature14129) and [licensed](https://creativecommons.org/licenses/by-nc-sa/3.0/) under the terms of the corresponding manuscript.
  * data/tcga_targetable_genes.xlsx -- Genes and pathways from Figure 3
  * data/putative_drivers_tcga_gistic_regions.xlsx -- Formed from the highlighted genes in Supplemental Data S2.1 
  * data/tcga_mut_freq.xlsx -- Significant (q < 0.1) genes and frequencies from Figure 2

Also need to supply the following files in the `external` folder:

* [Reactome FIs](https://reactome.org/download-data)
  * FIsInGene_070323_with_annotations.txt
* [Genomic Data Commons (TCGA-HNSC)](https://gdc.cancer.gov/about-data/publications/hnsc_2014):
  * classification_centroid_728genes.txt
  * all_thresholded.by_genes
* [TCGA-HNSC Supplementary Data](https://www.nature.com/articles/nature14129#Sec7)
  * nature14129-s2/7.2.xlsx
  * nature14129-s2/4.1.xlsx
* [Iorio et al 2016 Supplmental Information](https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html)
  * TableS4C.xlsx
* [Genomics of Drug Sensitivity in Cancer](https://www.cancerrxgene.org/downloads/bulk_download)
  * GDSC1_fitted_dose_response_27Oct23.xlsx
  * GDSC2_fitted_dose_response_27Oct23.xlsx
  * screened_compounds_rel_8.5.csv
* [Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads)
  * mutations_summary_20230202.csv
  * cnv_gistic_20191101.csv
  * model_list_20240110.csv
* [Hallmarks](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
  * h.all.v7.5.1.symbols.gmt
* [Entrez gene info file](https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/)
  * Homo_sapiens.gene_info.gz
* [PanCancer Pathways](https://doi.org/10.1016/j.cell.2018.03.035)
  * mmc3.xlsx


### Running the workflow

This workflow uses the `targets` R package.  To run simply use the command:

```{r}
targets::tar_make()
```

### R and Package Information

R version 4.4.2

| Package              |  Version |
-----------------------|-----------
| targets              |  1.9.0  |
| clusterProfiler      |  4.14.3 |
| ComplexHeatmap       |  2.22.0 |
| ConsensusClusterPlus |  1.70.0 |
| edgeR                |  4.4.0  |  
| ggplot2              |  3.5.1  | 
| ggrepel              |  0.9.6  |
| igraph               |  2.1.1  | 
| openxlsx             | 4.2.7.1 |
| patchwork            |  1.3.0  |  
| RColorBrewer         |  1.1-3  |
| tidyverse            |  2.0.0  |
| uwot                 |  0.2.2  |
| extrafont            |  0.19   |

### License

[MIT LICENSE](LICENSE.txt)

This code was developed by Daniel Bottomly, a member of the McWeeney Lab and is protected under copyright by the Oregon Health and Science University Knight Cancer Institute, 2024.