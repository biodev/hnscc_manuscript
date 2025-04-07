library(targets)
library(extrafont)
loadfonts()

source("R/inhibitor.R")
source("R/exprs.R")
source("R/integration.R")
source("R/alterations.R")
source("R/clinical.R")
source("R/gdsc_comparison.R")

tar_option_set(seed = 252551)

if (dir.exists("outputs") == F) {
  dir.create("outputs")
}

if (dir.exists("figures") == F) {
  dir.create("figures")
}

list(

  # data files

  ## data that can be downloaded from https://biodev.github.io/HNSCC/
  tar_target(
    clin_data_file,
    "data/hnscc_clinical.txt",
    format = "file"
  ),
  tar_target(
    inhib_file,
    "data/hnscc_inhibitor_data.txt",
    format = "file"
  ),
  tar_target(
    trgome_file,
    "data/hnscc_targetome_subset.txt",
    format = "file"
  ),
  tar_target(
    mut_file,
    "data/hnscc_wes_somatic_calls.txt",
    format = "file"
  ),
  tar_target(
    cna_file,
    "data/hnscc_gene_cna_calls.txt",
    format = "file"
  ),
  tar_target(
    exprs_file,
    "data/hnscc_rnaseq_counts.txt",
    format = "file"
  ),
  tar_target(
    wgcna_mm_file,
    "data/wgcna_mm.txt",
    format = "file"
  ),
  tar_target(
    hnscc_rna_annot_file,
    "data/hnscc_rna_annots.txt",
    format = "file"
  ),
  tar_target(
    cell_rppa_file,
    "data/hnscc_rppa.txt",
    format = "file"
  ),
  tar_target(
    rppa_annots_file,
    "data/hnscc_rppa_abs.txt",
    format = "file"
  ),

  ## derived files from TCGA-HNSC (https://doi.org:10.1038/nature14129) but distributed as part of this repo

  tar_target(
    tcga_actionable_gene_file,
    "data/tcga_targetable_genes.xlsx",
    format = "file"
  ),
  tar_target(
    tcga_driver_cnas_file,
    "data/putative_drivers_tcga_gistic_regions.xlsx"
  ),
  tar_target(
    tcga_mut_freq_file,
    "data/tcga_mut_freq.xlsx"
  ),

  ## combined and batch corrected version of TGCA-HNSC (from GDC) and HNSCC

  tar_target(
    tcga_hnscc_exprs_file,
    "data/tcga_plus_hnscc_exprs.txt.gz",
    format = "file"
  ),

  # preprocessing

  ## gene-lists

  tar_target(
    tcga_actionable,
    readxl::read_excel(tcga_actionable_gene_file)
  ),
  tar_target(
    tcga_driver_cnas,
    readxl::read_excel(tcga_driver_cnas_file, skip = 1)
  ),
  tar_target(
    tcga_driver_muts,
    readxl::read_excel("external/nature14129-s2/4.1.xlsx",
      range = "B4:B14",
      sheet = "Mutsig CV", col_names = "Gene"
    )
  ),

  # clinical

  tar_target(
    clin_data,
    suppressWarnings(readr::read_delim(clin_data_file, col_types = c("cncccccnccccccccc")))
  ),
  tar_target(
    clinical_summary,
    make.clin.summary.table(clin_data),
    packages = c("tidyverse", "openxlsx"),
    format = "file"
  ),

  # mutations

  tar_target(
    mut_data,
    readr::read_delim(mut_file)
  ),
  tar_target(
    mut_driver_oncoprint,
    make.mut.oncoprint(mut_data,
      tcga_mut_freq_file,                 
      clin = clin_data,
      gene.list = tcga_driver_muts$Gene
    ),
    packages = c("ComplexHeatmap", "tidyverse"),
    format = "file"
  ),
  tar_target(
    tcga_cna_data,
    process.gistic.calls("external/all_thresholded.by_genes"),
    packages = "tidyverse"
  ),
  tar_target(
    cna_data,
    readr::read_delim(cna_file)
  ),
  tar_target(
    cna_driver_oncoprint,
    make.cna.oncoprint(cna_data, clin_data, tcga_cna_data,
      gene.list = tcga_driver_cnas$Gene
    ),
    packages = c("ComplexHeatmap", "tidyverse"),
    format = "file"
  ),

  # expression analysis

  tar_target(
    hnscc_rna_annot,
    readr::read_delim(hnscc_rna_annot_file, col_types = "cccc")
  ),
  tar_target(
    tcga_hnscc_exprs,
    exprs.mat.from.file(tcga_hnscc_exprs_file),
    packages = "tidyverse"
  ),
  tar_target(
    centroids,
    as.matrix(read.table("external/classification_centroid_728genes.txt",
      stringsAsFactors = F
    ))
  ),
  tar_target(
    tcga_subtype_annot,
    readxl::read_excel("external/nature14129-s2/7.2.xlsx")
  ),
  tar_target(
    subtype_calls,
    classify.exprs.subtype(tcga_hnscc_exprs, centroids, tcga_subtype_annot),
    packages = c("tidyverse")
  ),
  tar_target(
    expression_umap,
    make.expression.plot(tcga_hnscc_exprs, subtype_calls, centroids),
    packages = c("tidyverse", "uwot"),
    format = "file"
  ),

  # wgcna analysis

  tar_target(
    wgcna_mm,
    readr::read_delim(wgcna_mm_file)
  ),
  tar_target(
    msigdb_hallmarks,
    read.msigdb("external/h.all.v7.5.1.symbols.gmt"),
    packages = c("clusterProfiler", "tidyverse")
  ),
  tar_target(
    tcga_hnscc_mes,
    form.mes.tcga.hnscc(wgcna_mm, tcga_hnscc_exprs, hnscc_rna_annot),
    packages = c("tidyverse")
  ),
  tar_target(
    hnscc_me_overlap_plot,
    plot.tcga.hnscc.me.overlap(tcga_hnscc_mes),
    packages = c("tidyverse")
  ),
  tar_target(
    mod_hallmark_plot,
    enrich.mods.hallmarks.plot(wgcna_mm, msigdb_hallmarks),
    packages = c("tidyverse", "clusterProfiler"),
    format = "file"
  ),
  tar_target(
    me_clin_cor_plot,
    correlate.me.with.clinical.plot(clin_data, tcga_hnscc_mes, subtype_calls, wgcna_mm),
    packages = c("tidyverse", "patchwork"),
    format = "file"
  ),

  # Inhibitor summaries

  tar_target(
    inhib,
    readr::read_delim(inhib_file, col_types = "ccccnnccnn")
  ),
  tar_target(
    inhib_list,
    make.inhib.table(inhib),
    packages = c("tidyverse", "openxlsx"),
    format = "file"
  ),
  tar_target(
    trgome,
    readr::read_delim(trgome_file)
  ),
  tar_target(
    pancan_pathways,
    form.pancan.pathways("external/mmc3.xlsx", "external/Homo_sapiens.gene_info.gz"),
    packages = c("tidyverse", "openxlsx")
  ),
  #https://reactome.org/download/current/ReactomePathways.gmt.zip
  tar_target(
    reactome_file,
    "external/ReactomePathways.gmt",
    format = "file"
  ),
  tar_target(
    reactome_pathways,
    read.reactome.gmt(reactome_file),
    packages = "tidyverse"
  ),
  tar_target(
    drug_targ_pathways_file,
    output.drug.target.pathways(trgome, pancan_pathways, reactome_pathways, reactome_hier),
    packages = c("tidyverse", "openxlsx", "igraph")
  ),
  tar_target(
    drug_pathways,
    annotate.drugs.pathways(inhib, trgome, pancan_pathways),
    packages = "tidyverse"
  ),
  tar_target(
    inhib_mat_list,
    drug.by.pt.mat(inhib),
    packages = "tidyverse"
  ),
  tar_target(
    trgome_mat,
    gene.by.drug.mat(trgome),
    packages = "tidyverse"
  ),
  tar_target(
    pathway_overlap_plot,
    overlap.tcga.path.plots(trgome, pancan_pathways, tcga_actionable, drug_pathways),
    packages = c("tidyverse", "RColorBrewer", "patchwork"),
    format = "file"
  ),
  tar_target(
    comb_ratio_plot,
    combination.ratio.plot(inhib, drug_pathways),
    packages = c("tidyverse"),
    format = "file"
  ),
  tar_target(
    gene_scores,
    compute.gene.scores(inhib_mat_list$`single-agent`, trgome_mat),
    packages = "tidyverse"
  ),

  # Inhibitor comparison with GDSC

  tar_target(
    gdsc_hnscc_drug_map,
    map.gdsc.hnscc.drugs(inhib, "external/screened_compounds_rel_8.5.csv"),
    packages = c("tidyverse")
  ),
  tar_target(
    gdsc_drug_data,
    get.gdsc.drug.data(gdsc_hnscc_drug_map, "external/GDSC1_fitted_dose_response_27Oct23.xlsx", "external/GDSC2_fitted_dose_response_27Oct23.xlsx"),
    packages = c("tidyverse")
  ),
  tar_target(
    hnsc_cls,
    limit.cell.lines.to.hnsc.drugs("external/model_list_20240110.csv", gdsc_drug_data),
    packages = c("tidyverse")
  ),
  tar_target(
    gdsc_assocs,
    proc.gdsc.assoc.genes("external/TableS4C.xlsx", gdsc_hnscc_drug_map),
    packages = "tidyverse"
  ),
  tar_target(
    mut_mat,
    make.mutation.mat(mut_data,
      gene.list = filter(gdsc_assocs, type == "Mut")$Gene
    ),
    packages = "tidyverse"
  ),
  tar_target(
    gdsc_mut_data,
    readr::read_delim("external/mutations_summary_20230202.csv") |>
      dplyr::rename(symbol = gene_symbol, lab_id = model_id) |>
      dplyr::filter(lab_id %in% hnsc_cls$model_id)
  ),
  tar_target(
    gdsc_mut_mat,
    make.mutation.mat(gdsc_mut_data,
      gene.list = dplyr::filter(gdsc_assocs, type == "Mut")$Gene
    ),
    packages = "tidyverse"
  ),
  tar_target(
    gdsc_cna_mat,
    make.gdsc.cna.mat(
      "external/cnv_gistic_20191101.csv",
      hnsc_cls$model_id,
      dplyr::filter(gdsc_assocs, type != "Mut")$Gene
    ),
    packages = "tidyverse"
  ),
  tar_target(
    cna_mat,
    make.cna.mat(cna_data,
      gene.list = dplyr::filter(gdsc_assocs, type != "Mut")$Gene
    ),
    packages = "tidyverse"
  ),
  tar_target(
    hnscc_cl_matches,
    find.cl.matches(mut_mat, gdsc_mut_mat, cna_mat, gdsc_cna_mat, gdsc_assocs),
    packages = c("tidyverse")
  ),
  tar_target(
    hnscc_cl_match_heatmap,
    plot.cl.matches(hnsc_cls, mut_mat, gdsc_mut_mat, cna_mat, gdsc_cna_mat, gdsc_assocs, inhib),
    packages = c("tidyverse", "ComplexHeatmap", "ConsensusClusterPlus"),
    format = "file"
  ),
  tar_target(
    hnscc_cl_inh_dist_summary_plot,
    plot.cor.dist.hnscc.cl(hnscc_cl_matches, inhib, gdsc_drug_data, hnsc_cls),
    packages = c("tidyverse", "ggrepel"),
    format = "file"
  ),
  tar_target(
    hnscc_cl_inh_plot,
    plot.hnscc.cl.inh(hnscc_cl_matches, inhib, gdsc_drug_data, hnsc_cls),
    packages = "tidyverse",
    format = "file"
  ),

  # Drug/alteration integrative analyses

  ## Network propagation

  ### Output for use in `run_rewiring.Rmd` which is run on HPC
  tar_target(
    gs_tbl_out,
    {
      readr::write_delim(gene_scores, file = "outputs/hnscc_inh_gene_scores.txt", delim = "\t", col_names = T)
      "outputs/hnscc_inh_gene_scores.txt"
    },
    format = "file"
  ),

  ### Make the original network for use below and with `run_rewiring.Rmd`
  tar_target(
    reactome_fi_lcc_graph,
    make.hc.reactome.w.lcc("external/FIsInGene_070323_with_annotations.txt"),
    packages = c("igraph")
  ),
  tar_target(
    reactome_fi_lcc_file,
    {
      orig.graph <- reactome_fi_lcc_graph
      save(orig.graph, file = "outputs/reactome_fi_graph.RData")
    },
    format = "file"
  ),

  ### Define and read back in results from run_rewiring.Rmd
  tar_target(
    netprop_results,
    "data/hnscc_netprop_results.txt",
    format = "file"
  ),
  tar_target(
    net_prop_results,
    readr::read_delim(netprop_results, col_types = "cclnnnnn")
  ),
  tar_target(
    net_prop_plus_gs,
    dplyr::full_join(gene_scores, net_prop_results, by = c("Gene", "lab_id"))
  ),

  ### Combine with alteration data

  tar_target(
    rppa_mat,
    rppa.mat.from.file(cell_rppa_file),
    packages = "tidyverse"
  ),
  tar_target(
    cell_rppa_z,
    mat.to.z.tbl(rppa_mat, "exprs_z"),
    packages = "tidyverse"
  ),
  tar_target(
    rppa_annots,
    readr::read_delim(rppa_annots_file)
  ),
  tar_target(
    cell_exprs,
    form.cell.exprs.mat(exprs_file, hnscc_rna_annot),
    packages = c("tidyverse", "edgeR")
  ),
  tar_target(
    cell_exprs_z,
    mat.to.z.tbl(cell_exprs, "exprs_z"),
    packages = "tidyverse"
  ),
  tar_target(
    net_prop_alt_tbl,
    combine.gs.mut.cna.exprs.rppa(
      net_prop_plus_gs, mut_data, cna_data,
      cell_exprs_z, cell_rppa_z, rppa_annots
    ),
    packages = "tidyverse"
  ),

  ## Network propagation barplots

  tar_target(
    net_prop_barplots,
    plot.network.prop.bars(net_prop_alt_tbl),
    packages = c("tidyverse", "ggplot2"),
    format = "file"
  ),
  
  tar_target(
    example_net_prop_barplots,
    plot.network.prop.bars(dplyr::filter(net_prop_alt_tbl, lab_id == "10058"), file.pref="pt10058"),
    packages = c("tidyverse", "ggplot2"),
    format = "file"
  ),

  ## Gene Score summary

  tar_target(
    gs_summary_plot,
    plot.gs.summary(net_prop_alt_tbl, use.genes = tcga_actionable$gene_name),
    packages = c("tidyverse", "ggplot2"),
    format = "file"
  ),

  ## Response Cards

  tar_target(
    example_response_cards_plot,
    plot.all.response.cards(
      dplyr::filter(net_prop_alt_tbl, lab_id == "10058"),
      trgome_mat, inhib_mat_list$`single-agent`,
      unique(dplyr::select(inhib, inhibitor, plot_name)),
      file.pref="pt10058", file.width=5, file.height=4
    ),
    packages = c("ComplexHeatmap", "tidyverse"),
    format = "file"
  ),
  
  tar_target(
    all_response_cards_plot,
    plot.all.response.cards(
      net_prop_alt_tbl,
      trgome_mat, inhib_mat_list$`single-agent`,
      unique(dplyr::select(inhib, inhibitor, plot_name))
    ),
    packages = c("ComplexHeatmap", "tidyverse"),
    format = "file"
  ),

  ## Cytoscape network plot

  tar_target(
    graph_for_cytoscape,
    make.graph.figure.by.pt(reactome_fi_lcc_graph, net_prop_alt_tbl,
      pt.id = "10058", top.n = 7
    ),
    packages = c("tidyverse", "igraph"),
    format = "file"
  )
)
