#' Read MSigDB GMT file and return a `tibble`
#'
#' @param mdb.file A GMT file
#' @returns A `tibble` containing columns for 'term' (Hallmark) and 'gene'.
read.msigdb <- function(mdb.file) {
  as_tibble(read.gmt(mdb.file))
}

#' Read a file containing expression data and convert to a matrix
#'
#' @param exprs.file A file containing rows for each gene and columns indicating:
#'   'gene_id' -- Gene identifier
#'   'gene_name' -- Gene symbol
#'   patient samples -- Expression values for one or more patient samples
#' @returns A gene x patient-sample matrix containing expression values
exprs.mat.from.file <- function(exprs.file) {
  exprs.tbl <- readr::read_delim(exprs.file)

  stopifnot(sum(duplicated(exprs.tbl$gene_name)) == 0)

  exprs.mat <- as.matrix(select(exprs.tbl, -gene_id, -gene_name))
  rownames(exprs.mat) <- exprs.tbl$gene_name

  exprs.mat
}

.process.filter.exprs.dge <- function(filt.dge, gene.info) {
  filt.dge <- calcNormFactors(filt.dge)

  logCPM <- cpm(filt.dge, log = TRUE, prior.count = 3)

  gene.avgs <- as_tibble(rowMeans(logCPM), rownames = "gene_id")

  gene.info <- inner_join(gene.info, gene.avgs, by = "gene_id") %>%
    arrange(desc(value)) %>%
    filter(!duplicated(gene_name))

  logCPM <- logCPM[gene.info$gene_id, ]
  rownames(logCPM) <- gene.info$gene_name

  logCPM
}

#' Read in a count file and process to a log2 CPM matrix limited to the cell culture data
#'
#' The count file is first read-in and converted to a matrix limited to the cell culture samples.
#' This matrix is then filtered by [filterByExpr()], TMM normalized and converted to log2 CPM using {edgeR}.
#' If there are multiple gene identifiers corresponding to a given gene symbol, the one with the highest
#' average log2 CPM is chosen (arbitrarily).
#'
#' @param exprs.file File containing the gene expression counts per sample.  The first
#'   two columns should be 'gene_id' and 'gene_name'
#' @param rna.annot A `tibble` of patient samples with expression data containing at least
#'   columns for:
#'   `lab_id` -- Patient sample identifier
#'   `sample_type` -- Either 'Tissue' or 'Cells'
#'   `sample_status` -- Either 'Tumor' or 'Normal'
#'   `barcode` -- Patient sample identifier present in the columns of `exprs.file`
#'
#' @returns A gene x patient-sample matrix of log2 CPM expression values
form.cell.exprs.mat <- function(exprs.file, rna.annot) {
  exprs.tbl <- readr::read_delim(exprs.file)

  gene.info <- exprs.tbl %>%
    select(gene_id, gene_name)

  exprs.mat <- as.matrix(select(exprs.tbl, -gene_id, -gene_name))
  rownames(exprs.mat) <- exprs.tbl$gene_id

  cell.samps <- filter(rna.annot, sample_type == "Cells" & sample_status == "Tumor")

  exprs.mat <- exprs.mat[, cell.samps$barcode]
  colnames(exprs.mat) <- cell.samps$lab_id

  dge <- DGEList(exprs.mat)

  # already known all samples are part of the same group
  keep.genes <- suppressWarnings(filterByExpr(dge))

  filt.dge <- dge[keep.genes, , keep.lib.sizes = FALSE]

  .process.filter.exprs.dge(filt.dge, gene.info)
}

#' Center/Scale an numeric matrix and convert to a `tibble`
#'
#' @param mat Gene x Patient sample numeric matrix
#' @param value.name Column name for the values in the resulting `tibble`
#' @returns A `tibble` with the following columns:
#'   `gene_name` -- Gene Symbol
#'   `lab_id` -- Patient sample identifier
#'   Column name of values as specified in value.name
mat.to.z.tbl <- function(mat, value.name = "exprs") {
  as_tibble(t(scale(t(mat))), rownames = "gene_name") %>%
    pivot_longer(cols = -gene_name, names_to = "lab_id", values_to = value.name)
}

#' Read in RPPA data and convert to a matrix
#'
#' Similar to [exprs.mat.from.file()] but involving RPPA data.  Instead of genes,
#' rows are named relative to antibodies.
#'
#' @param rppa.file File containing the RPPA data with antibodies in the rows and
#'   samples in the columns.  The first column however should be named 'Antibody Name'
#' @returns A antibody x patient sample matrix containing the RPPA values.
#' @note Adding RPPA here since there is relatively little of it that is different than exprs
rppa.mat.from.file <- function(rppa.file) {
  rppa.tbl <- readr::read_delim(rppa.file)

  rppa.mat <- as.matrix(select(rppa.tbl, -1))
  rownames(rppa.mat) <- rppa.tbl$`Antibody Name`

  rppa.mat
}

#' Assign subtypes based on correlation to centroids
#'
#' This approach is adapted from Ding et al 2019 (1).  Given a matrix of gene x subtype
#' centroids, calculate the specified correlation relative to the corresponding
#' expression values for each patient sample.  For each sample, compare the two highest correlation
#' values and assign the highest subtype if the difference between the correlations is
#' >= 0.1, other consider it an indeterminate call (IND).
#'
#' @param common.genes A character vector of genes in common between the expression
#'   and centroid data
#' @param exprs A gene x sample numeric expression matrix
#' @param centroids A gene x subtype numeric matrix of centroids
#' @param cor.method Correlation method, see [cor()] for available options.
#' @returns A `tibble` containing the following columns:
#'   `file_id` -- Sample identifier from `exprs`
#'   `subtype_best` -- Subtype with the best correlation
#'   `subtype_second_best` -- Subtype with the second best correlation
#'   `scor_best` -- Best correlation value
#'   `scor_second_best` -- Second best correlation value
#'   `call` -- Final call replacing the best subtypes with IND if the top two correlations
#'     were <= 0.1 apart.
#' @references (1) Ding YC, Steele L, Warden C, Wilczynski S, Mortimer J, Yuan Y, Neuhausen SL.
#'   Molecular subtypes of triple-negative breast cancer in women of different race and ethnicity.
#'   Oncotarget. 2019 Jan 4;10(2):198-208. doi: 10.18632/oncotarget.26559.
assign.subtype.cor <- function(common.genes, exprs, centroids, cor.method = "pearson") {
  cor.tbl <- cor(exprs[common.genes, ], centroids[common.genes, ], method = cor.method) %>%
    as_tibble(rownames = "file_id") %>%
    pivot_longer(cols = -file_id, names_to = "subtype", values_to = "scor")

  best.cor <- cor.tbl %>%
    group_by(file_id) %>%
    slice_max(order_by = scor, n = 2) %>%
    ungroup() %>%
    group_by(file_id) %>%
    arrange(desc(scor)) %>%
    mutate(idx = c("best", "second_best")) %>%
    pivot_wider(names_from = idx, values_from = c("subtype", "scor")) %>%
    ungroup()

  best.cor %>%
    mutate(call = ifelse(((scor_best - scor_second_best) >= .1) & (scor_best > 0), subtype_best, "IND"))
}

#' Call HNSCC subtypes
#'
#' Wrapper functions that matches the centroid and expression genes, median centers
#' the expression data and calls [assign.subtype.cor()] using pearson correlation.
#' This is an approach that roughly mirrors the methods described in Walter et al 2013 (1)
#' and performs the best in comparison to other potential approaches.  The resulting subtypes
#' are then combined with the supplied subtype annotation.
#'
#' @param logRPKM A gene x sample numeric matrix of expression values.
#' @param centroids A gene x subtype numeric matrix of centroids
#' @param subtype.annot A `tibble` with at least a column for 'Barcode' corresponding
#'   to the columns in `exprs` if '.' characters are converted to '-' and -01A' is added to it.
#' @references Walter, V. et al. Molecular Subtypes in Head and Neck Cancer
#'   Exhibit Distinct Patterns of Chromosomal Gain and Loss of Canonical Cancer Genes.
#'   PLOS ONE 8, e56823 (2013). https://doi.org:10.1371/journal.pone.0056823
#' @returns A `tibble` containing the output of [assign.subtype.cor()] merged with
#'   `subtype.annot`
classify.exprs.subtype <- function(logRPKM, centroids, subtype.annot) {
  # adjust labels to match expression matrix
  subtype.annot <- mutate(subtype.annot, file_id = paste0(gsub("\\.", "-", Barcode), "-01A"))

  common.genes <- intersect(rownames(centroids), rownames(logRPKM))
  stopifnot(length(common.genes) == 643)

  cent.exprs <- logRPKM - matrixStats::rowMedians(logRPKM)

  cent.pear <- assign.subtype.cor(common.genes, cent.exprs, centroids, cor.method = "pearson")

  subtype.calls <- left_join(cent.pear, subtype.annot, by = "file_id")

  subtype.calls
}

#' Create a `tibble` of module eigengene values
#'
#' Given a module map, first center and scale the data corresponding to TCGA
#' (not part of HNSCC) and compute their PCA scores saving the per-gene mean and standard
#' deviations.  Apply centering and scaling for the HNSCC samples using TCGA's means and
#' standard deviations and 'predict' the resulting PCA scores.  Combine the results
#' into a `tibble`
#'
#' @param mm A module map containg columns for 'module' and 'gene_name'
#' @param logRPKM A gene x sample numeric expression matrix
#' @param rna.annot A `tibble` of patient samples with expression data containing at least
#'   columns for:
#'   `lab_id` -- Patient sample identifier
#'   `sample_type` -- Either 'Tissue' or 'Cells'
#'   `sample_status` -- Either 'Tumor' or 'Normal'
#'   `barcode` -- Patient sample identifier present in the columns of `logRPKM`
#' @returns A `tibble` with the following columns:
#'   `lab_id` -- Patient sample identifier
#'   `PC1` -- PC1 score
#'   `PC2` -- PC2 score
#'   `module` -- Corresponding module
#'   `cohort` -- Cohort (TCGA or HNSCC) that the patient sample belongs to
form.mes.tcga.hnscc <- function(mm, logRPKM, rna.annot) {
  mm <- group_by(mm, module)

  hnscc.samples <- filter(rna.annot, sample_type == "Tissue" & sample_status == "Tumor")

  tcga.prs <- lapply(group_split(mm), function(x) {
    z.exprs <- scale(t(logRPKM[x$gene_name, setdiff(colnames(logRPKM), hnscc.samples$lab_id)]))

    av.expr <- rowMeans(z.exprs)

    tmp.prs <- prcomp(z.exprs, center = F, scale. = F)

    tmp.me <- tmp.prs$x[, 1:2]

    # per the typical use of PC1 as module eigengene
    if (sign(cor(tmp.me[names(av.expr), "PC1"], av.expr)) < 0) {
      tmp.me[, "PC1"] <- -tmp.me[, "PC1"]
    }

    me.tbl <- as_tibble(tmp.me, rownames = "lab_id") %>%
      mutate(module = x$module[1])

    list(mes = me.tbl, pca = tmp.prs, centers = attr(z.exprs, "scaled:center"), scales = attr(z.exprs, "scaled:scale"))
  })

  hnscc.prs <- bind_rows(lapply(tcga.prs, function(x) {
    z.exprs <- scale(t(logRPKM[, hnscc.samples$lab_id])[, names(x$centers)], center = x$centers, scale = x$scales)

    av.expr <- rowMeans(z.exprs)

    tmp.me <- z.exprs[, rownames(x$pca$rotation)] %*% x$pca$rotation

    if (sign(cor(tmp.me[names(av.expr), "PC1"], av.expr)) < 0) {
      tmp.me[, "PC1"] <- -tmp.me[, "PC1"]
    }

    me.tbl <- as_tibble(tmp.me[, 1:2], rownames = "lab_id") %>%
      mutate(module = x$mes$module[1])

    me.tbl
  }))


  tcga.prs <- bind_rows(lapply(tcga.prs, function(x) x$mes))

  all.eigs <- bind_rows(
    mutate(tcga.prs, cohort = "TCGA"),
    mutate(hnscc.prs, cohort = "HNSCC")
  )


  # remove grey module

  all.eigs <- filter(all.eigs, module != "Mod0 (grey)")

  all.eigs
}

#' Plot PC1 and PC2 for each module for TCGA overlaying HNSCC
#'
#' @param all.eigs A `tibble` of PC1 and 2 values for each module and cohort as
#'   produced by [form.mes.tcga.hnscc()]
#' @returns A PDF: 'figures/hnscc_tcga_me_comparison.pdf'
plot.tcga.hnscc.me.overlap <- function(all.eigs) {
  # visualize PC1 and 2 for validation and test
  p1 <- ggplot(data = all.eigs, mapping = aes(x = PC1, y = PC2, fill = cohort, alpha = cohort)) +
    geom_point(shape = 21) +
    scale_fill_manual(values = c("TCGA" = "grey", "HNSCC" = "purple")) +
    scale_alpha_manual(values = c("TCGA" = .25, "HNSCC" = 1)) +
    facet_wrap(~module, scales = "free") +
    theme_bw()

  ggsave(p1, file = "figures/hnscc_tcga_me_comparison.pdf", width = 11, height = 7)

  "figures/hnscc_tcga_me_comparison.pdf"
}

#' Compute and plot enrichment of Hallmarks in the modules
#'
#' First compute enrichment using the set of genes from any module as the universe.
#' Then generate a dotplot of the significant modules/Hallmarks
#'
#' @param mm -- A `tibble` containing the module map with at least columns for
#'   'module', 'gene_name' and 'label' a integer value which is used to order the
#'   modules.
#' @param gmt.tbl A `tibble` containing columns for 'term' (Hallmark) and 'gene' as
#'   returned by [read.msig]
#' @returns A PDF: 'figures/mm_hallmark_fig.pdf'
enrich.mods.hallmarks.plot <- function(mm, gmt.tbl) {
  # compute enrichment of hallmarks for each module

  mm <- group_by(mm, module)

  mod.enrich.list <- lapply(group_split(mm), function(x) {
    if (length(intersect(x$gene_name, gmt.tbl$gene)) == 0) {
      NULL
    } else {
      tmp <- clusterProfiler::enricher(
        gene = x$gene_name,
        pvalueCutoff = 1,
        pAdjustMethod = "none",
        universe = mm$gene_name,
        minGSSize = 5,
        maxGSSize = 1000,
        qvalueCutoff = 1,
        TERM2GENE = as.data.frame(gmt.tbl),
        TERM2NAME = NA
      )

      tmp.dt <- as_tibble(tmp)

      tmp.dt
    }
  })

  names(mod.enrich.list) <- group_keys(mm)$module

  mod.enrich <- bind_rows(mod.enrich.list, .id = "module")

  mod.enrich <- filter(mod.enrich, module != "Mod0 (grey)")

  sig.mods <- filter(mod.enrich, p.adjust < .05)

  mod.enrich <- filter(mod.enrich, module %in% sig.mods$module)

  mod.enrich <- mutate(mod.enrich, num_genes = as.numeric(sapply(strsplit(GeneRatio, "\\/"), function(x) x[1])))

  hall.sum <- summarize(mod.enrich, minp = min(pvalue), .by = ID) %>%
    arrange(desc(minp))

  mod.enrich <- mutate(mod.enrich, id_fac = factor(ID, levels = hall.sum$ID, ordered = T, labels = sub("HALLMARK_", "", hall.sum$ID)))

  mod.order <- mm %>%
    select(module, labels) %>%
    unique() %>%
    arrange(labels)

  mod.enrich <- mutate(mod.enrich, mod_fac = factor(module, levels = mod.order$module, ordered = T))

  mod.cols <- summarize(mod.enrich, n(), .by = module) %>%
    mutate(cols = str_match(module, "\\((\\w+)\\)")[, 2])

  mod.enrich <- mutate(mod.enrich, mod_fill = ifelse(p.adjust < .05, module, "none"))

  p1 <- ggplot(data = filter(mod.enrich, pvalue < .05), mapping = aes(x = mod_fac, y = id_fac, size = num_genes, fill = mod_fill)) +
    geom_point(shape = 21) +
    scale_fill_manual(values = c(setNames(mod.cols$cols, mod.cols$module), "none" = "white"), guide = "none") +
    scale_size_continuous(name = "Number of Genes") +
    theme_bw() +
    ylab("") +
    xlab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  ggsave(p1, file = "figures/mm_hallmark_fig.pdf", width = 7, height = 4.3)

  "figures/mm_hallmark_fig.pdf"
}

#' Create module association plots
#'
#' The top plot contains the median eigengenes for the HNSCC subtypes across modules.
#' in TCGA.  The bottom plot compares the T-statistics formed by comparing the
#' specified clinical groups.  White stars indicate significance.
#'
#' @param clin A `tibble` containing clinical information per patient.  At least
#'   it should contain: 'lab_id', Age', 'Gender' 'Disease_Site', 'Pack_Years', 'Alcohol_Use'
#'   and 'One_Year_RFS'
#' @param all.prs A `tibble` of module eigengene values as produced by [form.mes.tcga.hnscc()]
#' @param all.subtypes A `tibble` of called subtypes from TCGA as produced by [classify.exprs.subtype()]
#' @param mm A `tibble` containing the module map with at least columns for 'module' and
#'   numeric 'labels'.
#' @returns A PDF: 'figures/clinical_vs_mes_hnscc.pdf'
correlate.me.with.clinical.plot <- function(clin, all.prs, all.subtypes, mm) {
  mod.order <- mm %>%
    select(module, labels) %>%
    unique() %>%
    arrange(labels)

  # Make the top plot

  tcga.prs <- filter(all.prs, cohort == "TCGA") %>%
    inner_join(all.subtypes, by = c("lab_id" = "file_id"))

  med.pc1 <- tcga.prs %>%
    summarize(med_exprs = median(PC1), .by = c("module", "subtype_best")) %>%
    mutate(mod_fac = factor(module, levels = mod.order$module, ordered = T))

  p1 <- ggplot(data = med.pc1, mapping = aes(x = mod_fac, y = subtype_best, fill = med_exprs)) +
    geom_tile() +
    geom_hline(yintercept = seq(.5, length(unique(med.pc1$subtype_best)) + .5, by = 1)) +
    geom_vline(xintercept = seq(.5, length(unique(med.pc1$mod_fac)) + .5, by = 1)) +
    scale_fill_distiller(type = "div", palette = "PuOr", name = "Med. Eigengene") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      axis.line = element_blank()
    ) +
    ggtitle("TCGA")

  hnscc.prs <- filter(all.prs, cohort == "HNSCC")

  rna.clin <- filter(clin, lab_id %in% hnscc.prs$lab_id) %>%
    mutate(
      Age_fac = cut(Age, breaks = c(40, 62, 90)),
      Gender_fac = factor(Gender, levels = c("Female", "Male")),
      Larynx_vs_Oral_Cavity_fac = factor(Disease_Site, levels = c("Larynx", "Oral Cavity")),
      Oropharynx_vs_Oral_Cavity_fac = factor(Disease_Site, levels = c("Oropharynx", "Oral Cavity")),
      Smoker_Pack_Years_fac = droplevels(cut(Pack_Years, breaks = c(0, 15, 59, 100), include.lowest = T), "(15,59]"),
      Alcohol_Use_fac = factor(Alcohol_Use, levels = c("Heavy", "Minimal")),
      One_Year_RFS_fac = factor(One_Year_RFS, levels = c("Yes", "No"))
    )

  pr.list <- group_split(hnscc.prs, .by = module)

  sig.res <- lapply(grep("_fac", names(rna.clin), value = T), function(x) {
    tmp.rna.clin <- filter(rna.clin, is.na(.data[[x]]) == F)

    lapply(pr.list, function(mod) {
      tmp.mod <- inner_join(tmp.rna.clin, mod, by = "lab_id")
      tmp.form <- as.formula(paste("PC1~", x))
      broom::tidy(t.test(tmp.form, data = tmp.mod)) %>%
        mutate(module = mod$.by[1])
    }) %>%
      bind_rows() %>%
      mutate(
        comparison = sub("_fac", "", x),
        direction = paste(levels(tmp.rna.clin[[x]]), collapse = "-")
      )
  }) %>%
    bind_rows()

  sig.res <- mutate(sig.res,
    star_levels = case_when(
      p.value < .001 ~ "***",
      p.value < .01 ~ "**",
      p.value < .05 ~ "*",
      .default = ""
    ),
    comparison = gsub("_", " ", comparison),
    comp_fac = factor(comparison, levels = c("Age", "Gender", "Alcohol Use", "Smoker Pack Years", "Larynx vs Oral Cavity", "Oropharynx vs Oral Cavity", "One Year RFS"), ordered = T),
    mod_fac = factor(module, levels = mod.order$module, ordered = T)
  )

  p2 <- ggplot(data = sig.res, mapping = aes(x = mod_fac, y = comp_fac, fill = statistic)) +
    geom_tile() +
    geom_hline(yintercept = seq(.5, length(unique(sig.res$comparison)) + .5, by = 1)) +
    geom_vline(xintercept = seq(.5, length(unique(sig.res$module)) + .5, by = 1)) +
    geom_text(mapping = aes(label = star_levels), color = "white", size = 8, vjust = .75) +
    scale_fill_distiller(type = "div", palette = "RdGy", name = "T Statistic") +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      axis.line = element_blank()
    ) +
    ggtitle("HNSCC")


  full.p <- (p1 / p2) + plot_layout(
    design =
      "1
    2
    2"
  )

  full.p


  ggsave(full.p, file = "figures/clinical_vs_mes_hnscc.pdf", width = 6.3, height = 5)

  "figures/clinical_vs_mes_hnscc.pdf"
}

#' Produce a UMAP of TCGA subtypes based on expression with HNSCC overlaid
#'
#' The UMAP contains one dot for each TCGA patient sample colored by subtype
#' with HNSCC samples added and indicated by a solid dot.  Transparency of the
#' samples is increased if it was considered to be indeterminate (IND) per
#' [classify.exprs.subtype()].
#'
#' @param exprs A gene x patient sample matrix containing both TCGA and HNSCC
#' @param subtype A `tibble` of called subtypes from TCGA as produced by [classify.exprs.subtype()]
#' @param centroid.mat A gene x subtype numeric matrix of centroids (used for its rownames)
#' @returns A PDF: 'figures/expression_umap.pdf'
make.expression.plot <- function(exprs, subtypes, centroid.mat) {
  common.genes <- intersect(rownames(centroid.mat), rownames(exprs))

  tcga.z <- scale(t(exprs[common.genes, grepl("^\\d{5}$", colnames(exprs)) == F]))

  # center/scale HNSCC relative to mean/sd of TGCA
  hnscc.z <- scale(t(exprs[common.genes, grepl("^\\d{5}$", colnames(exprs)) == T])[, names(attr(tcga.z, "scaled:center"))],
    center = attr(tcga.z, "scaled:center"), scale = attr(tcga.z, "scaled:scale")
  )

  comb.mat <- rbind(
    tcga.z,
    hnscc.z
  )

  pc <- prcomp(comb.mat, center = F, scale. = F)

  # How many pcs does it take to explain 95% variation
  comp.idx <- 1:min(which(cumsum((pc$sdev^2) / sum(pc$sdev^2)) >= .95))

  tmp.lo.graph <- umap(X = pc$x[, comp.idx], n_neighbors = 30, min_dist = .5)

  colnames(tmp.lo.graph) <- c("Dim1", "Dim2")

  graph.dims <- as_tibble(tmp.lo.graph, rownames = "file_id")

  graph.dims <- mutate(graph.dims, is_hnscc = grepl("^\\d{5}$", file_id))

  graph.dims <- inner_join(graph.dims, subtypes, by = "file_id")

  p1 <- ggplot(data = graph.dims, mapping = aes(x = Dim1, y = Dim2)) +
    geom_point(size = 3, mapping = aes(fill = subtype_best, alpha = call != "IND"), shape=21) +
    geom_point(data = filter(graph.dims, is_hnscc == T), size = 3, color = "black") +
    scale_fill_discrete(name = "Subtype") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .25), guide = "none") +
    theme_bw()

  ggsave(p1, file = "figures/expression_umap.pdf", width = 4, height = 3)

  "figures/expression_umap.pdf"
}
