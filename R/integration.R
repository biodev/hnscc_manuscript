#' Combine Gene Score results with other alterations
#'
#' First produces a `tibble` with combined gene score and alterations.
#' The alterations are then further categorized taking into account degree of
#' concordance with the RPPA/expression.  Categories for genes with multiple
#' alterations are combined.
#'
#' Note some additional hard-coded filtering of the RPPA is performed.
#'
#' @param gs A `tibble` containing at least the following columns:
#' * `lab_id` Patient sample identifier
#' * `Gene` Gene name
#' @param muts A `tibble`containing information on gene mutations per sample
#'   containing at least the following columns:
#' * `lab_id` Patient sample identifier
#' * `Gene` Gene name
#' @param cnas A `tibble` containing gene-level copy number alterations containing
#'   at least the following columns:
#' *  `lab_id` Patient sample identifier
#' * `gene_name` Gene name
#' * `nTotal` Sum of the allele-specific copy number values (Can be NA)
#' * `cent_cna` Value indicating type and magnitude of the CNA, here +/-2.0 is
#'   used as the threshold for Amplification or Deletion
#' @param exprs A `tibble` of expression Zscore values derived from [mat.to.z.tbl()]
#' @param rppa A `tibble` of RPPA abundance Zscore values derived from [mat.to.z.tbl()]
#' @param rppa.annot A `tibble` containing the mapping from antibody name to
#'   gene name containing at least the following columns:
#' * `Antibody Name` RPPA antibody name
#' * `gene_name` Gene name
#' @results A `tibble` with at least the following columns:
#' * `Gene` Gene name
#' * `lab_id` Patient sample identifier
#' * `mut` Whether the patient has a mutation in the given gene (yes, no, NA)
#'      where NA indicates the patient did not have called mutations
#' * `nTotal` As carried over from the `cnas` param
#' * `cent_cna` As carried over from the `cnas` param
#' * `exprs_z` As carried over from the `exprs` param
#' * `rppa_z` As carried over from the `rppa` param
#' * `cna_class` One of the following classifications based on comparison of
#'   `cent_cna` with `exprs_z` and `rppa_z`:
#'   ** 'cA' Concordant Amplification
#'   ** 'dA' Discordant Amplification
#'   ** 'cD' Concordant Deletion
#'   ** 'dD' Discordant Deletion
#'   ** 'A' Amplification (No information on RPPA/exprs)
#'   ** 'D' Deletion (No information on RPPA/exprs)
#'   ** 'pA' Putative Amplification (No CNA data)
#'   ** 'pD' Putative Deletion (No CNA data)
#' * `combined_class` Combination of Mutation and CNA classifications delimited
#'   by ';'
combine.gs.mut.cna.exprs.rppa <- function(gs, muts, cnas, exprs, rppa, rppa.annots) {
  mut.genes <- select(muts, lab_id, Gene = symbol) %>%
    unique() %>%
    mutate(
      mut = "yes",
      lab_id = as.character(lab_id)
    )

  gs <- left_join(gs, mut.genes, by = c("lab_id", "Gene")) %>%
    mutate(mut = ifelse(is.na(mut), "no", mut))

  # if sample has no mutations, set to NA

  gs <- gs %>% mutate(mut = ifelse(lab_id %in% mut.genes$lab_id, mut, NA_character_))

  cna.genes <- select(cnas, lab_id, Gene = gene_name, nTotal, cent_cna) %>%
    mutate(lab_id = as.character(lab_id))

  gs <- left_join(gs, cna.genes,
    by = c("lab_id", "Gene")
  )

  gs <- left_join(gs, exprs, by = c(Gene = "gene_name", lab_id = "lab_id"))

  rppa <- inner_join(filter(rppa.annots, grepl("_p", `Antibody Name`) == F),
    rppa,
    by = c("Antibody Name" = "gene_name"), relationship = "many-to-many"
  )

  # remove remaining abs, that are not relevant to total-protein

  rppa <- filter(
    rppa,
    `Antibody Name` %in% c("Akt", "DM-Histone-H3", "DM-K9-Histone-H3", "Notch1-cleaved", "p38-MAPK") == F
  )

  summarize(rppa, n = n(), .by = c("lab_id", "gene_name")) %>%
    filter(n > 1) %>%
    {
      stopifnot(nrow(.) == 0)
    }

  gs <- left_join(gs, select(rppa, Gene = gene_name, lab_id, rppa_z = exprs_z), by = c("lab_id", "Gene"))

  gs <- mutate(gs,
    cna_class = case_when(
      cent_cna > 1 & ((exprs_z > 0 & is.na(rppa_z)) | (is.na(exprs_z) & rppa_z > 0) | (exprs_z > 0 & rppa_z > 0)) ~ "cA",
      cent_cna > 1 & ((exprs_z <= 0 & is.na(rppa_z)) | (is.na(exprs_z) & rppa_z <= 0) | (exprs_z <= 0 | rppa_z <= 0)) ~ "dA",
      cent_cna > 1 & is.na(exprs_z) & is.na(rppa_z) ~ "A",
      cent_cna < -1 & ((exprs_z < 0 & is.na(rppa_z)) | (is.na(exprs_z) & rppa_z < 0) | (exprs_z < 0 & rppa_z < 0)) ~ "cD",
      cent_cna < -1 & ((exprs_z >= 0 & is.na(rppa_z)) | (is.na(exprs_z) & rppa_z >= 0) | (exprs_z >= 0 | rppa_z >= 0)) ~ "dD",
      cent_cna < -1 & is.na(exprs_z) & is.na(rppa_z) ~ "D",
      is.na(cent_cna) & ((exprs_z > 2 & is.na(rppa_z)) | (is.na(exprs_z) & rppa_z > 2) | (exprs_z > 2 & rppa_z > 2)) ~ "pA",
      is.na(cent_cna) & ((exprs_z < -2 & is.na(rppa_z)) | (is.na(exprs_z) & rppa_z < -2) | (exprs_z < -2 & rppa_z < -2)) ~ "pD",
      .default = ""
    )
  )

  gs <- mutate(gs, combined_class = case_when(
    cna_class == "" & mut == "yes" ~ "Mut",
    cna_class == "" & (mut == "no" | is.na(mut)) ~ "",
    cna_class != "" & mut == "yes" ~ paste0("M;", cna_class),
    cna_class != "" & (mut == "no" | is.na(mut)) ~ cna_class,
  ))

  gs
}

#' Gene Score Significance Heatmap
#'
#' Produce a heatmap highlighting Gene Score significance for a given gene set.
#' Colors are relative to a P-value categorization.  Rows and columns are ordered
#' by hierarchical clustering of -log10(P-value).  Alteration classificatons
#' are also displayed.
#'
#' @param gs.tbl A `tibble` containing gene score and alteration data as generated
#'   by [combine.gs.mut.cna.exprs.rppa()]
#' @param use.genes A character vector of genes to plot
#' @returns A PDF: 'figures/tcga_actionable_gs.pdf'
plot.gs.summary <- function(gs.tbl, use.genes) {
  gs.tbl <- mutate(gs.tbl, level_fac = factor(case_when(
    level == 0 ~ "N.S.",
    level == 1 ~ "p<0.05",
    level == 2 ~ "p<0.01",
    level == 3 ~ "p<0.001"
  ), levels = c("N.S.", "p<0.05", "p<0.01", "p<0.001"), ordered = T))

  use.gs.tbl <- filter(gs.tbl, is.na(score) == F & Gene %in% use.genes)

  use.gs.tbl <- mutate(use.gs.tbl, nlp = -log10(pval))

  gs.tmat <- pivot_wider(use.gs.tbl, id_cols = Gene, names_from = lab_id, values_from = nlp)
  gs.mat <- as.matrix(select(gs.tmat, -1))
  rownames(gs.mat) <- gs.tmat[[1]]

  row.hc <- hclust(dist(gs.mat), method = "average")

  col.hc <- hclust(dist(t(gs.mat)), method = "average")

  use.gs.tbl <- mutate(use.gs.tbl,
    gene_fac = factor(Gene, levels = row.hc$labels[row.hc$order], ordered = T),
    sample_fac = factor(lab_id, levels = col.hc$labels[col.hc$order], ordered = T)
  )

  sig.pts <- filter(use.gs.tbl, level > 0)

  use.gs.tbl <- filter(use.gs.tbl, lab_id %in% sig.pts$lab_id)

  p1 <- ggplot(data = use.gs.tbl, mapping = aes(x = sample_fac, y = gene_fac, fill = level_fac)) +
    geom_tile(color = "black") +
    scale_fill_viridis_d(option = "G", direction = -1, name = "Sig. Category", begin = .25) +
    geom_text(mapping = aes(label = combined_class, color = -log10(pval) > 2), size=2) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    theme_classic() +
    theme(
      text = element_text(size=10, family = "Helvetica"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
      ) +
    xlab("") +
    ylab("")

  ggsave(p1, file = "figures/tcga_actionable_gs.pdf", width = 3.75, height = 2.25)

  "figures/tcga_actionable_gs.pdf"
}


#' Plot a PDF of multiple Response Cards
#'
#' Given one or more patient samples in a `tibble` from [combine.gs.mut.cna.exprs.rppa()]
#' generate a series of Response Cards, one-per-patient sample ordered by the
#' best significance level of the Gene Score (highest to lowest).
#'
#' @param gs.alt.tbl A `tibble` of combined Gene Score and alteration data as
#'   output from [combine.gs.mut.cna.exprs.rppa()].
#' @param htome.mat A gene x drug matrix as produced by [gene.by.drug.mat()]
#' @param score.mat A drug x patient matrix as produced by [drug.by.pt.mat()]
#' @param inhib.map A `tibble` containing a column for `inhibitor` and
#'   one for the corresponding `plot_name`
#' @param file.pref Prefix of the output pdf (default 'all')
#' @param file.width Width in inches of the output pdf
#' @param file.height Height in inches of the output pdf
#' @returns A PDF 'figures/hnscc_response_cards.pdf'
plot.all.response.cards <- function(gs.alt.tbl, htome.mat, score.mat, inhib.map, file.pref="all", file.width=7, file.height=6) {
  gs.alt.tbl <- filter(gs.alt.tbl, is.na(score) == F)

  # order patients from most significant to least
  pt.stars <- unique(select(gs.alt.tbl, lab_id, best_level)) %>%
    arrange(desc(best_level))

  gs.alt.tbl <- group_by(gs.alt.tbl, lab_id)
  pt.list <- setNames(group_split(gs.alt.tbl), group_keys(gs.alt.tbl)$lab_id)

  pt.list <- pt.list[pt.stars$lab_id]
  
  if (length(pt.list) > 26){
    stop("ERROR: Currently this function is limited to 26 patients due to labeling with the alphabet")
  }
  
  pdf(file = paste0("figures/",file.pref,"_response_cards.pdf"), width = file.width, height = file.height)

  for (pt_idx in seq_along(pt.list)) {
    message(names(pt.list)[pt_idx])
    plot.response.card(pt.list[[pt_idx]], htome.mat, score.mat, inhib.map, file.width, file.height, add.letter = letters[pt_idx])
  }

  dev.off()

  paste0("figures/",file.pref,"_response_cards.pdf")
}

#' Plot a Response Card
#'
#' Given a patient sample in a `tibble` derived from
#' [combine.gs.mut.cna.exprs.rppa()] generate a Response Cards.
#' Note the response card is plotted on the screen and so would need to be
#' captured into a PDF, PNG etc.
#'
#' @param result.tbl A `tibble` of combined Gene Score and alteration data as
#'   output from [combine.gs.mut.cna.exprs.rppa()].
#' @param htome.mat A gene x drug matrix as produced by [gene.by.drug.mat()]
#' @param score.mat A drug x patient matrix as produced by [drug.by.pt.mat()]
#' @param inhib.map A `tibble` containing a column for `inhibitor` and
#'   one for the corresponding `plot_name`
#' @param expected.height Expected plot width in inches (used for scaling the various plot components)
#' @param expected.height Expected plot height in inches (used for scaling the various plot components)
#' @param add.letter If present and non-NULL add the corresponding text to the top left of the plot (meant for multi-panel figures)
#' @returns Nothing is returned but a plot of the Response Card is produced as a side-effect
plot.response.card <- function(result.tbl, htome.mat, score.mat, inhib.map, expected.width, expected.height, add.letter=NULL) {
  cur.samp <- result.tbl$lab_id[1]

  # include target genes for drugs < -2

  cur.inhib <- score.mat[colnames(htome.mat), cur.samp]

  if (any(cur.inhib < -2, na.rm = T)) {
    lt2.inhib <- names(cur.inhib)[is.na(cur.inhib) == F & cur.inhib < -2]

    if (length(lt2.inhib) > 1) {
      lt2.targs <- rowSums(htome.mat[, lt2.inhib]) > 0
    } else {
      lt2.targs <- htome.mat[, lt2.inhib] > 0
    }

    addtl.targs <- names(lt2.targs)[lt2.targs == T]
  } else {
    addtl.targs <- character(0)
  }

  sig.results <- filter(result.tbl, (level == best_level) | (Gene %in% addtl.targs)) %>%
    arrange(pval)

  # provide a minimum of 10 targets

  if (nrow(sig.results) < 10) {
    top10.genes <- result.tbl %>%
      filter(Gene %in% addtl.targs == F) %>%
      arrange(pval) %>%
      head(n = 10 - length(addtl.targs))

    addtl.targs <- union(top10.genes$Gene, addtl.targs)

    sig.results <- filter(result.tbl, (level == best_level) | (Gene %in% addtl.targs)) %>%
      arrange(pval)
  }

  sig.results <- mutate(sig.results,
    rppa_z = ifelse(is.na(rppa_z), 0, rppa_z),
    exprs_z = ifelse(is.na(exprs_z), 0, exprs_z),
    cent_cna = ifelse(is.na(cent_cna), 0, cent_cna)
  )

  sig.results <- mutate(sig.results, sig_stars = sapply(level, function(x) {
    if (length(x) > 0) {
      paste(rep("*", x), collapse = "")
    } else {
      ""
    }
  }))

  score.ylim <- floor(min(sig.results$score, na.rm = T)) - nchar(sig.results$sig_stars[1])

  stopifnot(length(unique(sig.results$lab_id)) == 1)

  tmp.mat <- htome.mat[sig.results$Gene, ]

  # only keep inhibitors that are relevant to these genes

  inhib.pres <- colSums(tmp.mat) > 0

  keep.inhib <- names(inhib.pres)[inhib.pres == T]

  tmp.inhib <- sort(score.mat[keep.inhib, cur.samp], decreasing = F)

  inhib.vec <- setNames(inhib.map$plot_name, inhib.map$inhibitor)

  tmp.mat <- tmp.mat[, names(tmp.inhib)]
  colnames(tmp.mat) <- unname(inhib.vec[names(tmp.inhib)])

  stopifnot(all(rownames(tmp.mat) == sig.results$Gene))

  ha <- HeatmapAnnotation(
    `Inhibitor Zscore` = anno_barplot(tmp.inhib),
    height = unit((1/6)*expected.height, "inches")#was 3 cm
  )

  hl <- rowAnnotation(
    `Gene Score` = anno_barplot(sig.results$score,
      ylim = c(score.ylim, max(0, max(sig.results$score))),
      gp = gpar(fill = "grey")
    ),
    width = unit(.125*expected.height, "inches")#was 2 cm
  )

  hr.list <- list()
  has.rppa <- F
  has.exprs <- F

  # mutation
  # check if any of the genes has a mutation--surrogate for has been sequenced
  if (all(is.na(result.tbl$mut)) == F) {
    hr.list <- c(hr.list, rowAnnotation(
      Mutation = anno_simple(sig.results$mut, col = c(yes = "black", no = "white"), gp = gpar(col = "black")),
      annotation_name_rot = 90#was gap = unit((1/3)*expected.height, "inches") 10mm
    ))
  }

  # RPPA
  if (all(is.na(result.tbl$rppa_z)) == F) {
    rppa.ylim <- c(min(-2, min(sig.results$rppa_z)), c(max(2, max(sig.results$rppa_z))))

    hr.list <- c(hr.list, rowAnnotation(
      `RPPA Zscore` = anno_barplot(sig.results$rppa_z,
        ylim = rppa.ylim, width = unit(.05*expected.width, "inches"),#unit(1, "cm"),
        gp = gpar(fill = ifelse(sig.results$rppa_z > 0, "red", "blue"))
      ),
      annotation_name_rot = 90
    ))

    has.rppa <- T
  }

  # exprs

  if (all(is.na(result.tbl$exprs_z)) == F) {
    exprs.ylim <- c(min(-2, min(sig.results$exprs_z)), c(max(2, max(sig.results$exprs_z))))

    hr.list <- c(hr.list, rowAnnotation(
      `Exprs Zscore` = anno_barplot(sig.results$exprs_z,
        ylim = exprs.ylim, width = unit(.05*expected.width, "inches"),
        gp = gpar(fill = ifelse(sig.results$exprs_z > 0, "red", "blue"))
      ),
      annotation_name_rot = 90
    ))

    has.exprs <- T
  }


  # cnv
  if (all(is.na(result.tbl$cent_cna)) == F) {
    hr.list <- c(hr.list, rowAnnotation(
      CNV = anno_barplot(sig.results$cent_cna,
        gp = gpar(fill = ifelse(sig.results$cent_cna > 0, "red", "blue")),
        width = unit(.05*expected.width, "inches")
      ),
      annotation_name_rot = 90
    ))
  }

  hr.obj <- do.call(c, hr.list)

  if (length(hr.list) > 0) {
    hr.obj@gap <- rep(unit(.01*expected.width, "inches"), length(hr.obj) + 1)#was unit(2, "mm")
  }

  if (nrow(tmp.mat) > 15) {
    row.cex <- .75
  } else {
    row.cex <- 1
  }

  ht <- Heatmap(tmp.mat,
    rect_gp = gpar(col = "darkgrey"), col = circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "black")),
    cluster_columns = F, cluster_rows = F, top_annotation = ha, left_annotation = hl, right_annotation = hr.obj, column_names_gp = gpar(cex = .5),
    column_names_side = "top", show_heatmap_legend = F, column_title = paste("Sample:", cur.samp, "sig-level:", sig.results$best_level[1]),
    row_names_max_width = unit((1/6)*expected.width,"inches"), row_names_gp = gpar(cex = row.cex), row_names_side = "left" #was unit(30, "cm")
  )

  draw(ht, padding = unit(c(.01*expected.width, .1*expected.width, .01*expected.width, .1*expected.width), "inches"))
  #padding = bottom, left, top and right

  if (has.rppa == T) {
    decorate_annotation("RPPA Zscore", {
      grid.lines(unit(c(0, 0), "native"), c(0, 1), gp = gpar(col = "black", lty = 2))
    })
  }

  if (has.exprs == T) {
    decorate_annotation("Exprs Zscore", {
      grid.lines(unit(c(0, 0), "native"), c(0, 1), gp = gpar(col = "black", lty = 2))
    })
  }

  ro.results <- sig.results

  decorate_annotation("Gene Score", {
    grid.text(label = ro.results$sig_stars, x = unit(ro.results$score - .25, "native"), y = unit(nrow(ro.results) - seq_along(ro.results$score) + 1, "native"), just = c(1, .75))
  })
  
  if ((missing(add.letter) || is.null(add.letter))==F){
    grid.text(add.letter, x = 0.05, y = .975, gp = gpar(fontsize = 12))
  }
}

#' Perform Random Walk with Restarts
#'
#' Given a graph and a vector of seed proteins perform RWR to prioritize genes
#' relative to the seed proteins.  Note that symmetric normalization is hard-coded.
#'
#' @param i.graph An undirected {igraph} graph.
#' @param seed.prots Either a character vector of protein/gene names or a named
#'   numeric vector containing the (relative) values to be assigned to the seeds.
#'   The values will be normalized to sum to 1.
#' @param c.val Probability of restart at a seed node, default is 0.7.
#' @returns A named vector containing the RWR scores
run.rwr <- function(i.graph, seed.prots, c.val = .7) {
  threshold <- 1e-10
  maxit <- 100

  graph.mat <- as_adjacency_matrix(i.graph, type = "both", names = T, sparse = T)

  use.sums <- Matrix::colSums(graph.mat)

  stopifnot(all(rownames(graph.mat) == colnames(graph.mat)))

  d.graph.mat <- Diagonal(x = 1 / sqrt(use.sums))
  dimnames(d.graph.mat) <- list(rownames(graph.mat), colnames(graph.mat))

  # symmetric normalization as that is what RDPN seems best at
  graph.sp.mat <- d.graph.mat %*% graph.mat %*% d.graph.mat

  residue <- 1
  iter <- 1

  # probability of restart
  prox.vector <- rep(0, nrow(graph.sp.mat))
  names(prox.vector) <- rownames(graph.sp.mat)

  if (length(seed.prots) > 0) {
    if (is.numeric(seed.prots) && is.null(names(seed.prots)) == F) {
      common.genes <- intersect(names(seed.prots), names(prox.vector))
      prox.vector[common.genes] <- seed.prots[common.genes]
    } else if (is.character(seed.prots)) {
      prox.vector[seed.prots] <- 1
    } else {
      stop("ERROR: Unexpected type for 'seed.prots'")
    }

    prox.vector <- prox.vector / sum(prox.vector)

    restart.vector <- prox.vector

    while (residue > threshold && iter < maxit) {
      old.prox.vector <- prox.vector
      prox.vector <- as.numeric((1 - c.val) * (graph.sp.mat %*% prox.vector) + (c.val * restart.vector))
      residue <- as.numeric(dist(rbind(prox.vector, old.prox.vector), method = "euclidean"))
      iter <- iter + 1
    }
  }

  setNames(prox.vector, rownames(graph.sp.mat))
}

#' Produce a PDF of per-sample RWR score barplots for alterations
#'
#' For each patient sample, a barplot is produced which ranks alterations the top
#' 10 alterations by their highest RWR score.  There is one plot per page.
#' @param np.tbl A combined `tibble` of RWR scores and P-values and alterations
#'   as generated by [combine.gs.mut.cna.exprs.rppa()].  In addition it should
#'   have a column for RDPN P-value as `np_pval`.
#' @param file.pref Prefix for the output PDF (default 'all')
#' @results A PDF 'figures/net_prop_results.pdf'
plot.network.prop.bars <- function(np.tbl, file.pref="all") {
  np.tbl <- filter(
    np.tbl,
    combined_class %in% c("pA", "pD", "") == F &
      np_pval < .05
  )

  plot.list <- lapply(group_split(np.tbl, .by = lab_id), function(x) {
    tmp.x <- x %>%
      arrange(desc(np_score)) %>%
      head(n = 10) %>%
      mutate(
        gene_fac = factor(Gene, levels = Gene, ordered = T),
        mut_type = case_when(
          grepl("M", combined_class) ~ "Mut",
          grepl("A", combined_class) ~ "Amp",
          grepl("D", combined_class) ~ "Del",
        )
      )

    ggplot(data = tmp.x, mapping = aes(x = gene_fac, y = np_score, fill = mut_type)) +
      geom_col(color = "black") +
      scale_fill_manual(values = c(Mut = "purple", Amp = "red", Del = "blue"), guide="none") +
      geom_text(mapping = aes(label = combined_class, y = np_score + (max(np_score) * .1))) +
      theme_bw() +
      theme(
        axis.text = element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        plot.margin = margin(t = .25, unit="inches")
        ) +
      xlab("") +
      ylab("Network Propagation Score") +
      ggtitle(x$lab_id[1])
  })

  pdf(paste0("figures/", file.pref, "_net_prop_results.pdf"), width = 3.25, height = 2.5)

  for (p in seq_along(plot.list)) {
    show(plot.list[[p]])
    
    grid::grid.text(letters[p], x = 0.05, y = 1, just = "top", gp = grid::gpar(fontsize = 13))
  }

  dev.off()

  "figures/net_prop_results.pdf"
}

#' Produce a Steiner Tree graph
#'
#' This is an implementation of the merged steiner tree algorithm (termed STM) as
#' discussed by Sadeghi et al 2013 (1).  {SteinerNet} which is the accompanying
#' R package is no longer a part of CRAN, and hasn't been maintained for several
#' years.  More importantly it does not seem to scale to as large a network as we use.
#'
#' @param g An {igraph} graph object containing the full network
#' @param term.nodes A character vector containing the nodes to be included in
#'   the network
#' @param r Number of repetitions (A terminal node is randomly chosen to start)
#' @param summarize.by Whether to return the smallest graph (`min`) or the `union`
#'   of graphs across `r` iterations.
#' @param choose.single.sp TRUE or FALSE whether a single shortest path should be
#'   arbitrarily chosen or whether the union should be collected.
#' @param verbose TRUE or FALSE whether more verbose messages should be displayed
#' @returns An igraph object containing the merged Steiner tree
#' @reference (1) Sadeghi, A., Fröhlich, H. Steiner tree methods for optimal
#'   sub-network identification: an empirical study. BMC Bioinformatics 14, 144 (2013). https://doi.org/10.1186/1471-2105-14-144
steiner.tree <- function(g, term.nodes, r = 10, summarize.by = c("union", "min"), choose.single.sp = F, verbose = F) {
  summarize.by <- match.arg(summarize.by)

  missing.nodes <- setdiff(term.nodes, V(g)$name)

  if (length(missing.nodes) > 0) {
    message("Removing ", length(missing.nodes), " nodes not in graph")
  }

  # remove attributes

  for (i in edge_attr_names(g)) {
    g <- delete_edge_attr(g, i)
  }

  term.nodes <- intersect(term.nodes, V(g)$name)

  tree.list <- lapply(seq_len(r), function(cur_t) {
    if (verbose) {
      message("Starting iteration: ", cur_t)
    }
    .do.steiner.tree(g, term.nodes, choose.single.sp, verbose)
  })

  if (summarize.by == "min") {
    tree.size <- sapply(tree.list, function(x) {
      length(E(x))
    })

    use.tree <- which(tree.size == min(tree.size))
  } else if (summarize.by == "union") {
    use.tree <- seq_along(tree.list)
  }

  do.call(igraph::union, tree.list[use.tree])
}

.do.steiner.tree <- function(g, term.nodes, choose.single.sp = F, verbose = F) {
  if (is_connected(g) == F || is_directed(g) == T) {
    stop("ERROR: g needs to be a connected undirected graph")
  }

  if (length(term.nodes) <= 1 && all(term.nodes %in% V(g)$name) == F) {
    stop("ERROR: term.nodes needs to have more than 1 node in the graph")
  }

  init.term <- sample(term.nodes, 1)

  working.graph <- induced_subgraph(g, init.term)

  iter <- 1

  while (all(term.nodes %in% V(working.graph)$name) == F) {
    if (verbose) {
      message("    Starting terminal node sp iteration: ", iter)
    }

    rem.terms <- setdiff(term.nodes, V(working.graph)$name)

    # return a matrix of distances between the remaining terms (rows) and the nodes in the current working graph (columns)
    dist.mat <- distances(g, v = V(g)[name %in% rem.terms], to = V(g)[name %in% V(working.graph)$name], mode = "out", algorithm = "unweighted")

    min.dist <- apply(dist.mat, 1, min)

    min.terms <- names(min.dist)[min.dist == min(min.dist)]

    if (choose.single.sp) {
      min.terms <- min.terms[1]
    }

    sp.graphs <- lapply(min.terms, function(i) {
      # select the node(s) in the working graph that the current term is closest to
      # and compute the shortest path
      tmp.min <- dist.mat[i, , drop = F]
      tmp.min <- colnames(tmp.min)[tmp.min[1, ] == min(tmp.min[1, ])]
      cur.sps <- all_shortest_paths(g, from = V(g)[name == i], to = V(g)[name %in% tmp.min], mode = "out")

      if (choose.single.sp) {
        iter.sps <- cur.sps$res[1]
      } else {
        iter.sps <- cur.sps$res
      }
      
      tmp.sp.graph <- lapply(iter.sps, function(j) {
        # create a graph object from the shortest path results by ensuring
        # proper format: a->b b->c as opposed to a->b->c
        tmp.adj.l <- do.call(c, lapply(1:(length(j) - 1), function(k) {
          j[k:(k + 1)]
        }))

        subgraph_from_edges(g, eids = E(g, P = tmp.adj.l, directed = F), delete.vertices = TRUE)
      })
      
      if (length(tmp.sp.graph) == 1){
        tmp.sp.graph[[1]]
      }else{
        do.call(igraph::union, tmp.sp.graph)
      }
      
    })
    
    working.graph <- do.call(igraph::union, append(sp.graphs, list(working.graph)))

    iter <- iter + 1
  }

  if (verbose) {
    message("  Computing MST")
  }

  cand.mst <- mst(working.graph, algorithm = "unweighted")

  # remove non-terminal nodes of degree 1
  v.deg <- degree(cand.mst, mode = "out", loops = F)

  ret.mst <- delete_vertices(cand.mst, V(cand.mst)[name %in% setdiff(names(v.deg)[v.deg == 1], term.nodes)])

  ret.mst
}

#' Create an {igraph} graph object from Reactome FIs
#'
#' Reads in a Reactome FI file and thresholds the edges to only include those
#' with Scores > 0.95 and further restricts the graph to only the largest
#' connected component.
#'
#' @param fi.file Path to the Reactome FI graph file
#' @results An {igraph} graph object
make.hc.reactome.w.lcc <- function(fi.file) {
  graph.tbl <- readr::read_delim(fi.file) %>%
    dplyr::filter(Score > .95)

  fi.graph <- graph_from_data_frame(as.data.frame(graph.tbl), directed = F)

  # choose the largest connected component

  fi.comps <- groups(components(fi.graph))

  largest.comp <- which.max(sapply(fi.comps, length))

  orig.graph <- induced_subgraph(fi.graph, V(fi.graph)[name %in% fi.comps[[largest.comp]]])

  orig.graph
}

#' Output files to be use by Cytoscape to create a graph figure
#'
#' This function outputs a node table and edge table to be used with Cytoscape
#' based on the union of the 'seed' genes (top genes by Gene Score) and the closest
#' significant `top.n` alterations for a given patient sample `pt.id`.
#'
#' @param fi.graph An {igraph} graph object as used for the network propagation
#'   based prioritization
#' @param np.alts A `tibble` of merged Gene Score, network propagation and alterations
#'   as output by [combine.gs.mut.cna.exprs.rppa()]
#' @param pt.id The patient sample `lab_id` to use for the graph
#' @param top.n The number of top alterations to plot as ranked by network propagation.
#' @returns Two files: outputs/nodes_for_cytoscape_`pt.id`.txt and
#'   outputs/edges_for_cytoscape_`pt.id`.txt
make.graph.figure.by.pt <- function(fi.graph, np.alts, pt.id, top.n) {
  keep.targs <- filter(
    np.alts,
    lab_id %in% pt.id &
      combined_class %in% c("pA", "pD", "") == F &
      np_pval < .05
  ) %>%
    arrange(desc(np_score)) %>%
    head(n = top.n)

  seed.genes <- filter(np.alts, lab_id %in% pt.id & is_seed == T)
  
  sub.graph <- steiner.tree(fi.graph,
    term.nodes = union(keep.targs$Gene, seed.genes$Gene),
    r = 10,
    summarize.by = "union",
    choose.single.sp = F,
    verbose = T
  )

  all.genes <- bind_rows(keep.targs, seed.genes)

  all.genes <- mutate(all.genes,
    combined_class = ifelse(combined_class == "" & is_seed == T, "seed", combined_class)
  )

  # data frame representation to simplify annotating
  graph.list <- igraph::as_data_frame(sub.graph, what = c("both"))

  node.df <- inner_join(graph.list$vertices,
    select(all.genes, name = Gene, combined_class),
    by = c("name")
  )

  names(node.df)[1] <- "id"

  write.table(node.df,
    file = paste0("outputs/nodes_for_cytoscape_", pt.id, ".txt"),
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
  )

  write.table(setNames(graph.list$edges, c("source", "target")),
    file = paste0("outputs/edges_for_cytoscape_", pt.id, ".txt"),
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F
  )

  c(
    paste0("outputs/nodes_for_cytoscape_", pt.id, ".txt"),
    paste0("outputs/edges_for_cytoscape_", pt.id, ".txt")
  )
}
