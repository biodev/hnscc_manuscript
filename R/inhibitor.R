#' Create a drug x patient matrix given a `tibble`
#'
#' @param inhib A `tibble` that contains columns for:
#' * `type` Type of drug (e.g. single-agent vs combination)
#' * `inhibitor` Name of the drug
#' * `lab_id` Patient sample identifier
#' @param values The corresponding column in `inhib` to use as the values for
#'   the matrix.
drug.by.pt.mat <- function(inhib, values = "zscore") {
  inhib <- group_by(inhib, type)
  inhib.list <- setNames(group_split(inhib), group_keys(inhib)$type)

  lapply(inhib.list, function(tmp.inh) {
    tmp.tmat <- pivot_wider(tmp.inh,
      id_cols = inhibitor,
      names_from = lab_id,
      values_from = all_of(values)
    )

    tmp.mat <- as.matrix(select(tmp.tmat, -1))
    rownames(tmp.mat) <- tmp.tmat[[1]]

    tmp.mat
  })
}

#' Create a gene x drug binary matrix given a `tibble`
#'
#' @param tome A `tibble` that contains columns for:
#' * `Gene` Gene identifier
#' * `inhibitor` Name of a drug
gene.by.drug.mat <- function(tome) {
  tome <- mutate(tome, values = 1)

  tmp.tmat <- pivot_wider(tome,
    id_cols = Gene,
    names_from = inhibitor,
    values_from = "values",
    values_fill = 0
  )

  tmp.mat <- as.matrix(select(tmp.tmat, -1))
  rownames(tmp.mat) <- tmp.tmat[[1]]

  tmp.mat
}

#' Output a workbook containing a unique list of drugs
#'
#' @param inhib A `tibble` that contains columns for:
#' * `plot_name` Shortened/simplified name used in plots etcs.
#' * `type` Type of drug (e.g. single-agent vs combination)
#' * `inhibitor1` Name of the first drug (or single agent)
#' * `inhibitor2` Name of the second drug in a combination
#' @returns An Excel file: 'outputs/inhibitor_table.xlsx'
make.inhib.table <- function(inhib) {
  uniq.inhib <- select(
    inhib,
    display_name = plot_name, type, inhibitor1, inhibitor2
  ) %>%
    unique() %>%
    arrange(
      inhibitor1, inhibitor2
    )

  # Create a highlighted workbook, adapted from:
  # https://stackoverflow.com/questions/71600719/addstyle-function-in-openxlsx-does-not-fill-cells-in-excel-spreadsheet-with-the
  wb <- createWorkbook()

  # Add a worksheets
  addWorksheet(wb, "Inhibitor Table", gridLines = FALSE)

  # write data to worksheet 1
  writeData(wb, sheet = 1, uniq.inhib, rowNames = F)

  # Create style object in order to fill certain cells in the excel spreadsheet.

  Styling_object <- createStyle(fgFill = "lightgrey")

  cur.color <- T

  for (i in seq(5, nrow(uniq.inhib), by = 4)) {
    if (cur.color) {
      addStyle(wb, sheet = 1, style = Styling_object, rows = i:(i + 3), cols = 1:4, gridExpand = T)
      cur.color <- F
    } else {
      cur.color <- T
    }
  }

  # save the workbook
  saveWorkbook(wb, "outputs/inhibitor_table.xlsx", overwrite = TRUE)

  "outputs/inhibitor_table.xlsx"
}

#' Read and format the PanCan Oncogenic Signalling Pathways
#'
#' Reads the pathways from the Supplementary Data Excel file and compares them
#' to Entrez gene symbols only keeping those that are direct matches.
#'
#' @param pancan.pathway.excel Table S3 from https://doi.org/10.1016/j.cell.2018.03.035
#' @param gene.info.file Human gene info file from Entrez
#' @returns A `tibble` with three columns for 'pathway' and 'Gene'
form.pancan.pathways <- function(pancan.pathway.excel, gene.info.file) {
  gene.info <- readr::read_delim(gene.info.file) %>%
    filter(`#tax_id` == 9606)

  tabs <- getSheetNames(pancan.pathway.excel)
  tabs <- tabs[1:10]

  path.tbl <- bind_rows(lapply(setNames(tabs, tabs), function(x) {
    message(x)

    # the first tab is different than the others
    if (x == tabs[1]) {
      as_tibble(read.xlsx(pancan.pathway.excel, sheet = x, startRow = 3)) %>%
        select(Gene)
    } else {
      as_tibble(read.xlsx(pancan.pathway.excel, sheet = x)) %>%
        select(Gene)
    }
  }), .id = "pathway")

  # only keep those genes where there is a match (NOV is ambiguous)
  path.tbl.m <- inner_join(
    path.tbl,
    select(gene.info, GeneID, Gene = Symbol),
    by = "Gene"
  )

  # check to make sure there were no duplicates (genes or gene ids)
  select(path.tbl.m, Gene, GeneID) %>%
    unique() %>%
    summarize(n = n(), .by = c("Gene")) %>%
    {
      stopifnot(all(.$n == 1))
    }

  select(path.tbl.m, Gene, GeneID) %>%
    unique() %>%
    summarize(n = n(), .by = c("GeneID")) %>%
    {
      stopifnot(all(.$n == 1))
    }

  select(path.tbl.m, pathway, Gene)
}

#' Compute Gene Scores and Permutation P-values
#'
#' Given a drug x patient matrix of Zscores D and a gene x drug binary matrix T
#' compute the gene score as the sum of corresponding drug Zscores for each gene.
#' The permutation test permutes the rows of D and recomputes the gene scores.
#' P-values are produced by counting the number of permuted gene scores that are
#' greater than the original.  The P-values are then categorized into levels
#' per patient sample and the corresponding categories are added to the result.
#'
#' @param score.mat A drug x patient matrix as produced by `drug.by.pt.mat`
#' @param tome.mat A gene x drug matrix as produced by `gene.by.drug.mat`
#' @returns A `tibble` with the following columns:
#' * `Gene` Gene name
#' * `lab_id` Patient sample ID
#' * `score` Gene Score value
#' * `pval` Permutation P-value
#' * `level` 0-3 Where 0 indicates non-significance
#'   ** 1 P-value < .05
#'   ** 2 P-value < .01
#'   ** 3 P-value < .001
#' * `best_level` The highest level for any gene for the given patient.
compute.gene.scores <- function(score.mat, tome.mat) {
  # NAs are replaced by mean (zero as these are Zscores)

  score.mat[is.na(score.mat)] <- 0

  gene.scores <- tome.mat %*% score.mat[colnames(tome.mat), ]

  # compute permutation p-values

  #set.seed(4925225)

  perm.mat <- lapply(1:1000, function(perm) {
    tmp.inhibs <- sample(rownames(score.mat))
    tmp.score.mat <- score.mat[]
    rownames(tmp.score.mat) <- tmp.inhibs
    tmp.scores <- tome.mat %*% tmp.score.mat[colnames(tome.mat), ]
    # since lower values are better
    tmp.scores[rownames(gene.scores), colnames(gene.scores)] < gene.scores
  }) %>%
    Reduce(function(x, y) {
      x + y
    }, x = .)

  perm.tbl <- as_tibble(perm.mat, rownames = "Gene") %>%
    pivot_longer(cols = -Gene) %>%
    mutate(pval = (value + 1) / (1001)) %>%
    rename(lab_id = "name")

  score.tbl <- as_tibble(gene.scores, rownames = "Gene") %>%
    pivot_longer(cols = -Gene) %>%
    rename(lab_id = "name", score = "value")

  score.tbl <- inner_join(score.tbl, select(perm.tbl, -value), by = c("Gene", "lab_id"))

  score.tbl <- mutate(score.tbl,
    level = case_when(
      pval < .001 ~ 3,
      pval < .01 ~ 2,
      pval < .05 ~ 1,
      .default = 0
    )
  )

  best.level <- summarize(score.tbl, best_level = max(level), .by = lab_id)

  score.tbl <- inner_join(score.tbl, best.level, by = "lab_id")

  score.tbl
}

#' Create a mapping between PanCancer pathways and the HNSCC drugs based on Targetome
#'
#' In addition to the direct mapping for single-agent drugs, this also maps both
#' drugs in a given combination
#'
#' @param inhib A `tibble` that contains columns for:
#' * `type` Type of drug (e.g. single-agent vs combination)
#' * `inhibitor` Name of the drug
#' * `inhibitor1` Name of the first drug (or single agent)
#' * `inhibitor2` Name of the second drug in a combination
#' @param htome.summary A `tibble` containing an inhibitor -> Gene mapping with
#'   corresponding columns
#' @param pc.paths A `tibble` containing a Gene -> pathway mapping with
#'   corresponding columns
#' @returns A `tibble` containing the following columns:
#' * `inhibitor` Drug name
#' * `inhibitor1` Name of the first drug (or single agent)
#' * `inhibitor2` Name of the second drug in a combination
#' * `pathway1` Pathway targeted by `inhibitor1`
#' * `pathway2` Pathway targeted by `inhibitor2`
annotate.drugs.pathways <- function(inhib, htome.summary, pc.paths) {
  inhib.paths <- inner_join(
    htome.summary,
    pc.paths,
    by = "Gene"
  ) %>%
    select(inhibitor, pathway) %>%
    unique() %>%
    summarize(
      pathway = paste(sort(pathway), collapse = ";"),
      .by = inhibitor
    )

  combos <- filter(inhib, type == "combination")

  combo.paths <- select(
    combos,
    inhibitor, inhibitor1, inhibitor2
  ) %>%
    unique() %>%
    inner_join(
      select(inhib.paths, inhibitor1 = inhibitor, pathway1 = pathway),
      by = "inhibitor1"
    ) %>%
    inner_join(
      select(inhib.paths, inhibitor2 = inhibitor, pathway2 = pathway),
      by = "inhibitor2"
    )

  inhib.paths <- mutate(
    inhib.paths,
    inhibitor1 = inhibitor,
    inhibitor2 = NA_character_,
    pathway1 = pathway,
    pathway2 = NA_character_
  )

  inhib.paths <- bind_rows(
    inhib.paths[, names(combo.paths)],
    combo.paths
  )

  inhib.paths
}

#' Generate a combined set of plots displaying the relationship between drugs and pathways
#'
#' @param htome.summary A `tibble` containing an inhibitor -> Gene mapping with
#'   corresponding columns
#' @param pc.paths A `tibble` containing a Gene -> pathway mapping with
#'   corresponding columns
#' @param action.genes A `tibble` with a `gene_name` column indicating the genes
#'   reported to be therapeutic targets in the TCGA-HNSC.
#' @param drug.paths A `tibble` mapping of drugs to pathways as produced by
#'   [annotate.drugs.pathways()]
#' @returns A PDF figure in 'figures/summary_of_targets.pdf'
overlap.tcga.path.plots <- function(htome.summary, pc.paths, action.genes, drug.paths) {
  targ.paths <- inner_join(
    htome.summary,
    pc.paths,
    by = "Gene"
  ) %>%
    summarize(
      n_genes = length(unique(Gene)),
      .by = pathway
    )

  targ.paths <- arrange(targ.paths, n_genes) %>%
    mutate(path_fac = factor(pathway, levels = c(pathway, "None"), ordered = T))

  pathway.pal <- setNames(c(rev(brewer.pal(nrow(targ.paths), name = "Dark2")), "grey"), c(levels(targ.paths$path_fac)))

  p1 <- ggplot(data = targ.paths, mapping = aes(y = path_fac, x = n_genes, fill = pathway)) +
    geom_col(color = "black") +
    scale_fill_manual(values = pathway.pal, guide = "none") +
    theme_bw() +
    theme(text = element_text(family = "Helvetica", size = 12)) +
    xlab("Number of Targeted Genes") +
    ylab("PanCancer Pathways")


  targ.genes <- htome.summary %>%
    summarize(
      n = length(unique(inhibitor)),
      .by = Gene
    )

  targ.action.genes <- filter(targ.genes, Gene %in% action.genes$gene_name) %>%
    inner_join(
      pc.paths,
      by = "Gene"
    ) %>%
    arrange(
      desc(n), pathway
    ) %>%
    mutate(
      gene_fac = factor(Gene, levels = Gene, ordered = T),
      path_fac = factor(pathway, levels = levels(targ.paths$path_fac), ordered = T)
    )

  p2 <- ggplot(data = targ.action.genes, mapping = aes(x = gene_fac, y = n, fill = path_fac)) +
    geom_col(color = "black") +
    scale_fill_manual(values = pathway.pal, guide = "none") +
    theme_bw() +
    theme(
      text = element_text(family = "Helvetica", size = 12),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ylab("Number of Inhibitors") +
    xlab("TCGA Candidate therapeutic targets")

  # Summarize the drug combinations

  # key the pathway combinations by a unique value that covers both path1->path2 or path2->path1

  combo.drug.paths <- filter(drug.paths, is.na(inhibitor2) == F) %>%
    mutate(
      path_key = mapply(function(x, y) {
        paste(sort(c(x, y)), collapse = ";")
      }, pathway1, pathway2)
    )

  # create a count of pathways by drug

  combo.drug.paths.count <- summarize(
    combo.drug.paths,
    n = n(),
    .by = c("path_key")
  )

  # expand the pathway key into the different pathway combinations

  all.paths <- combo.drug.paths %>%
    {
      union(.$pathway1, .$pathway2)
    }

  expanded.paths <- expand.grid(
    list(pathway1_ro = all.paths, pathway2_ro = all.paths),
    stringsAsFactors = F
  ) %>%
    as_tibble() %>%
    mutate(
      path_key = mapply(function(x, y) {
        paste(sort(c(x, y)), collapse = ";")
      }, pathway1_ro, pathway2_ro)
    )

  # merge to create a symmetric matrix

  combo.drug.paths.count <- left_join(
    expanded.paths,
    combo.drug.paths.count,
    by = c("path_key")
  ) %>%
    mutate(
      n = ifelse(is.na(n), 0, n),
      path1_fac = factor(pathway1_ro, levels = c("PI3K", "RTK RAS", "NOTCH")),
      path2_fac = factor(pathway2_ro, levels = c("PI3K", "RTK RAS", "NOTCH"))
    )

  p3 <- ggplot(data = combo.drug.paths.count, mapping = aes(x = path1_fac, y = path2_fac, fill = n)) +
    scale_fill_gradient(low = "white", high = "purple", guide = "none") +
    geom_tile(color = "black") +
    geom_text(mapping = aes(label = n)) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica", size = 12)) +
    xlab("") +
    ylab("") +
    ggtitle("Pathways targeted by combinations")

  ggsave(p1 + p2 + p3, file = "figures/summary_of_targets.pdf", width = 13, height = 3.5)

  "figures/summary_of_targets.pdf"
}

#' Generate a plot displaying the combination ratio
#'
#' The combination ratio is defined as:
#' log(combination AUC / min(inhibitor1 AUC, inhibitor2 AUC)).  We test its
#' significance using a paired T-test comparing:
#' log(combination AUC) - log(min(inhibitor1 AUC, inhibitor2 AUC)).  Note that
#' the AUCs are rescaled to be between 0-1 to address any differences in
#' concentration. Note that only the combinations that map to pathways are
#' displayed.
#'
#' @param inhib A `tibble` that contains columns for:
#' * `lab_id` Patient sample identifier
#' * `plot_name` Shortened/simplified name used in plots etcs.
#' * `type` Type of drug (e.g. single-agent vs combination)
#' * `inhibitor1` Name of the first drug (or single agent)
#' * `inhibitor2` Name of the second drug in a combination
#' * `rescaled_auc` The AUC value divided by the log10 concentration range.
#' @param drug.paths A `tibble` of drug -> pathway mappings as produced by
#'   [annotate.drugs.pathways()]
#' @returns A PDF figure: 'figures/combination_ratio.pdf'
combination.ratio.plot <- function(inhib, drug.paths) {
  combos <- filter(inhib, type == "combination")
  sas <- filter(inhib, type == "single-agent")

  combos.w.sas <- left_join(
    combos,
    select(sas, lab_id, inhibitor1 = inhibitor, rescaled_auc1 = rescaled_auc),
    by = c("lab_id", "inhibitor1")
  )

  combos.w.sas <- left_join(
    combos.w.sas,
    select(sas, lab_id, inhibitor2 = inhibitor, rescaled_auc2 = rescaled_auc),
    by = c("lab_id", "inhibitor2")
  )

  combo.summary <- summarize(
    combos.w.sas,
    broom::tidy(
      t.test(
        x = log(rescaled_auc),
        y = log(pmin(rescaled_auc1, rescaled_auc2, na.rm = T)),
        paired = T
      )
    ),
    .by = c("inhibitor", "plot_name")
  )

  combo.summary <- left_join(combo.summary, drug.paths, by = "inhibitor")

  combo.summary <- mutate(
    combo.summary,
    pathway = mapply(function(x, y) {
      paste(sort(c(x, y)), collapse = " + ")
    }, pathway1, pathway2)
  )

  combo.summary <-
    filter(
      combo.summary,
      pathway != ""
    ) %>%
    arrange(
      pathway, estimate
    ) %>%
    mutate(
      name_fac = factor(plot_name, levels = plot_name, ordered = T),
      star_level = case_when(
        p.value < .001 ~ "***",
        p.value < .01 ~ "**",
        p.value < .05 ~ "*",
        .default = "N.S."
      )
    )

  comb.plot <- ggplot(data = combo.summary, mapping = aes(y = estimate, x = name_fac, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange(size = .25, color = "purple") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_text(mapping = aes(y = .5, label = star_level, vjust = ifelse(star_level == "N.S.", .5, .75)), angle = 90) +
    facet_grid(. ~ pathway, scales = "free_x", space = "free") +
    scale_color_discrete(name = "Significant") +
    theme_bw() +
    ylab("log Combination Ratio") +
    xlab("") +
    theme(
      text = element_text(family = "Helvetica", size = 12),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ylim(c(-1.25, .75))

  ggsave(comb.plot, file = "figures/combination_ratio.pdf", width = 9, height = 2.5)

  "figures/combination_ratio.pdf"
}
