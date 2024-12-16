#' Create a mapping `tibble` between the HNSCC and GDSC drugs
#'
#' Produce a mapping table between the HNSCC and GDSC drugs relative to the
#' single-agent drugs in HNSCC via synonym matching.
#'
#' @param inhib A `tibble` that contains columns for:
#' * `type` Type of drug (e.g. single-agent vs combination)
#' * `inhibitor` Name of the drug
#' * `plot_name` Shortened/simplified name used in plots etcs.
#' @param gdsc.drug.file A file containing info on the GDSC's drugs
#'   (e.g. 'screened_compounds_rel_8.5.csv')
#' @returns A `tibble` with the following columns:
#' * `inhibitor` Name of the drug in the HNSCC cohort
#' * `plot_name` Shortened/simplified name used in plots etc. for the HNSCC cohort
#' * `DRUG_NAME` Corresponding drug name in GDSC
map.gdsc.hnscc.drugs <- function(inhib, gdsc.drug.file) {
  gdsc.drugs <- readr::read_delim(gdsc.drug.file)

  gdsc.map <- separate_longer_delim(unique(select(gdsc.drugs, DRUG_NAME, SYNONYMS)), cols = SYNONYMS, delim = ", ")

  gdsc.map <- mutate(gdsc.map, plot_name = toupper(SYNONYMS))

  gdsc.map <- bind_rows(
    select(gdsc.map, DRUG_NAME, plot_name),
    select(mutate(gdsc.drugs, DRUG_NAME, plot_name = toupper(DRUG_NAME)), DRUG_NAME, plot_name)
  ) %>%
    unique()

  sa.inhib <- filter(inhib, type == "single-agent") %>%
    select(inhibitor, plot_name) %>%
    unique()

  sa.inhib <- left_join(sa.inhib, gdsc.map, by = "plot_name")

  sa.inhib
}

#' Collate drug response data from the GDSC prioritizing 'GDSC2'
#'
#' @param drug.map A `tibble` containing the matching drugs between HNSCC and
#'   GDSC as produced by [map.gdsc.hnscc.drugs()]
#' @param gdsc1.file An Excel file of drug responses from GDSC1
#'   (e.g. 'GDSC1_fitted_dose_response_27Oct23.xlsx')
#' @param gdsc2.file An Excel file of drug responses from GDSC2
#'   (e.g. 'GDSC2_fitted_dose_response_27Oct23.xlsx')
#' @returns A `tibble` with the following columns:
#' * `inhibitor` The HNSCC drug name
#' * `plot_name` The HNSCC drug display name for plots etc.
#' * `DRUG_NAME` The GDSC drug name
#' * `DATASET` Either GDSC1 or GDSC2 indicating where the data was derived from
#' * `SANGER_MODEL_ID` Cell line identifier from GDSC
#' * `DRUG_ID` GDSC drug identifier
#' * `AUC` GDSC area under the curve value
get.gdsc.drug.data <- function(drug.map, gdsc1.file, gdsc2.file) {
  drug.map <- drug.map %>%
    filter(is.na(DRUG_NAME) == F)

  gdsc.drugs <- bind_rows(
    select(
      readxl::read_excel(gdsc1.file),
      DATASET, SANGER_MODEL_ID, DRUG_ID, DRUG_NAME, AUC
    ),
    select(
      readxl::read_excel(gdsc2.file),
      DATASET, SANGER_MODEL_ID, DRUG_ID, DRUG_NAME, AUC
    )
  )

  # if data exists for both GDSC1 and 2, choose the highest drug_id for 2
  # assuming versions of drugs improved

  gdsc.drugs <- arrange(gdsc.drugs, desc(DATASET), desc(DRUG_ID)) %>%
    filter(!duplicated(paste(SANGER_MODEL_ID, DRUG_NAME)))

  gdsc.drugs <- inner_join(drug.map, gdsc.drugs, by = c("DRUG_NAME"))

  gdsc.drugs
}

#' Retrieve the HNSCC-relevant cell lines from GDSC
#'
#' @param cl.info.file File containing information on cell lines relative to the
#'   GDSC drug screen
#' @param gdsc.drugs A `tibble` containing information drug response for the GDSC
#'   cell lines as retrieved from [get.gdsc.drug.data()]
#' @returns A `tibble` containing the subset of `cl.info.file` with GDSC cell lines
limit.cell.lines.to.hnsc.drugs <- function(cl.info.file, gdsc.drugs) {
  # have to increase the file first pass due to inaccuracy of the guessing

  cl.info <- readr::read_csv(cl.info.file, guess_max = 2157)

  # limit to those that are head and neck derived

  hn.cls <- filter(cl.info, tissue == "Head and Neck")

  # of these limit to those with at least one overlapping drug

  hn.drug.cls <- filter(hn.cls, model_id %in% gdsc.drugs$SANGER_MODEL_ID)

  hn.drug.cls
}

.form.mut.tbl <- function(mut.mat, sample.name, cohort.name) {
  as_tibble(mut.mat, rownames = "Gene") %>%
    pivot_longer(cols = -Gene, names_to = sample.name, values_to = "num_var") %>%
    mutate(var = ifelse(num_var > 0, "Mut", "None")) %>%
    select(-num_var) %>%
    rename(any_of(setNames("var", cohort.name)))
}

.form.cna.tbl <- function(cna.mat, sample.name, cohort.name) {
  as_tibble(cna.mat, rownames = "Gene") %>%
    pivot_longer(cols = -Gene, names_to = sample.name, values_to = "num_var") %>%
    mutate(var = case_when(
      num_var > 1 ~ "Amp",
      num_var < -1 ~ "Del",
      .default = "None"
    )) %>%
    select(-num_var) %>%
    rename(any_of(setNames("var", cohort.name)))
}

#' Compute a weighted Jaccard distance based on gene alterations

#' Based on the approach from (1), compute a weighted Jaccard distance for the
#' specified gene set between the two sets of mutations and CNAs.
#'
#' @param mut.mat An gene x sample binary matrix of gene mutations from the cohort being studies
#' @param gdsc.mut.mat An gene x sample binary matrix of gene mutations from GDSC
#' @param cna.mat An gene x sample numeric matrix of cna alterations where negative values
#'   indicate deletions and positive values indicate amplifications.  By default abs(value) > 1 are
#'   used as the calls.
#' @param gdsc.cna.mat An gene x sample numeric matrix of cna alterations similar to `cna.mat` but
#'   derived from the GDSC data.
#' @param use.genes A `tibble` of genes to consider for the distance containing columns
#'   for Gene, type and weight where `type` indicates 'Mut', 'Amp' or 'Del'.
#' @returns A list with two elements:
#'   1.) `dist`: A `tibble` containing columns:
#'     lab_id--patient sample identifier
#'     gdsc_id--GDSC sample identifier
#'     mut_dist--resulting distance based on alterations
#'   2.) `data`: A `tibble` containing the alteration overlaps with the following columns:
#'     Gene--gene symbol
#'     lab_id--patient sample identifier
#'     hnscc--alteration type in HNSCC or 'None"
#'     gdsc_id--GDSC sample identifier
#'     gdsc--alteration type in GDSC or 'None'
#'     type--type of alteration, one of 'Mut', 'Amp' or 'Del'.
#'  @references (1) Sinha, R., Luna, A., Schultz, N. & Sander, C. A pan-cancer survey
#'    of cell line tumor similarity by feature-weighted molecular profiles.
#'    Cell Reports Methods 1 (2021). https://doi.org:10.1016/j.crmeth.2021.100039
find.cl.matches <- function(mut.mat, gdsc.mut.mat, cna.mat, gdsc.cna.mat, use.genes) {
  use.muts <- inner_join(
    .form.mut.tbl(mut.mat, "lab_id", "hnscc"),
    .form.mut.tbl(gdsc.mut.mat, "gdsc_id", "gdsc"),
    by = "Gene", relationship = "many-to-many"
  ) %>%
    filter(hnscc == "Mut" | gdsc == "Mut")

  use.muts <- inner_join(use.muts, filter(select(use.genes, Gene, type, weight), type == "Mut"), by = "Gene")

  use.cna <- inner_join(
    .form.cna.tbl(cna.mat, "lab_id", "hnscc"),
    .form.cna.tbl(gdsc.cna.mat, "gdsc_id", "gdsc"),
    by = "Gene", relationship = "many-to-many"
  )

  use.cna <- inner_join(use.cna, filter(select(use.genes, Gene, type, weight), type != "Mut"), by = "Gene")

  use.cna <- filter(
    use.cna,
    (type == "Amp" & (hnscc == "Amp" | gdsc == "Amp")) | (type == "Del" & (hnscc == "Del" | gdsc == "Del"))
  )

  comb.alts <- bind_rows(
    use.muts,
    use.cna
  )

  filter(comb.alts, hnscc == "None" & gdsc == "None") %>%
    {
      stopifnot(nrow(.) == 0)
    }

  mut.alt.dist <- summarize(comb.alts, mut_dist = sum(weight[hnscc != gdsc]) / sum(weight), .by = c("lab_id", "gdsc_id"))

  list(dist = mut.alt.dist, data = comb.alts)
}

.fill.in.missing <- function(use.mat, desired.samples) {
  missing.in.mat <- setdiff(desired.samples, colnames(use.mat))

  # have to set the missing samples to zero to not disrupt the clustering
  missing.mat <- matrix(0,
    nrow = nrow(use.mat),
    ncol = length(missing.in.mat),
    dimnames = list(rownames(use.mat), missing.in.mat)
  )

  use.mat <- cbind(use.mat, missing.mat)

  use.mat
}

#' Plot a heatmap of specified alterations for GDSC and the HNSCC cohort
#'
#' Samples without a given datatype are set to blanck (0).  HNSCC samples are limited
#'   to those with drug data and combined with the GDSC alterations.  Any samples or
#'   genes without at least one corresponding alteration were remove. {ConsensusClusterPlus}(1)
#'   was used to cluster the samples using the hierarchical clustering with Jaccard distance
#'   and split the samples into 8 clusters (currently hardcoded).  Rows were also clustered
#'   using hierarchical clustering with Jaccard distance.
#' @param gdsc.cls A `tibble` containing information on relevant cell lines as returned by
#'   [limit.cell.lines.to.hnsc.drugs()].
#' @param mut.mat An gene x sample binary matrix of gene mutations from the cohort being studies
#' @param gdsc.mut.mat An gene x sample binary matrix of gene mutations from GDSC
#' @param cna.mat An gene x sample numeric matrix of cna alterations where negative values
#'   indicate deletions and positive values indicate amplifications.  By default abs(value) > 1 are
#'   used as the calls.
#' @param gdsc.cna.mat An gene x sample numeric matrix of cna alterations similar to `cna.mat` but
#'   derived from the GDSC data.
#' @param gdsc.assocs A `tibble` of genes to consider for the distance containing columns
#'   for Gene, type and weight where `type` indicates 'Mut', 'Amp' or 'Del'.
#' @param pt.inhib A `tibble` containing a column `lab_id` indicating which patient have drug data.
#' @returns A PDF 'figures/hnscc_gdsc_alt_heatmap.pdf'
#' @references (1) Wilkerson, Matthew D., and D. Neil Hayes. "ConsensusClusterPlus:
#'   a class discovery tool with confidence assessments and item tracking." Bioinformatics 26.12 (2010): 1572-1573.
plot.cl.matches <- function(gdsc.cls, mut.mat, gdsc.mut.mat, cna.mat, gdsc.cna.mat, gdsc.assocs, pt.inhib) {
  gdsc.mut.mat <- gdsc.mut.mat[, gdsc.cls$model_id]
  colnames(gdsc.mut.mat) <- gdsc.cls$model_name

  cna.mat <- cna.mat %>%
    .fill.in.missing(desired.samples = colnames(mut.mat))

  # also limit to samples with drug data

  drug.samps <- intersect(pt.inhib$lab_id, colnames(mut.mat))

  mut.mat <- mut.mat[, drug.samps]
  cna.mat <- cna.mat[, drug.samps]

  gdsc.cna.mat <- gdsc.cna.mat %>%
    .fill.in.missing(desired.samples = gdsc.cls$model_id)
  gdsc.cna.mat <- gdsc.cna.mat[, gdsc.cls$model_id]
  colnames(gdsc.cna.mat) <- gdsc.cls$model_name

  comb.mut <- cbind(mut.mat, gdsc.mut.mat[rownames(mut.mat), ])
  rownames(comb.mut) <- paste0("mut.", rownames(comb.mut))

  amp.genes <- filter(gdsc.assocs, type == "Amp")

  comb.amp <- cbind(cna.mat[amp.genes$Gene, ], gdsc.cna.mat[amp.genes$Gene, ]) > 1
  rownames(comb.amp) <- paste0("amp.", rownames(comb.amp))

  del.genes <- filter(gdsc.assocs, type == "Del")

  comb.del <- cbind(cna.mat[del.genes$Gene, ], gdsc.cna.mat[del.genes$Gene, ]) < -1
  rownames(comb.del) <- paste0("del.", rownames(comb.del))

  comb.mat <- rbind(
    comb.mut,
    comb.amp[, colnames(comb.mut)],
    comb.del[, colnames(comb.mut)]
  )

  comb.mat <- comb.mat[rowSums(comb.mat) > 0, ]
  comb.mat <- comb.mat[, colSums(comb.mat) > 0]

  cc.res <- ConsensusClusterPlus(comb.mat,
    maxK = 15, reps = 1000, pItem = 0.8, pFeature = 1,
    distance = "binary", plot = "pdf", title = "outputs/mut_cc",
    seed = 135252
  )

  clust.res <- cc.res[[8]]$consensusClass[colnames(comb.mat)]

  ht.m <- Heatmap(comb.mat,
    rect_gp = gpar(col = "black"),
    col = circlize::colorRamp2(breaks = c(0, 1), colors = c("white", "purple3")),
    clustering_distance_columns = "binary",
    clustering_distance_rows = "binary",
    clustering_method_columns = "average",
    clustering_method_rows = "average",
    column_split = clust.res,
    column_title = NULL,
    show_heatmap_legend = F,
    row_names_gp = gpar(fontsize = 12, fontfamily = "Helvetica"),
    column_names_gp = gpar(
      fontsize = 12, fontfamily = "Helvetica",
      fontface = ifelse(grepl("\\d{5}", colnames(comb.mat)), "bold", "plain")
    )
  )

  pdf(file = "figures/hnscc_gdsc_alt_heatmap.pdf", width = 8.8, height = 4.5)

  draw(ht.m)

  dev.off()

  "figures/hnscc_gdsc_alt_heatmap.pdf"
}

#' Plot the relationship between alteration distance and correlation between GDSC and HNSCC AUC values
#'
#' @param cl.matches HNSCC-GDSC alteration distance information as returned by [find.cl.matches()]
#' @param pt.inhib Inhibitor data for each HNSCC patient sample with columns for:
#'   `type`--indicating single-agent or combination
#'   `lab_id`--patient sample identifier
#'   `inhibitor`--drug name
#'   `rescaled_auc`--rescaling of AUC values to between 0 and 1
#' @param gdsc.inhib Inhibitor data for each GDSC cell line with columns for:
#'   `SANGER_MODEL_ID`--cell line identifier
#'   `inhibitor`--drug name
#'   `AUC`--a rescaled AUC value that ranges between 0 and 1.
#' @param cl.info A `tibble` containing information on relevant cell lines as returned by
#'   [limit.cell.lines.to.hnsc.drugs()].
#' @returns A PDF 'figures/rel_mut_dist_inh.pdf'
plot.cor.dist.hnscc.cl <- function(cl.matches, pt.inhib, gdsc.inhib, cl.info) {
  sa.inhib <- filter(pt.inhib, type == "single-agent")

  drug.matches <- inner_join(cl.matches$dist, sa.inhib,
    by = "lab_id", relationship = "many-to-many"
  )

  drug.matches <- inner_join(drug.matches,
    select(gdsc.inhib, gdsc_id = SANGER_MODEL_ID, inhibitor, gdsc_auc = AUC),
    by = c("gdsc_id", "inhibitor")
  )

  drug.matches <- inner_join(
    select(cl.info, gdsc_id = model_id, model_name),
    drug.matches,
    by = "gdsc_id"
  )


  dm.cors <- summarize(drug.matches, broom::tidy(cor.test(x = rescaled_auc, y = gdsc_auc, method = "pearson")), .by = c("lab_id", "gdsc_id", "model_name", "mut_dist"))

  p1 <- ggplot(data = dm.cors, mapping = aes(x = mut_dist, y = estimate)) +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm", se = F) +
    facet_wrap(~lab_id, scales = "free") +
    theme_bw() +
    xlab("Alteration Distance") +
    ylab("Correlation")

  ggsave(p1, file = "figures/rel_mut_dist_inh.pdf", width = 10.5, height = 7.75)

  "figures/rel_mut_dist_inh.pdf"
}

#' Plot the relationship between  GDSC and HNSC AUC values
#'
#' Produces scatter plots of the relationship between GDSC and HNSC AUC values
#' faceted by patient samples and the matching cell lines.
#'
#' @param cl.matches HNSCC-GDSC alteration distance information as returned by [find.cl.matches()]
#' @param pt.inhib Inhibitor data for each HNSCC patient sample with columns for:
#'   `type`--indicating single-agent or combination
#'   `lab_id`--patient sample identifier
#'   `inhibitor`--drug name
#'   `rescaled_auc`--rescaling of AUC values to between 0 and 1
#' @param gdsc.inhib Inhibitor data for each GDSC cell line with columns for:
#'   `SANGER_MODEL_ID`--cell line identifier
#'   `inhibitor`--drug name
#'   `AUC`--a rescaled AUC value that ranges between 0 and 1.
#' @param cl.info A `tibble` containing information on relevant cell lines as returned by
#'   [limit.cell.lines.to.hnsc.drugs()].
#' @returns A PDF: 'figures/full_inhib_by_cl.pdf'
plot.hnscc.cl.inh <- function(cl.matches, pt.inhib, gdsc.inhib, cl.info) {
  sa.inhib <- filter(pt.inhib, type == "single-agent") %>%
    mutate(lab_id = as.character(lab_id))

  drug.matches <- inner_join(cl.matches$dist, sa.inhib,
    by = "lab_id", relationship = "many-to-many"
  )

  drug.matches <- inner_join(drug.matches,
    select(gdsc.inhib, gdsc_id = SANGER_MODEL_ID, inhibitor, gdsc_auc = AUC),
    by = c("gdsc_id", "inhibitor")
  )

  drug.matches <- inner_join(
    select(cl.info, gdsc_id = model_id, model_name),
    drug.matches,
    by = "gdsc_id"
  )

  dm.cors <- summarize(drug.matches, broom::tidy(cor.test(x = rescaled_auc, y = gdsc_auc, method = "pearson")), .by = c("lab_id", "model_name", "mut_dist"))

  # for each, choose the model with the best correlation

  dm.cors <- group_by(dm.cors, lab_id) %>%
    slice_min(order_by = mut_dist) %>%
    slice_max(order_by = estimate) %>%
    filter(mut_dist < 1) %>%
    mutate(star_levels = case_when(
      p.value < .001 ~ "***",
      p.value < .01 ~ "**",
      p.value < .05 ~ "**",
      .default = ""
    ))

  drug.matches <- inner_join(drug.matches, select(dm.cors, lab_id, model_name), by = c("lab_id", "model_name"))

  full.plot <- ggplot(data = drug.matches, mapping = aes(x = rescaled_auc, y = gdsc_auc)) +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm", se = F) +
    geom_text(data = dm.cors, mapping = aes(x = .5, y = .15, label = paste0("r: ", round(estimate, 3), star_levels))) +
    facet_wrap(~ lab_id + model_name, ncol = 2) +
    xlab("HNSCC AUC") +
    ylab("GDSC AUC") +
    theme_classic() +
    theme(text = element_text(family = "Helvetica", size = 12))

  ggsave(full.plot, file = "figures/full_inhib_by_cl.pdf", width = 4, height = 6.5)


  "figures/full_inhib_by_cl.pdf"
}
