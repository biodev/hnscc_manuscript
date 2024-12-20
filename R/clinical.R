.cat.summary <- function(clin.data, group.var, cat.var) {
  group_by(clin.data, {{ group.var }}, {{ cat.var }}) %>%
    summarize(n = n(), .groups = "drop") %>%
    pivot_wider(id_cols = {{ cat.var }}, names_from = {{ group.var }}, values_from = "n", values_fill = 0) %>%
    rename(name = {{ cat.var }})
}

.generate.clin.table <- function(clin.data, group.var) {
  n.samps <- summarize(clin.data, `Num. Patients` = n(), .by = {{ group.var }}) %>%
    pivot_longer(cols = c("Num. Patients")) %>%
    pivot_wider(id_cols = name, names_from = {{ group.var }}, values_from = value)

  age.var <- summarize(clin.data,
    `Age (Median)` = median(as.numeric(Age)),
    `Age (# < 40)` = sum(as.numeric(Age) < 40),
    .by = {{ group.var }}
  ) %>%
    pivot_longer(cols = c("Age (Median)", "Age (# < 40)")) %>%
    pivot_wider(id_cols = name, names_from = {{ group.var }}, values_from = value)

  gender.var <- .cat.summary(clin.data, {{ group.var }}, Gender)

  race.var <- .cat.summary(clin.data, {{ group.var }}, race_fac)

  alc.var <- .cat.summary(clin.data, {{ group.var }}, alc_fac)

  smok.var <- .cat.summary(clin.data, {{ group.var }}, smok_fac)

  pak.var <- .cat.summary(clin.data, {{ group.var }}, pac_fac)

  pak.mean <- summarize(clin.data,
    `Mean Packyears` = mean(`Pack_Years`, na.rm = T),
    .by = {{ group.var }}
  ) %>%
    pivot_longer(cols = c("Mean Packyears")) %>%
    pivot_wider(id_cols = name, names_from = {{ group.var }}, values_from = value)

  tis.var <- .cat.summary(clin.data, {{ group.var }}, site_fac)

  path.var <- .cat.summary(clin.data, {{ group.var }}, path_fac)

  path.t.var <- .cat.summary(clin.data, {{ group.var }}, path_t_fac)

  path.n.var <- .cat.summary(clin.data, {{ group.var }}, path_n_fac)

  one.year.var <- .cat.summary(clin.data, {{ group.var }}, one_year_fac)

  five.year.var <- .cat.summary(clin.data, {{ group.var }}, five_year_fac)

  # return in list for to facilitate highlighting in resulting excel spreadsheet
  list(
    n.samps, age.var, gender.var, race.var, alc.var, smok.var, pak.var, pak.mean,
    tis.var, path.var, path.t.var, path.n.var, one.year.var, five.year.var
  )
}

#' Create a summary clinical table and write to Excel file
#'
#' The summary table contains patient counts for (hardcoded) clinical categorizations
#' per HPV status (HPV+/-).
#'
#' @param clin.data A `tibble` containing clinical information for each patient
#'   sample including the following columns:
#'   Disease_Site, Race, Alcohol_Use, Smoking_Status, Pack_Years, Path_Final_Stage,
#'   Path_T_Stage, Path_N_Stage, One_Year_RFS, Five_Year_OS, HPV_Status.  Every
#'   other category is highlighted in grey.
#' @returns An Excel file 'outputs/clinical_summary.xlsx'
make.clin.summary.table <- function(clin.data) {
  # add in factors with corresponding levels

  clin.data <- mutate(clin.data,
    site_fac = factor(Disease_Site, levels = c("Oral Cavity", "Larynx", "Oropharynx", "Paranasal Sinus"), ordered = T),
    race_fac = factor(Race, levels = c("Black or African American", "White", "American Indian or Alaska Native"), ordered = T),
    alc_fac = factor(Alcohol_Use, levels = c("Minimal", "Moderate", "Heavy", "Unknown"), ordered = T),
    smok_fac = factor(Smoking_Status,
      levels = c("Lifelong Non-smoker", "Current reformed smoker for > 15 years", "Current reformed smoker for < or = 15 years", "Current smoker"),
      ordered = T
    ),
    pac_fac = factor(case_when(
      is.na(Pack_Years) ~ "Packyears NA",
      Pack_Years <= 10 ~ "Packyears <= 10",
      Pack_Years > 10 ~ "Packyears > 10"
    ), levels = c("Packyears <= 10", "Packyears > 10", "Packyears NA"), ordered = T),
    path_fac = factor(Path_Final_Stage,
      levels = c("1", "2", "3", "4a", "4b"),
      labels = c("Stage I", "Stage II", "Stage III", "Stage IVa", "Stage IVb"), ordered = T
    ),
    path_t_fac = factor(Path_T_Stage,
      levels = c("1", "2", "3", "4a"),
      labels = c("T0-2", "T0-2", "T3-4", "T3-4"),
      ordered = T
    ),
    path_n_fac = factor(Path_N_Stage,
      levels = c("0", "1", "2a", "2b", "2c", "3b"),
      labels = c("N0-1", "N0-1", "N2-3", "N2-3", "N2-3", "N2-3"),
      ordered = T
    ),
    one_year_fac = factor(One_Year_RFS,
      levels = c("Yes", "No", "Unknown"),
      labels = c("1 Year Alive w/o Relapse", "1 Year Deceased/Relapse", "1 Year NA"),
      ordered = T
    ),
    five_year_fac = factor(Five_Year_OS,
      levels = c("Yes", "No", "Unknown"),
      labels = c("Primary HNSCC 5 Years Survived", "Primary HNSCC 5 Years Deceased", "Primary HNSCC 5 Years NA"),
      ordered = T
    ),
    hpv_fac = factor(HPV_Status,
      levels = c("Negative", "Positive"),
      labels = c("HPV-", "HPV+"),
      ordered = T
    )
  )

  hpv.sum.list <- .generate.clin.table(clin.data, hpv_fac)

  # Create a highlighted workbook, adapted from:
  # https://stackoverflow.com/questions/71600719/addstyle-function-in-openxlsx-does-not-fill-cells-in-excel-spreadsheet-with-the
  wb <- createWorkbook()

  # Add a worksheets
  addWorksheet(wb, "Summary", gridLines = FALSE)

  # write data to worksheet 1
  writeData(wb, sheet = 1, bind_rows(hpv.sum.list), rowNames = F)

  # Create style object in order to fill certain cells in the excel spreadsheet.

  Styling_object <- createStyle(fgFill = "lightgrey")

  row.num <- 2

  for (i in seq_along(hpv.sum.list)) {
    cur.n.rows <- nrow(hpv.sum.list[[i]])

    if (i %% 2 == 1) {
      addStyle(wb, sheet = 1, style = Styling_object, rows = row.num:(row.num + cur.n.rows - 1), cols = 1:3, gridExpand = T)
    }

    row.num <- row.num + cur.n.rows
  }

  # save the workbook
  saveWorkbook(wb, "outputs/clinical_summary.xlsx", overwrite = TRUE)

  "outputs/clinical_summary.xlsx"
}
