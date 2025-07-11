#' DIA-NN version 2.0 Report Generator
#'
#' Creates a .tsv report from the 'report.tsv' from DIA-NN version 2.0 or
#'  greater. Has more customization options than the base reporter. Also creates
#'  more files.
#'
#' @importFrom arrow read_parquet
#' @import diann
#' @import dplyr
#' @importFrom janitor  clean_names
#' @import readr
#' @import stringr
#' @import tibble
#' @import tidyr
#' @import tibble
#' @importFrom magrittr  %>%
#'
#' @name diann_202_reporter
#'
#' @param report_in Name of the report.parquet from DIA-NN
#' @param pg_matrix_in Name of the report.pg_matrix.tsv from DIA-NN
#' @param report_out Name of the main output report
#' @param contam_out Name of the contaminant output report
#' @param filtered_out Name of the filtered ID output report
#' @param qc_out Name of the QC output report
#' @param qvalue_filter Value to filter on for the Q-Values
#' @param quality_filter Value to filter on for Quantity and PG MaxLFQ Qualities
#' @param pep_filter Value to filter on for the PEP
#' @param contaminants String to filter on for contaminants
#' @param protein_id_columns The name of the protein identifiers you want to
#'  keep in your output data
#'
#' @return Four tsv files with filtered protein abundances and number of
#'  peptides, contaminants, filtered precursor IDs, and QC metrics
#'
#' @export
diann_202_reporter <- function(report_in = 'report.parquet',
                               pg_matrix_in = 'report.pg_matrix.tsv',
                               report_out = 'diann_report.tsv',
                               contam_out = 'diann_contaminants.tsv',
                               filtered_out = 'diann_poor_quality.tsv',
                               qc_out = 'diann_qc.tsv',
                               qvalue_filter = 0.01,
                               quality_filter = 0,
                               pep_filter = 0.2,
                               contaminants = 'KRT|cRAP',
                               protein_id_columns = c(
                                 'Protein.Group',
                                 'Genes',
                                 'Protein.Names',
                                 'First.Protein.Description')
                               ) {

  contaminants = {{contaminants}}

  # Import the DIA-NN Outputs ----
  ## Read in the report file (the large parquet file).
  diann_report <- arrow::read_parquet(report_in)
  ## Read in the pg matrix file (a tsv file).
  ### This will be used to attach protein names to the final output.
  protein_names <- readr::read_tsv(pg_matrix_in,
                                   show_col_types = FALSE) %>%
    dplyr::select(dplyr::all_of(paste0(protein_id_columns)))

  # Create Secondary Data Frames ----
  ## These will be a series of data frames that will not be used in downstream
  ## analyses. But, we will still keep this data for future reference, or for
  ## data validation.
  ## The types of data that will be created are:
  ## 1. A data frame that contains only contaminants.
  ## 2. A data frame that contains low quality precursors that were removed from analysis.
  ## 3. A data frame that contains various QC metrics.

  ## Contaminant Data Frame ----
  contam_data <- diann_report %>%
    dplyr::filter(stringr::str_detect(.data$Genes,
                                      paste0(contaminants))) %>%
    dplyr::mutate(Run = sub("^[^_]*_([^_]*).*",
                            "\\1",
                            .data$Run))

  ## Low Quality Precursors Data Frame ----
  poor_quality_data <- diann_report %>%
    dplyr::filter(.data$Quantity.Quality <= quality_filter |
                    .data$PG.MaxLFQ.Quality <= quality_filter |
                    .data$Lib.Q.Value < qvalue_filter |
                    .data$Lib.PG.Q.Value < qvalue_filter |
                    .data$PEP < pep_filter |
                    .data$PG.Q.Value <= qvalue_filter) %>%
    dplyr::select(.data$Run,
                  .data$Protein.Group,
                  .data$Protein.Names,
                  .data$Genes,
                  .data$Precursor.Id,
                  .data$Modified.Sequence,
                  .data$Precursor.Charge,
                  .data$Proteotypic,
                  .data$Precursor.Mz,
                  .data$RT,
                  .data$RT.Start,
                  .data$RT.Stop,
                  .data$FWHM,
                  .data$Precursor.Quantity,
                  .data$ Precursor.Normalised,
                  .data$PG.MaxLFQ,
                  .data$Quantity.Quality,
                  .data$PG.MaxLFQ.Quality,
                  .data$PG.Q.Value,
                  .data$PEP,
                  .data$PG.PEP) %>%
  dplyr::mutate(Run = sub("^[^_]*_([^_]*).*",
                          "\\1",
                          .data$Run))

  # Prepare the Data For Downstream Analyses ----
  ## Filter the data ----
  ### Filter on the following variables:
  ###  1. Library Q-Value
  ###  2. Library Protein Group Q-Value
  ###  3. PEP
  ###  4. Additionally, this removes any contaminants listed in the arguments that are found in the "Gene" column.
  filtered_diann_report <- diann_report %>%
    dplyr::filter(.data$Lib.Q.Value <= qvalue_filter &
                    .data$Lib.PG.Q.Value <= qvalue_filter &
                    .data$PEP <= pep_filter &
                    .data$Quantity.Quality > quality_filter &
                    .data$PG.MaxLFQ.Quality > quality_filter &
                    .data$PG.Q.Value <= qvalue_filter &
                    !stringr::str_detect(.data$Genes,
                                         paste0(contaminants))) %>% # Searches the Genes column for any of the strings entered for the "contaminants" argument
    dplyr::mutate(Run = sub("^[^_]*_([^_]*).*",
                            "\\1",
                            .data$Run)) # Remove the researcher name and work order number from the run name.

  ## Quantify Number of Peptides ----
  ### Quantify the number of unique peptides quantified for each protein group within each sample.
  ### This will use unique peptides (not precursors) as the number of peptide IDs. This number also includes modified versions of the same peptide.
  n_peptide_data <- filtered_diann_report %>%
    dplyr::distinct(Run, ## Subset only the columns we will need for future steps
                    .data$Protein.Group,
                    .data$Genes,
                    .data$Precursor.Id,
                    .data$Modified.Sequence) %>%
    dplyr::select(-.data$Precursor.Id) %>%
    dplyr::group_by(Run, # We will count based upon the grouping of run name, protein group, and gene name
                    .data$Protein.Group,
                    .data$Genes) %>%
    dplyr::summarise(value = dplyr::n(), # Count the total number of unique peptides in each group.
                     .groups = 'drop' ) %>% # Remove the grouping that results from summarizing.
    dplyr::mutate(variable = 'Number.Peptides') %>% # Add a new column called "variable" with each row containing "Number.Peptides". This will become useful for when we combine MaxLFQ and number of peptides into one data frame.
    janitor::clean_names()

  ## Make a MaxLFQ Data Frame ----
  maxlfq_data <- filtered_diann_report %>%
    dplyr::distinct(.data$Run,
                    .data$Protein.Group,
                    .data$Genes,
                    .data$PG.MaxLFQ) %>%
    dplyr::rename('value' = .data$PG.MaxLFQ) %>%
    janitor::clean_names() %>%
    dplyr::mutate(variable = 'MaxLFQ') # Add a new column called "variable" with each row containing "MaxLFQ". This will become useful for when we combine MaxLFQ and number of peptides into one data frame.

  ## Combine Peptide Count and MaxLFQ Data Frames ----
  ### Combine the peptide count and MaxLFQ data frames together. Then make the data frame wide so that each row is a different protein group and each column is either the number of peptides or the MaxLFQ value of a given run.
  data_wide <- dplyr::bind_rows(x = n_peptide_data, # Append the MaxLFQ data frame to the end of the Peptide Count data frame.
                                y = maxlfq_data) %>%
    tidyr::pivot_wider(names_from = .data$variable, # Make the data frame wide by the "variable" column. This will create two new columns: Number.Peptides and MaxLFQ. But each row will still be identified by the protein group and run name.
                       values_from = .data$value) %>%
    tidyr::pivot_wider(names_from = Run, # Make the data frame wide by the "run" column. This will create a new column for each run with the values of Number.Peptides and MaxLFQ. Now each row is a different protein group and each column is either the number of peptides or the MaxLFQ value of a given run.
                       values_from = c(.data$Number.Peptides,
                                       .data$MaxLFQ)) %>%
    dplyr::rename('Protein.Group' = .data$protein_group,
                  'Genes' = .data$genes) %>%
    dplyr::left_join(y = protein_names, # Join the protein names data frame to the wide data frame by protein group and gene name. This will give us the protein name for each protein group.
                     by = c('Protein.Group',
                            'Genes')) %>%
    dplyr::rename('Protein.Name' = .data$First.Protein.Description,
                  'Gene.Name' = .data$Genes,
                  'Protein.Accession' = .data$Protein.Group) %>%
    dplyr::relocate(.data$Protein.Name,
                    .after = 'Gene.Name')

  # QC Data Frame ----
  ## We will use this later to create plots.
  qc_data <- filtered_diann_report %>%
    dplyr::distinct(.data$Run,
                    .data$Protein.Group,
                    .data$Protein.Names,
                    .data$Genes,
                    .data$Precursor.Id,
                    .data$Modified.Sequence,
                    .data$Precursor.Charge,
                    .data$Proteotypic,
                    .data$Precursor.Mz,
                    .data$RT,
                    .data$RT.Start,
                    .data$RT.Stop,
                    .data$FWHM,
                    .data$Precursor.Quantity,
                    .data$Precursor.Normalised,
                    .data$PG.MaxLFQ,
                    .data$Quantity.Quality,
                    .data$PG.MaxLFQ.Quality,
                    .data$PG.Q.Value,
                    .data$PEP,
                    .data$PG.PEP) %>%
    dplyr::mutate(Run = sub("^[^_]*_([^_]*).*",
                            "\\1",
                            .data$Run))

  # Write the Data to Files ----
  readr::write_tsv(data_wide,
                   report_out)
  readr::write_tsv(contam_data,
                   contam_out)
  readr::write_tsv(poor_quality_data,
                   filtered_out)
  readr::write_tsv(qc_data,
                   qc_out)

}
