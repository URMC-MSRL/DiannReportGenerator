#' DIA-NN version 1.8.1 Report Generator
#'
#' Creates a .tsv report from the 'report.tsv' from DIA-NN version 1.8.1. Has
#' more customization options than the base reporter.
#'
#' @import diann
#' @import dplyr
#' @importFrom janitor  clean_names
#' @import readr
#' @import stringr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr  %>%
#'
#' @name diann_181_reporter
#'
#' @param report_in Name of the report.tsv from DIA-NN
#' @param report_out Name of the output report
#' @param qvalue_filter Value to filter on for the Q-Values
#' @param pep_filter Value to filter on for the PEP
#' @param maxlfq_group Column name for the protein group in MaxLFQ calculation
#' @param maxlfq_id Column name for the Precursor ID to use in MaxLFQ calculation
#' @param maxlfq_quantity Column name to use for the quantity in MaxLFQ calculation
#'
#' @return A tsv file with filtered protein abundances and number of peptides
#'
#' @export
diann_181_reporter <- function(report_in = 'report.tsv',
                               report_out = 'diann_report.tsv',
                               qvalue_filter = 0.01,
                               pep_filter = 0.2,
                               maxlfq_group = 'Protein.Group',
                               maxlfq_id = 'Precursor.Id',
                               maxlfq_quantity = 'Precursor.Normalised') {
  # Import the DIA-NN Outputs ----
  ## Read in the report file (the large tsv file).
  diann_report <- readr::read_tsv(report_in,
                                  show_col_types = FALSE)

  # Prepare the Data For Downstream Analyses ----
  ## Filter the data ----
  ### Filter on the following variables:
  ###  1. Library Q-Value
  ###  2. Library Protein Group Q-Value
  ###  3. PEP
  ###  4. Additionally, this removes any Keratin proteins
  filtered_diann_report <- diann_report %>%
    dplyr::filter(.data$Lib.Q.Value <= qvalue_filter &
                    .data$Lib.PG.Q.Value <= qvalue_filter &
                    .data$PG.Q.Value <= qvalue_filter &
                    .data$PEP <= pep_filter &
                    !stringr::str_detect(.data$First.Protein.Description, # Searches the First.Protein.Description column for any Keratin.
                                         'Keratin')) #%>%

    #dplyr::mutate(Run = sub("^[^_]*_([^_]*).*",
    #                        "\\1",
    #                        .data$Run)) # Remove the researcher name and work order number from the run name.
  ## Make a data frame containing protein ID info (UniProt ID, Gene Name, and Protein Name).
  ### This will be used to attach protein names to the final output.
  protein_names <- filtered_diann_report %>%
    dplyr::distinct(.data$Protein.Group,
                    .data$Genes,
                    .data$First.Protein.Description) %>%
    janitor::clean_names()

  ## Quantify Number of Peptides ----
  ### Quantify the number of unique peptides quantified for each protein group within each sample.
  ### This will use unique peptides (not precursors) as the number of peptide IDs. This number also includes modified versions of the same peptide.
  n_peptide_data <- filtered_diann_report %>%
    dplyr::distinct(.data$Run, ## Subset only the columns we will need for future steps
                    .data$Protein.Group,
                    .data$Genes,
                    .data$Precursor.Id,
                    .data$Modified.Sequence) %>%
    dplyr::group_by(.data$Run,
                    .data$Protein.Group,
                    .data$Genes) %>% # We will count based upon the grouping of run name, protein group, and gene name.
    dplyr::summarise(value = dplyr::n(), # Count the total number of unique peptides in each group.
                     .groups = 'drop') %>% # Remove the grouping that results from summarizing.
    janitor::clean_names() %>%
    dplyr::mutate(variable = 'Number.Peptides') # Add a new column called "variable" with each row containing "Number.Peptides". This will become useful for when we combine MaxLFQ and number of peptides into one data frame.

  ## Make a data frame containing MaxLFQ Values ----
  maxlfq_data <- filtered_diann_report %>%
    dplyr::distinct(.data$Run,
                    .data$Protein.Group,
                    .data$Genes,
                    .data$PG.MaxLFQ) %>%
    dplyr::rename('value' = .data$PG.MaxLFQ) %>%
    janitor::clean_names() %>%
    dplyr::mutate(variable = 'MaxLFQ')

  ## Combine Peptide Count and MaxLFQ Data Frames ----
  ### Combine the peptide count and MaxLFQ data frames together. Then make the
  ### data frame wide so that each row is a different protein group and each
  ### column is either the number of peptides or the MaxLFQ value of a given run.
  data_wide <- dplyr::bind_rows(x = n_peptide_data, # Append the MaxLFQ data frame to the end of the Peptide Count data frame.
                                y = maxlfq_data) %>%
    tidyr::pivot_wider(names_from = .data$variable, # Make the data frame wide by the "variable" column. This will create two new columns: Number.Peptides and MaxLFQ. But each row will still be identified by the protein group and run name.
                       values_from = .data$value) %>%
    tidyr::pivot_wider(names_from = .data$run, # Make the data frame wide by the "run" column. This will create a new column for each run with the values of Number.Peptides and MaxLFQ. Now each row is a different protein group and each column is either the number of peptides or the MaxLFQ value of a given run.
                       values_from = c(.data$Number.Peptides,
                                       .data$MaxLFQ)) %>%
    dplyr::left_join(y = protein_names,
                     by = c('protein_group',
                            'genes')) %>%
    dplyr::relocate(.data$first_protein_description,
                    .after = 'genes') %>%
    dplyr::rename('Protein.Accession' = .data$protein_group,
                  'Protein.Name' = .data$first_protein_description,
                  'Genes' = .data$genes) %>%
    dplyr::rename_at(.vars = vars(dplyr::starts_with('MaxLFQ')),
                     .funs = ~sub('MaxLFQ_',
                                  '',
                                  .))

  # Write the Data to Files ----
  readr::write_tsv(data_wide,
                   report_out)
}
