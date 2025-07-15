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

  ## Quantify Number of Peptides ----
  ### Quantify the number of unique peptides quantified for each protein group within each sample.
  ### This will use unique peptides (not precursors) as the number of peptide IDs. This number also includes modified versions of the same peptide.
  n_peptide_data <- filtered_diann_report %>%
    dplyr::distinct(.data$Run, ## Subset only the columns we will need for future steps
                    .data$Protein.Group,
                    .data$Genes,
                    .data$First.Protein.Description,
                    .data$Precursor.Id,
                    .data$Modified.Sequence) %>%
    dplyr::group_by(.data$Run,
                    .data$Protein.Group,
                    .data$Genes,
                    .data$First.Protein.Description) %>% # We will count based upon the grouping of run name, protein group, and gene name.
    dplyr::summarise(value = dplyr::n(), # Count the total number of unique peptides in each group.
                     .groups = 'drop') %>% # Remove the grouping that results from summarizing.
    janitor::clean_names() %>%
    dplyr::mutate(variable = 'Number.Peptides') # Add a new column called "variable" with each row containing "Number.Peptides". This will become useful for when we combine MaxLFQ and number of peptides into one data frame.

  ## Calculate MaxLFQ ----
  ### Calculate MaxLFQ from our filtered precursors using variables defined in the function arguments:
  ###  Value: "maxlfq_quantity" (default is "Precursor.Normalised")
  ###  Group: "maxlfq_group" (default is "Protein.Group")
  ###  ID: "maxlfq_id" (default is "Precursor.Id")
  ### Prepare the data from MaxLFQ calculations ----
  #### Create a new data frame that contains the necessary columns for calculating MaxLFQ and put it into the format required for the MaxLFQ function from the "diann" package.
  maxlfq_input <- filtered_diann_report %>%
    dplyr::select(c(.data$Run,
                    .data$Precursor.Id,
                    .data$Genes,
                    .data$Protein.Group,
                    maxlfq_quantity,
                    .data$Modified.Sequence)) %>%
    dplyr::rename('File.Name' = .data$Run)

  ### Calculate MaxLFQ ----
  #### This is done using the MaxLFQ formula from the "diann" R package.
  #### This results in a matrix.
  maxlfq_matrix <- diann::diann_maxlfq(maxlfq_input,
                                       group.header = maxlfq_group,
                                       id.header = maxlfq_id,
                                       quantity.header = maxlfq_quantity)

  ### Convert the MaxLFQ matrix into a data frame ----
  #### Since the MaxLFQ functions results in a matrix, we need to convert it to
  #### a data frame to make it compatible with our n_peptide_data data frame.
  maxlfq_data <- maxlfq_matrix %>%
    as.data.frame() %>% # Converts the matrix to a data frame.
    tibble::rownames_to_column(var = 'Protein.Group') %>% # Converts the row names to a column called "Protein.Group"
    tibble::as_tibble() %>% # Converts the data frame to a more friendly tibble.
    tidyr::pivot_longer(cols = 2:tidyr::last_col(), # Pivots the data frame to a long format, where you have a column for each protein group in each sample.
                        names_to = 'File.Name',
                        values_to = 'MaxLFQ') %>%
    dplyr::inner_join(maxlfq_input %>% # Combine with the input file to obtain gene names.
                        dplyr::distinct(.data$File.Name,
                                        .data$Protein.Group,
                                        .data$Genes),
                      by = c('File.Name',
                             'Protein.Group'),
                      na_matches = 'na') %>%
    janitor::clean_names() %>%
    dplyr::rename('run' = .data$file_name,
                  'value' = .data$max_lfq) %>%
    dplyr::mutate(variable = 'MaxLFQ') # Add a new column called "variable" with each row containing "MaxLFQ". This will become useful for when we combine MaxLFQ and number of peptides into one data frame.

  ## Combine Peptide Count and MaxLFQ Data Frames ----
  ### Combine the peptide count and MaxLFQ data frames together. Then make the
  ### data frame wide so that each row is a different protein group and each
  ### column is either the number of peptides or the MaxLFQ value of a given run.
  data_wide <- dplyr::bind_rows(x = n_peptide_data %>%
                                  dplyr::select(-.data$first_protein_description), # Append the MaxLFQ data frame to the end of the Peptide Count data frame.
                                y = maxlfq_data) %>%
    tidyr::pivot_wider(names_from = .data$variable, # Make the data frame wide by the "variable" column. This will create two new columns: Number.Peptides and MaxLFQ. But each row will still be identified by the protein group and run name.
                       values_from = .data$value) %>%
    tidyr::pivot_wider(names_from = .data$run, # Make the data frame wide by the "run" column. This will create a new column for each run with the values of Number.Peptides and MaxLFQ. Now each row is a different protein group and each column is either the number of peptides or the MaxLFQ value of a given run.
                       values_from = c(.data$Number.Peptides,
                                       .data$MaxLFQ)) %>%
    dplyr::left_join(y = n_peptide_data %>%
                       dplyr::distinct(.data$protein_group,
                                       .data$genes,
                                       .data$first_protein_description),
                     by = c('protein_group',
                            'genes')) %>%
    dplyr::relocate(.data$first_protein_description,
                    .after = 'protein_group') %>%
    dplyr::rename('Protein.Accession' = .data$protein_group,
                  'Protein.Name' = .data$first_protein_description,
                  'Gene.Name' = .data$genes) %>%
    dplyr::rename_at(.vars = vars(dplyr::starts_with('MaxLFQ')),
                     .funs = ~sub('MaxLFQ-',
                                  '',
                                  .))

  # Write the Data to Files ----
  readr::write_tsv(data_wide,
                   report_out)
}
