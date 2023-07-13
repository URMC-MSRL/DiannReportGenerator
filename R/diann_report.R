#' URMC MSRL DIA-NN Report Generator
#'
#' Using the 'report.tsv' generating from a DIA-NN search, creates a new file that includes MaxLFQ and number of peptides in wide format
#' Requires the DiaNN R package to be installed from github (https://github.com/vdemichev/diann-rpackage) to run properly.
#' To install follow these steps:
#'     install.packages("devtools")
#'     library(devtools)
#'     isntall_github("https://github.com/vdemichev/diann-rpackage")
#'
#' @import diann
#' @import dplyr
#' @import purrr
#' @import readr
#' @import stringr
#' @import tibble
#' @import tidyr
#'
#' @name diann_reporter

#' @param df a data frame from an upstream step

peptide_count <- function(df) {
   as.data.frame(
      table(
         unique(df[,c(
            'Protein.Group',
            'Modified.Sequence')])$'Protein.Group'
      )
   )
}

#' @param df a data frame from an upstream step

eliminate_duplicate <- function(df) {
   dplyr::distinct(df,
            Protein.Group,
            .keep_all = T
   )
}

#' Reads the report.tsv as input and performs several functions including:
#'     Filters according to Library Q-Value Library Protein Group Q-Value (both 0.1), PEP (0.5), and removes keratin
#'     Performs MaxLFQ algorithm on the new data set
#'     Counts the number of peptides that were quantified for each protein group within each sample
#'     Writes a new report as a tab-separated file that contains in wide format:
#'         Protein Group Uniprot Accession
#'         Gene Name
#'         Protein Name
#'         MaxLFQ Intensity
#'         Number of Peptides Quantified
#'
#' @param report_in Name of report.tsv file generated from DIA-NN. Needs to be in quotes.
#' @param report_out Name of the output report. Can be whatever you want. Needs to be in quotes.
#'
#' @export

diann_reporter <- function(report_in, report_out
                           ) {

   filtered_diann_report <- diann::diann_load(report_in) %>%
      dplyr::filter(
         Lib.Q.Value <= 0.01 &
         Lib.PG.Q.Value <= 0.01 &
         PEP <= 0.5 &
         !stringr::str_detect(First.Protein.Description, 'Keratin')
      )

   maxlfq_matrix <- diann::diann_maxlfq(filtered_diann_report,
                                    group.header = 'Protein.Group',
                                    id.header = 'Precursor.Id',
                                    quantity.header = 'Precursor.Normalised')

   maxlfq_df <- maxlfq_matrix %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = 'Protein.Group') %>%
      tibble::as_tibble() %>%
      tidyr::pivot_longer(cols = 2:(ncol(maxlfq_matrix) + 1),
                          names_to = 'File.Name',
                          values_to = 'MaxLFQ') %>%
      dplyr::inner_join(
         filtered_diann_report %>%
            dplyr::select(-c(7:13, 15:38, 40:last_col()),
         .,
         by = c('File.Name', 'Protein.Group'),
         na_matches = 'na') %>%
      dplyr::select(-c(1, 4:5)) %>%
      dplyr::relocate(First.Protein.Description, .after = Genes) %>%
      tidyr::nest(data = c(Protein.Group:MaxLFQ))

   df_joined <- dplyr::left_join(
      maxlfq_df %>%
         tidyr::unnest(data) %>%
         dplyr::select(-c('Modified.Sequence')) %>%
         dplyr::rename('Sample.Name' = Run),
      maxlfq_df %>%
         dplyr::mutate(n.peptides = purrr::map(data, peptide_count)) %>%
         dplyr::select(c('Run', 'n.peptides')) %>%
         tidyr::unnest(n.peptides)  %>%
         purrr::set_names(
            'Sample.Name',
            'Protein.Group',
            'Number.Peptides'),
      by = c("Sample.Name","Protein.Group")) %>%
      dplyr::group_by(Sample.Name) %>%
      tidyr::nest() %>%
      dplyr::mutate(end.data = purrr::map(data, eliminate_duplicate)) %>%
      dplyr::select(-c(data)) %>%
      tidyr::unnest(end.data)

   fin_df <- dplyr::inner_join(
      df_joined %>%
         dplyr::select(-c(6)) %>%
         tidyr::pivot_wider(names_from = Sample.Name,
                            values_from = MaxLFQ),
      df_joined %>%
         dplyr::select(c(1:2, 6)) %>%
         tidyr::pivot_wider(names_from = Sample.Name,
                            values_from = Number.Peptides,
                            names_prefix = 'Number.Peptides'),
      by = 'Protein.Group') %>%
      dplyr::rename('Protein.Accession' = Protein.Group,
                    'Protein.Name' = First.Protein.Description)

   write_tsv(fin_df, report_out)
}
