#' DIA-NN Report Generator
#'
#' Creates a .tsv report from the 'report.tsv' from DIA-NN
#'
#' @import dplyr
#' @import readr
#' @import stringr
#' @import tidyr
#' @import tibble
#' @importFrom magrittr  %>%
#'
#' @name diann_reporter
#'
#' @param report_in Name of the report.tsv from DIA-NN
#' @param report_out Name of the output report
#'
#' @return A tsv file with filtered protein abundances and number of peptides
#'
#' @export
diann_reporter <- function(report_in,
                           report_out) {

  filtered_diann_report <- readr::read_tsv(report_in) %>%
    dplyr::filter(.data$Lib.Q.Value <= 0.01 &
                    .data$Lib.PG.Q.Value <= 0.01 &
                    .data$Quantity.Quality > 0 &
                    .data$PEP <= 0.2 &
                    !stringr::str_detect(.data$First.Protein.Description,
                                         'Keratin'))

  n_peptides <- filtered_diann_report %>%
    dplyr::distinct(.data$Run,
                    .data$Protein.Group,
                    .data$Genes,
                    .data$First.Protein.Description,
                    .data$Precursor.Id,
                    .data$Modified.Sequence) %>%
    dplyr::group_by(.data$Run,
                    .data$Protein.Group,
                    .data$Genes,
                    .data$First.Protein.Description) %>%
    dplyr::summarise(value = dplyr::n(),
                     .groups = 'drop') %>%
    dplyr::mutate(variable = 'Number.Peptides')

  maxlfq <- filtered_diann_report %>%
    dplyr::distinct(.data$Run,
                    .data$Protein.Group,
                    .data$Genes,
                    .data$First.Protein.Description,
                    .data$PG.MaxLFQ) %>%
    dplyr::rename('value' = .data$PG.MaxLFQ) %>%
    dplyr::mutate(variable = 'MaxLFQ')

  data_merge <- dplyr::bind_rows(x = n_peptides,
                                 y = maxlfq) %>%
    tidyr::pivot_wider(names_from = .data$variable,
                       values_from = .data$value)

  data_wide <- dplyr::full_join(
    x = data_merge %>%
      dplyr::select(-.data$MaxLFQ) %>%
      tidyr::pivot_wider(names_from = .data$Run,
                         values_from = .data$Number.Peptides,
                         names_prefix = 'Number.Peptides_'),
    y = data_merge %>%
      dplyr::select(-.data$Number.Peptides) %>%
      tidyr::pivot_wider(names_from = .data$Run,
                         values_from = .data$MaxLFQ),
    by = c('Protein.Group',
           'Genes',
           'First.Protein.Description')
  ) %>%
    dplyr::rename('Protein.Accession' = .data$Protein.Group,
                  'Protein.Name' = .data$First.Protein.Description)

  readr::write_tsv(data_wide,
                   report_out)
}
