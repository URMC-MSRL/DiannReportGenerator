#' DIA-NN Report Generator
#'
#' Creates a .tsv report from the 'report.tsv' from DIA-NN
#'
#' @import dplyr
#' @import readr
#' @import stringr
#' @import tidyr
#' @import tibble
#' @import magrittr
#'
#' @name diann_reporter
#'
#' @param report_in Name of the report.tsv from DIA-NN
#' @param report_out Name of the output report.
#'
#' @return A tsv file with filtered protein abundances and number of peptides
#'
#' @examples
#' diann_report(report_in = 'report.tsv',report_out = 'diann_report.tsv')
#'
#' @export
diann_reporter <- function(report_in,
                           report_out) {

  filtered_diann_report <- readr::read_tsv(report_in) %>%
    dpltr::filter(Lib.Q.Value <= 0.01 &
                    Lib.PG.Q.Value <= 0.01 &
                    Quantity.Quality > 0 &
                    PEP <= 0.2 &
                    !stringr::str_detect(First.Protein.Description,
                                         'Keratin'))

  n_peptides <- filtered_diann_report %>%
    dplyr::distinct(Run,
                    Protein.Group,
                    Genes,
                    First.Protein.Description,
                    Precursor.Id,
                    Modified.Sequence) %>%
    dplyr::select(-Precursor.Id) %>%
    dplyr::group_by(Run,
                    Protein.Group,
                    Genes,
                    First.Protein.Description,
                    Precursor.Id,
                    Modified.Sequence) %>%
    dplyr::select(-Precursor.Id) %>%
    dplyr::group_by(Run,
                    Protein.Group,
                    Genes,
                    First.Protein.Description) %>%
    dplyr::summarise(value = dplyr::n()) %>%
    dplyr::mutate(variable = 'Number.Peptides')

  maxlfq <- filtered_diann_report %>%
    dplyr::distinct(Run,
                    Protein.Group,
                    Genes,
                    First.Protein.Description,
                    PG.MaxLFQ) %>%
    dplyr::rename('value' = PG.MaxLFQ) %>%
    dplyr::mutate(variable = 'MaxLFQ')

  data_merge <- dplyr::bind_rows(x = n_peptides,
                                 y = maxlfq) %>%
    tidyr::pivot_wider(names_from = variable,
                       values_from = value)

  data_wide <- dplyr::full_join(
    x = data_merge %>%
      dplyr::select(-MaxLFQ) %>%
      tidyr::pivot_wider(names_from = Run,
                         values_from = Number.Peptides,
                         names_prefix = 'Number.Peptides_'),
    y = data_merge %>%
      dplyr::select(-Number.Peptides) %>%
      tidyr::pivot_wider(names_from = Run,
                         values_from = MaxLFQ),
    by = c('Protein.Group',
           'Genes')
  ) %>%
    dplyr::remame('Protein.Accession' = Protein.Group,
                  'Protein.Name' = First.Protein.Description)

  readr::write_tsv(data_wide,
                   report_out)
}
