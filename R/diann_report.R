#' URMC MSRL DIA-NN Report Generator
#'
#' Using the 'report.tsv' generated from a DIA-NN search, creates a new file that includes MaxLFQ and number of peptides in wide format
#'
#' @import dplyr
#' @import tidyr
#' @import readr
#' @import stringr
#' @import tibble
#'
#' @name diann_reporter

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

   filtered_diann_report <- readr::read_tsv('report.tsv') %>%
    dplyr::filter(Lib.Q.Value <= 0.01 &
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
                    First.Protein.Description) %>% 
    dplyr::summarise(value = n()) %>% 
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
  
  data_wide <- full_join(
    x = data_merge %>% 
      dplyr::select(-MaxLFQ) %>% 
      tidyr::pivot_wider(names_from = Run,
                         values_from = `Number.Peptides`,
                         names_prefix = 'Number.Peptides_'),
    y = data_merge %>% 
      dplyr::select(-Number.Peptides) %>% 
      tidyr::pivot_wider(names_from = Run,
                         values_from = MaxLFQ),
    by = c('Protein.Group',
           'Genes')
  ) %>% 
    dplyr::rename('Protein.Accession' = Protein.Group,
                  'Protein.Name' = First.Protein.Description)
  
  readr::write_tsv(data_wide, report_out)
}
