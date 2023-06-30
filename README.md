# DiannReportGenerator

The goal of DiannReportGenerator is to make a new document that contains, for each sample, the "MaxLFQ Intensity" and "Number of Peptides Quantified" for every Protein Group.

Precursors are filtered by both Library Q-Value and Protein Group Library Q-Value at 0.01 and PEP of 0.5, and keratin IDs are removed. After filtering, MaxLFQ-based quantification (<https://doi.org/10.1074/mcp.M113.031591>) is performed. The number of precursors/peptides that were quantified within a given protein group is counted and the data is then tidied to long format where every row corresponds to a protein group (reported as the FASTA accession number) and the rows contain:

-   Gene Name

-   Protein Name

-   MaxLFQ Intensity for each sample

-   Number of Peptides for each sample

## Installation

Prior to installing this package, you will need to install the latest version of RTools along with the "DiaNN R Package" developed by the Demichev lab.
You can download and install RTools from here:
<https://cran.r-project.org/bin/windows/Rtools/>

DiaNN R Package can be found here:
<https://github.com/vdemichev/diann-rpackage>

You can install DiannReportGenerator from:

<https://github.com/URMC-MSRL/DiannReportGenerator>

In short, instructions are as follows:

    install.packages("devtools")
    devtools::install_github("https://github.com/vdemichev/diann-rpackage")
    devtools::install_github("https://github.com/URMC-MSRL/DiannReportGenerator")
    library(DiannReportGenerator)


## Instructions

    library(DiannReportGenerator)

Put the name of the main DIA-NN report in quotes, followed by whatever you wish to call the new document in quotes.

    diann_reporter("report.tsv", 
                   "new_report.tsv")
