# InToxSa

Scripts and raw data for the publication: 

Hachani A, Giulieri SG, Guérillot R, et al. A novel cytotoxicity screening platform reveals agr-independent mutations in bacteraemia-associated Staphylococcus aureus that promote intracellular persistence. eLife2023;12:e84778 DOI: https://doi.org/10.7554/eLife.84778

Main folders:

0-Functions: contains R functions for standardisation, fitting, QC and outlier detection of raw PI uptake data.

1-plate_info: contains files mapping strain replicates to their position on 96 well plates

2-Raw_data: contains raw PI uptake data exported from the Clariostar

3-Data_parsing: contains R markdown scripts for each experiments 

4-Data_concatenation: scripts used to generate a concatenated dataframe of raw PI uptake data

5-Data_processing: R markdown scripts used to standardise and fit PI uptake data, exclude outlier and plates with low number of JE replicates and selection of genetic pairs with discordand PI uptake phenotype

6-InToxSa_paper: various scripts for GWAS, processing of convergent evolution and generation of figures for the paper

7-Operetta: scripts to process Operetta data

8-eLife_revision: supplementary analyses performed to address the comments by the eLife reviewers
