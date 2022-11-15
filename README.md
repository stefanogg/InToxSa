# InToxSa

Scripts and raw data for the publication: A novel cytotoxicity screening platform reveals agr-independent mutations in bacteraemia-associated Staphylococcus aureus that promote intracellular persistence (Hachani et al, 2022)

Main folders:

0-Functions: contains R functions for standardisation, fitting, QC and outlier detection of raw PI uptake data.

1-plate_info: contains files mapping strain replicates to their position on 96 well plates

2-Raw_data: contains raw PI uptake data exported from the Clariostar

3-Data_parsing: contains R markdown scripts for each experiments 

4-Data_concatenation: scripts used to generate a concatenated dataframe of raw PI uptake data

5-Genetic_pairs_analysis_Oct_2020: R markdown scripts used to standardise and fit PI uptake data, exclude outlier and plates with low number of JE replicates and selection of genetic pairs with discordand PI uptake phenotype

6-InToxSa_paper: various scripts for GWAS, processing of convergent evolution and generatio of figures for the paper

7-Operetta: scripts to process Operetta data
