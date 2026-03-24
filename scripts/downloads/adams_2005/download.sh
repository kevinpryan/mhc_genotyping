#!/bin/bash
# Download gold standard HLA calls NCI-60 (in HTML format, no TSV provided)
cd /hlamajority-paper/external/mhc_genotyping/downloads/pub/adams_2005/
# MHC-I
wget -O MHC_I_calls.html 'https://translational-medicine.biomedcentral.com/articles/10.1186/1479-5876-3-11/tables/2'
# MHC-II
wget -O MHC_II_calls.html 'https://translational-medicine.biomedcentral.com/articles/10.1186/1479-5876-3-11/tables/3'

# Then execute Rscript download_parse_tables.R to generate the RDS file
