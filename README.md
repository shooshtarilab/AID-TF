# AID-TF
This repository contains the codes for this paper:
"Integrative Data Analysis to Uncover Transcription Factors Involved in Gene Dysregulation of Nine Autoimmune and Inflammatory Diseases".<br>
In the following the description of folders and codes are provided:<br>
**Codes:**<br>
This directory contains all of the codes for the analysis and plotting the results. Each code is described in the following:<br>
**first_test.R**: This code performs the first step analysis which calculates the P values of enrichment of DHS samples for all samples in the dataset for a trait of interest. The inputs to the code are the trait name, the directory of DHS samples, the list of all DHS samples as a .txt file, address of the files having all immunochip data (reference) and Ci-SNPs (test) of all traits, and the output directory. This code calculates the P values of enrichment of all samples for the trait of interest and writes lists containing the list of processed samples, P values of enrichment of samples for the trait of interest, number of test SNPs (Ci-SNPs) and reference SNPs (all immunochip SNPs) matching the effect SNPs of samples. 

