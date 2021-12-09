# Parentage analysis pipeline

Isaac Miller-Crews, imillercrews@utexas.edu

The parentage analysis pipeline has been described in: https://www.fsigenetics.com/article/S1872-4973(21)00127-7/fulltext
Miller-Crews, I., Matz, M. V., & Hofmann, H. A. (2021). A 2b-RAD parentage analysis pipeline for complex and mixed DNA samples. Forensic Science International: Genetics, 55. https://doi.org/10.1016/j.fsigen.2021.102590

All sample data from the paper can be found on NCBI Bioproject RJNA754415 and the file 'renaming.bams' can be used to convert from SRA file name to file name used in scripts.
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA754415

## WOPI

* Script
* Example dataframe

## Supplemental Scripts

<details>
           <summary> 2bRAD mapping commands: from fastq to bams </summary>
           <p>Using command line script 'Bams_2bRAD_ParentageAnalysis.txt' goes from 2bRAD sequencing data as fastq files to output bams, with either a reference or reference-free de novo approach.  </p>
         </details>

<details>
           <summary> ANGSD commands: from bams to genotype probabilities and IBS matrix </summary>
           <p>Using command line script 'ANGSD_2bRAD_ParentageAnalysis.txt' goes from bam files to output genotype probabilities and IBS matrix using ANGSD.  </p>
         </details>
         
<details>
           <summary> Converting genotype liklihood file for use with WPOI </summary>
           <p>In the folder 'Filenames_2bRAD_ParentageAnalysis' are the R scripts and all files needed to generate appropriate sample files. This is mainly important for creating the 'data.geno.filter.trans.RData' from the genotype probabilities output from ANGSD in beagle file format. The first R script 'Create ID file for samples.R' is needed to create the 'bams.ind.csv' and 'bams.ind.geno.csv'. These files along with output from the ANGSD commands, 'filtered2.pos'and 'ddB.geno.likelihood.beagle', are used by the R script 'Pat geno likelihoods filtered.R' to generate the 'data.geno.filter.trans.RData' used in the WOPI example. To create the data used in the WOPI example from 'data.geno.filter.trans.RData' use the R script 'Create files for WOPI.R'.</p>
         </details>

<details>
           <summary> IBS known triads Script </summary>
           <p>In the folder 'IBS_2bRAD_ParentageAnalysis' are the R scripts and all files, including IBS matrices, needed to generate clustered IBS matrices for known triad samples. </p>
         </details>
         
<details>
           <summary> IBS community Script </summary>
           <p>In the folder 'IBS_community_scripts' are the R scripts and all files, including IBS matrices, needed to generate clustered IBS matrices for each naturalistic community. </p>
         </details>
         
<details>
           <summary> Weighted CPI Script </summary>
           <p>An attempt at a more traditional CPI paternity test utilizing the function PaternityIndex can be found in the R scripts 'Phase 1 Tank Examples Pat geno likelihoods filtered.R' and 'Phase 2 Tank Examples Pat geno likelihoods filtered.R'</p>
         </details>      
<details>
           <summary> Heterozygosity Script </summary>
           <p>To calculate heterozygosity statistics from bam files use 'README_2bRAD_ParentageAnalysis_Heterozygosity.txt'. </p>
         </details>      
