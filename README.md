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
           <summary> ANGSD Commands </summary>
           <p>Using command line script 'README_2bRAD_ParentageAnalysis_ANGSD.txt' take raw 2bRAD fastq files to bams and then use ANGSD to output genotype probabilities and IBS matrix.  </p>
         </details>
         
<details>
           <summary> IBS known triads Script </summary>
           <p>In the folder 'IBS_2bRAD_ParentageAnalysis' are the rscripts and all files, including IBS matrices, needed to generate clustered IBS matrices for known triad samples. </p>
         </details>
         
<details>
           <summary> IBS community Script </summary>
           <p>In the folder 'IBS_community_scripts' are the rscripts and all files, including IBS matrices, needed to generate clustered IBS matrices for each naturalistic community. </p>
         </details>
         
<details>
           <summary> Weighted CPI Script </summary>
           <p>Content 1 Content 1 Content 1 Content 1 Content 1</p>
         </details>      
<details>
           <summary> Heterozygosity Script </summary>
           <p>To calculate heterozygosity statistics from bam files use 'README_2bRAD_ParentageAnalysis_Heterozygosity.txt'. </p>
         </details>      
