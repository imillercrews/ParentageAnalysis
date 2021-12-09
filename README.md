# Parentage analysis pipeline

Isaac Miller-Crews, imillercrews@utexas.edu

The parentage analysis pipeline has been described in: https://www.fsigenetics.com/article/S1872-4973(21)00127-7/fulltext
Miller-Crews, I., Matz, M. V., & Hofmann, H. A. (2021). A 2b-RAD parentage analysis pipeline for complex and mixed DNA samples. Forensic Science International: Genetics, 55. https://doi.org/10.1016/j.fsigen.2021.102590

All sample data from the paper can be found on NCBI Bioproject RJNA754415 and the file 'renaming.bams' can be used to convert from SRA file name to file name used in scripts.
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA754415

## WOPI

The function for Weighted Outlier Paternity Index, along with example running through the samples, is in the R script 'WOPI_2bRAD_ParentageAnalysis.R'. TO run in the example, load in included 'WOPI.example.data.RData'. This command requires genotyping probability data to be in long format across both sites and genotypes (ex. AA, AB, BB). See below for an example data frame. It is recommended to remove sites with 'low info scores', that have equal probability to be any genotype at a specific site. For example, 'site.1' for 'Male1' in the below example.

**Variable descriptions**

Sample.name: Column with unique names for each sample

Sex: Column used to identify which type of sample is included: KnownParent, AllegedParent, Child. In the above example this would be 'F', 'M', 'B', respectively.

Marker: Unique site marker name

Allele: One of the genotypes for each site that a genotype probability can be assigned and need to be in the format: 'AA', 'AB', and 'BB'. 

GenotypeProbability: The probability that a specific individual is a certain genotype at a given site. 

*Optional
Only needed if user wishes to supply their own allele frequency data instead of calculating from supplied data. 

AlleleMarker: Denotes the site markers for supplied allele frequency 

Afreq: A allele frequency at a specific site

Bfreq: B allele frequency at a specific site

**Example data**

| Sample.name  | Sex | Marker  | Allele | GenotypeProbability |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| Male1  | M | site.1 | AA  | 0.33  |
| Male1  | M | site.1 | AB  | 0.33  |
| Male1  | M | site.1 | BB  | 0.33  |
| Male2  | M | site.1 | AA  | 0.25  |
| Male2  | M | site.1 | AB  | 0.50  |
| Male2  | M | site.1 | BB  | 0.25  |
| Female1  | F | site.1 | AA  | 0  |
| Female1  | F | site.1 | AB  | 0.25  |
| Female1  | F | site.1 | BB  | 0.75  |
| Child  | B | site.1 | AA  | 0  |
| Child  | B | site.1 | AB  | 0.25  |
| Child  | B | site.1 | BB  | 0.75  |
| Male1  | M | site.2 | AA  | 0.90  |
| Male1  | M | site.2 | AB  | 0.10  |
| Male1  | M | site.2 | BB  | 0  |
| Male2  | M | site.2 | AA  | 0.90  |
| Male2  | M | site.2 | AB  | 0.10  |
| Male2  | M | site.2 | BB  | 0  |
| Female1  | F | site.2 | AA  | 0  |
| Female1  | F | site.2 | AB  | 0.05  |
| Female1  | F | site.2 | BB  | 0.95  |
| Child  | B | site.2 | AA  | 0  |
| Child  | B | site.2 | AB  | 0.10  |
| Child  | B | site.2 | BB  | 0.90  |

*Missing sites and NA will be ignored 

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
           <summary> Converting genotype likelihood file for use with WPOI </summary>
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
           <p>An attempt at a more traditional CPI paternity test utilizing the function PaternityIndex can be found in the R scripts 'Phase 1 Tank Examples Pat geno likelihoods filtered.R' and 'Phase 2 Tank Examples Pat geno likelihoods filtered.R' found in the folder 'CPI.test'.</p>
         </details>      
<details>
           <summary> Heterozygosity Script </summary>
           <p>To calculate heterozygosity statistics from bam files use 'README_2bRAD_ParentageAnalysis_Heterozygosity.txt'. </p>
         </details>      
