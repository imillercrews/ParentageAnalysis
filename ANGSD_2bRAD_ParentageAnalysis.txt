#### Parentage Analysis 2bRAD
### From raw fastq files to genotype probabilities, IBS matrix, and Heterozygosity
### Following 2bRAD directions from: https://github.com/z0on/2bRAD_GATK/blob/master/2bRAD_README.txt

# =============================
### Installations


#need to load before bowtie
module load intel/17.0.4
#load in bowtie2
module load bowtie/2.3.2
#load in samtools
module load samtools

#others
module load perl
module load picard-tools

------- cd-hit:

git clone https://github.com/weizhongli/cdhit.git
cd cd-hit
make

##ANGSD: 

# install xz first from https://tukaani.org/xz/
cd
wget https://tukaani.org/xz/xz-5.2.3.tar.gz --no-check-certificate
tar vxf xz-5.2.3.tar.gz 
cd xz-5.2.3/
./configure --prefix=$HOME/xz-5.2.3/
make
make install

# edit .bashrc:
nano .bashrc
   export LD_LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LD_LIBRARY_PATH
   export LIBRARY_PATH=$HOME/xz-5.2.3/lib:$LIBRARY_PATH
   export C_INCLUDE_PATH=$HOME/xz-5.2.3/include:$C_INCLUDE_PATH
#logout

# re-login

# now, install htslib:
cd
git clone https://github.com/samtools/htslib.git
cd htslib
make CFLAGS=" -g -Wall -O2 -D_GNU_SOURCE -I$HOME/xz-5.2.3/include"
#if error, try running below code before re-running make command
# git submodule update --init --recursive

#install angsd
cd
git clone https://github.com/ANGSD/angsd.git 
cd angsd
make HTSSRC=../htslib

# now adding ANGSD to $PATH
cd
nano .bashrc
# section 2:
   export PATH=$HOME/angsd:$PATH
   export PATH=$HOME/angsd/misc:$PATH
# save (Ctl-O, Ctl-X)
#logout and login afterwards

export PATH=$HOME/trinityrnaseq-v2.13.2:$PATH

## downloading and installing all 2bRAD scripts in $HOME/bin (or change to whatever directory you want)
cd
mkdir bin 
cd ~/bin 
# cloning github repositories
git clone https://github.com/z0on/2bRAD_denovo.git
git clone https://github.com/z0on/2bRAD_GATK.git
# move scripts to ~/bin from sub-directories
mv 2bRAD_denovo/* . 
mv 2bRAD_GATK/* . 
# remove now-empty directory
rm -rf 2bRAD_denovo 
rm -rf 2bRAD_GATK 

# designating all .pl and .py files (perl and python scripts) as executable
chmod +x *.pl 
chmod +x *.py
chmod +x *.R

# adding ~/bin to your $PATH 
cd
nano .bashrc
# paste this where appropriate (note: .bashrc configuration might be specific to your cluster, consult your sysadmin if in doubt)
   export PATH=$HOME/bin:$PATH

# Ctl-o, Ctl-x  (to save and exit in nano)
# log out and re-login to make sure .bashrc changes took effect

# does it work?
# try running a script from $HOME:
cd
2bRAD_trim_launch.pl
# if you get "command not found" something is wrong

### checking fastq files
#Count reads in fastq
echo $(cat *.fastq|wc -l)/4|bc
#count files with .fastq ending
ls -dq *fastq | wc -l

# =============================
### demultiplexing pooled fastq files
## concatenating pairs of raw files according to the value given as Index in the file name:
ngs_concat.pl fastq "Pool-(.+)_L00(..)"

##Fastq quality control
module load fastqc
fastqc *.fq
export PATH="/work/projects/BioITeam/ls5/opt/multiqc-1.0:$PATH"
export PYTHONPATH="/work/projects/BioITeam/ls5/lib/python2.7/annab-packages:$PYTHONPATH"
multiqc .
#Open multiqc_report.html by copying to folder and then open with chrome

### Splitting by in-read barcode, deduplicating and quality-filtering the reads
2bRAD_trim_launch_dedup.pl fq > trims
#>>> make cc a launcher job and run it
launcher_creator.py -j trims -n trims -t 5:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch trims.slurm

## do we have expected number of *.tr0 files created?
ls -l *.tr0 | wc -l
#144
# quality filtering using fastx_toolkit (install fastx_toolkit if you don't have this module)
module load fastx_toolkit
ls *.tr0 | perl -pe 's/^(\S+)\.tr0$/cat $1\.tr0 \| fastq_quality_filter -q 20 -p 100 >$1\.trim/' >filt0

# NOTE: run the next line ONLY if your qualities are 33-based 
# (if you don't know just try to see if it works eventually, if you get errors from fastx_toolkit, try the other one):
	cat filt0 | perl -pe 's/filter /filter -Q33 /' > filt
#if you did NOT run the line above, run this one (after removing # symbol):
#	mv filt0 filt

# execute all commands in filt file (serial or parallel using Launcher, if your system allows) 
launcher_creator.py -j filt -n filt -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch filt.slurm
# do we have expected number of *.trim files created?
ls -l *.trim | wc -l
#144

## trimmed and filtered files added to BioProject

#====================
### Mapping to denovo genome
# denovo RAD business (skip to next "#========" if you have genome reference):

# 'uniquing' ('stacking') individual trimmed fastq reads:
ls *.trim | perl -pe 's/^(.+)$/uniquerOne.pl $1 >$1\.uni/' >unii

# execute all commands written to unii...
launcher_creator.py -j unii -n unii -t 1:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch unii.slurm

# Done! do you have .uni for all your samples?... 
ls -l *.uni | wc -l  

# collecting common tags (= major alleles)
# merging uniqued files (set minInd to >10, or >10% of total number of samples, whichever is greater)
#14
mergeUniq.pl uni minInd=14 >all.uniq

# discarding tags that have more than 7 observations without reverse-complement
awk '!($3>7 && $4==0) && $2!="seq"' all.uniq >all.tab

# creating fasta file out of merged and filtered tags:
awk '{print ">"$1"\n"$2}' all.tab > all.fasta

# clustering allowing for up to 3 mismatches (-c 0.91); the most abundant sequence becomes reference
/home1/04315/imc/cdhit/cd-hit-est -i all.fasta -o cdh_alltags.fas -aL 1 -aS 1 -g 1 -c 0.91 -M 0 -T 0  

#------------
# making fake reference genome (of 30 chromosomes) out of major-allele tags
# need bowtie2 and samtools for indexing

concatFasta.pl fasta=cdh_alltags.fas num=20

# formatting fake genome
export GENOME_FASTA=cdh_alltags_cc.fasta
export GENOME_DICT=cdh_alltags_cc.dict 
bowtie2-build $GENOME_FASTA $GENOME_FASTA
samtools faidx $GENOME_FASTA

#-------------
# Mapping reads to de novo reference and formatting bam files 

# for denovo: map reads to fake genome: 
GENOME_FASTA=cdh_alltags_cc.fasta

# mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > maps

launcher_creator.py -j maps -n maps -t 8:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch maps.slurm

#check alignment rates
#two samples have low alignment rates
>alignmentRates
for F in `ls *trim`; do 
M=`grep -E '^[ATGCN]+$' $F | wc -l | grep -f - maps.e* -A 4 | tail -1 | perl -pe 's/maps\.e\d+-|% overall alignment rate//g'` ;
echo "$F.sam $M">>alignmentRates;
done

ls *.sam > sams
cat sams | wc -l  # number should match number of trim files

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
module load samtools
>s2b
for file in *.sam; do
echo "samtools sort -O bam -o ${file/.sam/}.bam $file && samtools index ${file/.sam/}.bam">>s2b;
done

launcher_creator.py -j s2b -n s2b -t 8:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch s2b.slurm
	 
# run all commands listed in s2b file
ls *bam | wc -l  # should be the same number as number of trim files

# BAM files are the input into various genotype calling / popgen programs, this is the main interim result of the analysis. Archive them.


# =============================
### Mapping to reference genomes
### getting reference genome
##Haplochromis burtoni genome
wget ftp://ftp.ncbi.nih.gov/genomes/Haplochromis_burtoni/CHR_Un/hbu_ref_AstBur1.0_chrUn.fa.gz
# assuming we have a fasta file mygenome.fasta and it lives in the directory $WORK/db
export GENOME_FASTA=$WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa
export GENOME_DICT=$WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.dict 
# indexing genome for bowtie2 mapper
bowtie2-build hbu_ref_AstBur1.0_chrUn.fa $GENOME_FASTA
samtools faidx $GENOME_FASTA
export GENOME_DICT=$WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.dict

##Oreochromis niloticus genome
wget -r -l1 --no-parent -A.fa.gz ftp://ftp.ncbi.nih.gov/genomes/Oreochromis_niloticus/Assembled_chromosomes/seq/ 
#concatenate genome
cat *.fa > O_niloticus_genome.fasta
# assuming we have a fasta file mygenome.fasta and it lives in the directory $WORK/db
export GENOME_FASTA=$WORK/db/Oreochromis/O_niloticus_genome.fasta
export GENOME_DICT=$WORK/db/Oreochromis/O_niloticus_genome.dict 
# indexing genome for bowtie2 mapper
bowtie2-build O_niloticus_genome.fasta $GENOME_FASTA
samtools faidx $GENOME_FASTA
export GENOME_DICT=$WORK/db/O_niloticus_genome.dict 
java -jar $TACC_PICARD_T_DIR/picard.jar CreateSequenceDictionary R=$GENOME_FASTA  O=$GENOME_DICT

###Mapping reads to reference (created Oreochromis and burtoni reference folders)
# for reference-based: mapping with --local option, enables clipping of mismatching ends (guards against deletions near ends of RAD tags)
##Oreochromis
export GENOME_FASTA=$WORK/db/Oreochromis/O_niloticus_genome.fasta
export GENOME_DICT=$WORK/db/Oreochromis/O_niloticus_genome.dict 
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > bt2
launcher_creator.py -j bt2 -n bt2 -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch bt2.slurm

##Burtoni
export GENOME_FASTA=$WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa
export GENOME_DICT=$WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.dict 
2bRAD_bowtie2_launch.pl '\.trim$' $GENOME_FASTA > bt2
launcher_creator.py -j bt2 -n bt2B -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch bt2A.slurm

# what are mapping efficiencies? 
cat maps.e*

# find out mapping efficiency for a particular input file (O9.fastq in this case)
# (assuming all input files have different numbers of reads)
grep -E '^[ATGCN]+$' O9.*trim | wc -l | grep -f - maps.e* -A 4 


## Move forward with Burtoni files from this point as they had better mapping efficienry 

ls *.bt2.sam > sams
cat sams | wc -l  # number should match number of trim files

# next stage is compressing, sorting and indexing the SAM files, so they become BAM files:
cat sams | perl -pe 's/(\S+)\.sam/samtools import \$GENOME_FASTA $1\.sam $1\.unsorted\.bam && samtools sort -o $1\.sorted\.bam $1\.unsorted\.bam && java -Xmx5g -jar \$TACC_PICARD_T_DIR\/picard\.jar AddOrReplaceReadGroups INPUT=$1\.sorted\.bam OUTPUT=$1\.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$1 && samtools index $1\.bam/' >s2b
launcher_creator.py -j s2b -n s2b -t 4:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 12

rm *sorted*
ls *bam | wc -l  # should be the same number as number of trim files


# =============================
### ANGSD
# "FUZZY genotyping" with ANGSD - without calling actual genotypes but working with genotype likelihoods at each SNP. Optimal for low-coverage data (<10x).
# install ANGSD first (see Installations section above)
# listing all bam filenames 
ls *bam >bams

###creating list of bams for adults and broods (removing duplicates) [NOTE: all tags file names listed as reverse complement]
##Brood duplicates
# TGTTAG_S59_GTGA.trim.bt2.bam,TGTTAG_S59_TCAG.trim.bt2.bam,TGTTAG_S59_GCTT.trim.bt2.bam,TGTTAG_S59_CTAC.trim.bt2.bam,TGTTAG_S59_TGTC.trim.bt2.bam,TGTTAG_S59_TCAC.trim.bt2.bam,TGTTAG_S59_GACT.trim.bt2.bam
# TTCTAT_S60_GTGA*.bam, TTCTAT_S60_TCAG*.bam, TTCTAT_S60_GCTT*.bam, TTCTAT_S60_CTAC*.bam, TTCTAT_S60_TGTC*.bam, TTCTAT_S60_TCAC*.bam, TTCTAT_S60_GACT*.bam
ls {TGCAAA*.bam,TGGCAC*.bam,TGTTAG*.bam,TTCTAT*.bam} >bams.broods
grep -v -e TGTTAG_S59_GTGA.trim.bt2.bam -e TGTTAG_S59_TCAG.trim.bt2.bam -e TGTTAG_S59_GCTT.trim.bt2.bam -e TGTTAG_S59_CTAC.trim.bt2.bam -e TGTTAG_S59_TGTC.trim.bt2.bam -e TGTTAG_S59_TCAC.trim.bt2.bam -e TGTTAG_S59_GACT.trim.bt2.bam bams.broods > bams.broods.removed.duplicates

##Adult duplicates
# TCACAT_S55_GCTT*.bam, TCACAT_S55_CTAC*.bam, TCACAT_S55_TGTC*.bam, TCACAT_S55_TCAC*.bam, TCACAT_S55_GACT*.bam
# TCTATA_S56_GCTT*.bam, TCTATA_S56_CTAC*.bam, TCTATA_S56_TGTC*.bam, TCTATA_S56_TCAC*.bam, TCTATA_S56_GACT*.bam
ls {GCACCC*.bam,GCAGGA*.bam,GCCGCG*.bam,GGCGGT*.bam,GTATTA*.bam,TACGTG*.bam,TCACAT*.bam,TCTATA*.bam} >bams.adults
grep -v -e TCACAT_S55_GCTT*.bam -e TCACAT_S55_CTAC*.bam -e TCACAT_S55_TGTC*.bam -e TCACAT_S55_TCAC*.bam -e TCACAT_S55_GACT*.bam bams.adults > bams.adults.removed.duplicates
#combine removed duplicates
cat bams.adults.removed.duplicates bams.broods.removed.duplicates > bams.removed.duplicates

##assessing base qualities and coverage depth
# (optional if not going to slurm) entering interactive session, giving all node's memory to one process:
# idev -tpn 1 -N 1

#For burtoni reference
export GENOME_FASTA=$WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa

#for Oreochromis reference
export GENOME_FASTA=$WORK/db/Oreochromis/O_niloticus_genome.fasta

# angsd settings:
# -uniqueOnly 1 :
# -remove_bads 1 :
# -minMapQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa : set reference genome 
# -maxDepth 1000 : highest total depth (sum over all samples) to assess; set to 10x number of samples
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -baq 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa -maxDepth 1000"

# T O   D O : 
# -doQsDist 1 :
# -doDepth 1 :
# -doCounts 1 :
# -dumpCounts 2 : prints depth of each individual (each line corresponds to same line in position file)
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

# in the following line, -r argument is one chromosome or contig to work with (no need to do this for whole genome as long as the chosen chromosome or contig is long enough)
# (look up lengths of your contigs in the header of *.sam files)
#RUN:
#angsd -b bams -GL 1 $FILTERS $TODO -P 12 -out ddB

#Could run with just broods as well
#angsd -b bams.broods.removed.duplicates  -GL 1 $FILTERS $TODO -P 1 -out ddB.broods

##create angsd command
echo "angsd -b bams -GL 1 -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -baq 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa -maxDepth 1000 -doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -P 12 -out ddB_QC" > ddB_QC
launcher_creator.py -j ddB_QC -n ddB_QC -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch ddB_QC.slurm

##check results for adults and broods 
# summarizing results (using modified script by Matteo Fumagalli)
Rscript ~/bin/plotQC.R ddB_QC
cat ddB_QC.info 
# scp dd.pdf to laptop to look at distribution of base quality scores, fraction of sites in each sample passing coverage thresholds, and fraction of sites passing genotyping rates cutoffs. Use these to guide choices of -minQ,  -minIndDepth and -minInd filters in subsequent ANGSD runs

# =============================
###Allele counts
# (optional if not going to slurm) entering interactive session, giving all node's memory to one process:
# idev -tpn 1 -N 1 -t 00:90:00

# angsd settings:
# -minInd 50 
# -uniqueOnly 1 
# -remove_bads 1  
# -minQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa : set reference genome 
# -setMaxDepth 1440 : highest total depth (sum over all samples) to assess; set to 10x number of samples
#-minInd filter the number of individuals
# -setMinDepth 2 : If the total depth is below this value, the site is discarded
FILTERS="-minInd 50 -uniqueOnly 1 -remove_bads 1 -minQ 20 -baq 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa -setMaxDepth 1440 -setMinDepth 2"

# T O   D O : 
# -dumpCounts 2 : prints depth of each individual (each line corresponds to same line in position file)
# -doQsDist 1 :
# -doDepth 1 :
# -doCounts 1 :
TODO="-dumpCounts 2 -doDepth 1 -doQsDist 1 -doCounts 1"

angsd -b bams  -GL 1 $FILTERS $TODO -P 1 -out ddB.ind

# T O   D O AGAIN : 
# -dumpCounts 4 :prints depth for each of four bases for each individual at each site (each line corresponds to same line in position file) [Should have 4 times as many columns as individuals]
FILTERS="-minInd 50 -uniqueOnly 1 -remove_bads 1 -minQ 20 -baq 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa -setMaxDepth 1440 -setMinDepth 2"
TODO="-dumpCounts 4 -doDepth -doQsDist 1 -doCounts 1"
angsd -b bams  -GL 1 $FILTERS $TODO -P 1 -out ddB.ind.bases


# =============================
### Genotype likelihood
# (optional if not going to slurm) entering interactive session, giving all node's memory to one process:
# idev -tpn 1 -N 1 -t 00:90:00

# angsd settings:
# -minQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -setMaxDepth : highest total depth (sum over all samples) to assess; set to 10x number of samples
# -minInd 50 : filter sites that are not present in 50 individuals
# -setMinDepth 2 : If the total depth is below this value, the site is discarded

FILTERS="-minInd 50 -uniqueOnly 1 -remove_bads 1 -minQ 20 -baq 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa -setMaxDepth 1440 -setMinDepth 2"

# T O   D O : 
# -doGlf 2 : beagle genotype likelihood format. Column 1 (Marker), Col 2 (Major allele: codes as 0=A, 1=C, 2=G, 3=T), Col 3 (Minor allele), Col 4 (Genotype likelihood for the major/major genotype for the first individual), Col 5 (Major/minor genotype likelihood for first individual), Col 6 (minor/minor genotype likelihood for first individual)  
# -doMaf 8 : Allele frequency based on base counts. This method does not rely on genotype likelihood or probabilities but instead infers the allele frequency directly on the base counts.
# -doMajorMinor 2 : Pre specified Major using a reference. You can force the major allele according to the reference states if you have defined those -ref. The minor allele will be inferred based on the genotype likelihood (see do major minor 1). This is the approach used by both GATK and Samtools
TODO="-doGlf 2 -doMajorMinor 4 -doCounts 1"
#ANGSD command
angsd -b bams  -GL 1 $FILTERS $TODO -P 1 -out ddB.geno.likelihood

# =============================
##Reduce reads to variable SNPS
#Could also try -minMaf set to some value…
#Run ANGSD to get counts for eachsite
FILTERS="-minInd 50 -uniqueOnly 1 -remove_bads 1 -minQ 20 -baq 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa -setMaxDepth 1440 -setMinDepth 2"
TODO="-dumpCounts 3 -doCounts 1"
echo "angsd -b bams  -GL 1 $FILTERS $TODO -P 1 -out ddB.geno.likelihood.count" > geno.count
launcher_creator.py -j geno.count -n geno.count -t 4:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch geno.count.slurm

#simply post-filter the output text file made by -dumpCounts 3 using this program, which should select rows in which at least two BAMs have alternative alleles with counts >1: 
#Copy-paste this in nano, save as filterCounts.awk :
{   
    nz=0
    nbams=1+NF/4
    mincount=2 # minimal counts to call a new allele
    for(i=1;i<=NF;i++) 
        if($i>=mincount) nz++
    if (nz>nbams) 
	print 
}

#Run filterCounts.awk
cat ddB.geno.likelihood.count.counts | awk -f filterCounts.awk > filtered.txt


#join count and pos files?
paste ddB.geno.likelihood.count.pos ddB.geno.likelihood.count.counts > ddB.geno.likelihood.count.count.out
#Copy-paste this in nano, save as filterCounts2.awk :
{   
    nz=0
    nbams=1+NF/4
    mincount=2 # minimal counts to call a new allele
    for(i=4;i<=NF;i++) 
        if($i>=mincount) nz++
    if (nz>nbams) 
	print 
}

#Run filterCounts2.awk
cat ddB.geno.likelihood.count.count.out | awk -f filterCounts2.awk > filtered2.txt

# =============================
### Create IBS and covariance matrix
# (optional if not going to slurm) entering interactive session, giving all node's memory to one process:
# idev -tpn 1 -N 1 -t 00:90:00

## angsd settings:
## Filter options
# -GL 1 : samtools model for genotype likelihoods from mapped reads
# -P 1 : number of cores to use
# -minQ 20 : only highly unique mappings (prob of erroneous mapping = 1%)
# -minMapQ 30 : mapping score
# -SNP_pval 2e-6 : only significant snps
# -baq 1 : realign around indels (not terribly relevant for 2bRAD reads mapped with --local option) 
# -setMaxDepth 1440 : highest total depth (sum over all samples) to assess; set to 10x number of samples
#-minInd 50 : filter the number of individuals
# -setMinDepth 2 : If the total depth is below this value, the site is discarded
# -uniqueOnly 1 : remove reads that have multiple best hits
# -remove_bads 1 : remove read with flag above 255
# -minMaf 0.01 : only uses sites with allele freq above number
# -makeMatrix 1 : pairwise IBS matrix, avg. distance between pairs of individuals (*.ibsMat file)
# -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa : uses a. burtoni as reference genome

## To Do options
# -doIBS 1 : random sampled read create IBS matrix (sample single base)
# -doCounts 1 : counts different bases at each position
# -doCov 1 : covariance matrix for PCA (*.covMat file) 
# -doMajorMinor 4 : Pre specified Major using a reference. You can force the major allele according to the reference states if you have defined those -ref. The minor allele will be inferred based on the genotype likelihood (see do major minor 1). This is the approach used by both GATK and Samtools
# -doMaf 8 : Allele frequency based on base counts. This method does not rely on genotype likelihood or probabilities but instead infers the allele frequency directly on the base counts

##run IBS with only adults
#Use this to create file with sites variable across adults
#filters
FILTERS="-minQ 20 -minMapQ 30 -SNP_pval 2e-6 -baq 1 -setMaxDepth 1440 -minInd 10 -setMinDepth 2 -uniqueOnly 1 -remove_bads 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa"
# T O   D O : 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 1 -doMajorMinor 4 -doMaf 1"
#ANGSD
echo "angsd -b bams.adults.removed.duplicates -GL 1 $FILTERS $TODO -P 1 -out Adults.filter/adults.filter " > Adults.filter.cmd
launcher_creator.py -j Adults.filter.cmd -n Adults.filter.cmd -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch Adults.filter.cmd.slurm

# use R to convert .pos file to csv and remove depth row and remove col names to make adults.filter.pos
#remove col headers
#index file
angsd sites index adults.filter.pos

#move files into folder first? 
xargs -a bams.phase1.males.broods cp -t Matrix_phase1_males_broods/

##run IBS for known triads
FILTERS="-makeMatrix 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa"
TODO="-doIBS 1 -doCov 1 -doCounts 1 -doMajorMinor 4"
echo "angsd -sites adults.filter.pos -b bamsphase1 -GL 1 $FILTERS $TODO -P 1 -out Matrix_phase1/ddB.IBS.phase1all" > IBS.phase1all
launcher_creator.py -j IBS.phase1all -n IBS.phase1all -t 0:10:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch IBS.phase1all.slurm


## Run IBS for naturalistic communities
#run IBS for Phase 2 males and then females
#replace tank G2
#males and broods
FILTERS="-makeMatrix 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa"
TODO="-doIBS 1 -doCov 1 -doCounts 1 -doMajorMinor 4"
echo "angsd -sites adults.filter.pos -b bamsphase2g2males -GL 1 $FILTERS $TODO -P 1 -out Matrix_phase2/Tank_G2/IBS.TankG2.males" > IBS.phase2tankG2males
launcher_creator.py -j IBS.phase2tankG2males -n IBS.phase2tankG2males -t 0:10:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch IBS.phase2tankG2males.slurm

#females and broods
FILTERS="-makeMatrix 1 -ref $WORK/db/Burtoni/hbu_ref_AstBur1.0_chrUn.fa"
TODO="-doIBS 1 -doCov 1 -doCounts 1 -doMajorMinor 4"
echo "angsd -sites adults.filter.pos -b bamsphase2g2females -GL 1 $FILTERS $TODO -P 1 -out Matrix_phase2/Tank_G2/IBS.TankG2.females" > IBS.phase2tankG2females
launcher_creator.py -j IBS.phase2tankG2females -n IBS.phase2tankG2females -t 0:10:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1
sbatch IBS.phase2tankG2females.slurm

# =============================
### rename files for SRA upload
## need a file 'names.txt' that is tab seperated with <OldName> <NewName>
## Create new folder with copy of all bams
#may need to run 'dos2unix' on tab seperated naming file
while read a b; do mv "$a" "$b"; done < renaming.bams.txt








