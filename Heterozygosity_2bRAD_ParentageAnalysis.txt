#### Inbreeding estimates 
### running HWE test and individual heterozygosity on bam files

###HWE test
##http://www.popgen.dk/angsd/index.php/HWE_test
##original bams
#all sites
echo "angsd -bam bams.adults.removed.duplicates -doHWE 1 -domajorminor 1 -GL 1 -out adults.all.sites.HWE" > adults.all.sites.HWE

launcher_creator.py -j adults.all.sites.HWE -n adults.all.sites.HWE -t 8:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch adults.all.sites.HWE.slurm

#variable sites only
echo "angsd -bam bams.adults.removed.duplicates  -doHWE 1 -domajorminor 1 -GL 1 -doMaf 1 -SNP_pval 1e-6 -out adults.variable.sites.HWE" > adults.variable.sites.HWE

launcher_creator.py -j adults.variable.sites.HWE -n adults.variable.sites.HWE -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch adults.variable.sites.HWE.slurm

#variable sites only
#filter minInd
echo "angsd -bam bams.adults.removed.duplicates  -doHWE 1 -domajorminor 1 -GL 1 -doMaf 1 -SNP_pval 1e-6 -minInd 90 -out adults.variable.sites.HWE.minind" > adults.variable.sites.HWE.minind

launcher_creator.py -j adults.variable.sites.HWE.minind -n adults.variable.sites.HWE.minind -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch adults.variable.sites.HWE.minind.slurm

##de novo bams
#variable sites only
echo "angsd -bam bams.adults.removed.duplicates  -doHWE 1 -domajorminor 1 -GL 1 -doMaf 1 -SNP_pval 1e-6 -out adults.variable.sites.HWE" > adults.variable.sites.HWE

launcher_creator.py -j adults.variable.sites.HWE -n adults.variable.sites.HWE -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch adults.variable.sites.HWE.slurm

#variable sites only
#filter by sites
#add in p-value 
echo "angsd -bam bams.adults.removed.duplicates  -doHWE 1 -domajorminor 1 -GL 1 -doMaf 1 -sites filter.adults.all.sites -out adults.variable.sites.HWE.filter" > adults.variable.sites.HWE.filter

launcher_creator.py -j adults.variable.sites.HWE.filter -n adults.variable.sites.HWE.filter -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch adults.variable.sites.HWE.filter.slurm

#variable sites only
#filter by sites
echo "angsd -bam bams.adults.removed.duplicates  -doHWE 1 -domajorminor 1 -GL 1 -doMaf 1 -sites filter.adults.all.sites -SNP_pval 1e-6 -out adults.variable.sites.HWE.filter.snp" > adults.variable.sites.HWE.filter.snp

launcher_creator.py -j adults.variable.sites.HWE.filter.snp -n adults.variable.sites.HWE.filter.snp -t 1:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch adults.variable.sites.HWE.filter.snp.slurm

## move files into R
#copy .hwe file (will need to gunzip) for analysis

###Individual heterozygosity
##https://github.com/z0on/2bRAD_denovo/blob/master/2bRAD_README.txt
## note: did not use the 'regions' and 'regions.cut' files, and no '-rf' flag, as mentioned as this led to absurdly long run times
# ------individual heterozygosity, by running angsd -dosaf on individual bams
# using all sites (variable and invariable) - as it should be done
## similar to http://www.popgen.dk/angsd/index.php/Heterozygosity
##using de novo reference
export GENOME_REF=/work2/04315/imc/stampede2/Burtoni_paternity_2brad/2bRAD_data/fastq/cdh_alltags_cc.fasta

### create list of all sites 
# common name prefix for all files to be generated
export NAM=filter.adults.all

# collecting filter-passing sites to work on (all sites, not just variable ones)
# most important filter here is -minInd:  the site must be genotyped in at least that many individuals (guards against DNA contaminations)
FILTERS='-minInd 90 -uniqueOnly 1 -skipTriallelic 1 -doHWE 1 -hetbias_pval 1e-5'
TODO="-doSaf 1 -anc $GENOME_REF -ref $GENOME_REF -doMajorMinor 1 -doMaf 1 -dosnpstat 1 -doPost 1 -doGlf 2"
# list of all sites passing minInd and hetbias filters
echo "angsd -b bams.adults.removed.duplicates -GL 1 -P 4 $FILTERS $TODO -out $NAM" > filter.adults.sites

launcher_creator.py -j filter.adults.sites -n filter.adults.sites -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 1

sbatch filter.adults.sites.slurm

# writing down list of sites and indexing
zcat ${NAM}.mafs.gz | cut -f 1,2 | tail -n +2 >${NAM}.sites
angsd sites index ${NAM}.sites

### run on all sites
# estimating individual heterozygosity for each bam, for filter-passing sites with coverage 10+
export NAM=filter.adults.all
export GENOME_FASTA=/work2/04315/imc/stampede2/Burtoni_paternity_2brad/2bRAD_data/fastq/cdh_alltags_cc.fasta
FILTERS="-sites ${NAM}.sites"

>adults.hets.ind
>bams.adults.removed.duplicates.ind.het
for F in `cat bams.adults.removed.duplicates`; do
echo "angsd -i $F -anc $GENOME_FASTA $FILTERS -GL 1 -doSaf 1 -out ind.${F/.bam/} && realSFS ind.${F/.bam/}.saf.idx >ind.${F/.bam/}.ml && awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ind.${F/.bam/}.ml >>bams.adults.removed.duplicates.ind.het">>adults.hets.ind;
done
# execute all commands in hets; it will produce a tab-delimited table my.hets:

# check if it works (ctl-C the process if you see no errors immediately)
head -1 adults.hets | bash

##bams.adults.removed.duplicates.het
# [bam filename]   [total number of sites passing filters]   [heterozygosity]


#launch command on TACC
launcher_creator.py -j adults.hets.ind -n adults.hets.ind -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 4
sbatch adults.hets.ind.slurm
#took 1 min on 4 nodes...

### run on all sites with additional depth filter
# estimating individual heterozygosity for each bam, for filter-passing sites with coverage 10+
export NAM=filter.adults.all
export GENOME_FASTA=/work2/04315/imc/stampede2/Burtoni_paternity_2brad/2bRAD_data/fastq/cdh_alltags_cc.fasta
FILTERS="-sites ${NAM}.sites -setMinDepth 10 -setMaxDepth 920"

>adults.hets.ind.depth
>bams.adults.removed.duplicates.ind.depth.het
for F in `cat bams.adults.removed.duplicates`; do
echo "angsd -i $F -anc $GENOME_FASTA $FILTERS -GL 1 -doSaf 1 -docounts 1 -out ind.depth.${F/.bam/} && realSFS ind.depth.${F/.bam/}.saf.idx >ind.depth.${F/.bam/}.ml && awk -v file=$F '{print file\"\t\"(\$1+\$2+\$3)\"\t\"\$2/(\$1+\$2+\$3)}' ind.depth.${F/.bam/}.ml >>bams.adults.removed.duplicates.ind.depth.het">>adults.hets.ind.depth;
done
# execute all commands in hets; it will produce a tab-delimited table my.hets:

# check if it works (ctl-C the process if you see no errors immediately)
head -1 adults.hets.ind.depth | bash

##bams.adults.removed.duplicates.het
# [bam filename]   [total number of sites passing filters]   [heterozygosity]


#launch command on TACC
launcher_creator.py -j adults.hets.ind.depth -n adults.hets.ind.depth -t 2:00:00 -e imillercrews@utexas.edu -a NeuroEthoEvoDevo -q normal -N 4
sbatch adults.hets.ind.depth.slurm
#took 1 min on 4 nodes...

##extract files and combine 






























