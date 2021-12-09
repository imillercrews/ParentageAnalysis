#### Converting genotype likelihood beagle file into long format
#Set working directory
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity Testing/2bRAD_data/Data")
#load library
library(tidyverse)
library(paternity)
library(parallel)
library(reshape2)

#### convert genotype beagle file into long format ####
#load genotype likelihood data
data.geno=read.delim('ddB.geno.likelihood.beagle')

data.geno[c(1:10),c(1:10)]
#load name data
data.pat.names=read.csv('bams.ind.csv')
data.pat.names.phase.1=data.pat.names
data.pat.names=data.pat.names[,-c(1:3,6:12,15:18)]
data.pat.names[]=lapply(data.pat.names, function(x) as.character(x))
data.pat.names.phase.1=subset(data.pat.names.phase.1,Phase=='1')
data.pat.names.phase.1=data.pat.names.phase.1[,-c(1:3,7:12,15:18)]
data.pat.names.phase.1[]=lapply(data.pat.names.phase.1, function(x) as.character(x))

#Load pos data
data.pos=read.table('filtered2.pos')
colnames(data.pos)=c('chr','pos')
data.pos=data.pos[-1,]
data.pos$marker=paste(data.pos$chr,data.pos$pos,sep='_')

#Filter data.geno to sites in data.pos
data.geno.filter=merge(data.geno,data.pos,by='marker',all=F)
data.geno.filter=subset(data.geno.filter,select=-c(chr,pos))

#Load names
data.bam.names.geno=read.csv('bams.ind.geno.csv')

#Change row names to marker
rownames(data.geno.filter)=data.geno.filter$marker
data.geno.filter=data.geno.filter[,-c(1:3)]
data.geno.filter.trans=t(data.geno.filter)

#Merge dataframes with names
rownames(data.bam.names.geno)=data.bam.names.geno$File.name.ind.geno
data.bam.names.geno=data.bam.names.geno[,-c(1,8)]
data.geno.filter.trans=merge(data.bam.names.geno,data.geno.filter.trans,by='row.names')

#save file out
save(data.geno.filter.trans,file='data.geno.filter.trans.RData')

#### check filter options ####
#load data file
load('data.geno.filter.trans.RData')
     
#check hist to see distribution and what cutoff number to pick 
data.geno.filter.trans.melt=melt(data.geno.filter.trans,id=c(colnames(data.geno.filter.trans)[1:8]))
ggplot(subset(data.geno.filter.trans.melt,value>=0.00001))+geom_histogram(binwidth=0.01,aes(value))+scale_x_continuous(name='Genotype filter options',breaks = seq(0,1, by = 0.1))
ggsave('Genotype filter options histogram.jpeg')

#create matrix of filter options
Filters.type=rep(seq(0.6:1,by=0.05),times=5)
length(seq(0.6:1,by=0.05))
Age.type=c(rep(c('Adult','Brood'),each=9,times=2),rep('Total',times=9))
Phase.type=c(rep(c(1,2),each=9,times=1),rep(c(2,1),each=9,times=1),rep('Total',times=9))
filter.options=data.frame(Filters.type=Filters.type,Phase.type=Phase.type,Age.type=Age.type)
length(which(subset(data.geno.filter.trans.melt,Sex!='B'&Phase==1)$value>0.95))

filter.options$value=mapply(function(x,y,z) ifelse(z=='Total',length(which(data.geno.filter.trans.melt$value>x)),
                                                   ifelse(z=='Adult'&y==1,length(which(subset(data.geno.filter.trans.melt,Sex!='B'&Phase==y)$value>x)),
                                                          ifelse(z=='Adult'&y==2,length(which(subset(data.geno.filter.trans.melt,Sex!='B'&Phase==y)$value>x)),ifelse(z=='Brood'&y==1,length(which(subset(data.geno.filter.trans.melt,Sex=='B'&Phase==y)$value>x)),ifelse(z=='Brood'&y==2,length(which(subset(data.geno.filter.trans.melt,Sex=='B'&Phase==y)$value>x)),NA)))))
  ,x=filter.options$Filters.type
  ,y=filter.options$Phase.type
  ,z=filter.options$Age.type
                            )
filter.options$ID=paste(filter.options$Age.type,filter.options$Phase.type,sep=".")
#graph filter options
ggplot(filter.options,aes(x=Filters.type,y=value,group=ID,color=ID))+geom_line()+geom_point()+ggtitle('Comparison of genotype filter threshold count')+ylab('counts')+xlab('Genotype filter options')
ggsave('Comparison of genotype filter threshold count.jpeg')

#filter option percentages
filter.options.2=filter.options
filter.options.2$Percent=ifelse(filter.options$ID=='Adult.1',filter.options$value/max(subset(filter.options,ID=='Adult.1')$value),ifelse(filter.options$ID=='Adult.2',filter.options$value/max(subset(filter.options,ID=='Adult.2')$value),ifelse(filter.options$ID=='Brood.1',filter.options$value/max(subset(filter.options,ID=='Brood.1')$value),ifelse(filter.options$ID=='Brood.2',filter.options$value/max(subset(filter.options,ID=='Brood.2')$value),ifelse(filter.options$ID=='Total.Total',filter.options$value/max(subset(filter.options,ID=='Total.Total')$value),NA)))))

#graph filter option percentages
ggplot(filter.options.2,aes(x=Filters.type,y=Percent,group=ID,color=ID))+geom_line()+geom_point()+ggtitle('Comparison of genotype filter threshold percentage')+ylab('Percentage')+xlab('Genotype filter options')+theme( panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour ="black"),panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "black"))
ggsave('Comparison of genotype filter threshold percentage.jpeg')


#### Convert likelihoods into major/minor allele calls ####
#Can change ifelse(x>=#... to a higher value to be more stringent
#use above code to graph out filters
load('data.geno.filter.trans.RData')

data.geno.filter.trans$MajorMinor=as.character(data.geno.filter.trans$MajorMinor)
data.geno.filter.trans[9:ncol(data.geno.filter.trans)]=mclapply(data.geno.filter.trans[9:ncol(data.geno.filter.trans)], function(x) ifelse(x>=0.6,data.geno.filter.trans$MajorMinor,""),mc.cores = 16)

#Convert genotype likelihoods to character
data.geno.filter.trans[9:ncol(data.geno.filter.trans)]=lapply(data.geno.filter.trans[9:ncol(data.geno.filter.trans)], function(x) as.character(x))
data.geno.filter.trans=subset(data.geno.filter.trans,select=-c(Row.names,MajorMinor))

#Collapse into one row
#could you do this with lapply?
data.geno.filter.trans.ind=aggregate(data.geno.filter.trans[9:ncol(data.geno.filter.trans)], list(data.geno.filter.trans$Sample.name), function(x) paste(unique(x),collapse = ''))

#merge with names
data.bam.names=read.csv('bams.ind.csv')
data.bam.names.ind=data.bam.names[,c(4:6,12:14)]
colnames(data.geno.filter.trans.ind)[1]="Sample.name"
data.geno.filter.trans.ind.merge=merge(data.bam.names.ind,data.geno.filter.trans.ind,by='Sample.name')

#save data
save(data.geno.filter.trans.ind.merge,file='data.geno.filter.trans.ind.merge.RData')
load('data.geno.filter.trans.ind.merge.RData')
