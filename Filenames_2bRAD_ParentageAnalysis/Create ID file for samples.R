####Creating ID file for samples
#Set working directory
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity Testing/2bRAD_data/Data")
#load library
library(tidyverse)
library(Biostrings)
library(seqRFLP)
#load sample data
data=read.csv("2bRAD_IMC_samples.csv")
head(data)
#create column to separate
data$Brood.date=NA
data$Sample.name=as.character(data$Sample.name)
data$Sample.name.2=data$Sample.name
#subset adults and separate sample.name
data.adults=subset(data,Plate=='Adult')
head(data.adults)
data.adults.2=separate(data.adults,col = Sample.name.2,into = c('Sample','Tank'),sep = '[.]')
#subset brood and separate sample.name
data.brood=subset(data,Plate=='Brood')
head(data.brood)
data.brood.2=separate(data.brood,col = Sample.name.2,into = c('Tank','Sample','Brood.date'),sep = '[.]',extra = 'merge')
#recombine into adults and broods
data.2=rbind(data.adults.2,data.brood.2)
#Reverse complement of primer codes
data.2$ILLRAD.barcode.sequence.RC=sapply(data.2$ILLRAD.barcode.sequence,revComp)
data.2$X3ill.barcode.sequence.RC=substr(data.2$X3ill.barcode.sequence,1,4)
data.2$X3ill.barcode.sequence.RC=sapply(data.2$X3ill.barcode.sequence.RC,revComp)
#create file name for each sample
data.2$ILLRAD.barcode.lane.number=rep(49:60,each=12)
data.2$ILLRAD.barcode.lane.number=paste('S',data.2$ILLRAD.barcode.lane.number,sep='')
data.2$File.name=paste(data.2$ILLRAD.barcode.sequence.RC,data.2$ILLRAD.barcode.lane.number,data.2$X3ill.barcode.sequence.RC,sep = '_')
data.2$File.name=paste(data.2$File.name,'.trim.bt2.bam',sep = '')
#create names file for total counts
bam.names=read.delim('bams',header=FALSE)
colnames(bam.names)='File.name'
bam.names$File.name.ind=0:(nrow(bam.names)-1)
bam.names$File.name.ind=paste('ind',bam.names$File.name.ind,sep='')
bam.names$File.name.ind=paste(bam.names$File.name.ind,'TotDepth',sep='')
data.bam.names=merge(data.2,bam.names,by='File.name')
write.csv(data.bam.names,file = 'bams.ind.csv')
#create names file for bases counts
Bases=rep(c('A','C','G','T'),times=144)
Number.position=rep(0:143,each=4)
ind.df=data.frame(Bases=Bases,Number.position=Number.position)
ind.df$File.name.ind.bases=paste(ind.df$Number.position,ind.df$Bases,sep='_')
ind.df$File.name.ind.bases=paste('ind',ind.df$File.name.ind.bases,sep='')
ind.df$File.name.ind=rep(bam.names$File.name.ind,each=4)
data.bam.names.bases=merge(data.bam.names,ind.df,by='File.name.ind')
write.csv(data.bam.names.bases,file = 'bams.ind.bases.csv')
#create names file for genotype
data.bam.names=read.csv('bams.ind.csv')
data.bam.names.geno=data.bam.names[,c(4:6,12:14,18)]
data.bam.names.geno$File.name.ind.geno=0:(nrow(data.bam.names.geno)-1)
data.bam.names.geno$File.name.ind.geno=paste('Ind',data.bam.names.geno$File.name.ind.geno,sep='')
Geno.df=data.frame(rep(c("",1,2),times=144),rep(data.bam.names.geno$Sample.name,each=3))
colnames(Geno.df)=c("MajorMinor","Sample.name")
data.bam.names.geno.2=merge(data.bam.names.geno,Geno.df,by='Sample.name')
data.bam.names.geno.2$File.name.ind.geno=ifelse(data.bam.names.geno.2$MajorMinor=="",data.bam.names.geno.2$File.name.ind.geno,paste(data.bam.names.geno.2$File.name.ind.geno,data.bam.names.geno.2$MajorMinor,sep='.'))
data.bam.names.geno.2$MajorMinor=ifelse(data.bam.names.geno.2$MajorMinor=="",'AA',data.bam.names.geno.2$MajorMinor)
data.bam.names.geno.2$MajorMinor=ifelse(data.bam.names.geno.2$MajorMinor=="2",'AB',data.bam.names.geno.2$MajorMinor)
data.bam.names.geno.2$MajorMinor=ifelse(data.bam.names.geno.2$MajorMinor=="3",'BB',data.bam.names.geno.2$MajorMinor)
data.bam.names.geno.2=data.bam.names.geno.2[,-c(7)]
write.csv(data.bam.names.geno.2,file = 'bams.ind.geno.csv')









