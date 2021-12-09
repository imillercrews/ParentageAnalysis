#Tank Examples Phase 1 to go with 'Pat geno likelihoods filtered.R'
#Set working directory
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity_Testing/2bRAD_data/Data")
#load library
library(tidyverse)
library(paternity)
library(parallel)
library(Relatedness)

load('data.geno.filter.trans.ind.merge.RData')

###########
#load name data
data.pat.names=read.csv('bams.ind.csv')
data.pat.names.phase.1=data.pat.names
data.pat.names=data.pat.names[,-c(1:3,6:12,15:18)]
data.pat.names[]=lapply(data.pat.names, function(x) as.character(x))
data.pat.names.phase.1=subset(data.pat.names.phase.1,Phase=='1')
data.pat.names.phase.1=data.pat.names.phase.1[,-c(1:3,7:12,15:18)]
data.pat.names.phase.1[]=lapply(data.pat.names.phase.1, function(x) as.character(x))

data.phase.1=subset(data.geno.filter.trans.ind.merge,Phase==1)
data.phase.1=droplevels(data.phase.1)
levels(data.phase.1$Tank)

#transpose so col for markers and then each sample
data.phase.1.sample=data.phase.1[,-c(2:6)]
rownames(data.phase.1.sample)=data.phase.1.sample$Sample.name
data.phase.1.sample=data.phase.1.sample[,-c(1)]
data.phase.1.sample.trans=t(data.phase.1.sample)
data.phase.1.sample.trans=as.data.frame(data.phase.1.sample.trans)
data.phase.1.sample.trans=rownames_to_column(data.phase.1.sample.trans,var="Marker")

#subset into males
data.phase.1.males=subset(data.pat.names.phase.1,Sex=='M')
#subset broods with mothers
data.pat.names.phase.1$Mother=ifelse(data.pat.names.phase.1$Sex=='B',paste(data.pat.names.phase.1$Sample,data.pat.names.phase.1$Tank,sep = '.'),NA)
data.phase.1.brood=subset(data.pat.names.phase.1,Sex=='B')

#create two B5A.BLUE.7.31 with two mothers BLUE_B.B5A and BLUE_A.B5A
data.phase.1.brood[6,6]='BLUE_A.B5A'
#duplicate B5A.BLUE.7.31 row
B5A.BLUE.7.31.row=data.phase.1.brood[6,]
B5A.BLUE.7.31.row[1,6]='BLUE_B.B5A'
data.phase.1.brood=rbind(data.phase.1.brood,B5A.BLUE.7.31.row)

#Remove B5A.BLUE.7.31
#data.phase.1.brood=subset(data.phase.1.brood,Sample.name!='B5A.BLUE.7.31')

########
#Check relatedness
#need 

data(Genotype)
data(Frequencies)
data(Cross)

par(mar=c(1,1,1,1))

RelatednessCoefficient <- RelCoef(IndividualGenom=matrix(0,ncol=0,nrow=0),
                                  ParentalLineGenom=Genotype,
                                  Freq=Frequencies,Crossing=Cross,
                                  ParentPop=rep(1,8),Phased=TRUE,NbCores=2,Details = T)
print(RelatednessCoefficient$Delta3)


RelatednessCoefficient <- RelCoef(IndividualGenom=matrix(0,ncol=0,nrow=0),
                                  Freq=Frequencies,Details = T)
print(RelatednessCoefficient$Delta3)


data.phase.1.cor.test=data.phase.1.sample.trans.2[,4:ncol(data.phase.1.sample.trans.2)]

data.phase.1.cor.test[]=lapply(data.phase.1.cor.test, function(x) ifelse(x=="",NA,x))
#data.phase.1.cor.test[]=lapply(data.phase.1.cor.test, function(x) as.numeric(x))

View(cor(data.phase.1.cor.test,use = "complete.obs"))

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = cor(data.phase.1.cor.test,use = "complete.obs"),col=col, symm = TRUE)

data.phase.1.cor.test=data.phase.1.sample.trans.2[,c(data.phase.1.males$Sample.name)]
data.phase.1.cor.test[]=lapply(data.phase.1.cor.test, function(x) ifelse(x=="",NA,x))
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = cor(data.phase.1.cor.test,use = "complete.obs",method = "spearman"),col=col, symm = TRUE)
View(cor(data.phase.1.cor.test,use = "complete.obs",method = 'spearman'))









##########
#check number of calls for each sample
data.phase.1.sample.counts=mclapply(colnames(data.phase.1.sample.trans), function(x) sum(data.phase.1.sample.trans[[x]]!=""),mc.cores=16)
data.phase.1.sample.counts.2=data.frame(matrix(unlist(data.phase.1.sample.counts),ncol=1))
colnames(data.phase.1.sample.counts.2)=c('Counts')
data.phase.1.sample.counts.2$Sample.name=colnames(data.phase.1.sample.trans)

data.phase.1.sample.counts.3=merge(data.pat.names.phase.1,data.phase.1.sample.counts.2,by='Sample.name')
data.phase.1.sample.counts.3$Percent=data.phase.1.sample.counts.3$Counts/data.phase.1.sample.counts.2$Counts[data.phase.1.sample.counts.2$Sample.name=='Marker']

ggplot(data.phase.1.sample.counts.3)+geom_histogram(aes(Percent,fill=Sex))
#check number of calls for each sample after genotype filter
data.phase.2.sample.counts=mclapply(colnames(data.phase.2.sample.trans.2), function(x) sum(data.phase.2.sample.trans.2[[x]]!=""),mc.cores=16)

data.phase.2.sample.counts.2=data.frame(matrix(unlist(data.phase.2.sample.counts),ncol=1))
colnames(data.phase.2.sample.counts.2)=c('Counts')
data.phase.2.sample.counts.2$Sample.name=colnames(data.phase.2.sample.trans.2)

data.phase.2.sample.counts.3=merge(data.pat.names.phase.2,data.phase.2.sample.counts.2,by='Sample.name')
data.phase.2.sample.counts.3$Percent=data.phase.2.sample.counts.3$Counts/data.phase.2.sample.counts.2$Counts[data.phase.2.sample.counts.2$Sample.name=='Marker']

ggplot(data.phase.2.sample.counts.3)+geom_histogram(aes(Percent,fill=Sex))


###############
###Randomly sample to test for mixed broods
##calculate genotype frequency beforehand with only adults!
#subset data.phase.1.sample.trans to just adults
data.phase.1.sample.trans.adults=data.phase.1.sample.trans[,c(1,7:20,36:45)]

#create length for each marker
phase.1.AA_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='AA')),mc.cores = 16)
phase.1.BB_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='BB')),mc.cores = 16)
phase.1.AB_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='AB')),mc.cores = 16)
#create dataframe
phase.1.AA_count.df=data.frame(matrix(unlist(phase.1.AA_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.BB_count.df=data.frame(matrix(unlist(phase.1.BB_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.AB_count.df=data.frame(matrix(unlist(phase.1.AB_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.genotypes=cbind(phase.1.AA_count.df,phase.1.BB_count.df,phase.1.AB_count.df)
colnames(phase.1.genotypes)=c('AA_count','BB_count','AB_count')
phase.1.genotypes$Total_count=phase.1.genotypes$AA_count+phase.1.genotypes$BB_count+phase.1.genotypes$AB_count
phase.1.genotypes$Afreq=(phase.1.genotypes$AA_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes$Bfreq=(phase.1.genotypes$BB_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes.2=data.frame(phase.1.genotypes$Afreq,phase.1.genotypes$Bfreq)
colnames(phase.1.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.1.sample.trans.2=cbind(phase.1.genotypes.2,data.phase.1.sample.trans)

#remove rows with 100% Major allele freq
data.phase.1.sample.trans.2=subset(data.phase.1.sample.trans.2, Afreq!=1 & Bfreq!=1)

#add create columns of random numbers
#50:50 ratio of 0 to 1 or B4C.PINK.6.7 to B5A.YELLOW.7.6
data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50=rbinom(nrow(data.phase.1.sample.trans.2),1,0.5)
data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50=ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50==0,data.phase.1.sample.trans.2$B4C.PINK.6.7,data.phase.1.sample.trans.2$B5A.YELLOW.7.6)
data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50=ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50==1,"",
                                                                 ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50==2,"AA",
                                                                        ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_50==3,"AB","BB")))

#75:25 ratio of 0 to 1 or B4C.PINK.6.7 to B5A.YELLOW.7.6
data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25=rbinom(nrow(data.phase.1.sample.trans.2),1,0.25)
data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25=ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25==0,data.phase.1.sample.trans.2$B4C.PINK.6.7,data.phase.1.sample.trans.2$B5A.YELLOW.7.6)
data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25=ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25==1,"",
                                                                 ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25==2,"AA",
                                                                        ifelse(data.phase.1.sample.trans.2$B4C.PINK.6.7_B5A.YELLOW.7.6_25==3,"AB","BB")))

#50:50 ratio of 0 to 1 or B5A.BLUE.7.31 to B5A.YELLOW.7.6 
data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50=rbinom(nrow(data.phase.1.sample.trans.2),1,0.5)
data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50=ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50==0,data.phase.1.sample.trans.2$B5A.BLUE.7.31,data.phase.1.sample.trans.2$B5A.YELLOW.7.6)
data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50=ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50==1,"",
                                                                  ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50==2,"AA",
                                                                         ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_50==3,"AB","BB")))

#75:25 ratio of 0 to 1 or B5A.BLUE.7.31 to B5A.YELLOW.7.6
data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25=rbinom(nrow(data.phase.1.sample.trans.2),1,0.25)
data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25=ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25==0,data.phase.1.sample.trans.2$B5A.BLUE.7.31,data.phase.1.sample.trans.2$B5A.YELLOW.7.6)
data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25=ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25==1,"",
                                                                  ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25==2,"AA",
                                                                         ifelse(data.phase.1.sample.trans.2$B5A.BLUE.7.31_B5A.YELLOW.7.6_25==3,"AB","BB")))

#subset data.phase.1.brood to include just B4C.PINK.6.7 & B5A.YELLOW.7.6 & B5A.BLUE.7.31 and random mixed broods
data.phase.1.brood.random=data.phase.1.brood
data.phase.1.brood.random=subset(data.phase.1.brood.random,Sample.name=='B4C.PINK.6.7'|Sample.name=='B5A.YELLOW.7.6'|Mother=='BLUE_A.B5A')
#Create new rows with random brood names for both mothers
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B4C.PINK.6.7_B5A.YELLOW.7.6_50',Sex='B',Phase=1,Sample='PINK',Tank='B4C',Mother='PINK.B4C'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B4C.PINK.6.7_B5A.YELLOW.7.6_50',Sex='B',Phase=1,Sample='YELLOW',Tank='B5A',Mother='YELLOW.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B4C.PINK.6.7_B5A.YELLOW.7.6_25',Sex='B',Phase=1,Sample='PINK',Tank='B4C',Mother='PINK.B4C'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B4C.PINK.6.7_B5A.YELLOW.7.6_25',Sex='B',Phase=1,Sample='YELLOW',Tank='B5A',Mother='YELLOW.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31_B5A.YELLOW.7.6_50',Sex='B',Phase=1,Sample='BLUE',Tank='B5A',Mother='BLUE_A.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31_B5A.YELLOW.7.6_50',Sex='B',Phase=1,Sample='YELLOW',Tank='B5A',Mother='YELLOW.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31_B5A.YELLOW.7.6_25',Sex='B',Phase=1,Sample='BLUE',Tank='B5A',Mother='BLUE_A.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31_B5A.YELLOW.7.6_25',Sex='B',Phase=1,Sample='YELLOW',Tank='B5A',Mother='YELLOW.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31_B5A.YELLOW.7.6_50',Sex='B',Phase=1,Sample='BLACK',Tank='B5A',Mother='BLACK.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31_B5A.YELLOW.7.6_25',Sex='B',Phase=1,Sample='BLACK',Tank='B5A',Mother='BLACK.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.YELLOW.7.6',Sex='B',Phase=1,Sample='BLACK',Tank='B5A',Mother='BLACK.B5A'))
data.phase.1.brood.random=rbind(data.phase.1.brood.random,data.frame(Sample.name='B5A.BLUE.7.31',Sex='B',Phase=1,Sample='BLACK',Tank='B5A',Mother='BLACK.B5A'))

#create new col to identify which tank it is associated with
data.phase.1.brood.random$Sample.name.2=c("B5A.YELLOW.7.6",'B5A.BLUE.7.31',"B4C.PINK.6.7","B4C.B4C.PINK.6.7_B5A.YELLOW.7.6_50", "B5A.B4C.PINK.6.7_B5A.YELLOW.7.6_50","B4C.B4C.PINK.6.7_B5A.YELLOW.7.6_25", "B5A.B4C.PINK.6.7_B5A.YELLOW.7.6_25",'B5A.BLUE.7.31_B5A.YELLOW.7.6_50.B','B5A.BLUE.7.31_B5A.YELLOW.7.6_50.Y','B5A.BLUE.7.31_B5A.YELLOW.7.6_25.B','B5A.BLUE.7.31_B5A.YELLOW.7.6_25.Y','B5A.BLUE.7.31_B5A.YELLOW.7.6_50.BL','B5A.BLUE.7.31_B5A.YELLOW.7.6_25.BL','B5A.YELLOW.7.6.BL','B5A.BLUE.7.31.BL')
#create new list of alleged fathers
data.phase.1.males.random=data.phase.1.males
data.phase.1.males.random=rbind(data.phase.1.males.random,data.frame(Sample.name=c('YELLOW.B5A','BLUE_A.B5A'),Sex='F',Phase=1,Sample=c('YELLOW','BLUE'),Tank='B5A'))

#run paternity test and graph on just random broods
#Loop through all broods and males in a tank to pro6duce table of PI for each brood for each male
data.phase.1.males.2=mapply(function(x, y) data.phase.1.males.random[[x]]=mclapply(data.phase.1.males.random$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                            data.phase.1.brood.random$Sample.name,  #brood names as x
                            data.phase.1.brood.random$Mother) #mother names as y


#convert into matrix
data.phase.1.males.3=data.frame(matrix(unlist(data.phase.1.males.2), nrow=nrow(data.phase.1.males.2)))
#brood names
phase.1.brood.names=c(data.phase.1.brood.random$Sample.name.2)
#rename columns
colnames(data.phase.1.males.3)=phase.1.brood.names
#combine paternity index with male data
data.phase.1.males.PI=cbind(data.phase.1.males.random,data.phase.1.males.3)

#Name father
phase.1.father=lapply(phase.1.brood.names, function(x) ifelse(data.phase.1.males.PI$Tank==substr(x,1,3) &data.phase.1.males.PI$Sex=='M', 1,0))
phase.1.father=data.frame(matrix(unlist(phase.1.father), ncol=length(phase.1.father)))
colnames(phase.1.father)=c(paste(phase.1.brood.names,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI=cbind(data.phase.1.males.PI,phase.1.father)

#Graph
mclapply(phase.1.brood.names, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_random/",x,".jpeg",sep=""),plot=ggplot(data.phase.1.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI[[x]]),y=data.phase.1.males.PI[[x]],fill=data.phase.1.males.PI[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)


#melt males PI data
data.phase.1.males.PI.melt=cbind(data.phase.1.males.random,data.phase.1.males.3)
data.phase.1.males.PI.melt=melt(data.phase.1.males.PI.melt,id=c(colnames(data.phase.1.males.PI.melt)[1:5]))
#create father coloumn
phase.1.father.melt=cbind(data.phase.1.males.PI[,1:5],phase.1.father)
phase.1.father.melt=melt(phase.1.father.melt,id=c(colnames(phase.1.father.melt)[1:5]))
phase.1.father.melt=phase.1.father.melt[,7]
data.phase.1.males.PI.melt=cbind(data.phase.1.males.PI.melt,phase.1.father.melt)

length(unique(subset(data.phase.1.males.PI.melt,value>=0&phase.1.father.melt==0)$variable))
length(subset(data.phase.1.males.PI.melt,phase.1.father.melt==1&value>=7)$variable)


#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 1 PI Distribution Across Fathers')
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_random/Phase 1 PI Distribution Across Fathers.jpeg')

#Create Pat score based off the PI count for brood
data.phase.1.males.PI.melt$value.exp=exp(data.phase.1.males.PI.melt$value)

data.phase.1.brood.total.PI=aggregate(data.phase.1.males.PI.melt$value.exp,list(variable=data.phase.1.males.PI.melt$variable),sum)

colnames(data.phase.1.brood.total.PI)=c('variable','brood.total')

data.phase.1.males.PI.melt.2=merge(data.phase.1.males.PI.melt,data.phase.1.brood.total.PI)

data.phase.1.males.PI.melt.2$Percent.PI=data.phase.1.males.PI.melt.2$value.exp/data.phase.1.males.PI.melt.2$brood.total

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers')

#Filter pat score by percentage and PI count for brood
data.phase.1.males.PI.melt.2$Percent.PI.filter.27=ifelse(data.phase.1.males.PI.melt.2$value>=27,data.phase.1.males.PI.melt.2$Percent.PI,NA)

data.phase.1.males.PI.melt.2$Percent.PI.filter.20=ifelse(data.phase.1.males.PI.melt.2$value>=20,data.phase.1.males.PI.melt.2$Percent.PI,NA)

#Graphing success of paternity metric
ggplot(subset(data.phase.1.males.PI.melt.2))+geom_bar(color='black',stat='identity',aes(x=reorder(variable,Tank), y=Percent.PI.filter.20,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(subset(data.phase.1.males.PI.melt.2))+geom_bar(color='black',stat='identity',aes(x=reorder(variable,Tank), y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_random/Paternity metric success.jpeg')
###############
##calculate genotype frequency beforehand with only adults!
##DOES NOT WORK
#NEW paternityIndex method
#change paternityIndex so every site gets a min score of 0.00000000001 
#change code to:  alldata$PI = with(alldata, COEFA * (1/a) + COEFB * (1/b) + COEFAB * (1/(a + b)) + 0.00000000001)
paternityIndex2=edit(paternityIndex)
#subset data.phase.1.sample.trans to just adults
data.phase.1.sample.trans.adults=data.phase.1.sample.trans[,c(1,7:20,36:45)]

#create length for each marker
phase.1.AA_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='AA')),mc.cores = 16)
phase.1.BB_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='BB')),mc.cores = 16)
phase.1.AB_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='AB')),mc.cores = 16)
#create dataframe
phase.1.AA_count.df=data.frame(matrix(unlist(phase.1.AA_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.BB_count.df=data.frame(matrix(unlist(phase.1.BB_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.AB_count.df=data.frame(matrix(unlist(phase.1.AB_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.genotypes=cbind(phase.1.AA_count.df,phase.1.BB_count.df,phase.1.AB_count.df)
colnames(phase.1.genotypes)=c('AA_count','BB_count','AB_count')
phase.1.genotypes$Total_count=phase.1.genotypes$AA_count+phase.1.genotypes$BB_count+phase.1.genotypes$AB_count
phase.1.genotypes$Afreq=(phase.1.genotypes$AA_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes$Bfreq=(phase.1.genotypes$BB_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes.2=data.frame(phase.1.genotypes$Afreq,phase.1.genotypes$Bfreq)
colnames(phase.1.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.1.sample.trans.2=cbind(phase.1.genotypes.2,data.phase.1.sample.trans)

#remove rows with 100% Major allele freq
data.phase.1.sample.trans.2=subset(data.phase.1.sample.trans.2, Afreq!=1 & Bfreq!=1)

#run paternity test and graph
#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.phase.1.males.2=mapply(function(x, y) data.phase.1.males[[x]]=mclapply(data.phase.1.males$Sample.name,function(z) log(prod(subset(paternityIndex2(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex2(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                            data.phase.1.brood$Sample.name,  #brood names as x
                            data.phase.1.brood$Mother) #mother names as y


#convert into matrix
data.phase.1.males.3=data.frame(matrix(unlist(data.phase.1.males.2), nrow=nrow(data.phase.1.males.2)))
#rename second B5A.BLUE.7.31 sample to B5A.BLUE.7.31_2
colnames(data.phase.1.males.2)[21]='B5A.BLUE.7.31_2'
colnames(data.phase.1.males.3)=colnames(data.phase.1.males.2)

#combine paternity index with male data
data.phase.1.males.PI=cbind(data.phase.1.males,data.phase.1.males.3)

#Tank phase.1 brood names
phase.1.brood.names=names(data.phase.1.males.PI[6:ncol(data.phase.1.males.PI)])

#Name father
phase.1.father=lapply(phase.1.brood.names, function(x) ifelse(data.phase.1.males.PI$Tank==substr(x,1,3) &data.phase.1.males.PI$Sex=='M', 1,0))
phase.1.father=data.frame(matrix(unlist(phase.1.father), ncol=length(phase.1.father)))
colnames(phase.1.father)=c(paste(phase.1.brood.names,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI=cbind(data.phase.1.males.PI,phase.1.father)

#Graph
mclapply(phase.1.brood.names, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_no_patexclusion/",x,".jpeg",sep=""),plot=ggplot(data.phase.1.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI[[x]]),y=data.phase.1.males.PI[[x]],fill=data.phase.1.males.PI[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)


#melt males PI data
data.phase.1.males.PI.melt=cbind(data.phase.1.males,data.phase.1.males.3)
data.phase.1.males.PI.melt=melt(data.phase.1.males.PI.melt,id=c(colnames(data.phase.1.males.PI.melt)[1:5]))
#create father coloumn
phase.1.father.melt=cbind(data.phase.1.males.PI[,1:5],phase.1.father)
phase.1.father.melt=melt(phase.1.father.melt,id=c(colnames(phase.1.father.melt)[1:5]))
phase.1.father.melt=phase.1.father.melt[,7]
data.phase.1.males.PI.melt=cbind(data.phase.1.males.PI.melt,phase.1.father.melt)

length(unique(subset(data.phase.1.males.PI.melt,value>=0&phase.1.father.melt==0)$variable))
length(subset(data.phase.1.males.PI.melt,phase.1.father.melt==1&value>=7)$variable)


#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 1 PI Distribution Across Fathers')
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_no_patexclusion/Phase 1 PI Distribution Across Fathers.jpeg')

#Create Pat score based off the PI count for brood
data.phase.1.males.PI.melt$value.exp=exp(data.phase.1.males.PI.melt$value)

data.phase.1.brood.total.PI=aggregate(data.phase.1.males.PI.melt$value.exp,list(variable=data.phase.1.males.PI.melt$variable),sum)

colnames(data.phase.1.brood.total.PI)=c('variable','brood.total')

data.phase.1.males.PI.melt.2=merge(data.phase.1.males.PI.melt,data.phase.1.brood.total.PI)

data.phase.1.males.PI.melt.2$Percent.PI=data.phase.1.males.PI.melt.2$value.exp/data.phase.1.males.PI.melt.2$brood.total

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers')

#Filter pat score by percentage and PI count for brood
data.phase.1.males.PI.melt.2$Percent.PI.filter=ifelse(data.phase.1.males.PI.melt.2$value>=27,data.phase.1.males.PI.melt.2$Percent.PI,NA)

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI.filter,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers filtered')


#Graphing success of paternity metric
ggplot(subset(data.phase.1.males.PI.melt.2))+geom_bar(color='black',stat='identity',aes(x=variable, y=Percent.PI.filter,fill=phase.1.father.melt))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_no_patexclusion/Paternity metric success.jpeg')

#Creating figure for outcome of every brood
data.phase.1.males.PI.melt.3=subset(data.phase.1.males.PI.melt.2,phase.1.father.melt==1)
#remove replicated fathers
data.phase.1.males.PI.melt.3$Sample.name.2=paste(data.phase.1.males.PI.melt.3$Sample.name,data.phase.1.males.PI.melt.3$variable,sep = '.')
data.phase.1.males.PI.melt.3 = data.phase.1.males.PI.melt.3[order(data.phase.1.males.PI.melt.3[,'variable'],-data.phase.1.males.PI.melt.3[,'Percent.PI.filter']),]
data.phase.1.males.PI.melt.3 = data.phase.1.males.PI.melt.3[!duplicated(data.phase.1.males.PI.melt.3$variable),]

#give pat success
data.phase.1.males.PI.melt.3$Pat.identity=1
data.phase.1.males.PI.melt.3$Pat.success=0
data.phase.1.males.PI.melt.3$Pat.success=ifelse(data.phase.1.males.PI.melt.3$Percent.PI>=0.90&data.phase.1.males.PI.melt.3$value>=26,2,data.phase.1.males.PI.melt.3$Pat.success)
data.phase.1.males.PI.melt.3$Pat.success=ifelse(data.phase.1.males.PI.melt.3$Percent.PI>=0.90&data.phase.1.males.PI.melt.3$value<=26,1,data.phase.1.males.PI.melt.3$Pat.success)
#rename pat success
data.phase.1.males.PI.melt.3$Pat.success.2=ifelse(data.phase.1.males.PI.melt.3$Pat.success==2,'Correct father above threshold',ifelse(data.phase.1.males.PI.melt.3$Pat.success==1,'Correct father below threshold','Incorrect father below threshold'))

ggplot(data.phase.1.males.PI.melt.3)+geom_histogram(color='black',stat='identity',aes(x=Pat.success.2, y=Pat.identity,fill=Pat.success.2))+ theme(axis.text.x = element_blank())+ylab('Number of Broods')+xlab('Paternity testing success')+labs(fill='Known triad paternity')+ theme(legend.position = c(0.4, 0.7))
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log_no_patexclusion/Paternity testing success.jpeg')

###############
##calculate genotype frequency beforehand with only adults!
#subset data.phase.1.sample.trans to just adults
data.phase.1.sample.trans.adults=data.phase.1.sample.trans[,c(1,7:20,36:45)]

#create length for each marker
phase.1.AA_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='AA')),mc.cores = 16)
phase.1.BB_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='BB')),mc.cores = 16)
phase.1.AB_count=mclapply(data.phase.1.sample.trans.adults$Marker, function(x) length(which(data.phase.1.sample.trans.adults[data.phase.1.sample.trans.adults$Marker==x,2:ncol(data.phase.1.sample.trans.adults)]=='AB')),mc.cores = 16)
#create dataframe
phase.1.AA_count.df=data.frame(matrix(unlist(phase.1.AA_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.BB_count.df=data.frame(matrix(unlist(phase.1.BB_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.AB_count.df=data.frame(matrix(unlist(phase.1.AB_count), nrow=nrow(data.phase.1.sample.trans.adults)))
phase.1.genotypes=cbind(phase.1.AA_count.df,phase.1.BB_count.df,phase.1.AB_count.df)
colnames(phase.1.genotypes)=c('AA_count','BB_count','AB_count')
phase.1.genotypes$Total_count=phase.1.genotypes$AA_count+phase.1.genotypes$BB_count+phase.1.genotypes$AB_count
phase.1.genotypes$Afreq=(phase.1.genotypes$AA_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes$Bfreq=(phase.1.genotypes$BB_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes.2=data.frame(phase.1.genotypes$Afreq,phase.1.genotypes$Bfreq)
colnames(phase.1.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.1.sample.trans.2=cbind(phase.1.genotypes.2,data.phase.1.sample.trans)

#remove rows with 100% Major allele freq
#cuts from 23392 markers down to 3349
data.phase.1.sample.trans.2=subset(data.phase.1.sample.trans.2, Afreq!=1 & Bfreq!=1)

#run paternity test and graph
#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.phase.1.males.2=mapply(function(x, y) data.phase.1.males[[x]]=mclapply(data.phase.1.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                            data.phase.1.brood$Sample.name,  #brood names as x
                            data.phase.1.brood$Mother) #mother names as y


#convert into matrix
data.phase.1.males.3=data.frame(matrix(unlist(data.phase.1.males.2), nrow=nrow(data.phase.1.males.2)))
#rename second B5A.BLUE.7.31 sample to B5A.BLUE.7.31_2
colnames(data.phase.1.males.2)[21]='B5A.BLUE.7.31_2'
colnames(data.phase.1.males.3)=colnames(data.phase.1.males.2)

#combine paternity index with male data
data.phase.1.males.PI=cbind(data.phase.1.males,data.phase.1.males.3)

#Tank phase.1 brood names
phase.1.brood.names=names(data.phase.1.males.PI[6:ncol(data.phase.1.males.PI)])

#Name father
phase.1.father=lapply(phase.1.brood.names, function(x) ifelse(data.phase.1.males.PI$Tank==substr(x,1,3) &data.phase.1.males.PI$Sex=='M', 1,0))
phase.1.father=data.frame(matrix(unlist(phase.1.father), ncol=length(phase.1.father)))
colnames(phase.1.father)=c(paste(phase.1.brood.names,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI=cbind(data.phase.1.males.PI,phase.1.father)

#Graph
mclapply(phase.1.brood.names, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log/",x,".jpeg",sep=""),plot=ggplot(data.phase.1.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI[[x]]),y=data.phase.1.males.PI[[x]],fill=data.phase.1.males.PI[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)

####
#create sample graph with B5A.YELLOW.7.6
data.phase.1.males.PI.B5A.YELLOW.7.6=data.phase.1.males.PI[,c(1:5,10,31)]
data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.Father=as.character(data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.Father)
data.phase.1.males.PI.B5A.YELLOW.7.6=arrange(data.phase.1.males.PI.B5A.YELLOW.7.6,-B5A.YELLOW.7.6)
#remove duplicate fathers
data.phase.1.males.PI.B5A.YELLOW.7.6 = data.phase.1.males.PI.B5A.YELLOW.7.6[-c(5,8),]

data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.Father.list=c('Father','Male 1','Male 2','Male 3','Male 4','Male 5','Male 6','Male 7')
data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.log=data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6
data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6=exp(data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6)
data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.log.10=log10(data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6)
 
ggplot(data.phase.1.males.PI.B5A.YELLOW.7.6)+geom_bar(stat='identity',aes(x=reorder(B5A.YELLOW.7.6.Father.list,-B5A.YELLOW.7.6),y=B5A.YELLOW.7.6,fill=B5A.YELLOW.7.6.Father))+ylab('CPI (log scale)')+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5,size = 28),axis.text.y=element_text(size=28),axis.line.x = element_line(colour = 'black', size = 2),axis.line.y = element_line(colour = 'black', size = 2),text = element_text(size=30),axis.ticks = element_line(size=2))+
  guides(fill=FALSE)+geom_hline(yintercept = exp(27), linetype = 'dashed',size=2)+scale_y_log10()+ scale_fill_manual(values=c( "#BDBDBD","#3C51A3"))
ggsave('Poster_sample/CPI.B5A.YELLOW.7.6.pdf',width=12,height = 7.5, dpi=320, units = 'in')


#percent PI
data.phase.1.males.PI.B5A.YELLOW.7.6$Percent.PI=data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6/sum(data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6)

data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.Father.2=ifelse(data.phase.1.males.PI.B5A.YELLOW.7.6$B5A.YELLOW.7.6.Father==1,'Father', 'Other Males')

ggplot(data.phase.1.males.PI.B5A.YELLOW.7.6)+geom_bar(stat='identity',aes(x=reorder(B5A.YELLOW.7.6.Father.2,-B5A.YELLOW.7.6),y=Percent.PI,fill=B5A.YELLOW.7.6.Father))+ylab('Relative CPI')+theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5,size = 28),axis.text.y=element_text(size=28),axis.line.x = element_line(colour = 'black', size = 2),axis.line.y = element_line(colour = 'black', size = 2),text = element_text(size=30),axis.ticks = element_line(size=2))+guides(fill=FALSE)+ scale_y_continuous(labels = scales::percent)+ scale_fill_manual(values=c( "#BDBDBD","#3C51A3"))
ggsave('Poster_sample/PercentPI.B5A.YELLOW.7.6.pdf',width=5.75,height = 5, dpi=320, units = 'in')

###
#melt males PI data
data.phase.1.males.PI.melt=cbind(data.phase.1.males,data.phase.1.males.3)
data.phase.1.males.PI.melt=melt(data.phase.1.males.PI.melt,id=c(colnames(data.phase.1.males.PI.melt)[1:5]))
#create father coloumn
phase.1.father.melt=cbind(data.phase.1.males.PI[,1:5],phase.1.father)
phase.1.father.melt=melt(phase.1.father.melt,id=c(colnames(phase.1.father.melt)[1:5]))
phase.1.father.melt=phase.1.father.melt[,7]
data.phase.1.males.PI.melt=cbind(data.phase.1.males.PI.melt,phase.1.father.melt)

length(unique(subset(data.phase.1.males.PI.melt,value>=0&phase.1.father.melt==0)$variable))
length(subset(data.phase.1.males.PI.melt,phase.1.father.melt==1&value>=7)$variable)


#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 1 PI Distribution Across Fathers')
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log/Phase 1 PI Distribution Across Fathers.jpeg')


#Create Pat score based off the PI count for brood
data.phase.1.males.PI.melt$value.exp=exp(data.phase.1.males.PI.melt$value)

data.phase.1.brood.total.PI=aggregate(data.phase.1.males.PI.melt$value.exp,list(variable=data.phase.1.males.PI.melt$variable),sum)

colnames(data.phase.1.brood.total.PI)=c('variable','brood.total')

data.phase.1.males.PI.melt.2=merge(data.phase.1.males.PI.melt,data.phase.1.brood.total.PI)

data.phase.1.males.PI.melt.2$Percent.PI=data.phase.1.males.PI.melt.2$value.exp/data.phase.1.males.PI.melt.2$brood.total

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers')

#Filter pat score by percentage and PI count for brood
data.phase.1.males.PI.melt.2$Percent.PI.filter=ifelse(data.phase.1.males.PI.melt.2$value>=27,data.phase.1.males.PI.melt.2$Percent.PI,NA)

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI.filter,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers filtered')


#Graphing success of paternity metric
ggplot(subset(data.phase.1.males.PI.melt.2))+geom_bar(color='black',stat='identity',aes(x=variable, y=Percent.PI.filter,fill=phase.1.father.melt))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log/Paternity metric success.jpeg')

##create boxplot of fathers
data.phase.1.males.PI.melt.father=data.phase.1.males.PI.melt.2

data.phase.1.males.PI.melt.father %>% 
  
ggplot(data.phase.1.males.PI.melt.father,aes(y=Percent.PI,x=phase.1.father.melt,fill=phase.1.father.melt==1))+geom_boxplot(outlier.shape = NA)+geom_jitter()+guides(fill=guide_legend(title="Father"))+xlab('Male')+ggtitle('Phase 1 CPI Distribution Across Fathers')+ylab('Percent PI')
  
ggplot(data.phase.1.males.PI.melt.father,aes(y=value.exp,x=phase.1.father.melt,fill=phase.1.father.melt==1))+geom_boxplot(outlier.shape = NA)+geom_jitter()+guides(fill=guide_legend(title="Father"))+xlab('Male')+ggtitle('Phase 1 CPI Distribution Across Fathers')+ylab('CPI')

#Creating figure for outcome of every brood
data.phase.1.males.PI.melt.3=subset(data.phase.1.males.PI.melt.2,phase.1.father.melt==1)

#remove replicated fathers
data.phase.1.males.PI.melt.3$Sample.name.2=paste(data.phase.1.males.PI.melt.3$Sample.name,data.phase.1.males.PI.melt.3$variable,sep = '.')
data.phase.1.males.PI.melt.3 = data.phase.1.males.PI.melt.3[order(data.phase.1.males.PI.melt.3[,'variable'],-data.phase.1.males.PI.melt.3[,'Percent.PI.filter']),]
data.phase.1.males.PI.melt.3 = data.phase.1.males.PI.melt.3[!duplicated(data.phase.1.males.PI.melt.3$variable),]

#give pat success
data.phase.1.males.PI.melt.3$Pat.identity=1
data.phase.1.males.PI.melt.3$Pat.success=0
data.phase.1.males.PI.melt.3$Pat.success=ifelse(data.phase.1.males.PI.melt.3$Percent.PI>=0.90&data.phase.1.males.PI.melt.3$value>=26,2,data.phase.1.males.PI.melt.3$Pat.success)
data.phase.1.males.PI.melt.3$Pat.success=ifelse(data.phase.1.males.PI.melt.3$Percent.PI>=0.90&data.phase.1.males.PI.melt.3$value<=26,1,data.phase.1.males.PI.melt.3$Pat.success)
#rename pat success
data.phase.1.males.PI.melt.3$Pat.success.2=ifelse(data.phase.1.males.PI.melt.3$Pat.success==2,'Correct father above threshold',ifelse(data.phase.1.males.PI.melt.3$Pat.success==1,'Correct father below threshold','Incorrect father below threshold'))

ggplot(data.phase.1.males.PI.melt.3)+geom_histogram(color='black',stat='identity',aes(x=Pat.success.2, y=Pat.identity,fill=Pat.success.2))+ theme(axis.text.x = element_blank())+ylab('Number of Broods')+xlab('Paternity testing success')+labs(fill='Known triad paternity')+ theme(legend.position = c(0.4, 0.7))
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_adults_no_log/Paternity testing success.jpeg')



###############
##calculate genotype frequency beforehand! 
phase.1.AA_count=mclapply(data.phase.1.sample.trans$Marker, function(x) length(which(data.phase.1.sample.trans[data.phase.1.sample.trans$Marker==x,2:ncol(data.phase.1.sample.trans)]=='AA')),mc.cores = 16)
phase.1.BB_count=mclapply(data.phase.1.sample.trans$Marker, function(x) length(which(data.phase.1.sample.trans[data.phase.1.sample.trans$Marker==x,2:ncol(data.phase.1.sample.trans)]=='BB')),mc.cores = 16)
phase.1.AB_count=mclapply(data.phase.1.sample.trans$Marker, function(x) length(which(data.phase.1.sample.trans[data.phase.1.sample.trans$Marker==x,2:ncol(data.phase.1.sample.trans)]=='AB')),mc.cores = 16)
#create dataframe
phase.1.AA_count.df=data.frame(matrix(unlist(phase.1.AA_count), nrow=nrow(data.phase.1.sample.trans)))
phase.1.BB_count.df=data.frame(matrix(unlist(phase.1.BB_count), nrow=nrow(data.phase.1.sample.trans)))
phase.1.AB_count.df=data.frame(matrix(unlist(phase.1.AB_count), nrow=nrow(data.phase.1.sample.trans)))
phase.1.genotypes=cbind(phase.1.AA_count.df,phase.1.BB_count.df,phase.1.AB_count.df)
colnames(phase.1.genotypes)=c('AA_count','BB_count','AB_count')
phase.1.genotypes$Total_count=phase.1.genotypes$AA_count+phase.1.genotypes$BB_count+phase.1.genotypes$AB_count
phase.1.genotypes$Afreq=(phase.1.genotypes$AA_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes$Bfreq=(phase.1.genotypes$BB_count+0.5*phase.1.genotypes$AB_count)/phase.1.genotypes$Total_count
phase.1.genotypes.2=data.frame(phase.1.genotypes$Afreq,phase.1.genotypes$Bfreq)
colnames(phase.1.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.1.sample.trans.2=cbind(phase.1.genotypes.2,data.phase.1.sample.trans)

#remove rows with 100% Major allele freq
data.phase.1.sample.trans.2=subset(data.phase.1.sample.trans.2, Afreq!=1 & Bfreq!=1)

#run paternity test and graph
#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.phase.1.males.2=mapply(function(x, y) data.phase.1.males[[x]]=mclapply(data.phase.1.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex(dataframe =data.phase.1.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                            data.phase.1.brood$Sample.name,  #brood names as x
                            data.phase.1.brood$Mother) #mother names as y


#convert into matrix
data.phase.1.males.3=data.frame(matrix(unlist(data.phase.1.males.2), nrow=nrow(data.phase.1.males.2)))
#rename second B5A.BLUE.7.31 sample to B5A.BLUE.7.31_2
colnames(data.phase.1.males.2)[21]='B5A.BLUE.7.31_2'
colnames(data.phase.1.males.3)=colnames(data.phase.1.males.2)

#combine paternity index with male data
data.phase.1.males.PI=cbind(data.phase.1.males,data.phase.1.males.3)

#Tank phase.1 brood names
phase.1.brood.names=names(data.phase.1.males.PI[6:ncol(data.phase.1.males.PI)])

#Name father
phase.1.father=lapply(phase.1.brood.names, function(x) ifelse(data.phase.1.males.PI$Tank==substr(x,1,3) &data.phase.1.males.PI$Sex=='M', 1,0))
phase.1.father=data.frame(matrix(unlist(phase.1.father), ncol=length(phase.1.father)))
colnames(phase.1.father)=c(paste(phase.1.brood.names,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI=cbind(data.phase.1.males.PI,phase.1.father)

#Graph
mclapply(phase.1.brood.names, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_no_log/",x,".jpeg",sep=""),plot=ggplot(data.phase.1.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI[[x]]),y=data.phase.1.males.PI[[x]],fill=data.phase.1.males.PI[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)


#melt males PI data
data.phase.1.males.PI.melt=cbind(data.phase.1.males,data.phase.1.males.3)
data.phase.1.males.PI.melt=melt(data.phase.1.males.PI.melt,id=c(colnames(data.phase.1.males.PI.melt)[1:5]))
#create father coloumn
phase.1.father.melt=cbind(data.phase.1.males.PI[,1:5],phase.1.father)
phase.1.father.melt=melt(phase.1.father.melt,id=c(colnames(phase.1.father.melt)[1:5]))
phase.1.father.melt=phase.1.father.melt[,7]
data.phase.1.males.PI.melt=cbind(data.phase.1.males.PI.melt,phase.1.father.melt)

length(unique(subset(data.phase.1.males.PI.melt,value>=0&phase.1.father.melt==0)$variable))
length(subset(data.phase.1.males.PI.melt,phase.1.father.melt==1&value>=7)$variable)


#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 1 PI Distribution Across Fathers')
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_no_log/Phase 1 PI Distribution Across Fathers.jpeg')

#Create Pat score based off the PI count for brood
data.phase.1.males.PI.melt$value.exp=exp(data.phase.1.males.PI.melt$value)

data.phase.1.brood.total.PI=aggregate(data.phase.1.males.PI.melt$value.exp,list(variable=data.phase.1.males.PI.melt$variable),sum)

colnames(data.phase.1.brood.total.PI)=c('variable','brood.total')

data.phase.1.males.PI.melt.2=merge(data.phase.1.males.PI.melt,data.phase.1.brood.total.PI)

data.phase.1.males.PI.melt.2$Percent.PI=data.phase.1.males.PI.melt.2$value.exp/data.phase.1.males.PI.melt.2$brood.total

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers')

#Filter pat score by percentage and PI count for brood
data.phase.1.males.PI.melt.2$Percent.PI.filter=ifelse(data.phase.1.males.PI.melt.2$value>=7,data.phase.1.males.PI.melt.2$Percent.PI,NA)

ggplot(data.phase.1.males.PI.melt.2)+geom_histogram(binwidth=0.05,aes(Percent.PI.filter,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI Percentage')+ggtitle('Phase 1 CPI Percent Distribution Across Fathers filtered')


#Graphing success of paternity metric
ggplot(subset(data.phase.1.males.PI.melt.2))+geom_bar(color='black',stat='identity',aes(x=variable, y=Percent.PI.filter,fill=phase.1.father.melt))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('Phase_1_Pat_Brood/filter_0.6_filter_genotypeFreq_no_log/Paternity metric success.jpeg')

#############
#run paternity test and graph
#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.phase.1.males.2=mapply(function(x, y) data.phase.1.males[[x]]=mclapply(data.phase.1.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.1.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.1.sample.trans)),paternityIndex(dataframe =data.phase.1.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.1.sample.trans))$PI>0)$PI)),mc.cores = 16),
                            data.phase.1.brood$Sample.name,  #brood names as x
                            data.phase.1.brood$Mother) #mother names as y

#convert into matrix
data.phase.1.males.3=data.frame(matrix(unlist(data.phase.1.males.2), nrow=nrow(data.phase.1.males.2)))
#rename second B5A.BLUE.7.31 sample to B5A.BLUE.7.31_2
colnames(data.phase.1.males.2)[21]='B5A.BLUE.7.31_2'
colnames(data.phase.1.males.3)=colnames(data.phase.1.males.2)

#combine paternity index with male data
data.phase.1.males.PI=cbind(data.phase.1.males,data.phase.1.males.3)

#Tank phase.1 brood names
phase.1.brood.names=names(data.phase.1.males.PI[6:ncol(data.phase.1.males.PI)])

#Name father
phase.1.father=lapply(phase.1.brood.names, function(x) ifelse(data.phase.1.males.PI$Tank==substr(x,1,3) &data.phase.1.males.PI$Sex=='M', 1,0))
phase.1.father=data.frame(matrix(unlist(phase.1.father), ncol=length(phase.1.father)))
colnames(phase.1.father)=c(paste(phase.1.brood.names,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI=cbind(data.phase.1.males.PI,phase.1.father)

#Graph
mclapply(phase.1.brood.names, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/filter_0.6/",x,".jpeg",sep=""),plot=ggplot(data.phase.1.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI[[x]]),y=data.phase.1.males.PI[[x]],fill=data.phase.1.males.PI[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)


#melt males PI data
data.phase.1.males.PI.melt=cbind(data.phase.1.males,data.phase.1.males.3)
data.phase.1.males.PI.melt=melt(data.phase.1.males.PI.melt,id=c(colnames(data.phase.1.males.PI.melt)[1:6]))
#create father coloumn
phase.1.father.melt=cbind(data.phase.1.males.PI[,1:6],phase.1.father)
phase.1.father.melt=melt(phase.1.father.melt,id=c(colnames(phase.1.father.melt)[1:6]))
phase.1.father.melt=phase.1.father.melt[,8]
data.phase.1.males.PI.melt=cbind(data.phase.1.males.PI.melt,phase.1.father.melt)

#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 1 PI Distribution Across Fathers')
ggsave('Phase_1_Pat_Brood/filter_0.6/Phase 1 PI Distribution Across Fathers.jpeg')

#############
#run paternity test with only maternal snps and graph
#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.phase.1.males.2=mapply(function(x, y) data.phase.1.males[[x]]=mclapply(data.phase.1.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =subset(data.phase.1.sample.trans,y!=""), marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.1.sample.trans)),paternityIndex(dataframe =subset(data.phase.1.sample.trans,y!=""), marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.1.sample.trans))$PI>0)$PI)),mc.cores = 16),
                            data.phase.1.brood$Sample.name,  #brood names as x
                            data.phase.1.brood$Mother) #mother names as y


#convert into matrix
data.phase.1.males.3=data.frame(matrix(unlist(data.phase.1.males.2), nrow=nrow(data.phase.1.males.2)))
#rename second B5A.BLUE.7.31 sample to B5A.BLUE.7.31_2
colnames(data.phase.1.males.2)[21]='B5A.BLUE.7.31_2'
colnames(data.phase.1.males.3)=colnames(data.phase.1.males.2)

#combine paternity index with male data
data.phase.1.males.PI=cbind(data.phase.1.males,data.phase.1.males.3)

#Tank phase.1 brood names
phase.1.brood.names=names(data.phase.1.males.PI[6:ncol(data.phase.1.males.PI)])

#Name father
phase.1.father=lapply(phase.1.brood.names, function(x) ifelse(data.phase.1.males.PI$Tank==substr(x,1,3) &data.phase.1.males.PI$Sex=='M', 1,0))
phase.1.father=data.frame(matrix(unlist(phase.1.father), ncol=length(phase.1.father)))
colnames(phase.1.father)=c(paste(phase.1.brood.names,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI=cbind(data.phase.1.males.PI,phase.1.father)

#Graph
mclapply(phase.1.brood.names, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/MaternalSnps_filter_0.6/",x,".jpeg",sep=""),plot=ggplot(data.phase.1.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI[[x]]),y=data.phase.1.males.PI[[x]],fill=data.phase.1.males.PI[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)


#melt males PI data
data.phase.1.males.PI.melt=cbind(data.phase.1.males,data.phase.1.males.3)
data.phase.1.males.PI.melt=melt(data.phase.1.males.PI.melt,id=c(colnames(data.phase.1.males.PI.melt)[1:5]))
#create father coloumn
phase.1.father.melt=cbind(data.phase.1.males.PI[,1:5],phase.1.father)
phase.1.father.melt=melt(phase.1.father.melt,id=c(colnames(phase.1.father.melt)[1:5]))
phase.1.father.melt=phase.1.father.melt[,7]
data.phase.1.males.PI.melt=cbind(data.phase.1.males.PI.melt,phase.1.father.melt)

#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 1 PI Distribution Across Fathers')
ggsave('Phase_1_Pat_Brood/MaternalSnps_filter_0.6/Phase 1 PI Distribution Across Fathers.jpeg')
###############
#test if reduce genotype to just males in tank
#Loop through all broods and males in a tank to produce table of PI for each brood for each male
#use genofreq of just males, brood,mother
data.phase.1.males.2.test=mapply(function(x, y) data.phase.1.males[[x]]=mclapply(data.phase.1.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.1.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = grep(paste(c(x,y,data.phase.1.males$Sample.name),collapse = '|'),colnames(data.phase.1.sample.trans))),paternityIndex(dataframe =data.phase.1.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = grep(paste(c(x,y,data.phase.1.males$Sample.name),collapse = '|'),colnames(data.phase.1.sample.trans)))$PI>0)$PI)),mc.cores = 16),
                                 data.phase.1.brood$Sample.name,  #brood names as x
                                 data.phase.1.brood$Mother) #mother names as y
#convert into matrix
data.phase.1.males.3.test=data.frame(matrix(unlist(data.phase.1.males.2.test), nrow=nrow(data.phase.1.males.2.test)))
#rename second B5A.BLUE.7.31 sample to B5A.BLUE.7.31_2
colnames(data.phase.1.males.2.test)[21]='B5A.BLUE.7.31_2'
colnames(data.phase.1.males.3.test)=colnames(data.phase.1.males.2.test)

#combine paternity index with male data
data.phase.1.males.PI.test=cbind(data.phase.1.males,data.phase.1.males.3.test)

#Tank phase.1 brood names
phase.1.brood.names.test=names(data.phase.1.males.PI.test[7:ncol(data.phase.1.males.PI.test)])

#Name father
phase.1.father.test=lapply(phase.1.brood.names.test, function(x) ifelse(data.phase.1.males.PI.test$Tank==substr(x,1,3) &data.phase.1.males.PI.test$Sex=='M', 1,0))
phase.1.father.test=data.frame(matrix(unlist(phase.1.father.test), ncol=length(phase.1.father.test)))
colnames(phase.1.father.test)=c(paste(phase.1.brood.names.test,'Father',sep='.'))

#combine paternity index with father
data.phase.1.males.PI.test=cbind(data.phase.1.males.PI.test,phase.1.father.test)

#Graph
mclapply(phase.1.brood.names.test, function(x) ggsave(filename=paste("Phase_1_Pat_Brood/Lapply/Test/",x,"Test.jpeg",sep=""),plot=ggplot(data.phase.1.males.PI.test)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.phase.1.males.PI.test[[x]]),y=data.phase.1.males.PI.test[[x]],fill=data.phase.1.males.PI.test[paste(x,'Father',sep='.')]==1))+ggtitle(paste("Brood",x,"Paternity Results Test",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores=16)

#melt males PI data
data.phase.1.males.PI.melt.test=cbind(data.phase.1.males,data.phase.1.males.3.test)
data.phase.1.males.PI.melt.test=melt(data.phase.1.males.PI.melt.test,id=c(colnames(data.phase.1.males.PI.melt.test)[1:6]))
#create father coloumn
phase.1.father.melt.test=cbind(data.phase.1.males.PI.test[,1:6],phase.1.father)
phase.1.father.melt.test=melt(phase.1.father.melt.test,id=c(colnames(phase.1.father.melt.test)[1:6]))
phase.1.father.melt.test=phase.1.father.melt.test[,8]
data.phase.1.males.PI.melt.test=cbind(data.phase.1.males.PI.melt.test,phase.1.father.melt.test)

#graph histogram with fathers
ggplot(data.phase.1.males.PI.melt.test)+geom_histogram(binwidth=1,aes(value,fill=phase.1.father.melt==1))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index Test')+ggtitle('Phase 1 PI Distribution Across Fathers Test')
ggsave('Phase_1_Pat_Brood/Lapply/Test/Phase 1 PI Distribution Across Fathers Test.jpeg')









############
#Replace code for single Brood
#replace B2.ORANGE.8.24 with brood name B2.ORANGE.8.24
#replace ORANGE.B2 with mother name
#replace B2 with tank name


#subset into B2 and remove other broods
data.phase.2.B2=subset(data.phase.2,Tank=='B2')
data.phase.2.B2=subset(data.phase.2.B2,Sample.name=='B2.ORANGE.8.24'|Sex!='B')
#transpose so col for markers and then each sample
data.phase.2.B2.sample=data.phase.2.B2[,-c(2:6)]
rownames(data.phase.2.B2.sample)=data.phase.2.B2.sample$Sample.name
data.phase.2.B2.sample=data.phase.2.B2.sample[,-c(1)]
data.phase.2.B2.sample.trans=t(data.phase.2.B2.sample)
data.phase.2.B2.sample.trans=as.data.frame(data.phase.2.B2.sample.trans)
data.phase.2.B2.sample.trans=rownames_to_column(data.phase.2.B2.sample.trans,var="Marker")

#Tank dataframe
Males=colnames(data.phase.2.B2.sample.trans)[2:ncol(data.phase.2.B2.sample.trans)]

#Change to brood name
PI.B2.ORANGE.8.24=NA
data.B2.pat=data.frame(Sample.name=Males,PI.B2.ORANGE.8.24=PI.B2.ORANGE.8.24)
data.B2.pat=merge(data.B2.pat,data.pat.names,by='Sample.name')
data.B2.pat$Relation=ifelse(data.B2.pat$Sample.name=='ORANGE.B2','Mother',ifelse(
  data.B2.pat$Sex=='B','Brood',ifelse(
    data.B2.pat$Tank=='B2'&data.B2.pat$Sex=='M','Male','Female')))

##Paternity Test B2.ORANGE.8.24
data.B2.pat$PI.B2.ORANGE.8.24=mclapply(colnames(data.phase.2.B2.sample.trans)[2:ncol(data.phase.2.B2.sample.trans)],function(x) log(prod(subset(paternityIndex(dataframe =data.phase.2.B2.sample.trans, marker = "Marker", mother ="ORANGE.B2", child = "B2.ORANGE.8.24", AF = x, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.B2.sample.trans)),paternityIndex(dataframe =data.phase.2.B2.sample.trans, marker = "Marker", mother ="ORANGE.B2", child = "B2.ORANGE.8.24", AF = x, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.B2.sample.trans))$PI>0)$PI)),mc.cores = 16)

#Graph
data.B2.pat$PI.B2.ORANGE.8.24=as.numeric(data.B2.pat$PI.B2.ORANGE.8.24)
ggplot(data.B2.pat)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-PI.B2.ORANGE.8.24),y=PI.B2.ORANGE.8.24,fill=Relation))+ggtitle('Brood B2.ORANGE.8.24 Paternity Results')+xlab('Sample')+ylab('Log composite paternity index')+theme(axis.text.x = element_text(angle=90,vjust=0.5))
ggsave('Brood B2.ORANGE.8.24 Paternity Results.png')