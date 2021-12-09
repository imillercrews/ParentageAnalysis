#Tank Examples Phase 2 to go with 'Pat geno likelihoods filtered.R'
#Set working directory
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity Testing/2bRAD_data/Data")
#load library
library(tidyverse)
library(paternity)
library(parallel)

load('data.geno.filter.trans.ind.merge.RData')

#load name data
data.pat.names=read.csv('bams.ind.csv')
data.pat.names.phase.1=data.pat.names
data.pat.names.phase.2=data.pat.names
data.pat.names=data.pat.names[,-c(1:3,6:12,15:18)]
data.pat.names[]=lapply(data.pat.names, function(x) as.character(x))
data.pat.names.phase.1=subset(data.pat.names.phase.1,Phase=='1')
data.pat.names.phase.1=data.pat.names.phase.1[,-c(1:3,7:12,15:18)]
data.pat.names.phase.1[]=lapply(data.pat.names.phase.1, function(x) as.character(x))
data.pat.names.phase.2=subset(data.pat.names.phase.2,Phase=='2')
data.pat.names.phase.2=data.pat.names.phase.2[,-c(1:3,6:12,15:18)]
data.pat.names.phase.2[]=lapply(data.pat.names.phase.2, function(x) as.character(x))


data.phase.2=subset(data.geno.filter.trans.ind.merge,Phase==2)
data.phase.2=droplevels(data.phase.2)
levels(data.phase.2$Tank)

##########
#check number of calls for each sample
#transpose so col for markers and then each sample
data.phase.2.sample=data.phase.2[,-c(2:6)]
rownames(data.phase.2.sample)=data.phase.2.sample$Sample.name
data.phase.2.sample=data.phase.2.sample[,-c(1)]
data.phase.2.sample.trans=t(data.phase.2.sample)
data.phase.2.sample.trans=as.data.frame(data.phase.2.sample.trans)
data.phase.2.sample.trans=rownames_to_column(data.phase.2.sample.trans,var="Marker")

data.phase.2.sample.counts=mclapply(colnames(data.phase.2.sample.trans), function(x) sum(data.phase.2.sample.trans[[x]]!=""),mc.cores=16)
data.phase.2.sample.counts.2=data.frame(matrix(unlist(data.phase.2.sample.counts),ncol=1))
colnames(data.phase.2.sample.counts.2)=c('Counts')
data.phase.2.sample.counts.2$Sample.name=colnames(data.phase.2.sample.trans)

data.phase.2.sample.counts.3=merge(data.pat.names.phase.2,data.phase.2.sample.counts.2,by='Sample.name')
data.phase.2.sample.counts.3$Percent=data.phase.2.sample.counts.3$Counts/data.phase.2.sample.counts.2$Counts[data.phase.2.sample.counts.2$Sample.name=='Marker']

ggplot(data.phase.2.sample.counts.3)+geom_histogram(aes(Percent,fill=Sex))
#check number of calls for each sample after genotype filter
data.phase.2.sample.counts=mclapply(colnames(data.phase.2.sample.trans.2), function(x) sum(data.phase.2.sample.trans.2[[x]]!=""),mc.cores=16)

data.phase.2.sample.counts.2=data.frame(matrix(unlist(data.phase.2.sample.counts),ncol=1))
colnames(data.phase.2.sample.counts.2)=c('Counts')
data.phase.2.sample.counts.2$Sample.name=colnames(data.phase.2.sample.trans.2)

data.phase.2.sample.counts.3=merge(data.pat.names.phase.2,data.phase.2.sample.counts.2,by='Sample.name')
data.phase.2.sample.counts.3$Percent=data.phase.2.sample.counts.3$Counts/data.phase.2.sample.counts.2$Counts[data.phase.2.sample.counts.2$Sample.name=='Marker']

ggplot(data.phase.2.sample.counts.3)+geom_histogram(aes(Percent,fill=Sex))
ggsave("Phase_2_Pat_Brood/Filter_0.6_filter_genotype_all_pop/Histogram counts.jpeg")

#Need to remove BLUE.G2, GREEN.G2 [Should keep Green.G2?]

##############
#create correlation matrix for males
data.phase.2.cor.test=data.phase.2.sample.trans.2[,4:ncol(data.phase.2.sample.trans.2)]
data.phase.2.cor.test[]=lapply(data.phase.2.cor.test, function(x) ifelse(x=="",NA,x))
unique(data.pat.names.phase.2.males$Tank)

jpeg('Phase_2_Pat_Brood/Heatmap/Tank G2 Male Heatmap.jpeg')
heatmap(x = cor(data.phase.2.cor.test[,c(subset(data.pat.names.phase.2.males, Tank=='G2')$Sample.name[c(1:2,4:8)])],use = "complete.obs",method = "spearman"),col=colorRampPalette(c("blue", "white", "red"))(20), symm = TRUE)
dev.off()


##############
#Check for missing mothers
#remove D2.LIGHTBLUE.6.23
data.pat.names.2=subset(data.pat.names,Sex!='M')
data.pat.names.2$Mother=ifelse(data.pat.names.2$Sex=='B',paste(data.pat.names.2$Sample,data.pat.names.2$Tank,sep = '.'),NA)

data.pat.names.2$Check=ifelse(data.pat.names.2$Sex=='B',data.pat.names.2$Mother,data.pat.names.2$Sample.name)

nrow(data.pat.names.2)
nrow(data.pat.names.2[data.pat.names.2$Sex=='B',])
nrow(data.pat.names.2[data.pat.names.2$Sex=='F',])

nrow(data.pat.names.2[!duplicated(data.pat.names.2$Check), ])

nrow(data.pat.names.2[duplicated(data.pat.names.2$Check), ])

##########
#Paternity test and graph using ONLY adults for genotype freq
#transpose so col for markers and then each sample
data.phase.2.sample=data.phase.2[,-c(2:6)]
rownames(data.phase.2.sample)=data.phase.2.sample$Sample.name
data.phase.2.sample=data.phase.2.sample[,-c(1)]
data.phase.2.sample.trans=t(data.phase.2.sample)
data.phase.2.sample.trans=as.data.frame(data.phase.2.sample.trans)
data.phase.2.sample.trans=rownames_to_column(data.phase.2.sample.trans,var="Marker")

#subset into males
data.pat.names.phase.2.males=subset(data.pat.names.phase.2,Sex=='M')
#subset broods with mothers
data.pat.names.phase.2$Mother=ifelse(data.pat.names.phase.2$Sex=='B',paste(data.pat.names.phase.2$Sample,data.pat.names.phase.2$Tank,sep = '.'),NA)
data.pat.names.phase.2.brood=subset(data.pat.names.phase.2,Sex=='B')
#remove brood with missing mother D2.LIGHTBLUE.6.23
data.pat.names.phase.2.brood=subset(data.pat.names.phase.2.brood,Sample.name!='D2.LIGHTBLUE.6.23')

##calculate genotype frequency
#Subset into only adults
data.phase.2.sample.trans.adults=data.phase.2.sample.trans[,c(1,5:22,30:33,51:101)]

#create length data
phase.2.AA_count=mclapply(data.phase.2.sample.trans.adults$Marker, function(x) length(which(data.phase.2.sample.trans.adults[data.phase.2.sample.trans.adults$Marker==x,2:ncol(data.phase.2.sample.trans.adults)]=='AA')),mc.cores = 16)
phase.2.BB_count=mclapply(data.phase.2.sample.trans.adults$Marker, function(x) length(which(data.phase.2.sample.trans.adults[data.phase.2.sample.trans.adults$Marker==x,2:ncol(data.phase.2.sample.trans.adults)]=='BB')),mc.cores = 16)
phase.2.AB_count=mclapply(data.phase.2.sample.trans.adults$Marker, function(x) length(which(data.phase.2.sample.trans.adults[data.phase.2.sample.trans.adults$Marker==x,2:ncol(data.phase.2.sample.trans.adults)]=='AB')),mc.cores = 16)
#create dataframe
phase.2.AA_count.df=data.frame(matrix(unlist(phase.2.AA_count), nrow=nrow(data.phase.2.sample.trans.adults)))
phase.2.BB_count.df=data.frame(matrix(unlist(phase.2.BB_count), nrow=nrow(data.phase.2.sample.trans.adults)))
phase.2.AB_count.df=data.frame(matrix(unlist(phase.2.AB_count), nrow=nrow(data.phase.2.sample.trans.adults)))
phase.2.genotypes=cbind(phase.2.AA_count.df,phase.2.BB_count.df,phase.2.AB_count.df)
colnames(phase.2.genotypes)=c('AA_count','BB_count','AB_count')
phase.2.genotypes$Total_count=phase.2.genotypes$AA_count+phase.2.genotypes$BB_count+phase.2.genotypes$AB_count
phase.2.genotypes$Afreq=(phase.2.genotypes$AA_count+0.5*phase.2.genotypes$AB_count)/phase.2.genotypes$Total_count
phase.2.genotypes$Bfreq=(phase.2.genotypes$BB_count+0.5*phase.2.genotypes$AB_count)/phase.2.genotypes$Total_count
phase.2.genotypes.2=data.frame(phase.2.genotypes$Afreq,phase.2.genotypes$Bfreq)
colnames(phase.2.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.2.sample.trans.2=cbind(phase.2.genotypes.2,data.phase.2.sample.trans)

#remove rows with 100% Major allele freq
data.phase.2.sample.trans.2=subset(data.phase.2.sample.trans.2, Afreq!=1 & Bfreq!=1)

#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.pat.names.phase.2.males.2=mapply(function(x, y) data.pat.names.phase.2.males[[x]]=mclapply(data.pat.names.phase.2.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.2.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex(dataframe =data.phase.2.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                                      data.pat.names.phase.2.brood$Sample.name,  #brood names as x
                                      data.pat.names.phase.2.brood$Mother) #mother names as y


#convert into matrix
data.pat.names.phase.2.males.3=data.frame(matrix(unlist(data.pat.names.phase.2.males.2), nrow=nrow(data.pat.names.phase.2.males.2)))
colnames(data.pat.names.phase.2.males.3)=colnames(data.pat.names.phase.2.males.2)
#combine paternity index with male data
data.pat.names.phase.2.males.PI=cbind(data.pat.names.phase.2.males,data.pat.names.phase.2.males.3)
#Tank Phase 2 brood names
Phase.2.brood.names=names(data.pat.names.phase.2.males.PI[6:ncol(data.pat.names.phase.2.males.PI)])

#Remove BLUE.G2
data.pat.names.phase.2.males.PI=subset(data.pat.names.phase.2.males.PI,Sample.name!='BLUE.G2')

#Graph
mclapply(Phase.2.brood.names, function(x) ggsave(filename=paste("Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_G2removed/",x,".jpeg",sep=""),plot=ggplot(subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2)))+geom_bar(stat='identity',aes(x=reorder(Sample.name,-subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]]),y=subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]],fill=subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]]>=0))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores = 16)

#melt males PI data
data.phase.2.males.PI.melt=cbind(data.pat.names.phase.2.males,data.pat.names.phase.2.males.3)
data.phase.2.males.PI.melt=melt(data.phase.2.males.PI.melt,id=c(colnames(data.phase.2.males.PI.melt)[1:5]))

#graph histogram with fathers
ggplot(data.phase.2.males.PI.melt)+geom_histogram(binwidth=1,aes(value))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 2 PI Distribution Across Fathers')
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_G2removed/Phase 2 PI Distribution Across Fathers.jpeg')


#Create Pat score based off the PI count for brood
data.phase.2.males.PI.melt$value.exp=exp(data.phase.2.males.PI.melt$value)

#For unknown tanks will want to keep tank as part of list
data.phase.2.brood.total.PI=aggregate(data.phase.2.males.PI.melt$value.exp,list(Tank=data.phase.2.males.PI.melt$Tank,variable=data.phase.2.males.PI.melt$variable),sum)

data.phase.2.brood.total.PI$Match=ifelse(data.phase.2.brood.total.PI$Tank==substr(data.phase.2.brood.total.PI$variable,1,2),1,0)

data.phase.2.brood.total.PI=subset(data.phase.2.brood.total.PI,Match==1)

data.phase.2.brood.total.PI=data.phase.2.brood.total.PI[,1:3]

colnames(data.phase.2.brood.total.PI)=c('Tank','variable','brood.total')

data.phase.2.males.PI.melt.2=merge(data.phase.2.males.PI.melt,data.phase.2.brood.total.PI,all=T)

data.phase.2.males.PI.melt.2$Percent.PI=data.phase.2.males.PI.melt.2$value.exp/data.phase.2.males.PI.melt.2$brood.total

data.phase.2.males.PI.melt.3=subset(data.phase.2.males.PI.melt.2,Percent.PI>=0)

ggplot(data.phase.2.males.PI.melt.3)+geom_histogram(aes(Percent.PI,fill=value>=7))


ggplot(subset(data.phase.2.males.PI.melt.3))+geom_histogram(aes(Percent.PI,fill=value>=7))


#Select the top likely father candidate
data.phase.2.males.PI.melt.3$Percent.PI.filter=ifelse(data.phase.2.males.PI.melt.3$value>=27,data.phase.2.males.PI.melt.3$Percent.PI,NA)

Fish.Cols<-c("BLACK"="black","BLUE"="blue","BROWN"="brown","GREEN"="green","PURPLE"="purple","RED"="red","YELLOW"="yellow","WHITE"="grey","YELLOW_A"="yellow","YELLOW_B"="yellow")


ggplot(subset(data.phase.2.males.PI.melt.3,Tank=='D2'))+geom_bar(color='black',stat='identity',aes(x=variable, y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=Fish.Cols,aesthetics = c("colour", "fill"))+xlab('Brood name')


data.phase.2.males.PI.melt.3$Tank.2=as.factor(data.phase.2.males.PI.melt.3$Tank)
data.phase.2.males.PI.melt.3$Tank.2=as.numeric(data.phase.2.males.PI.melt.3$Tank.2)

ggplot(subset(data.phase.2.males.PI.melt.3))+geom_bar(color='black',stat='identity',aes(x=reorder(variable,Tank.2), y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=Fish.Cols,aesthetics = c("colour", "fill"))+xlab('Brood name')+geom_vline(xintercept = c(4.5,11.5,13.5,18.5,20.5),size=2)
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_G2removed/Percent.PI.jpg',width = 10, height = 10)

#create paternity success dataframe
data.phase.2.paternity=data.phase.2.males.PI.melt.3
data.phase.2.paternity=subset(data.phase.2.paternity,Percent.PI>=0.01)
data.phase.2.paternity$Pat.success='Unknown father'
data.phase.2.paternity$Pat.success=ifelse(data.phase.2.paternity$value>=27&data.phase.2.paternity$Percent.PI>=0,'Partial father',data.phase.2.paternity$Pat.success)
data.phase.2.paternity$Pat.success=ifelse(data.phase.2.paternity$value>=27&data.phase.2.paternity$Percent.PI>=0.99,'Full father',data.phase.2.paternity$Pat.success)
data.phase.2.paternity$Pat.identity=1
data.phase.2.paternity$Pat.identity=ifelse(data.phase.2.paternity$Pat.success=='Unknown father'&data.phase.2.paternity$Percent.PI<=0.9,2,data.phase.2.paternity$Pat.identity)
data.phase.2.paternity=subset(data.phase.2.paternity,Pat.identity==1)


ggplot(data.phase.2.paternity)+geom_tile(color='black',stat='identity',aes(x=reorder(variable,Tank.2), y=Pat.identity,fill=Pat.success))+ theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank())+xlab('Brood by tank')+ylab('')+geom_vline(xintercept = c(4.5,11.5,13.5,18.5,20.5),size=2)
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_G2removed/Paternity success.jpg',width = 10, height = 10)

data.phase.2.paternity=data.phase.2.paternity[,c(1:3,5,11,13)]
write.csv(data.phase.2.paternity,file='data.phase.2.paternity.csv')

##########
#Paternity test and graph using all samples for genotype freq
#transpose so col for markers and then each sample
data.phase.2.sample=data.phase.2[,-c(2:6)]
rownames(data.phase.2.sample)=data.phase.2.sample$Sample.name
data.phase.2.sample=data.phase.2.sample[,-c(1)]
data.phase.2.sample.trans=t(data.phase.2.sample)
data.phase.2.sample.trans=as.data.frame(data.phase.2.sample.trans)
data.phase.2.sample.trans=rownames_to_column(data.phase.2.sample.trans,var="Marker")

#subset into males
data.pat.names.phase.2.males=subset(data.pat.names.phase.2,Sex=='M')
#subset broods with mothers
data.pat.names.phase.2$Mother=ifelse(data.pat.names.phase.2$Sex=='B',paste(data.pat.names.phase.2$Sample,data.pat.names.phase.2$Tank,sep = '.'),NA)
data.pat.names.phase.2.brood=subset(data.pat.names.phase.2,Sex=='B')
#remove brood with missing mother D2.LIGHTBLUE.6.23
data.pat.names.phase.2.brood=subset(data.pat.names.phase.2.brood,Sample.name!='D2.LIGHTBLUE.6.23')

##calculate genotype frequency
phase.2.AA_count=mclapply(data.phase.2.sample.trans$Marker, function(x) length(which(data.phase.2.sample.trans[data.phase.2.sample.trans$Marker==x,2:ncol(data.phase.2.sample.trans)]=='AA')),mc.cores = 16)
phase.2.BB_count=mclapply(data.phase.2.sample.trans$Marker, function(x) length(which(data.phase.2.sample.trans[data.phase.2.sample.trans$Marker==x,2:ncol(data.phase.2.sample.trans)]=='BB')),mc.cores = 16)
phase.2.AB_count=mclapply(data.phase.2.sample.trans$Marker, function(x) length(which(data.phase.2.sample.trans[data.phase.2.sample.trans$Marker==x,2:ncol(data.phase.2.sample.trans)]=='AB')),mc.cores = 16)
#create dataframe
phase.2.AA_count.df=data.frame(matrix(unlist(phase.2.AA_count), nrow=nrow(data.phase.2.sample.trans)))
phase.2.BB_count.df=data.frame(matrix(unlist(phase.2.BB_count), nrow=nrow(data.phase.2.sample.trans)))
phase.2.AB_count.df=data.frame(matrix(unlist(phase.2.AB_count), nrow=nrow(data.phase.2.sample.trans)))
phase.2.genotypes=cbind(phase.2.AA_count.df,phase.2.BB_count.df,phase.2.AB_count.df)
colnames(phase.2.genotypes)=c('AA_count','BB_count','AB_count')
phase.2.genotypes$Total_count=phase.2.genotypes$AA_count+phase.2.genotypes$BB_count+phase.2.genotypes$AB_count
phase.2.genotypes$Afreq=(phase.2.genotypes$AA_count+0.5*phase.2.genotypes$AB_count)/phase.2.genotypes$Total_count
phase.2.genotypes$Bfreq=(phase.2.genotypes$BB_count+0.5*phase.2.genotypes$AB_count)/phase.2.genotypes$Total_count
phase.2.genotypes.2=data.frame(phase.2.genotypes$Afreq,phase.2.genotypes$Bfreq)
colnames(phase.2.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.2.sample.trans.2=cbind(phase.2.genotypes.2,data.phase.2.sample.trans)

#remove rows with 100% Major allele freq
data.phase.2.sample.trans.2=subset(data.phase.2.sample.trans.2, Afreq!=1 & Bfreq!=1)

#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.pat.names.phase.2.males.2=mapply(function(x, y) data.pat.names.phase.2.males[[x]]=mclapply(data.pat.names.phase.2.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.2.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex(dataframe =data.phase.2.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                                    data.pat.names.phase.2.brood$Sample.name,  #brood names as x
                                    data.pat.names.phase.2.brood$Mother) #mother names as y


#convert into matrix
data.pat.names.phase.2.males.3=data.frame(matrix(unlist(data.pat.names.phase.2.males.2), nrow=nrow(data.pat.names.phase.2.males.2)))
colnames(data.pat.names.phase.2.males.3)=colnames(data.pat.names.phase.2.males.2)
#combine paternity index with male data
data.pat.names.phase.2.males.PI=cbind(data.pat.names.phase.2.males,data.pat.names.phase.2.males.3)
#Tank Phase 2 brood names
Phase.2.brood.names=names(data.pat.names.phase.2.males.PI[6:ncol(data.pat.names.phase.2.males.PI)])

#Remove BLUE.G2
data.pat.names.phase.2.males.PI=subset(data.pat.names.phase.2.males.PI,Sample.name!='BLUE.G2')

#Graph
mclapply(Phase.2.brood.names, function(x) ggsave(filename=paste("Phase_2_Pat_Brood/Filter_0.6_filter_genotype_all_pop_G2removed/",x,".jpeg",sep=""),plot=ggplot(subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2)))+geom_bar(stat='identity',aes(x=reorder(Sample.name,-subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]]),y=subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]],fill=subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]]>=0))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores = 16)

#melt males PI data
data.phase.2.males.PI.melt=cbind(data.pat.names.phase.2.males,data.pat.names.phase.2.males.3)
data.phase.2.males.PI.melt=melt(data.phase.2.males.PI.melt,id=c(colnames(data.phase.2.males.PI.melt)[1:5]))

#graph histogram with fathers
ggplot(data.phase.2.males.PI.melt)+geom_histogram(binwidth=1,aes(value))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 2 PI Distribution Across Fathers')
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_all_pop_G2removed/Phase 2 PI Distribution Across Fathers.jpeg')


#Create Pat score based off the PI count for brood
data.phase.2.males.PI.melt$value.exp=exp(data.phase.2.males.PI.melt$value)

#For unknown tanks will want to keep tank as part of list
data.phase.2.brood.total.PI=aggregate(data.phase.2.males.PI.melt$value.exp,list(Tank=data.phase.2.males.PI.melt$Tank,variable=data.phase.2.males.PI.melt$variable),sum)

data.phase.2.brood.total.PI$Match=ifelse(data.phase.2.brood.total.PI$Tank==substr(data.phase.2.brood.total.PI$variable,1,2),1,0)

data.phase.2.brood.total.PI=subset(data.phase.2.brood.total.PI,Match==1)

data.phase.2.brood.total.PI=data.phase.2.brood.total.PI[,1:3]

colnames(data.phase.2.brood.total.PI)=c('Tank','variable','brood.total')

data.phase.2.males.PI.melt.2=merge(data.phase.2.males.PI.melt,data.phase.2.brood.total.PI,all=T)

data.phase.2.males.PI.melt.2$Percent.PI=data.phase.2.males.PI.melt.2$value.exp/data.phase.2.males.PI.melt.2$brood.total

data.phase.2.males.PI.melt.3=subset(data.phase.2.males.PI.melt.2,Percent.PI>=0)

ggplot(data.phase.2.males.PI.melt.3)+geom_histogram(aes(Percent.PI,fill=value>=7))


ggplot(subset(data.phase.2.males.PI.melt.3))+geom_histogram(aes(Percent.PI,fill=value>=7))


#Select the top likely father candidate
data.phase.2.males.PI.melt.3$Percent.PI.filter=ifelse(data.phase.2.males.PI.melt.3$value>=7,data.phase.2.males.PI.melt.3$Percent.PI,NA)

Fish.Cols<-c("BLACK"="black","BLUE"="blue","BROWN"="brown","GREEN"="green","PURPLE"="purple","RED"="red","YELLOW"="yellow","WHITE"="grey")


ggplot(subset(data.phase.2.males.PI.melt.3,Tank=='G2'))+geom_bar(color='black',stat='identity',aes(x=variable, y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=Fish.Cols,aesthetics = c("colour", "fill"))+xlab('Brood name')


data.phase.2.males.PI.melt.3$Tank.2=as.factor(data.phase.2.males.PI.melt.3$Tank)
data.phase.2.males.PI.melt.3$Tank.2=as.numeric(data.phase.2.males.PI.melt.3$Tank.2)

ggplot(subset(data.phase.2.males.PI.melt.3))+geom_bar(color='black',stat='identity',aes(x=reorder(variable,Tank.2), y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=Fish.Cols,aesthetics = c("colour", "fill"))+xlab('Brood name')+geom_vline(xintercept = c(4.5,11.5,13.5,18.5,20.5),size=2)
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_all_pop_G2removed/Percent.PI.jpg',width = 10, height = 10)


##########
#Replace code for tank
#replace D2 with tank name


#subset into D2 and remove other broods
data.phase.D2=subset(data.phase.2,Tank=='D2')
#data.phase.2.D2=subset(data.phase.2.D2,Sample.name=='D2.ORANGE.8.24'|Sex!='B')

#transpose so col for markers and then each sample
data.phase.2.D2.sample=data.phase.D2[,-c(2:6)]
rownames(data.phase.2.D2.sample)=data.phase.2.D2.sample$Sample.name
data.phase.2.D2.sample=data.phase.2.D2.sample[,-c(1)]
data.phase.2.D2.sample.trans=t(data.phase.2.D2.sample)
data.phase.2.D2.sample.trans=as.data.frame(data.phase.2.D2.sample.trans)
data.phase.2.D2.sample.trans=rownames_to_column(data.phase.2.D2.sample.trans,var="Marker")

#Tank dataframe
data.D2.pat=subset(data.pat.names,Tank=='D2')
#subset into males
data.D2.pat.males=subset(data.D2.pat,Sex=='M')
#subset broods with mothers
data.D2.pat$Mother=ifelse(data.D2.pat$Sex=='B',paste(data.D2.pat$Sample,data.D2.pat$Tank,sep = '.'),NA)
data.D2.pat.brood=subset(data.D2.pat,Sex=='B')

#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.D2.pat.males.2=mapply(function(x, y) data.D2.pat.males[[x]]=mclapply(data.D2.pat.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.2.D2.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.D2.sample.trans)),paternityIndex(dataframe =data.phase.2.D2.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.D2.sample.trans))$PI>0)$PI)),mc.cores = 16),
       data.D2.pat.brood$Sample.name,  #brood names as x
       data.D2.pat.brood$Mother) #mother names as y
#convert into matrix
data.D2.pat.males.3=data.frame(matrix(unlist(data.D2.pat.males.2), nrow=nrow(data.D2.pat.males.2)))
colnames(data.D2.pat.males.3)=colnames(data.D2.pat.males.2)
#combine paternity index with male data
data.D2.pat.males.PI=cbind(data.D2.pat.males,data.D2.pat.males.3)
#Tank D2 brood names
D2.brood.names=names(data.D2.pat.males.PI[5:ncol(data.D2.pat.males.PI)])

#Graph
lapply(D2.brood.names, function(x) ggsave(filename=paste("Phase_2_Pat_Brood/Tank_D2/",x,".jpeg",sep=""),plot=ggplot(data.D2.pat.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.D2.pat.males.PI[[x]]),y=data.D2.pat.males.PI[[x]],fill=data.D2.pat.males.PI[[x]]>=0))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)))



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

##########
#Paternity test and graph using ONLY adults for genotype freq NEW paternityIndex method
##DOES NOT WORK
#change paternityIndex so every site gets a min score of 0.00000000001 
#change code to:  alldata$PI = with(alldata, COEFA * (1/a) + COEFB * (1/b) + COEFAB * (1/(a + b)) + 0.00000000001)
paternityIndex2=edit(paternityIndex)

#transpose so col for markers and then each sample
data.phase.2.sample=data.phase.2[,-c(2:6)]
rownames(data.phase.2.sample)=data.phase.2.sample$Sample.name
data.phase.2.sample=data.phase.2.sample[,-c(1)]
data.phase.2.sample.trans=t(data.phase.2.sample)
data.phase.2.sample.trans=as.data.frame(data.phase.2.sample.trans)
data.phase.2.sample.trans=rownames_to_column(data.phase.2.sample.trans,var="Marker")

#subset into males
data.pat.names.phase.2.males=subset(data.pat.names.phase.2,Sex=='M')
#subset broods with mothers
data.pat.names.phase.2$Mother=ifelse(data.pat.names.phase.2$Sex=='B',paste(data.pat.names.phase.2$Sample,data.pat.names.phase.2$Tank,sep = '.'),NA)
data.pat.names.phase.2.brood=subset(data.pat.names.phase.2,Sex=='B')
#remove brood with missing mother D2.LIGHTBLUE.6.23
data.pat.names.phase.2.brood=subset(data.pat.names.phase.2.brood,Sample.name!='D2.LIGHTBLUE.6.23')

##calculate genotype frequency
#Subset into only adults
data.phase.2.sample.trans.adults=data.phase.2.sample.trans[,c(1,5:22,30:33,51:101)]

#create length data
phase.2.AA_count=mclapply(data.phase.2.sample.trans.adults$Marker, function(x) length(which(data.phase.2.sample.trans.adults[data.phase.2.sample.trans.adults$Marker==x,2:ncol(data.phase.2.sample.trans.adults)]=='AA')),mc.cores = 16)
phase.2.BB_count=mclapply(data.phase.2.sample.trans.adults$Marker, function(x) length(which(data.phase.2.sample.trans.adults[data.phase.2.sample.trans.adults$Marker==x,2:ncol(data.phase.2.sample.trans.adults)]=='BB')),mc.cores = 16)
phase.2.AB_count=mclapply(data.phase.2.sample.trans.adults$Marker, function(x) length(which(data.phase.2.sample.trans.adults[data.phase.2.sample.trans.adults$Marker==x,2:ncol(data.phase.2.sample.trans.adults)]=='AB')),mc.cores = 16)
#create dataframe
phase.2.AA_count.df=data.frame(matrix(unlist(phase.2.AA_count), nrow=nrow(data.phase.2.sample.trans.adults)))
phase.2.BB_count.df=data.frame(matrix(unlist(phase.2.BB_count), nrow=nrow(data.phase.2.sample.trans.adults)))
phase.2.AB_count.df=data.frame(matrix(unlist(phase.2.AB_count), nrow=nrow(data.phase.2.sample.trans.adults)))
phase.2.genotypes=cbind(phase.2.AA_count.df,phase.2.BB_count.df,phase.2.AB_count.df)
colnames(phase.2.genotypes)=c('AA_count','BB_count','AB_count')
phase.2.genotypes$Total_count=phase.2.genotypes$AA_count+phase.2.genotypes$BB_count+phase.2.genotypes$AB_count
phase.2.genotypes$Afreq=(phase.2.genotypes$AA_count+0.5*phase.2.genotypes$AB_count)/phase.2.genotypes$Total_count
phase.2.genotypes$Bfreq=(phase.2.genotypes$BB_count+0.5*phase.2.genotypes$AB_count)/phase.2.genotypes$Total_count
phase.2.genotypes.2=data.frame(phase.2.genotypes$Afreq,phase.2.genotypes$Bfreq)
colnames(phase.2.genotypes.2)=c('Afreq','Bfreq')

#combine genotype data with sample data
data.phase.2.sample.trans.2=cbind(phase.2.genotypes.2,data.phase.2.sample.trans)

#remove rows with 100% Major allele freq
data.phase.2.sample.trans.2=subset(data.phase.2.sample.trans.2, Afreq!=1 & Bfreq!=1)

#Loop through all broods and males in a tank to produce table of PI for each brood for each male
data.pat.names.phase.2.males.2=mapply(function(x, y) data.pat.names.phase.2.males[[x]]=mclapply(data.pat.names.phase.2.males$Sample.name,function(z) log(prod(subset(paternityIndex2(dataframe =data.phase.2.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA),paternityIndex2(dataframe =data.phase.2.sample.trans.2, marker = "Marker", mother =y, child = x, AF = z, Afreq = "Afreq", Bfreq = "Bfreq", afcalcol = NA)$PI>0)$PI)),mc.cores = 16),
                                      data.pat.names.phase.2.brood$Sample.name,  #brood names as x
                                      data.pat.names.phase.2.brood$Mother) #mother names as y


#convert into matrix
data.pat.names.phase.2.males.3=data.frame(matrix(unlist(data.pat.names.phase.2.males.2), nrow=nrow(data.pat.names.phase.2.males.2)))
colnames(data.pat.names.phase.2.males.3)=colnames(data.pat.names.phase.2.males.2)
#combine paternity index with male data
data.pat.names.phase.2.males.PI=cbind(data.pat.names.phase.2.males,data.pat.names.phase.2.males.3)
#Tank Phase 2 brood names
Phase.2.brood.names=names(data.pat.names.phase.2.males.PI[6:ncol(data.pat.names.phase.2.males.PI)])

#Remove BLUE.G2
data.pat.names.phase.2.males.PI=subset(data.pat.names.phase.2.males.PI,Sample.name!='BLUE.G2')

#Graph (Takes a long time so runn when needed!)
#mclapply(Phase.2.brood.names, function(x) ggsave(filename=paste("Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_G2removed/",x,".jpeg",sep=""),plot=ggplot(subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2)))+geom_bar(stat='identity',aes(x=reorder(Sample.name,-subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]]),y=subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]],fill=subset(data.pat.names.phase.2.males.PI,Tank==substr(x,1,2))[[x]]>=0))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE)),mc.cores = 16)

#melt males PI data
data.phase.2.males.PI.melt=cbind(data.pat.names.phase.2.males,data.pat.names.phase.2.males.3)
data.phase.2.males.PI.melt=melt(data.phase.2.males.PI.melt,id=c(colnames(data.phase.2.males.PI.melt)[1:5]))

#graph histogram with fathers
ggplot(data.phase.2.males.PI.melt)+geom_histogram(binwidth=1,aes(value))+guides(fill=guide_legend(title="Father"))+xlab('Composite PI index')+ggtitle('Phase 2 PI Distribution Across Fathers')
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_no_pat_exclusion/Phase 2 PI Distribution Across Fathers.jpeg')


#Create Pat score based off the PI count for brood
data.phase.2.males.PI.melt$value.exp=exp(data.phase.2.males.PI.melt$value)

#For unknown tanks will want to keep tank as part of list
data.phase.2.brood.total.PI=aggregate(data.phase.2.males.PI.melt$value.exp,list(Tank=data.phase.2.males.PI.melt$Tank,variable=data.phase.2.males.PI.melt$variable),sum)

data.phase.2.brood.total.PI$Match=ifelse(data.phase.2.brood.total.PI$Tank==substr(data.phase.2.brood.total.PI$variable,1,2),1,0)

data.phase.2.brood.total.PI=subset(data.phase.2.brood.total.PI,Match==1)

data.phase.2.brood.total.PI=data.phase.2.brood.total.PI[,1:3]

colnames(data.phase.2.brood.total.PI)=c('Tank','variable','brood.total')

data.phase.2.males.PI.melt.2=merge(data.phase.2.males.PI.melt,data.phase.2.brood.total.PI,all=T)

data.phase.2.males.PI.melt.2$Percent.PI=data.phase.2.males.PI.melt.2$value.exp/data.phase.2.males.PI.melt.2$brood.total

data.phase.2.males.PI.melt.3=subset(data.phase.2.males.PI.melt.2,Percent.PI>=0)

ggplot(data.phase.2.males.PI.melt.3)+geom_histogram(aes(Percent.PI,fill=value>=7))


ggplot(subset(data.phase.2.males.PI.melt.3))+geom_histogram(aes(Percent.PI,fill=value>=7))


#Select the top likely father candidate
data.phase.2.males.PI.melt.3$Percent.PI.filter=ifelse(data.phase.2.males.PI.melt.3$value>=27,data.phase.2.males.PI.melt.3$Percent.PI,NA)

Fish.Cols<-c("BLACK"="black","BLUE"="blue","BROWN"="brown","GREEN"="green","PURPLE"="purple","RED"="red","YELLOW"="yellow","WHITE"="grey","YELLOW_A"="yellow","YELLOW_B"="yellow")


ggplot(subset(data.phase.2.males.PI.melt.3,Tank=='D2'))+geom_bar(color='black',stat='identity',aes(x=variable, y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=Fish.Cols,aesthetics = c("colour", "fill"))+xlab('Brood name')


data.phase.2.males.PI.melt.3$Tank.2=as.factor(data.phase.2.males.PI.melt.3$Tank)
data.phase.2.males.PI.melt.3$Tank.2=as.numeric(data.phase.2.males.PI.melt.3$Tank.2)

ggplot(subset(data.phase.2.males.PI.melt.3))+geom_bar(color='black',stat='identity',aes(x=reorder(variable,Tank.2), y=Percent.PI.filter,fill=Sample))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+scale_color_manual(values=Fish.Cols,aesthetics = c("colour", "fill"))+xlab('Brood name')+geom_vline(xintercept = c(4.5,11.5,13.5,18.5,20.5),size=2)
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_no_pat_exclusion/Percent.PI.jpg',width = 10, height = 10)

#create paternity success dataframe
data.phase.2.paternity=data.phase.2.males.PI.melt.3
data.phase.2.paternity=subset(data.phase.2.paternity,Percent.PI>=0.01)
data.phase.2.paternity$Pat.success='Unknown father'
data.phase.2.paternity$Pat.success=ifelse(data.phase.2.paternity$value>=27&data.phase.2.paternity$Percent.PI>=0,'Partial father',data.phase.2.paternity$Pat.success)
data.phase.2.paternity$Pat.success=ifelse(data.phase.2.paternity$value>=27&data.phase.2.paternity$Percent.PI>=0.99,'Full father',data.phase.2.paternity$Pat.success)
data.phase.2.paternity$Pat.identity=1
data.phase.2.paternity$Pat.identity=ifelse(data.phase.2.paternity$Pat.success=='Unknown father'&data.phase.2.paternity$Percent.PI<=0.9,2,data.phase.2.paternity$Pat.identity)
data.phase.2.paternity=subset(data.phase.2.paternity,Pat.identity==1)


ggplot(data.phase.2.paternity)+geom_tile(color='black',stat='identity',aes(x=reorder(variable,Tank.2), y=Pat.identity,fill=Pat.success))+ theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank())+xlab('Brood by tank')+ylab('')+geom_vline(xintercept = c(4.5,11.5,13.5,18.5,20.5),size=2)
ggsave('Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_no_pat_exclusion/Paternity success.jpg',width = 10, height = 10)

data.phase.2.paternity=data.phase.2.paternity[,c(1:3,5,11,13)]
write.csv(data.phase.2.paternity,file='Phase_2_Pat_Brood/Filter_0.6_filter_genotype_adults_all_pop_no_pat_exclusion/data.phase.2.paternity.csv')
#############

########Trash


##Paternity Test B2.ORANGE.8.24
#data.B2.pat.males$PI.B2.ORANGE.8.24=mclapply(data.B2.pat.males$Sample.name,function(x) log(prod(subset(paternityIndex(dataframe =data.phase.2.B2.sample.trans, marker = "Marker", mother ="ORANGE.B2", child = "B2.ORANGE.8.24", AF = x, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.B2.sample.trans)),paternityIndex(dataframe =data.phase.2.B2.sample.trans, marker = "Marker", mother ="ORANGE.B2", child = "B2.ORANGE.8.24", AF = x, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.B2.sample.trans))$PI>0)$PI)),mc.cores = 16)

#data.B2.pat.males.2=mapply(function(x, y) cbind(data.B2.pat.males, setNames( mclapply(data.B2.pat.males$Sample.name,function(z) log(prod(subset(paternityIndex(dataframe =data.phase.2.B2.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.B2.sample.trans)),paternityIndex(dataframe =data.phase.2.B2.sample.trans, marker = "Marker", mother =y, child = x, AF = z, Afreq = "calculate", Bfreq = "calculate", afcalcol = 2:ncol(data.phase.2.B2.sample.trans))$PI>0)$PI)),mc.cores = 16), data.B2.pat.brood$Sample.name)),data.B2.pat.brood$Sample.name,  #brood names as x data.B2.pat.brood$Mother) #mother names as y

#cbind(data.B2.pat.males, setNames( lapply(data.B2.pat.brood$Sample.name, function(x) x=NA), data.B2.pat.brood$Sample.name) )

#Graph
#plots.B2=lapply(B2.brood.names, function(x) ggplot(data.B2.pat.males.PI)+geom_bar(stat='identity',aes(x=reorder(Sample.name,-data.B2.pat.males[x]),y=x,fill=data.B2.pat.males$data.B2.pat.males[x]>=0))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE))
  
#plots.B2=lapply(B2.brood.names, function(x) ggplot(data.B2.pat.males)+geom_bar(stat='identity',aes(x=Sample.name,y=B2.brood.names[x]))+ggtitle(paste("Brood",x,"Paternity Results",sep=" "))+xlab('Males')+ylab('Log composite PI')+theme(axis.text.x = element_text(angle=90,vjust=0.5))+guides(fill=FALSE))


#lapply(names(plots.B2), function(x) ggsave(filename=paste("Phase_2_Pat_Brood/",x,".jpeg",sep=""), plot=plots.B2[[x]]))


#plotserieslines <- function(yvar){
#  ggplot(data.B2.pat.males, aes(x=Sample.name,y=as.name(yvar))) +
#    geom_bar(stat='identity')
#}
#lapply(names(data.B2.pat.males[5:ncol(data.B2.pat.males)]), plotserieslines)


#######convert to data.frame
#data.B2.pat.males.PI.colnames=data.frame(matrix(unlist(data.B2.pat.males.PI), nrow=nrow(data.B2.pat.males.2)))
#colnames(data.B2.pat.males.PI.colnames)=colnames(data.B2.pat.males.PI)



