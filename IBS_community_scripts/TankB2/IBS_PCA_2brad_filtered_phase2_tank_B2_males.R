##IBS and covariance matrix of 2brad data
#IBS matrix used fitered sites and run on phase 2 tank B2 males

#set working directory
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity_Testing/2bRAD_data/Data/Matrix/Phase2/TankB2")

#load libraries
library(tidyverse)
library(ggfortify)
library(pvclust)
library(gplots)
library(igraph)
library(GGally)


#read in bam/names order file
#change to list of bams name used
samples=read.delim('bamsphase2b2males',header=F)
colnames(samples)='File.name'
samples$order=seq(1, nrow(samples), by=1)

#change to sample names
data.bam.names=read.csv('bams.ind.csv')
rownames(data.bam.names)=data.bam.names$File.name.ind
data.bam.names.reduced=data.bam.names[,c(2,4:6,12:14)]
#merge sample names and counts data
sample.data=merge(data.bam.names.reduced,samples)

#load in ibs matrix
m.ibs <- as.matrix(read.table("IBS.TankB2.males.ibsMat"))
colnames(m.ibs)=sample.data$Sample.name
rownames(m.ibs)=sample.data$Sample.name

#load in covariance matrix
m.cov <- as.matrix(read.table("IBS.TankB2.males.covMat"))
colnames(m.cov)=sample.data$Sample.name
rownames(m.cov)=sample.data$Sample.name

######
##subsetting matrices
#subset names by sex
male.names=subset(sample.data$Sample.name,sample.data$Sex=='M')
female.names=subset(sample.data$Sample.name,sample.data$Sex=='F')
brood.names=subset(sample.data$Sample.name,sample.data$Sex=='B')
adult.names=brood.names=subset(sample.data$Sample.name,sample.data$Sex!='B')

#subset names by phase
phase.1.names=subset(sample.data$Sample.name,sample.data$Phase==1)
phase.2.names=subset(sample.data$Sample.name,sample.data$Phase==2)

#subset phase 1 by no females
phase.1.names.no.females=subset(sample.data$Sample.name,sample.data$Phase==1&sample.data$Sex!='F')

#subset phase 1 by no males
phase.1.names.no.males=subset(sample.data$Sample.name,sample.data$Phase==1&sample.data$Sex!='M')


#subset phase 2 by tank
phase.2.B2.names=subset(sample.data$Sample.name,sample.data$Tank=='B2')
phase.2.C2.names=subset(sample.data$Sample.name,sample.data$Tank=='C2')
phase.2.D2.names=subset(sample.data$Sample.name,sample.data$Tank=='D2')
phase.2.E2.names=subset(sample.data$Sample.name,sample.data$Tank=='E2')
phase.2.F2.names=subset(sample.data$Sample.name,sample.data$Tank=='F2')
phase.2.G2.names=subset(sample.data$Sample.name,sample.data$Tank=='G2')

#subset phase 2 by tank and no females
phase.2.B2.names.no.females=subset(sample.data$Sample.name,sample.data$Tank=='B2'&sample.data$Sex!='F')
phase.2.C2.names.no.females=subset(sample.data$Sample.name,sample.data$Tank=='C2'&sample.data$Sex!='F')
phase.2.D2.names.no.females=subset(sample.data$Sample.name,sample.data$Tank=='D2'&sample.data$Sex!='F')
phase.2.E2.names.no.females=subset(sample.data$Sample.name,sample.data$Tank=='E2'&sample.data$Sex!='F')
phase.2.F2.names.no.females=subset(sample.data$Sample.name,sample.data$Tank=='F2'&sample.data$Sex!='F')
phase.2.G2.names.no.females=subset(sample.data$Sample.name,sample.data$Tank=='G2'&sample.data$Sex!='F')

#subset phase 2 by tank and no males
phase.2.B2.names.no.males=subset(sample.data$Sample.name,sample.data$Tank=='B2'&sample.data$Sex!='M')
phase.2.C2.names.no.males=subset(sample.data$Sample.name,sample.data$Tank=='C2'&sample.data$Sex!='M')
phase.2.D2.names.no.males=subset(sample.data$Sample.name,sample.data$Tank=='D2'&sample.data$Sex!='M')
phase.2.E2.names.no.males=subset(sample.data$Sample.name,sample.data$Tank=='E2'&sample.data$Sex!='M')
phase.2.F2.names.no.males=subset(sample.data$Sample.name,sample.data$Tank=='F2'&sample.data$Sex!='M')
phase.2.G2.names.no.males=subset(sample.data$Sample.name,sample.data$Tank=='G2'&sample.data$Sex!='M')

#subset phase 2 by tank and males
phase.2.B2.names.males=subset(sample.data$Sample.name,sample.data$Tank=='B2'&sample.data$Sex=='M')
phase.2.C2.names.males=subset(sample.data$Sample.name,sample.data$Tank=='C2'&sample.data$Sex=='M')
phase.2.D2.names.males=subset(sample.data$Sample.name,sample.data$Tank=='D2'&sample.data$Sex=='M')
phase.2.E2.names.males=subset(sample.data$Sample.name,sample.data$Tank=='E2'&sample.data$Sex=='M')
phase.2.F2.names.males=subset(sample.data$Sample.name,sample.data$Tank=='F2'&sample.data$Sex=='M')
phase.2.G2.names.males=subset(sample.data$Sample.name,sample.data$Tank=='G2'&sample.data$Sex=='M')

#subset phase 2 by tank and males
phase.2.B2.names.females=subset(sample.data$Sample.name,sample.data$Tank=='B2'&sample.data$Sex=='F')
phase.2.C2.names.females=subset(sample.data$Sample.name,sample.data$Tank=='C2'&sample.data$Sex=='F')
phase.2.D2.names.females=subset(sample.data$Sample.name,sample.data$Tank=='D2'&sample.data$Sex=='F')
phase.2.E2.names.females=subset(sample.data$Sample.name,sample.data$Tank=='E2'&sample.data$Sex=='F')
phase.2.F2.names.females=subset(sample.data$Sample.name,sample.data$Tank=='F2'&sample.data$Sex=='F')
phase.2.G2.names.females=subset(sample.data$Sample.name,sample.data$Tank=='G2'&sample.data$Sex=='F')

#########
## clustering  of tank B2 in phase 2
#subset IBS matrix by phase 2 Tank B2
#m.ibs.B2=as.data.frame(m.ibs)
#m.ibs.B2=m.ibs.B2[rownames(m.ibs.B2) %in% phase.2.B2.names,colnames(m.ibs.B2) %in% phase.2.B2.names]
#m.ibs.B2=data.matrix(m.ibs.B2)

#IBS matrix whole tank
#png("pvclut.phase.2.tank.B2.png", width = 1800, height = 700)
#plot(pvclust(m.ibs.B2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2')
#dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.B2.no.females=as.data.frame(m.ibs)
m.ibs.B2.no.females=m.ibs.B2.no.females[rownames(m.ibs.B2.no.females) %in% phase.2.B2.names.no.females,colnames(m.ibs.B2.no.females) %in% phase.2.B2.names.no.females]
m.ibs.B2.no.females=data.matrix(m.ibs.B2.no.females)

png("pvclut.phase.2.tank.B2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.B2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2 No Females')
dev.off()

#network undirected matrix with edge's labeled
g=graph_from_adjacency_matrix(m.ibs.B2.no.females,mode='undirected',weighted = TRUE)
layout <- layout.mds(g, dist = as.matrix(m.ibs.B2.no.females))

#color nodes
g=set_vertex_attr(g, "Sex", value=sample.data$Sex)
V(g)$color=V(g)$Sex 
V(g)$color=gsub("2","red",V(g)$color) 
V(g)$color=gsub("1","blue",V(g)$color)

png("mds.network.phase.2.tank.B2.no.females.png", width = 1800, height = 700)
plot(g, layout = layout, vertex.size = 3, edge.label=round(E(g)$weight, 3), main='Tank B2 Males')
dev.off()

#just edges under 0.2
m.ibs.B2.no.females.2=ifelse(m.ibs.B2.no.females>=0.2,0,m.ibs.B2.no.females)
h=graph_from_adjacency_matrix(m.ibs.B2.no.females.2,mode='undirected',weighted = TRUE)

#color nodes
h=set_vertex_attr(h, "Sex", value=sample.data$Sex)
V(h)$color=V(h)$Sex 
V(h)$color=gsub("2","red",V(h)$color) 
V(h)$color=gsub("1","blue",V(h)$color)

png("cut.mds.network.phase.2.tank.B2.no.females.png", width = 1800, height = 700)
plot(h, layout = layout, vertex.size = 3, edge.label=round(E(h)$weight, 3), main= 'Tank B2 Males Cut')
dev.off()

##subset IBS matrix by phase 1 males
m.ibs.B2.males=as.data.frame(m.ibs)
m.ibs.B2.males=m.ibs.B2.males[rownames(m.ibs.B2.males) %in% phase.2.B2.names.males,colnames(m.ibs.B2.males) %in% phase.2.B2.names.males]
m.ibs.B2.males=data.matrix(m.ibs.B2.males)

png("pvclut.phase.2.tank.B2.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.B2.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2 Males')
dev.off()


#network undirected matrix with edge's labeled
g=graph_from_adjacency_matrix(m.ibs.B2.males,mode='undirected',weighted = TRUE)
layout <- layout.mds(g, dist = as.matrix(m.ibs.B2.males))

png("mds.network.phase.2.tank.B2.males.png", width = 1800, height = 700)
plot(g, layout = layout, vertex.size = 3, edge.label=round(E(g)$weight, 3), main='Tank B2 Males')
dev.off()

#just edges under 0.2
m.ibs.B2.males.2=ifelse(m.ibs.B2.males>=0.2,0,m.ibs.B2.males)
h=graph_from_adjacency_matrix(m.ibs.B2.males.2,mode='undirected',weighted = TRUE)

png("cut.mds.network.phase.2.tank.B2.males.png", width = 1800, height = 700)
plot(h, layout = layout, vertex.size = 3, edge.label=round(E(h)$weight, 3), main= 'Tank B2 Males Cut')
dev.off()

##ggnet2

ggnet2(h, mode= layout,label=T, edge.label=round(E(h)$weight, 3))+ggtitle('Cut Tank B2 Males')+scale_y_continuous(limits = c(-1,1), breaks = NULL)+scale_x_continuous(limits = c(-1,1), breaks = NULL)





#subset IBS matrix by phase 1 no males
m.ibs.B2.no.males=as.data.frame(m.ibs)
m.ibs.B2.no.males=m.ibs.B2.no.males[rownames(m.ibs.B2.no.males) %in% phase.2.B2.names.no.males,colnames(m.ibs.B2.no.males) %in% phase.2.B2.names.no.males]
m.ibs.B2.no.males=data.matrix(m.ibs.B2.no.males)

png("pvclut.phase.2.tank.B2.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.B2.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2 No males')
dev.off()

#subset IBS matrix by phase 1 females
#m.ibs.B2.females=as.data.frame(m.ibs)
#m.ibs.B2.females=m.ibs.B2.females[rownames(m.ibs.B2.females) %in% phase.2.B2.names.females,colnames(m.ibs.B2.females) %in% phase.2.B2.names.females]
#m.ibs.B2.females=data.matrix(m.ibs.B2.females)

#png("pvclut.phase.2.tank.B2.females.png", width = 1800, height = 700)
#plot(pvclust(m.ibs.B2.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2 females')
#dev.off()

