##IBS and covariance matrix of 2brad data
#IBS matrix used filtered sites and only run on broods and males from phase 1

#set working directory
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity_Testing/2bRAD_data/Data/Matrix/Matrix_filter_sites_phase1")

#load libraries
library(tidyverse)
library(ggfortify)
library(pvclust)
library(gplots)
library(ggdendro)
library(dendextend)


#####
#read in bam/names order file
#change to list of bams name used
samples=read.delim('bamsphase1',header=F)
colnames(samples)='File.name'
samples$order=seq(1, nrow(samples), by=1)

#change to sample names
data.bam.names=read.csv('bams.ind.csv')
rownames(data.bam.names)=data.bam.names$File.name.ind
data.bam.names.reduced=data.bam.names[,c(2,4:6,12:14)]
#merge sample names and counts data
sample.data=merge(data.bam.names.reduced,samples)

#load in ibs matrix
m.ibs <- as.matrix(read.table("ddB.IBS.phase1all.ibsMat"))
colnames(m.ibs)=sample.data$Sample.name
rownames(m.ibs)=sample.data$Sample.name

#load in covariance matrix
m.cov <- as.matrix(read.table("ddB.IBS.phase1all.covMat"))
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

#########
## heatmap / clustering / trees of all
#heat map ibs matrix
heatmap(m.ibs)
#neighbour joining
plot(ape::nj(m.ibs))
plot(hclust(dist(m.ibs), "ave"))
plot(pvclust(m.ibs, method.hclust = "average",method.dist = 'euclidian'), main='All Fish')

## MDS with ibs matrix
mds <- cmdscale(as.dist(m.ibs))
plot(mds,lwd=2,ylab="Dist",xlab="Dist",main="multidimensional scaling",col=rep(1:3,each=10))
#autoplot(mds)

##PCA with covariance
e <- eigen(m.cov)
plot(e$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 1",main="Principal components",col=rep(1:3,each=10),pch=16)
#autoplot(prcomp(m.cov))

#########
## heatmap / clustering / trees of phase 1
#subset IBS matrix by phase 1
m.ibs.phase.1=as.data.frame(m.ibs)
m.ibs.phase.1=m.ibs.phase.1[rownames(m.ibs.phase.1) %in% phase.1.names,colnames(m.ibs.phase.1) %in% phase.1.names]
m.ibs.phase.1=data.matrix(m.ibs.phase.1)
#subset cov matrix by phase 1
m.cov.phase.1=as.data.frame(m.cov)
m.cov.phase.1=m.cov.phase.1[rownames(m.cov.phase.1) %in% phase.1.names,colnames(m.cov.phase.1) %in% phase.1.names]
m.cov.phase.1=data.matrix(m.cov.phase.1)
#heat map ibs matrix
#heatmap(m.ibs.phase.1)
#heatmap.2(m.ibs.phase.1,dendrogram="none")
#neighbour joining
#plot(ape::nj(m.ibs.phase.1))
#plot(hclust(dist(m.ibs.phase.1), "ave"))

png("pvclut.phase.1.png", width = 1800, height = 700)
plot(pvclust(m.ibs.phase.1, method.hclust = "average",method.dist = 'euclidian'), main='Pvclust Phase 1')
dev.off()

#subset IBS matrix by phase 1 males and broods
m.ibs.phase.1.no.females=as.data.frame(m.ibs)
m.ibs.phase.1.no.females=m.ibs.phase.1.no.females[rownames(m.ibs.phase.1.no.females) %in% phase.1.names.no.females,colnames(m.ibs.phase.1.no.females) %in% phase.1.names.no.females]
m.ibs.phase.1.no.females=data.matrix(m.ibs.phase.1.no.females)
#neighbour joining
png("pvclut.phase.1.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.phase.1.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Pvclust Phase 1 No Females')
dev.off()

#subset IBS matrix by phase 1 females and broods
m.ibs.phase.1.no.males=as.data.frame(m.ibs)
m.ibs.phase.1.no.males=m.ibs.phase.1.no.males[rownames(m.ibs.phase.1.no.males) %in% phase.1.names.no.males,colnames(m.ibs.phase.1.no.males) %in% phase.1.names.no.males]
m.ibs.phase.1.no.males=data.matrix(m.ibs.phase.1.no.males)
#neighbour joining
png("pvclut.phase.1.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.phase.1.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Pvclust Phase 1 No males')
dev.off()

#subset phase 1 by no broods
phase.1.names.no.broods=subset(sample.data$Sample.name,sample.data$Phase==1&sample.data$Sex!='B')
#subset IBS matrix by phase 1 no broods
m.ibs.phase.1.no.broods=as.data.frame(m.ibs)
m.ibs.phase.1.no.broods=m.ibs.phase.1.no.broods[rownames(m.ibs.phase.1.no.broods) %in% phase.1.names.no.broods,colnames(m.ibs.phase.1.no.broods) %in% phase.1.names.no.broods]
m.ibs.phase.1.no.broods=data.matrix(m.ibs.phase.1.no.broods)
#neighbour joining
png("pvclut.phase.1.no.broods.png", width = 1800, height = 700)
plot(pvclust(m.ibs.phase.1.no.broods, method.hclust = "average",method.dist = 'euclidian'), main='Pvclust Phase 1 No broods')
dev.off()

## MDS with ibs matrix
#mds.phase.1 <- cmdscale(as.dist(m.ibs.phase.1))
#plot(mds.phase.1,lwd=2,ylab="Dist",xlab="Dist",main="multidimensional scaling",col=rep(1:3,each=10))
#autoplot(mds.phase.1)
##PCA with covariance
#e.phase.1 <- eigen(m.cov.phase.1)
#plot(e.phase.1$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 1",main="Principal components",col=rep(1:3,each=10),pch=16)
#autoplot(prcomp(m.cov.phase.1))

#########
###Sample IBS from known paternity for poster
# sample: B5A.YELLOW.7.6
## heatmap / clustering / trees of phase 1

#create list of fathers and sample brood
phase.1.names.males.B5A.YELLOW.7.6=rownames(m.ibs.phase.1.no.females)[c(1:8,15)]

#subset IBS matrix by phase 1 males and broods
m.ibs.phase.1.no.females.B5A.YELLOW.7.6=as.data.frame(m.ibs)
m.ibs.phase.1.no.females.B5A.YELLOW.7.6=m.ibs.phase.1.no.females.B5A.YELLOW.7.6[rownames(m.ibs.phase.1.no.females.B5A.YELLOW.7.6) %in% phase.1.names.males.B5A.YELLOW.7.6,colnames(m.ibs.phase.1.no.females.B5A.YELLOW.7.6) %in% phase.1.names.males.B5A.YELLOW.7.6]
phase.1.names.males.B5A.YELLOW.7.6=data.matrix(phase.1.names.males.B5A.YELLOW.7.6)
#neighbour joining
clust.phase.1.no.females.B5A.YELLOW.7.6=pvclust(m.ibs.phase.1.no.females.B5A.YELLOW.7.6, method.hclust = "average",method.dist = 'euclidian')
# #graph
# pdf("pvclut.phase.1.no.females.B5A.YELLOW.7.6.pdf", width = 1200, height = 550)
# plot(clust.phase.1.no.females.B5A.YELLOW.7.6, main='Pvclust Phase 1 No Females B5A.YELLOW.7.6')
# dev.off()

#rename to males
m.ibs.phase.1.no.females.B5A.YELLOW.7.6.names=m.ibs.phase.1.no.females.B5A.YELLOW.7.6
males.brood.name=c('Father','Male 1','Male 2','Male 3','Male 4','Male 5','Male 6','Male 7','Offspring')

rownames(m.ibs.phase.1.no.females.B5A.YELLOW.7.6.names) =males.brood.name
colnames(m.ibs.phase.1.no.females.B5A.YELLOW.7.6.names) =males.brood.name
#neighbour joining
clust.phase.1.no.females.B5A.YELLOW.7.6.names=pvclust(m.ibs.phase.1.no.females.B5A.YELLOW.7.6.names, method.hclust = "average",method.dist = 'euclidian')

#graph
pdf("pvclut.phase.1.no.females.B5A.YELLOW.7.6.pdf", width = 12, height = 6,pointsize = 28)
par(mar=c(0,0,0,0))
plot(clust.phase.1.no.females.B5A.YELLOW.7.6.names, main='',lwd=5,ylab='',axes=F)
dev.off()

#create list of mothers and sample brood
phase.1.names.no.males.B5A.YELLOW.7.6=droplevels(phase.1.names.no.males[c(1:13,19)])

#subset IBS matrix by phase 1 females and broods
m.ibs.phase.1.no.males.B5A.YELLOW.7.6=as.data.frame(m.ibs)
m.ibs.phase.1.no.males.B5A.YELLOW.7.6=m.ibs.phase.1.no.males.B5A.YELLOW.7.6[rownames(m.ibs.phase.1.no.males.B5A.YELLOW.7.6) %in% phase.1.names.no.males.B5A.YELLOW.7.6,colnames(m.ibs.phase.1.no.males.B5A.YELLOW.7.6) %in% phase.1.names.no.males.B5A.YELLOW.7.6]
m.ibs.phase.1.no.males.B5A.YELLOW.7.6=data.matrix(m.ibs.phase.1.no.males.B5A.YELLOW.7.6)
#neighbour joining
clust.phase.1.no.males.B5A.YELLOW.7.6=pvclust(m.ibs.phase.1.no.males.B5A.YELLOW.7.6, method.hclust = "average",method.dist = 'euclidian')

# #graph
# pdf("pvclut.phase.1.no.males.B5A.YELLOW.7.6.pdf", width = 12, height = 5.5)
# plot(clust.phase.1.no.males.B5A.YELLOW.7.6, main='Pvclust Phase 1 No males B5A.YELLOW.7.6')
# dev.off()

#rename to females
m.ibs.phase.1.no.males.B5A.YELLOW.7.6.names=m.ibs.phase.1.no.males.B5A.YELLOW.7.6
females.brood.name=c('Female 1','Female 2','Female 3','Female 4','Female 5','Female 6','Mother','Female 7','Female 8','Female 9','Female 10','Female 11','Female 12','Offspring')

rownames(m.ibs.phase.1.no.males.B5A.YELLOW.7.6.names) =females.brood.name
colnames(m.ibs.phase.1.no.males.B5A.YELLOW.7.6.names) =females.brood.name
#neighbour joining
clust.phase.1.no.males.B5A.YELLOW.7.6.names=pvclust(m.ibs.phase.1.no.males.B5A.YELLOW.7.6.names, method.hclust = "average",method.dist = 'euclidian')

#graph
pdf("pvclut.phase.1.no.males.B5A.YELLOW.7.6.pdf", width = 12, height = 6,pointsize = 28)
par(mar=c(0,0,0,0))
plot(clust.phase.1.no.males.B5A.YELLOW.7.6.names, main='',lwd=5,ylab='',axes=F)
dev.off()

# #graph using hclust
# dd.no.males.B5A.YELLOW.7.6 <- dist(m.ibs.phase.1.no.males.B5A.YELLOW.7.6.names, method = "euclidean")
# hc.no.males.B5A.YELLOW.7.6 <- hclust(dd.no.males.B5A.YELLOW.7.6, method = "average")
# # Build dendrogram object from hclust results
# dend.no.males.B5A.YELLOW.7.6 <- as.dendrogram(hc.no.males.B5A.YELLOW.7.6)
# 
# # Extract the data (for rectangular lines)
# # Type can be "rectangle" or "triangle"
# dend_data.no.males.B5A.YELLOW.7.6 <- dendro_data(dend.no.males.B5A.YELLOW.7.6, type = "rectangle")
# ggplot(dend_data.no.males.B5A.YELLOW.7.6$segments) +
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
#   theme(legend.position = 'none',axis.title.x = element_blank(), axis.line.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank())+
#   geom_text(data = dend_data.no.males.B5A.YELLOW.7.6$labels, aes(x, y, label = label), hjust =1, angle = 90, size = 28)
# ggsave('pvclut.phase.1.no.males.B5A.YELLOW.7.6.pdf',width=6,height = 7, dpi=320, units = 'in')



#########
## heatmap / clustering / trees of adults
#subset IBS matrix by phase 1
m.ibs.sex=as.data.frame(m.ibs)
m.ibs.sex=m.ibs.sex[rownames(m.ibs.sex) %in% adult.names,colnames(m.ibs.sex) %in% adult.names]
m.ibs.sex=data.matrix(m.ibs.sex)
#subset cov matrix by phase 1
m.cov.sex=as.data.frame(m.cov)
m.cov.sex=m.cov.sex[rownames(m.cov.sex) %in% adult.names,colnames(m.cov.sex) %in% adult.names]
m.cov.sex=data.matrix(m.cov.sex)
#heat map ibs matrix
heatmap(m.ibs.sex)
#neighbour joining
plot(ape::nj(m.ibs.sex))
plot(hclust(dist(m.ibs.sex), "ave"))

png("pvclut.adults.png", width = 1800, height = 700)
plot(pvclust(m.ibs.sex, method.hclust = "average",method.dist = 'euclidian'), main='Pvclust Adults')
dev.off()

## MDS with ibs matrix
mds.sex <- cmdscale(as.dist(m.ibs.sex))
plot(mds.sex,lwd=2,ylab="Dist",xlab="Dist",main="multidimensional scaling",col=rep(1:3,each=10))
#autoplot(mds.sex, data = subset(sample.data,sample.data$Sex!='B'), colour = 'Sex')
##PCA with covariance
e.sex <- eigen(m.cov.sex)
plot(e.sex$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 1",main="Principal components",col=rep(1:3,each=10),pch=16)
autoplot(prcomp(m.cov.sex), data = subset(sample.data,sample.data$Sex!='B'), colour = 'Sex')

#########
## clustering  of tank B2 in phase 2
#subset IBS matrix by phase 1
m.ibs.B2=as.data.frame(m.ibs)
m.ibs.B2=m.ibs.B2[rownames(m.ibs.B2) %in% phase.2.B2.names,colnames(m.ibs.B2) %in% phase.2.B2.names]
m.ibs.B2=data.matrix(m.ibs.B2)

#IBS matrix whole tank
png("pvclut.phase.2.tank.B2.png", width = 1800, height = 700)
plot(pvclust(m.ibs.B2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2')
dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.B2.no.females=as.data.frame(m.ibs)
m.ibs.B2.no.females=m.ibs.B2.no.females[rownames(m.ibs.B2.no.females) %in% phase.2.B2.names.no.females,colnames(m.ibs.B2.no.females) %in% phase.2.B2.names.no.females]
m.ibs.B2.no.females=data.matrix(m.ibs.B2.no.females)

png("pvclut.phase.2.tank.B2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.B2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank B2 No Females')
dev.off()

#########
## clustering  of tank C2 in phase 2
#subset IBS matrix by phase 1
m.ibs.C2=as.data.frame(m.ibs)
m.ibs.C2=m.ibs.C2[rownames(m.ibs.C2) %in% phase.2.C2.names,colnames(m.ibs.C2) %in% phase.2.C2.names]
m.ibs.C2=data.matrix(m.ibs.C2)

#IBS matrix whole tank
png("pvclut.phase.2.tank.C2.png", width = 1800, height = 700)
plot(pvclust(m.ibs.C2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank C2')
dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.C2.no.females=as.data.frame(m.ibs)
m.ibs.C2.no.females=m.ibs.C2.no.females[rownames(m.ibs.C2.no.females) %in% phase.2.C2.names.no.females,colnames(m.ibs.C2.no.females) %in% phase.2.C2.names.no.females]
m.ibs.C2.no.females=data.matrix(m.ibs.C2.no.females)

png("pvclut.phase.2.tank.C2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.C2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank C2 No Females')
dev.off()

#subset IBS matrix by phase 1 no males
m.ibs.C2.no.males=as.data.frame(m.ibs)
m.ibs.C2.no.males=m.ibs.C2.no.males[rownames(m.ibs.C2.no.males) %in% phase.2.C2.names.no.males,colnames(m.ibs.C2.no.males) %in% phase.2.C2.names.no.males]
m.ibs.C2.no.males=data.matrix(m.ibs.C2.no.males)

png("pvclut.phase.2.tank.C2.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.C2.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank C2 No males')
dev.off()
#########
## clustering  of tank D2 in phase 2
#subset IBS matrix by phase 1
m.ibs.D2=as.data.frame(m.ibs)
m.ibs.D2=m.ibs.D2[rownames(m.ibs.D2) %in% phase.2.D2.names,colnames(m.ibs.D2) %in% phase.2.D2.names]
m.ibs.D2=data.matrix(m.ibs.D2)

#IBS matrix whole tank
png("pvclut.phase.2.tank.D2.png", width = 1800, height = 700)
plot(pvclust(m.ibs.D2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank D2')
dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.D2.no.females=as.data.frame(m.ibs)
m.ibs.D2.no.females=m.ibs.D2.no.females[rownames(m.ibs.D2.no.females) %in% phase.2.D2.names.no.females,colnames(m.ibs.D2.no.females) %in% phase.2.D2.names.no.females]
m.ibs.D2.no.females=data.matrix(m.ibs.D2.no.females)

png("pvclut.phase.2.tank.D2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.D2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank D2 No Females')
dev.off()

#subset IBS matrix by phase 1 no males
m.ibs.D2.no.males=as.data.frame(m.ibs)
m.ibs.D2.no.males=m.ibs.D2.no.males[rownames(m.ibs.D2.no.males) %in% phase.2.D2.names.no.males,colnames(m.ibs.D2.no.males) %in% phase.2.D2.names.no.males]
m.ibs.D2.no.males=data.matrix(m.ibs.D2.no.males)

png("pvclut.phase.2.tank.D2.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.D2.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank D2 No males')
dev.off()
#########
## clustering  of tank E2 in phase 2
#subset IBS matrix by phase 1
m.ibs.E2=as.data.frame(m.ibs)
m.ibs.E2=m.ibs.E2[rownames(m.ibs.E2) %in% phase.2.E2.names,colnames(m.ibs.E2) %in% phase.2.E2.names]
m.ibs.E2=data.matrix(m.ibs.E2)

#IBS matrix whole tank
png("pvclut.phase.2.tank.E2.png", width = 1800, height = 700)
plot(pvclust(m.ibs.E2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank E2')
dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.E2.no.females=as.data.frame(m.ibs)
m.ibs.E2.no.females=m.ibs.E2.no.females[rownames(m.ibs.E2.no.females) %in% phase.2.E2.names.no.females,colnames(m.ibs.E2.no.females) %in% phase.2.E2.names.no.females]
m.ibs.E2.no.females=data.matrix(m.ibs.E2.no.females)

png("pvclut.phase.2.tank.E2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.E2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank E2 No Females')
dev.off()

#subset IBS matrix by phase 1 no males
m.ibs.E2.no.males=as.data.frame(m.ibs)
m.ibs.E2.no.males=m.ibs.E2.no.males[rownames(m.ibs.E2.no.males) %in% phase.2.E2.names.no.males,colnames(m.ibs.E2.no.males) %in% phase.2.E2.names.no.males]
m.ibs.E2.no.males=data.matrix(m.ibs.E2.no.males)

png("pvclut.phase.2.tank.E2.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.E2.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank E2 No males')
dev.off()
#########
## clustering  of tank F2 in phase 2
#subset IBS matrix by phase 1
m.ibs.F2=as.data.frame(m.ibs)
m.ibs.F2=m.ibs.F2[rownames(m.ibs.F2) %in% phase.2.F2.names,colnames(m.ibs.F2) %in% phase.2.F2.names]
m.ibs.F2=data.matrix(m.ibs.F2)

#IBS matrix whole tank
png("pvclut.phase.2.tank.F2.png", width = 1800, height = 700)
plot(pvclust(m.ibs.F2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank F2')
dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.F2.no.females=as.data.frame(m.ibs)
m.ibs.F2.no.females=m.ibs.F2.no.females[rownames(m.ibs.F2.no.females) %in% phase.2.F2.names.no.females,colnames(m.ibs.F2.no.females) %in% phase.2.F2.names.no.females]
m.ibs.F2.no.females=data.matrix(m.ibs.F2.no.females)

png("pvclut.phase.2.tank.F2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.F2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank F2 No Females')
dev.off()

#subset IBS matrix by phase 1 no males
m.ibs.F2.no.males=as.data.frame(m.ibs)
m.ibs.F2.no.males=m.ibs.F2.no.males[rownames(m.ibs.F2.no.males) %in% phase.2.F2.names.no.males,colnames(m.ibs.F2.no.males) %in% phase.2.F2.names.no.males]
m.ibs.F2.no.males=data.matrix(m.ibs.F2.no.males)

png("pvclut.phase.2.tank.F2.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.F2.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank F2 No males')
dev.off()

#########
## clustering  of tank G2 in phase 2
#subset IBS matrix by phase 1
m.ibs.G2=as.data.frame(m.ibs)
m.ibs.G2=m.ibs.G2[rownames(m.ibs.G2) %in% phase.2.G2.names,colnames(m.ibs.G2) %in% phase.2.G2.names]
m.ibs.G2=data.matrix(m.ibs.G2)

#IBS matrix whole tank
png("pvclut.phase.2.tank.G2.png", width = 1800, height = 700)
plot(pvclust(m.ibs.G2, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank G2')
dev.off()

#subset IBS matrix by phase 1 no females
m.ibs.G2.no.females=as.data.frame(m.ibs)
m.ibs.G2.no.females=m.ibs.G2.no.females[rownames(m.ibs.G2.no.females) %in% phase.2.G2.names.no.females,colnames(m.ibs.G2.no.females) %in% phase.2.G2.names.no.females]
m.ibs.G2.no.females=data.matrix(m.ibs.G2.no.females)

png("pvclut.phase.2.tank.G2.no.females.png", width = 1800, height = 700)
plot(pvclust(m.ibs.G2.no.females, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank G2 No Females')
dev.off()
#subset IBS matrix by phase 1 no males
m.ibs.G2.no.males=as.data.frame(m.ibs)
m.ibs.G2.no.males=m.ibs.G2.no.males[rownames(m.ibs.G2.no.males) %in% phase.2.G2.names.no.males,colnames(m.ibs.G2.no.males) %in% phase.2.G2.names.no.males]
m.ibs.G2.no.males=data.matrix(m.ibs.G2.no.males)

png("pvclut.phase.2.tank.G2.no.males.png", width = 1800, height = 700)
plot(pvclust(m.ibs.G2.no.males, method.hclust = "average",method.dist = 'euclidian'), main='Phase 2 Tank G2 No males')
dev.off()