##Create files needed for WOPI function
##remove info sites and compare prod of outlier tails 

#set working directory to source file location
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity_Testing/2bRAD_data/Data")

#load library
library(reshape2)
library(tidyverse) #load last

##### Load in sample data ####
#Tank Examples Phase 1 to go with 'Pat geno likelihoods filtered.R'
#Filtered to variable sites from filtered2.pos
#Each individual has three rows and sample data
#Column for each marker with genotype probabilities
load('data.geno.filter.trans.RData')

#convert into long format for marker
data.geno.filter.trans.melt=melt(data.geno.filter.trans,id=c(colnames(data.geno.filter.trans)[1:8]))

#Keep just phase 1 samples
data.geno.filter.trans.melt.phase.1 = data.geno.filter.trans.melt %>% 
  filter(Phase==1)

###Remove zero info sites
data.geno.filter.trans.melt.phase.1.zero = data.geno.filter.trans.melt.phase.1 %>% 
  group_by(variable, Sample.name) %>% 
  mutate(Info.score = sum(max(value)-min(value))) %>% 
  filter(Info.score != 0)

#convert to dataframe
data.geno.filter.trans.melt.phase.1.zero = 
  as.data.frame(data.geno.filter.trans.melt.phase.1.zero)

###create list of brood and mothers
##create list of broods
brood.list = data.geno.filter.trans.melt.phase.1.zero %>%
  filter(Sex == 'B') %>%
  distinct(Sample.name, .keep_all = T )
#remove extra columns
brood.list = brood.list[,c(2,6,7)]
#rename to brood.name
colnames(brood.list)[1]='Brood.name'
##create list of mothers
mother.list = data.geno.filter.trans.melt.phase.1.zero %>%
  filter(Sex == 'F') %>%
  distinct(Sample.name, .keep_all = T )
#remove extra columns
mother.list = mother.list[,c(2,6,7)]
#rename to mother.name
colnames(mother.list)[1]='Mother.name'
#remove mother Blue_B.B5A
mother.list = mother.list %>% 
  filter(Mother.name!='BLUE_B.B5A')
#change blue_A to blue
mother.list[13,2]='BLUE'
##merge brood and mother names
#note this list is shorter because no Blue.B5A mother
brood.mother.list = merge(brood.list,mother.list)

#create Child and Tank dataframe
child.tank = brood.list[,c(1,3)]
#rename col to Child
colnames(child.tank)[1] = 'Child'

# ###Get A and B allele freq
# #only use adults
# data.geno.filter.trans.melt.adult=data.geno.filter.trans.melt %>% 
#   filter(Sex!='B' & Phase == 1)
# #count number of A alleles weighted by genotype probability
# data.geno.filter.trans.melt.adult$Acount=ifelse(data.geno.filter.trans.melt.adult$MajorMinor=='AA', 2*data.geno.filter.trans.melt.adult$value,
#                                                 ifelse(data.geno.filter.trans.melt.adult$MajorMinor=='AB', 1*data.geno.filter.trans.melt.adult$value,0))
# #count number of B alleles weighted by genotype probability
# data.geno.filter.trans.melt.adult$Bcount=ifelse(data.geno.filter.trans.melt.adult$MajorMinor=='BB', 2*data.geno.filter.trans.melt.adult$value,
#                                                 ifelse(data.geno.filter.trans.melt.adult$MajorMinor=='AB', 1*data.geno.filter.trans.melt.adult$value,0))
# 
# #create new dataframe of allele sums per marker
# allelefreq.phase1=data.geno.filter.trans.melt.adult %>% 
#   group_by(variable) %>% 
#   summarize(Asum=sum(Acount,na.rm=T),Bsum=sum(Bcount,na.rm=T))
# #calculate total number of alleles per marker
# allelefreq.phase1$Total=allelefreq.phase1$Asum+allelefreq.phase1$Bsum
# #calc A allele freq
# allelefreq.phase1$Afreq=allelefreq.phase1$Asum/allelefreq.phase1$Total
# #calc B allele freq
# allelefreq.phase1$Bfreq=allelefreq.phase1$Bsum/allelefreq.phase1$Total
# #calc major allele freq
# allelefreq.phase1$MajorAllele=pmax(allelefreq.phase1$Afreq,allelefreq.phase1$Bfreq)



##### Phase 2  sample data ####
#Tank Examples Phase 2 to go with 'Pat geno likelihoods filtered.R'
#Filtered to variable sites from filtered2.pos
#Each individual has three rows and sample data
#Column for each marker with genotype probabilities
#if not loaded before for phase 1 run:
#make sure to load and convert data
# load('data.geno.filter.trans.RData')
# 
# #convert into long format for marker
# data.geno.filter.trans.melt=melt(data.geno.filter.trans,id=c(colnames(data.geno.filter.trans)[1:8]))

#Keep just phase 2 samples
data.geno.filter.trans.melt.phase.2 = data.geno.filter.trans.melt %>% 
  filter(Phase==2) 

###Remove zero info sites
data.geno.filter.trans.melt.phase.2.zero = data.geno.filter.trans.melt.phase.2 %>% 
  group_by(variable, Sample.name) %>% 
  mutate(Info.score = sum(max(value)-min(value))) %>% 
  filter(Info.score != 0)

#convert to dataframe
data.geno.filter.trans.melt.phase.2.zero = 
  as.data.frame(data.geno.filter.trans.melt.phase.2.zero %>% 
                  droplevels())

#count number of males per tank
# data.geno.filter.trans.melt.phase.2.zero %>%
#   filter(Sex=='M') %>%
#   select(Tank,Sample.name) %>%
#   filter(!str_detect(Sample.name, '_B.')) %>%
#   distinct() %>%
#   mutate(Present = 1) %>%
#   group_by(Tank) %>%
#   summarize(total = sum(Present))

###create list of brood and mothers
##create list of broods
brood.list.phase2 = data.geno.filter.trans.melt.phase.2.zero %>%
  filter(Sex == 'B') %>%
  distinct(Sample.name, .keep_all = T )
#remove extra columns
brood.list.phase2 = brood.list.phase2[,c(2,6,7)]
#rename to brood.name
colnames(brood.list.phase2)[1]='Brood.name'
##create list of mothers
mother.list.phase2 = data.geno.filter.trans.melt.phase.2.zero %>%
  filter(Sex == 'F') %>%
  distinct(Sample.name, .keep_all = T )
#remove extra columns
mother.list.phase2 = mother.list.phase2[,c(2,6,7)]
#rename to mother.name
colnames(mother.list.phase2)[1]='Mother.name'
##merge brood and mother names
#note this list is shorter because no Blue.B5A mother
brood.mother.list.phase2 = merge(brood.list.phase2,mother.list.phase2)

#create tank list
tank.list = data.geno.filter.trans.melt.phase.2.zero %>% 
  pull(Tank) %>% 
  unique() %>% 
  droplevels() %>% 
  as.character()

##### remove objects and save required ones ####
rm(data.geno.filter.trans)
rm(data.geno.filter.trans.melt)
rm(data.geno.filter.trans.melt.phase.1)
rm(data.geno.filter.trans.melt.phase.2)

save.image(file = 'WOPI.example.data.RData')







