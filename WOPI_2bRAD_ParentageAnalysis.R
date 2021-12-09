##WOPI function
##remove info sites and compare prod of outlier tails 

#set working directory to source file location
setwd("/stor/work/Hofmann/projects/aburtoni_paternity_testing/Paternity_Testing/2bRAD_data/Data")

#load library
library(parallel)
library(pbmcapply)
library(pbapply)
library(reshape2)
library(tidyverse) #load last


##### Function for distribution of weighted PI ####
#Every marker needs an allele freq
#Works for missing, NA, and 0 genotype probabilities
totalWPI =
  function (dataframe = samples, Sample = "Sample",  Marker = "Markers", Allele = "Allele", GenotypeProbability = "Value", Type = "Type", Child = "C", KnownParent = 'M', AllegedParent = 'AF', allelefreq = 'calculate', Afreq = 'Afreq', Bfreq = 'Bfreq', AlleleMarker = 'Markers') {
    ##Create sample dataframe
    data = data.frame(Sample = dataframe[, Sample],
                      Marker = dataframe[, Marker],
                      Type = dataframe[, Type],
                      Allele = dataframe[, Allele],
                      GenotypeProbability = dataframe[, GenotypeProbability],
                      stringsAsFactors = FALSE)
    ##create allelefreq data
    if(allelefreq == 'calculate'){
      #only use adults
      data.adults=data %>% 
        filter(Type!=Child)
      #count number of A alleles weighted by genotype probability
      data.adults$Acount=ifelse(data.adults$Allele=='AA', 2*data.adults$GenotypeProbability,
                                ifelse(data.adults$Allele=='AB', 1*data.adults$GenotypeProbability,0))
      #count number of A alleles weighted by genotype probability
      data.adults$Bcount=ifelse(data.adults$Allele=='BB', 2*data.adults$GenotypeProbability,
                                ifelse(data.adults$Allele=='AB', 1*data.adults$GenotypeProbability,0))
      
      #create new dataframe of allele sums per marker
      allelefreq=data.adults %>% 
        group_by(Marker) %>% 
        summarize(Asum=sum(Acount,na.rm=T),Bsum=sum(Bcount,na.rm=T))
      #calculate total number of alleles per marker
      allelefreq$Total=allelefreq$Asum+allelefreq$Bsum
      #calc A allele freq
      allelefreq$Afreq=allelefreq$Asum/allelefreq$Total
      #calc B allele freq
      allelefreq$Bfreq=allelefreq$Bsum/allelefreq$Total
    }
    else{
      allelefreq = data.frame(Marker = allelefreq[ , AlleleMarker],
                              Afreq = allelefreq[, Afreq],
                              Bfreq = allelefreq[, Bfreq],
                              stringsAsFactors = FALSE)
    }
    
    #add type ID variable
    data$ID=paste(data$Type,data$Marker,data$Allele, sep='~')
    
    ##create dataframe of PI
    KnownParent.a <- c("AA", "AA", "AB", "AB", "BB", "BB", "AB", "AA", 
                       "AA", "AB", "AB", "BB", "BB", "AB", "AA", "AA", "AB", 
                       "AB", "BB", "BB", "AB")
    Child.a <- c("AA", "AB", "AA", "AB", "BB", "AB", "BB", "AA", 
                 "AB", "AA", "AB", "BB", "AB", "BB", "AA", "AB", "AA", 
                 "AB", "BB", "AB", "BB")
    AllegedParent.a <- c("AA", "AA", "AA", "AA", "AA", "AA", "AA", "AB", "AB", 
                         "AB", "AB", "AB", "AB", "AB", "BB", "BB", "BB", "BB", 
                         "BB", "BB", "BB")
    Formula <- c("1/a", 0, "1/a", "1/(a+b)", 0, "1/a", 0, "0.5/a", 
                 "0.5/a", "0.5/a", "1/(a+b)", "0.5 / b", "0.5/a", "0.5/b", 
                 0, "1/b", 0, "1/(a+b)", "1/b", 0, "1/b")
    COEFA <- c(1, 0, 1, 0, 0, 1, 0, 0.5, 0.5, 0.5, 0, 0, 0.5, 
               0, 0, 0, 0, 0, 0, 0, 0)
    COEFAB <- c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
                0, 0, 1, 0, 0, 0)
    COEFB <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0.5, 
               0, 1, 0, 0, 1, 0, 1)
    coefft <- data.frame(KnownParent.a, Child.a, AllegedParent.a, Formula,
                         COEFA, COEFAB, COEFB)
    
    ##merge allelefreq and coefft dataframes  
    WPI=merge(allelefreq,coefft)
    
    ##combine coefft with sample data
    #add type ID variable 
    WPI$MID=paste(KnownParent,WPI$Marker,WPI$KnownParent, sep='~')
    WPI$CID=paste(Child,WPI$Marker,WPI$Child, sep='~')
    WPI$AFID=paste(AllegedParent,WPI$Marker,WPI$AllegedParent, sep='~')
    
    #Use ID to add prob for KnownParent and Child
    WPI$M.prob=data$GenotypeProbability[match(WPI$MID,data$ID)]
    WPI$C.prob=data$GenotypeProbability[match(WPI$CID,data$ID)]
    
    #List of AllegedParent
    AllegedParent.samples= data %>% 
      filter(Type==AllegedParent) %>% 
      pull(Sample) %>% 
      unique() %>% 
      droplevels() %>% 
      as.character()
    
    #Use ID to add prob for each AllegedParent
    for(i in AllegedParent.samples){
      WPI[, ncol(WPI) + 1] = subset(data, Sample==i)$GenotypeProbability[match(WPI$AFID,subset(data, Sample==i)$ID)]
      names(WPI)[ncol(WPI)] = paste0("AF.prob.", i)
    }
    
    #calc PI
    WPI$PI= WPI$COEFA * (1/WPI$Afreq)+ WPI$COEFB * (1/WPI$Bfreq) + WPI$COEFAB * (1/(WPI$Afreq + WPI$Bfreq))
    
    #calc total prob for each AllegedParent
    for(i in AllegedParent.samples){
      WPI[, ncol(WPI) + 1] = WPI$M.prob*WPI$C.prob*WPI[[paste0("AF.prob.", i)]]
      names(WPI)[ncol(WPI)] = paste0("Prob.", i)
    }
    
    #PI by prob for each AllegedParent
    for(i in AllegedParent.samples){
      WPI[, ncol(WPI) + 1] = WPI$PI * WPI[[paste0('Prob.', i)]]
      names(WPI)[ncol(WPI)] = paste0("WPI.", i)
    }
    
    #list of AllegedParent column names
    AllegedParent.samples.col = paste0("WPI.",AllegedParent.samples)
    
    ##aggregate WPI score
    #by marker
    #for each male
    WPI.sum = WPI %>% 
      group_by(Marker) %>% 
      summarize_at(AllegedParent.samples.col, sum) %>% 
      pivot_longer(
        cols = starts_with("WPI."),
        names_to = "AllegedParent",
        names_prefix = "WPI.",
        values_to = "WPI.sum"
      )
    
    #Add known parent name
    WPI.sum$KnownParent = data %>% 
      filter(Type==KnownParent) %>% 
      select(Sample) %>% 
      unique() %>% 
      pull()
    #Add offspring name
    WPI.sum$Child = data %>% 
      filter(Type==Child) %>% 
      select(Sample) %>% 
      unique() %>% 
      pull()
    
    return(WPI.sum)
  }
##### Load in sample data ####
load('WOPI.example.data.RData')

##### Phase 1 run across all males and samples use all sites ####
#convert allelefreq.phase1 to dataframe
# allelefreq.phase1.df=as.data.frame(allelefreq.phase1)

#took ~2 min w/ calculate
WPI.sum.total = pbmapply(function(x, y)
  #run paternity test   
  totalWPI(dataframe = 
             data.geno.filter.trans.melt.phase.1.zero  %>% 
             filter(Sex == 'M' | Sample.name == x | Sample.name == y), 
           Sample = "Sample.name",  
           Marker = "variable", 
           Allele = "MajorMinor", 
           GenotypeProbability = "value", 
           Type = "Sex", 
           Child = 'B', 
           KnownParent = 'F', 
           AllegedParent = 'M', 
           allelefreq = 'calculate',
           Afreq = 'Afreq', 
           Bfreq = 'Bfreq', 
           AlleleMarker = 'variable'),
  brood.mother.list$Brood.name,  #brood names as x
  brood.mother.list$Mother.name, #mother names as y
  SIMPLIFY = F
) %>% 
  bind_rows() %>% 
  merge(child.tank) %>% 
  mutate(APTank = str_sub(AllegedParent,
                          start = -3)) %>% 
  mutate(Correct = ifelse(APTank == Tank,
                          1,
                          0
  )) 

##### Phase 1 calculate tail prod ####
##calculate product tails
WPI.sum.total.tail.prod= WPI.sum.total %>% 
  filter(WPI.sum>=1.5 | WPI.sum<=0.5) %>% #set outlier range
  group_by(Child, AllegedParent, Correct) %>% 
  summarize(tails.prod = prod(WPI.sum)) 

#graph tail product
ggplot(WPI.sum.total.tail.prod,
       aes(x = Child,
           y= log(tails.prod),
           shape= as.character(Correct))) +
  geom_point()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle('Triad log WCPI tail products') + 
  labs(shape = "Correct Father")
ggsave('WOPI/Triad/Triad log WCPI tail products.png')

###compare to sample without correct father and duplicates removed
#check distribution
WPI.sum.total.tail.prod.no.father =
  WPI.sum.total.tail.prod %>%
  filter(!str_detect(AllegedParent, '_B.')) %>% 
  filter(Correct!=1) %>%
  group_by(Child) %>%
  summarize(Mean = mean(tails.prod),
            Sd = sd(tails.prod))

#merge
WPI.sum.total.tail.prod.merge = 
  merge(WPI.sum.total.tail.prod %>% 
          filter(!str_detect(AllegedParent, '_B.')),
        WPI.sum.total.tail.prod.no.father)

#zscore
WPI.sum.total.tail.prod.merge = 
  WPI.sum.total.tail.prod.merge %>% 
  mutate(Zscore = (tails.prod - Mean)/Sd,
         Log.zscore = log(Zscore))

#graph child zscore
ggplot(WPI.sum.total.tail.prod.merge,
       aes(x=Child, 
           y= Log.zscore,
           shape=as.character(Correct))) +
  geom_point() +
  geom_hline(yintercept = 20) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle('Triad WOPI per child') + 
  labs(shape = "Correct Father")
ggsave('WOPI/Triad/Triad WOPI per child.png')

#graph boxplot
ggplot(WPI.sum.total.tail.prod.merge,
       aes(x=as.factor(Correct), 
           y= Log.zscore,
           shape=as.character(Correct))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point()+
  geom_hline(yintercept = 20) +
  theme_classic() +
  ggtitle('Triad WOPI results') +
  xlab('Correct Father')+
  labs(shape = "Correct Father")
ggsave('WOPI/Triad/Triad WOPI results.png')

###compare to sample without correct father and duplicates removed
##remove top scoring male (i.e. father)
#check distribution
WPI.sum.total.tail.prod.no.father.no.top =
  WPI.sum.total.tail.prod %>%
  filter(!str_detect(AllegedParent, '_B.')) %>% 
  filter(Correct!=1) %>%
  group_by(Child) %>%
  mutate(Top.score= max(tails.prod),
         Top.NoFather = ifelse(Top.score==tails.prod,
                               1,0)) %>% 
  filter(Top.NoFather==0) %>% 
  summarize(Mean = mean(tails.prod),
            Sd = sd(tails.prod))

#merge
WPI.sum.total.tail.prod.merge.no.top = 
  merge(WPI.sum.total.tail.prod %>% 
          filter(!str_detect(AllegedParent, '_B.')),
        WPI.sum.total.tail.prod.no.father.no.top)

#zscore
WPI.sum.total.tail.prod.merge.no.top = 
  WPI.sum.total.tail.prod.merge.no.top %>% 
  mutate(Zscore = (tails.prod - Mean)/Sd,
         Log.zscore = log(Zscore))

#graph child zscore
ggplot(WPI.sum.total.tail.prod.merge.no.top,
       aes(x=Child, 
           y= Log.zscore,
           shape=as.character(Correct))) +
  geom_point() +
  geom_hline(yintercept = 20) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle('Triad WOPI per child father removed') +
  xlab('Correct Father') +
  labs(shape = "Correct Father")
ggsave('WOPI/Triad/Triad WOPI per child father removed.png')

#graph boxplot
ggplot(WPI.sum.total.tail.prod.merge.no.top,
       aes(x=as.factor(Correct), 
           y= Log.zscore,
           shape=as.character(Correct))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point()+
  geom_hline(yintercept = 20) +
  theme_classic() +
  ggtitle('Triad WOPI results father removed') +
  xlab('Correct Father') +
  labs(shape = "Correct Father")
ggsave('WOPI/Triad/Triad WOPI results father removed.png')

# ##relative tails prod score
# #remove duplicates
# WPI.sum.total.tail.prod = WPI.sum.total.tail.prod %>% 
#   filter(!str_detect(AllegedParent, '_B.')) %>% 
#   group_by(Child) %>% 
#   mutate(Child.tails.prod = sum(tails.prod))
# 
# #calculate relative 
# WPI.sum.total.tail.prod = WPI.sum.total.tail.prod %>% 
#   filter(!str_detect(AllegedParent, '_B.')) %>% 
#   mutate(Relative.tails.prod = tails.prod/Child.tails.prod)
# 
# #graph relative tail product
# ggplot(WPI.sum.total.tail.prod,
#        aes(x = Child,
#            y= Relative.tails.prod,
#            color=Correct)) +
#   geom_point()+
#   theme(axis.text.x = element_text(angle = 90))

##### Phase 2 run across all males and samples use all sites ####
#took ~2 min w/ calculate
WPI.sum.total.phase2 = pblapply(setNames(tank.list,tank.list), function(z)
pbmapply(function(x, y)
  #run paternity test   
  totalWPI(dataframe = 
             data.geno.filter.trans.melt.phase.2.zero  %>% 
             filter(Tank == z) %>% 
             filter(Sex == 'M' | Sample.name == x | Sample.name == y), 
           Sample = "Sample.name",  
           Marker = "variable", 
           Allele = "MajorMinor", 
           GenotypeProbability = "value", 
           Type = "Sex", 
           Child = 'B', 
           KnownParent = 'F', 
           AllegedParent = 'M', 
           allelefreq = 'calculate',
           Afreq = 'Afreq', 
           Bfreq = 'Bfreq', 
           AlleleMarker = 'variable'),
  brood.mother.list.phase2 %>% 
  filter(Tank == z) %>% 
    pull(Brood.name),  #brood names as x
  brood.mother.list.phase2 %>% 
    filter(Tank == z) %>% 
    pull(Mother.name), #mother names as y
  SIMPLIFY = F
) %>% 
  bind_rows()
)


#create long format dataframe with tank
WPI.sum.total.phase2.all = WPI.sum.total.phase2 %>% 
  bind_rows(.id='tank.list') 

##### Phase 2 calculate tail prod ####
WPI.sum.total.phase2.all.prod= WPI.sum.total.phase2.all %>% 
  filter(WPI.sum>=1.5 | WPI.sum<=0.5) %>% #set outlier range
  group_by(Child, AllegedParent,tank.list) %>% 
  summarize(tails.prod = prod(WPI.sum)) 

#graph tail product
ggplot(WPI.sum.total.phase2.all.prod,
       aes(x = Child,
           y= log(tails.prod))) +
  geom_point()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle('Community log WCPI tail products')
ggsave('WOPI/Community/Community log WCPI tail products.png')

###zscore compare to sample without duplicates removed
#check distribution with top males removed
WPI.sum.total.phase2.all.prod.NoTop =
  WPI.sum.total.phase2.all.prod %>%
  filter(!str_detect(AllegedParent, '_B.')) %>% 
  group_by(Child) %>%
  mutate(Max = max(tails.prod),
         Top.male = ifelse(Max==tails.prod,
                           1,
                           0)) %>% 
  filter(Top.male!=1) %>% 
  summarize(Mean = mean(tails.prod),
            Sd = sd(tails.prod))

#merge
WPI.sum.total.phase2.all.prod.merge = 
  left_join(WPI.sum.total.phase2.all.prod %>% 
          filter(!str_detect(AllegedParent, '_B.')),
        WPI.sum.total.phase2.all.prod.NoTop)

#zscore
WPI.sum.total.phase2.all.prod.merge = 
  WPI.sum.total.phase2.all.prod.merge %>% 
  mutate(Zscore = (tails.prod - Mean)/Sd,
         LogZscore = log(Zscore),
         Above = ifelse(LogZscore>=20,
                        'Above 20',
                        'Below 20'))%>%
  group_by(Child) %>% 
  mutate(Max = max(tails.prod),
         Top.male = ifelse(tails.prod==Max,
                           'Top','Other'))


#graph child zscore
ggplot(WPI.sum.total.phase2.all.prod.merge %>% 
         filter(!is.na(LogZscore)),
       aes(x=Child, 
           y= LogZscore,
           shape=reorder(Above, desc(Above)))) +
  geom_point() +
  geom_hline(yintercept = 20) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))+
  ggtitle('Community WOPI per child father') +
  labs(shape = "Threshold")
ggsave('WOPI/Community/Community WOPI per child.png')

#graph boxplot of top scoring male
ggplot(WPI.sum.total.phase2.all.prod.merge %>% 
         filter(!is.na(Above)),
       aes(x=Top.male, 
           y= LogZscore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape=reorder(Above, desc(Above))))+
  geom_hline(yintercept = 20) +
  theme_classic() +
  ggtitle('Community WOPI results father') +
  labs(shape = "Community WPOI top scoring male")
ggsave('WOPI/Community/Community WOPI results top male.png')


# ##relative tails prod score
# #remove duplicates
# WPI.sum.total.phase2.all.prod = WPI.sum.total.phase2.all.prod %>% 
#   filter(!str_detect(AllegedParent, '_B.')) %>% 
#   group_by(Child) %>% 
#   mutate(Child.tails.prod = sum(tails.prod))
# 
# #calculate relative 
# WPI.sum.total.phase2.all.prod = WPI.sum.total.phase2.all.prod %>% 
#   filter(!str_detect(AllegedParent, '_B.')) %>% 
#   mutate(Relative.tails.prod = tails.prod/Child.tails.prod)
# 
# #graph relative tail product
# ggplot(WPI.sum.total.phase2.all.prod,
#        aes(x = Child,
#            y= Relative.tails.prod)) +
#   geom_point()+
#   theme(axis.text.x = element_text(angle = 90))