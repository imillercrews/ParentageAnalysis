##WOPI function
##compare prod of paternity index outlier tails 
##Isaac Miller-Crews, imillercrews@utexas.edu
#https://github.com/imillercrews/ParentageAnalysis
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
