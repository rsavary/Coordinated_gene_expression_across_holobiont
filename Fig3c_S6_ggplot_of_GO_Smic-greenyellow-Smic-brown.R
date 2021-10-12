### Fig 3c, Barplot of GO enrichement (BP) of the Smic-greenyellow module and for all modules.
### Fig S6, Barplot of GO enrichement (BP) of the Smic-brown module


#install.packages("tidyverse")
#install.packages("egg")
library(tidyverse)
library(egg)


###### Smic
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/02_Go_enrichment_of_module/")

# Smic

filesT1 <- list.files(pattern =  "\\Smicro")

list.data<-list()

for (i in 1:length(filesT1))
{
  list.data[[i]]<-read.csv2(filesT1[i],header=T)
}


list.data_t1<-list.data


##### only the 20th first GO all plot all category of GO (BP, MF, CC)
for(i in 1:length(list.data_t1)){
  df1<-as.data.frame(list.data_t1[i])
  df1<-df1[-c((dim(df1)[1]-9):dim(df1)[1]),]
  df2<-subset(df1,df1$padj<=0.05)
  if(length(df2$ontology)>0){
    a<-paste(df2$term,paste(df2$numDEInCat,"/",df2$numInCat, sep=""))
    b<- -log10(df2$padj)
    b2<-as.numeric(unique(b))
    b<-as.numeric(gsub("Inf",(sort(b2,decreasing = T)[2]*1.5),b))
    d<-data.frame(a=factor(a),b=as.numeric(b),onto=df2$ontology)
    d<-na.omit(d)
    d<-d[order(d$onto,-d$b),]
    #d<-d[,1:2]
    d$a <- factor(d$a, levels = d$a[1:dim(d)[1]])
    p<-ggplot(data=d, aes(x=a, y=b, fill=onto, group = onto)) + geom_bar(stat="identity", fill="green1") + 
      xlab("") + ylab("-log10 (P value)") + ggtitle(filesT1[i]) +
      facet_wrap(~ onto)
    p + coord_flip()
    setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/03_Plot_GO/")
    ggsave(paste(filesT1[i],".pdf",sep=""))#limitsize = FALSE)# height = dim(d)[1]/3
  }
}


### Fig 3c #####
#### Smic greenyellow GO enrichement barplot
i=17
df1<-as.data.frame(list.data_t1[i])
df1<-df1[-c((dim(df1)[1]-9):dim(df1)[1]),]
df1<-subset(df1, df1$ontology=="BP")
df2<-subset(df1,df1$padj<=0.05)
if(length(df2$ontology)>0){
  a<-paste(df2$term,paste(df2$numDEInCat,"/",df2$numInCat, sep=""))
  b<- -log10(df2$padj)
  b2<-as.numeric(unique(b))
  b<-as.numeric(gsub("Inf",(sort(b2,decreasing = T)[2]*1.5),b))
  d<-data.frame(a=factor(a),b=as.numeric(b),onto=df2$ontology)
  d<-na.omit(d)
  d<-d[order(d$onto,-d$b),]
  #d<-d[,1:2]
  d$a <- factor(d$a, levels = d$a[1:dim(d)[1]])
  p<-ggplot(data=d, aes(x=a, y=b, fill=onto, group = onto)) + geom_bar(stat="identity",width = 0.8, fill="grey40") + 
    xlab("") + ylab("-log10 (P value)") + ggtitle(paste("A",filesT1[i],sep=" ")) +
    facet_wrap(~ onto)
  p1<-p + coord_flip()
  #setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/07_GOplot_Smicro/")
  #ggsave(paste(filesT1[i],".pdf",sep=""))#limitsize = FALSE)# height = dim(d)[1]/3
}


### Fig S6 #####
#### Smic brown GO enrichement barplot
i=14
df1<-as.data.frame(list.data_t1[i])
df1<-df1[-c((dim(df1)[1]-9):dim(df1)[1]),]
#df1<-subset(df1, df1$ontology=="BP")
df2<-subset(df1,df1$padj<=0.05)
if(length(df2$ontology)>0){
  a<-paste(df2$term,paste(df2$numDEInCat,"/",df2$numInCat, sep=""))
  b<- -log10(df2$padj)
  b2<-as.numeric(unique(b))
  b<-as.numeric(gsub("Inf",(sort(b2,decreasing = T)[2]*1.5),b))
  d<-data.frame(a=factor(a),b=as.numeric(b),onto=df2$ontology)
  d<-na.omit(d)
  d<-d[order(d$onto,-d$b),]
  #d<-d[,1:2]
  d$a <- factor(d$a, levels = d$a[1:dim(d)[1]])
  p<-ggplot(data=d, aes(x=a, y=b, fill=onto, group = onto)) + geom_bar(stat="identity",width = 0.8, fill="grey40") + 
    xlab("") + ylab("-log10 (P value)") + ggtitle(paste("A",filesT1[i],sep=" ")) +
    facet_wrap(~ onto)
  p1<-p + coord_flip()
  #setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/07_GOplot_Smicro/")
  #ggsave(paste(filesT1[i],".pdf",sep=""))#limitsize = FALSE)# height = dim(d)[1]/3
}



####### Spis barplot of enrichment of all modules


setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Spist_modulesize30/02_Go_enrichment_of_module/")

# Smicro

filesT1 <- list.files(pattern =  "\\Spist")

list.data<-list()

for (i in 1:length(filesT1))
{
  list.data[[i]]<-read.csv2(filesT1[i],header=T)
}


list.data_t1<-list.data


##### only the 20th first GO
for(i in 1:length(list.data_t1)){
  df1<-as.data.frame(list.data_t1[i])
  df1<-df1[-c((dim(df1)[1]-9):dim(df1)[1]),]
  df2<-subset(df1,df1$padj<=0.05)
  if(length(df2$ontology)>0){
    a<-paste(df2$term,paste(df2$numDEInCat,"/",df2$numInCat, sep=""))
    b<- -log10(df2$padj)
    b2<-as.numeric(unique(b))
    b<-as.numeric(gsub("Inf",(sort(b2,decreasing = T)[2]*1.5),b))
    d<-data.frame(a=factor(a),b=as.numeric(b),onto=df2$ontology)
    d<-na.omit(d)
    d<-d[order(d$onto,-d$b),]
    #d<-d[,1:2]
    d$a <- factor(d$a, levels = d$a[1:dim(d)[1]])
    p<-ggplot(data=d, aes(x=a, y=b, fill=onto, group = onto)) + geom_bar(stat="identity", fill="brown1") + 
      xlab("") + ylab("-log10 (P value)") + ggtitle(filesT1[i]) +
      facet_wrap(~ onto)
    p + coord_flip()
    setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Spist_modulesize30/03_GOplot/")
    ggsave(paste(filesT1[i],".pdf",sep=""),limitsize = FALSE)# height = dim(d)[1]/3
  }
}

