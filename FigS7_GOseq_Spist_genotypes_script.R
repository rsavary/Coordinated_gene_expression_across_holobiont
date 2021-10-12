######## GOseq enrichement of files for Spist

################# GOseq necessary elements for Spist
#GO enrichement with correction for gene length first for Spist
#source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicFeatures")
#install.packages("reshape")
library("GenomicFeatures")
library("reshape")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("goseq")
library("goseq")

###########
# Creat 4 files, file 1, length of sequences of gene or? sum of exon
###########

###########


###########
########### file 1 length
###########

annotgff3M<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/Spis.genome.annotation.gff3",header=F)
annotgff3M2<-subset(annotgff3M,annotgff3M$V3=="exon")
head(annotgff3M2)

rm(annotgff3M)
dim(annotgff3M2)
annotgff3M3= transform(annotgff3M2, V9 = colsplit(V9, split = ";", names = c('ID','Parent')))
annotgff3M4= transform(annotgff3M3$V9, Names = colsplit(Parent, split = "=", names = c('V9', 'Names')))
LengthM<-as.data.frame(annotgff3M3$V5-annotgff3M3$V4)
LengthM[,2]<-annotgff3M4$Names.Names
colnames(LengthM)<-c("exonlength","gene")
LengthM2<-aggregate(exonlength ~ gene,LengthM,sum)

hist(LengthM2$exonlength, n=100)

######
###### file 2  get GO annotation and creat a file with gene name and one GO per gene name
######
library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/spis_tabulated_annots_csv.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, split = " ", names = c('Name','Namelong')))
annotGOSm3<-as.data.frame(cbind(as.character(annotGOSm2$Query[,1]),as.character(annotGOSm2$GO.terms)))
#annotGOSmbis<-as.data.frame(cbind(as.character(annotGOSm2$Query[,1]),as.character(annotGOSm2$Hit.description), as.character(annotGOSm2$GO.terms)))
dim(annotGOSm3)



########
######## merge file 1 and file 2 =f1f2 with length, GO, gene name old and gene name new
########
#newannotGOSM

dim(annotGOSm3)
head(annotGOSm3)
colnames(annotGOSm3)<-c("gene","GO")
f1f2<-merge(annotGOSm3,LengthM2,by="gene")
dim(f1f2)
##### new file with new name and the length
f1f2new<-f1f2[,c(1,3)]
genesM<-f1f2new$gene
LengthMnew<-f1f2new$exonlength


#### creat file with a single gene name and associated GO with the following loop then same the file and reload it
#write.csv(annotGOSm3,"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/Spist_Annotation_with_gene_name_GO.csv")

annotyGO<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/Spist_one_genename_one_GO.csv")
annotyGO2<-transform(annotyGO$gene.GO, gene.GO=colsplit(annotyGO$gene.GO, split = " ", names = c('Gene', 'GO')))
head(annotyGO2)
annotyGO3<-annotyGO2[,c(2:3)]
colnames(annotyGO3)<-c("Gene","GO")
annotyGO4<-annotyGO3[grepl("GO",annotyGO3$GO),]


#####files
annotyGO4 # annotation ONE GENE PER ONE GO
dim(annotyGO4)

### to adjust i think
LengthMnew # length of each gene 25769
genesM # list of gene 25769
length(genesM)
length(LengthMnew)










########### for list of gene DEGs 
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/02_DEG_list_without_RSS_34-5C_samples/")

# load file with pattern Spist

filesT1 <- list.files(pattern =  "\\Spist")

list.data<-list()

for (i in 1:length(filesT1))
{
  list.data[[i]]<-read.csv(filesT1[i],header=T)
}

list.data_t1<-list.data

#list.data_t1
#filesT1

for(i in 1:length(list.data_t1)){
  if(dim(list.data_t1[[i]])[1]==0){
    next
  }
  DETS2<-list.data_t1[[i]][1] #have a list 
  DETS2<-as.character(DETS2$X)
  inter<-as.data.frame(cbind(as.character(genesM),rep(0,times = 25769)))
  inter$V2<-as.numeric(as.character(inter$V2))
  rownamesDETS<-rownames(inter[inter$V1 %in% DETS2, ])
  for(k in 1:length(rownamesDETS)){
  j<-rownamesDETS[k]
  inter[j,2]<-"1"
  }
  fileone<-inter
  pwf<-nullp(as.numeric(as.vector(fileone$V2)),genesM,bias.data=LengthMnew)
  rownames(pwf)<-genesM
  GO.wall=goseq(pwf,id = genesM,gene2cat =annotyGO4,use_genes_without_cat = F )
  newfile<-as.data.frame(cbind(GO.wall,p.adjust(GO.wall$over_represented_pvalue,method="BH")))
  colnames(newfile)[8]<-"padj"
  newfile2<-subset(newfile,newfile$padj<0.05)
  newfile3<-head(newfile,n=10)
  newfile4<-rbind(newfile2,newfile3)
  name<-gsub(".csv","",filesT1[i])
  dimen<-dim(newfile2)[1]
  write.csv2(newfile4,paste("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/02_Go_enrichment/Go_enrichment_RSS_witout_34-5C_samples/",name,"GO_enrichement_",dimen,"_terms.csv"))
}
  












########### for list of gene uniq or common 
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/03_Save_common_uniq_genes/")

# load file with pattern Spist

filesT1 <- list.files(pattern =  "\\Spist")

list.data<-list()

for (i in 1:length(filesT1))
{
  list.data[[i]]<-read.csv(filesT1[i],header=T)
}

list.data_t1<-list.data

#list.data_t1
#filesT1

for(i in 1:length(list.data_t1)){
  if(dim(list.data_t1[[i]])[1]==0){
    next
  }
  DETS2<-list.data_t1[[i]][2] #have a list 
  DETS2<-as.character(DETS2$X)
  inter<-as.data.frame(cbind(as.character(genesM),rep(0,times = 25769)))
  inter$V2<-as.numeric(as.character(inter$V2))
  rownamesDETS<-rownames(inter[inter$V1 %in% DETS2, ])
  for(k in 1:length(rownamesDETS)){
    j<-rownamesDETS[k]
    inter[j,2]<-"1"
  }
  fileone<-inter
  pwf<-nullp(as.numeric(as.vector(fileone$V2)),genesM,bias.data=LengthMnew)
  rownames(pwf)<-genesM
  GO.wall=goseq(pwf,id = genesM,gene2cat =annotyGO4,use_genes_without_cat = F )
  newfile<-as.data.frame(cbind(GO.wall,p.adjust(GO.wall$over_represented_pvalue,method="BH")))
  colnames(newfile)[8]<-"padj"
  newfile2<-subset(newfile,newfile$padj<0.05)
  newfile3<-head(newfile,n=10)
  newfile4<-rbind(newfile2,newfile3)
  name<-gsub(".csv","",filesT1[i])
  dimen<-dim(newfile2)[1]
  write.csv2(newfile4,paste("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/04_GO_enrichement_Spist/",name,"GO_enrichement_",dimen,"_terms.csv"))
}






for(i in 1:length(list.data_t1)){
  if(dim(list.data_t1[[i]])[1]==0){
    next
  }
  print(i)
  print(filesT1[i])
  print(dim(list.data_t1[[i]])[1])
  
}
  





