###### Enrichment of DEG list between genotypes Smicro

######################
########## Go enrichment of these modules ####
######################


library("GenomicFeatures")
library("reshape")
library("goseq")

###########
# Creat 4 files, file 1, length of sequences of gene or? sum of exon
###########

###########
########### file 1 length
###########

annotgff3M<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/Smic.genome.annotation.gff3",header=F)
annotgff3M2<-subset(annotgff3M,annotgff3M$V3=="exon")


rm(annotgff3M)
dim(annotgff3M2)
#annotgff3M2<-head(annotgff3M2,n=300)

annotgff3M3= transform(annotgff3M2, V9 = colsplit(V9, split = ";", names = c('ID','Parent')))
annotgff3M4= transform(annotgff3M3$V9, Names = colsplit(Parent, split = "=", names = c('V9', 'Names')))
LengthM<-as.data.frame(annotgff3M3$V5-annotgff3M3$V4)
LengthM[,2]<-annotgff3M4$Names.Names
colnames(LengthM)<-c("exonlength","gene")
LengthM2<-aggregate(exonlength ~ gene,LengthM,sum)


######
###### file 2  get GO annotation and creat a file with gene name and one GO per gene name
######
library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/smic_tabulated_annots-NEW.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, split = " ", names = c('Name','Namelong')))
annotGOSm3<-as.data.frame(cbind(as.character(annotGOSm2$Query[,1]),as.character(annotGOSm2$GO.terms)))
#annotGOSmbis<-as.data.frame(cbind(as.character(annotGOSm2$Query[,1]),as.character(annotGOSm2$Hit.description), as.character(annotGOSm2$GO.terms)))
dim(annotGOSm3)


########
######## merge file 1 and file 2 =f1f2 with length, GO, gene name old and gene name new
########
#
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

#######
### This work but python way waaaaay faster so download file from my python script /Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/Python_first_script/Python_script_creat_annotation_file_one_gene_one_GO.py New script with only Smic
######

annotyGO<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/Python_first_script/Smicro_genename_Smicro_GO.csv")
annotyGO2<-transform(annotyGO$gene.V2, gene.V2=colsplit(annotyGO$gene.V2, split = " ", names = c('Gene', 'GO')))
head(annotyGO2)
annotyGO3<-annotyGO2[,c(2:3)]
colnames(annotyGO3)<-c("Gene","GO")
annotyGO4<-annotyGO3[grepl("GO",annotyGO3$GO),]

####### files
dim(annotyGO4)
length(LengthMnew)
length(genesM)

########### for list of gene DEGs 
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/")

# load file with pattern Spist

filesT1 <- list.files(pattern =  "\\Smicro")

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
  DETS2<-as.character(DETS2$Query.Name)
  inter<-as.data.frame(cbind(as.character(genesM),rep(0,times = 49109)))
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
  write.csv2(newfile4,paste("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/02_Go_enrichment/",name,"GO_enrichement_",dimen,"_terms.csv"))
}




###### list all the GO term in all genotypes comparisons of genotypes and do a barplot

setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/02_Go_enrichment/Go_enrichment_RSS_all_samples/")

filesT1 <- list.files(pattern =  "\\Smicro")

#a<-read.csv2("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/02_Go_enrichment/ Smicro_CBASS_G13A_vs_14A_temp_and_gen_model GO_enrichement_ 4 _terms.csv",header=T)
b<-0
for (i in 1:length(filesT1))
{
  a<-read.csv2(filesT1[i],header=T)
  a<-head(a,-10)
  b<-c(b,as.character(a$term))
}
sort(table(b))

#barplot(sort(table(b)),las=2)
par(mar=c(17,3,3,3),mfrow=c(1,1))
#barplot(tail(sort(table(b)),50),las=2)
barplot(head(sort(table(b),decreasing = TRUE),50),las=2,ylim=c(0,35),ylab="nb")
abline(h=30)
text(45,31,"Nb genotype comparisons for GO enrichment")
abline(h=22)
text(42,23,"Nb of comparison with at least one GO term enriched")

par(mar=c(3,25,3,3))
barplot(tail(sort(table(b)),50),las=2,xlim=c(0,35),ylab="nb",horiz = T)
abline(v =30)
text(31,40,"Nb genotype comparisons for GO enrichment", srt=270)
abline(v=20)
text(21,40,"Nb of comparison with at least one GO term enriched",srt=270)


#####Only biological process
#Fig 4

b<-0
for (i in 1:length(filesT1))
{
  a<-read.csv2(filesT1[i],header=T)
  a<-head(a,-10)
  a<-subset(a,a$ontology=="BP")
  b<-c(b,as.character(a$term))
}
sort(table(b))

#barplot(sort(table(b)),las=2)
par(mar=c(17,3,3,3))
#barplot(tail(sort(table(b)),50),las=2)
barplot(head(sort(table(b),decreasing = TRUE),50),las=2,ylim=c(0,35),ylab="nb")
abline(h=30)
text(45,31,"Nb genotype comparisons for GO enrichment")
abline(h=20)
text(42,21,"Nb of comparison with at least one GO term enriched")


