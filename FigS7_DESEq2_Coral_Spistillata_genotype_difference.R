# Differential expression between genotypes of S pistillata (Spis)

###Library needed
library("DESeq2")
library("limma")
library("edgeR")
library("ggplot2")
library("tximport")
library("adegenet")
library("pheatmap")
library("ape")





###Kallisto for Stylophora pistillata ####

setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Abundance_Kallisto_RSS_vs_CBASS_without_field/")

filesToProcessSK <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcessSK)<-gsub('_abundance.tsv$','',filesToProcessSK,perl=TRUE)
samplesSK<-read.table('/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/02_Kallisto_abundance_Smicro-Spist/Abundance_Kallisto_RSS_vs_CBASS_without_field/CB-T1-27-G13A_abundance.tsv',h=T)

tx2geneSK<-cbind.data.frame(samplesSK$target_id,gsub('.[0-9].v6.1$','',samplesSK$target_id,perl=TRUE))
colnames(tx2geneSK)<-c('TXNAME','GENEID')
#add tximport function manualy, as the tximport is not made for this version. 
txiSK <- tximport(filesToProcessSK, type="kallisto", tx2gene=tx2geneSK, countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length # countsFromAbundance="lengthScaledTPM" or countsFromAbundance="no"
summary(txiSK)
#######following protocol manual deseq2 page 8
#labels and design
countsSK<-txiSK[[2]]
header<-colnames(countsSK) # header<-gsub("\\d", "", colnames(counts) )


# get temp from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

# get T1/T3 from header
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]


#exp
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]


#genotype
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

gen[c(8,13,23,28,38,48,63)]<-"G15A2" # genotype G15-A2 of the merged colonie identified in Savary et al., 2021 => see SNPs calling Savary et al., supp Fig S3
gen[c(3,18,33,43,53,58,68,73,78)]<-"G15A1" # genotype G15-A1 of the merged colonie identified in Savary et al., 2021 => see SNPs calling Savary et al., supp Fig S3

#
####### 1)  design only short-term heat stress (CBASS) Spis ####
#
# get temp from header

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) #
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(41:80)] # keep only Spis libraries from the short-term heat stress (CBASS)
design2<-design[-c(41:80),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene

ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")
plotPCA(ddsSKtempvst, intgroup="gen")

#Make a counts table that is scaled by the size factors

resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )


temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 5)/8)*100   #### keep only site that have at least 5reads in 88%
}
new2<-subset(new,new[,1]>=80)
dim(new2)

#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Smicro_CBASS_table_rowMeans_per_treatments")


resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
plotPCA(ddsSKtempvst, intgroup="gen")



library("gtools")
library("reshape2")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/spis_tabulated_annots_csv.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, " ", names = c('Name','Namelong')))
colnames(annotGOSm2)[1]



#####results contrast genotypes "G8A","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
#plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G8A_vs_G13A_temp_and_gen_model.csv")

#

#####results contrast genotypes "G8A","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
#plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G8A_vs_G14A_temp_and_gen_model.csv")

#

#####results contrast genotypes "G8A","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
#plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G8A_vs_G9A_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G8A_vs_G15A1_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
#plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A2"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G8A_vs_G15A2_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A2_vs_G15A1_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A2_vs_G9A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A2_vs_G13A_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A2_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A1_vs_G13A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A1","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A1_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G15A1_vs_G9A_temp_and_gen_model.csv")




#####results contrast genotypes "G13A","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G13A_vs_G9A_temp_and_gen_model.csv")



#####results contrast genotypes "G13A","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G13A_vs_14A_temp_and_gen_model.csv")


#####results contrast genotypes "G9A","14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G9A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G9A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_CBASS_G9A_vs_14A_temp_and_gen_model.csv")








#
#### 2)  design only long-term heat stees (RSS) Spis ####
#
# get temp from header

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]


temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # 
#####
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(1:40)] #keep only Spis and long-term heat stress (RSS)
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")
plotPCA(ddsSKtempvst, intgroup="gen")


#### reorder factor and model

resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )

temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)

table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 5)/8)*100   #### keep only site that have at least 10reads in 90%
}
new2<-subset(new,new[,1]>=80)
dim(new2)

#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Smicro_RSS_table_rowMeans_per_treatments")

########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
############
plotPCA(ddsSKtempvst, intgroup="gen")



#####results contrast genotypes "G8A","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G8A_vs_G13A_temp_and_gen_model.csv")

#

#####results contrast genotypes "G8A","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G8A_vs_G14A_temp_and_gen_model.csv")

#

#####results contrast genotypes "G8A","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G8A_vs_G9A_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G8A_vs_G15A1_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A2"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G8A_vs_G15A2_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A2_vs_G15A1_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A2_vs_G9A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A2_vs_G13A_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A2_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A1_vs_G13A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A1","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A1_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G15A1_vs_G9A_temp_and_gen_model.csv")

#####results contrast genotypes "G13A","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G13A_vs_G9A_temp_and_gen_model.csv")



#####results contrast genotypes "G13A","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G13A_vs_14A_temp_and_gen_model.csv")


#####results contrast genotypes "G9A","14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G9A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G9A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/01_DEG_list/Spist_RSS_G9A_vs_14A_temp_and_gen_model.csv")












