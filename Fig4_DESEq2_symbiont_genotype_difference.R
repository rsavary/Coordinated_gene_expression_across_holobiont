# Differential expression between genotypes

#
#Script of RNA-seq comparison Kallisto Stylophora Pistillata 
###Library needed
library("DESeq2")
library("limma")
library("edgeR")
library("ggplot2")
library("tximport")
library("adegenet")
library("pheatmap")
library("ape")




########################################
##########Kallisto for Stylophora pistillata ####
########################################
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
#temp=c(23,23,23,32,34,29,29,34,34,34,27,27,27,27,29,29,29,29,32,32,27,27,27,27,27,29,29,29,32,32,32,34,34) #temp
#temp=c(4,4,4,4,4,4,4,2,2,1,1,1,2,1,1,2,1,1,4,4,1,2,1,1,2,2,1,2,2,2,2,4,4) # lane
#temp=c("T0","T0","T0","T1","T1","T3","T3","T3","T3","T3","T1","T1","T1","T1","T1","T1","T1","T1","T1","T1","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3","T3") # time
#temp<-c(1,1,1,1,34,1,1,34,1,1,1,1,1,1,1,1,1,1,1,1,34,34)
#temp<-c("Field","Field","Field","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","CBASS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS","RSS")

######## modification to the data after first analyze:
# 1) exchange the two field samples that were exchanged.
# 2) remove low counts libraries.




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

gen[c(8,13,23,28,38,48,63)]<-"G15A2"
gen[c(3,18,33,43,53,58,68,73,78)]<-"G15A1"
###### summary matrix of mapped reads to Mesculenta after Kallisto ####

#
#### 1  design only CBASS Smicro ####
#
# get temp from header

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library CTL8,Q11, G4, G12 for CTL and AMF
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,-c(41:80)] #to remove lib
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


#to do DEG which one with model and follow up on the simple model deseq2 and then more difficult model limma
#### reorder factor and model
#ddsCTLAMF$condition<-factor(ddsCTLAMF$condition, levels=c("CTL","AMF"))
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
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/smic_tabulated_annots-NEW.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, split = " ", names = c('Name','Namelong')))
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
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G8A_vs_G13A_temp_and_gen_model.csv")

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
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G8A_vs_G14A_temp_and_gen_model.csv")

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
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G8A_vs_G9A_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G8A_vs_G15A1_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A2"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G8A_vs_G15A2_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A2_vs_G15A1_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A2_vs_G9A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A2_vs_G13A_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A2_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A1_vs_G13A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A1","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A1_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G15A1_vs_G9A_temp_and_gen_model.csv")




#####results contrast genotypes "G13A","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G13A_vs_G9A_temp_and_gen_model.csv")



#####results contrast genotypes "G13A","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G13A_vs_14A_temp_and_gen_model.csv")


#####results contrast genotypes "G9A","14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G9A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G9A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_CBASS_G9A_vs_14A_temp_and_gen_model.csv")











#
#### 1  design only RSS Smicro ####
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
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library CTL8,Q11, G4, G12 for CTL and AMF
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,-c(1:40)] #to remove lib
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")
plotPCA(ddsSKtempvst, intgroup="gen")


#to do DEG which one with model and follow up on the simple model deseq2 and then more difficult model limma
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
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G8A_vs_G13A_temp_and_gen_model.csv")

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
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G8A_vs_G14A_temp_and_gen_model.csv")

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
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G8A_vs_G9A_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G8A_vs_G15A1_temp_and_gen_model.csv")


#####results contrast genotypes "G8A","G15A2" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G8A","G15A2"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G8A","G15A2"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G8A_vs_G15A2_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G15A1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G15A1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G15A1"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A2_vs_G15A1_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A2_vs_G9A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A2","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A2_vs_G13A_temp_and_gen_model.csv")



#####results contrast genotypes "G15A2","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A2","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A2","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A2_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G13A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G13A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G13A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A1_vs_G13A_temp_and_gen_model.csv")

#####results contrast genotypes "G15A1","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A1_vs_G14A_temp_and_gen_model.csv")


#####results contrast genotypes "G15A1","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G15A1","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G15A1","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G15A1_vs_G9A_temp_and_gen_model.csv")

#####results contrast genotypes "G13A","G9A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G9A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G9A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G13A_vs_G9A_temp_and_gen_model.csv")



#####results contrast genotypes "G13A","G14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G13A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G13A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G13A_vs_14A_temp_and_gen_model.csv")


#####results contrast genotypes "G9A","14A" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('gen',"G9A","G14A"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Smic18218", intgroup="gen",main="Nitrate reductase")
plotCounts(ddsSKtemp, gene="Smic4625", intgroup="gen",main="High affinity nitrate transporter")
plotCounts(ddsSKtemp, gene="Smic42197", intgroup="gen",main="High affinity nitrate transporter")
which.min(resutlsSKtemp$padj)

### save results 0.05 sign "G9A","G14A"

dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/01_DEG_list/Smicro_RSS_G9A_vs_14A_temp_and_gen_model.csv")

















#
#### 2  design only CBASS spistillata ####
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
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(41:80)] #to remove lib
design2<-design[-c(41:80),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")



#to do DEG which one with model and follow up on the simple model deseq2 and then more difficult model limma
#### reorder factor and model
#ddsCTLAMF$condition<-factor(ddsCTLAMF$condition, levels=c("CTL","AMF"))
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
  new[i,1]<-(sum(table[i,] > 1)/8)*100   #### keep only site that have at least 10reads in 90%
}
new2<-subset(new,new[,1]>=80)
dim(new2)

resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Spist_CBASS_table_rowMeans_per_treatments")

########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
counts(resSKtemp)
############


library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/spis_tabulated_annots_csv.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, split = " ", names = c('Name','Namelong')))
annotGOSm3<-as.data.frame(cbind(as.character(annotGOSm2$Query[,1]),as.character(annotGOSm2$GO.terms)))


par(mar=c(3,3,3,3))
#####results contrast "34T1","27T1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis4494", intgroup="temp",main="Spist CBASS Spis4494 heat shock protein Hsp-16.2")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_34T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "32T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_32T1_vs_27T1_temp_and_gen_model.csv")

#####results contrast "29T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_29T1_vs_27T1_temp_and_gen_model.csv")





#####results contrast "34T3","27T3" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T3","27T3"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis1000.t1", intgroup="temp",main="Spist CBASS Spis1000.t1 Common gene in recovery for CBASS and RSS for 29,32,34")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_34T3_vs_27T3_temp_and_gen_model.csv")


#####results contrast "32T3","27T3" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T3","27T3"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_32T3_vs_27T3_temp_and_gen_model.csv")

#####results contrast "29T3","27T3" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T3","27T3"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_29T3_vs_27T3_temp_and_gen_model.csv")




#######
####### Normal vs recovery
#######
#####results contrast "34T3","34T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T3","34T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_34T3_vs_34T1_temp_and_gen_model.csv")

#####results contrast "32T3","32T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T3","32T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis18245", intgroup="temp",main="Oxidative stress-responsive serine-rich protein 1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T3 vs 32T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_32T3_vs_32T1_temp_and_gen_model.csv")

#####results contrast "29T3","29T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T3","29T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 29T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_29T3_vs_29T1_temp_and_gen_model.csv")

#####results contrast "27T3","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"27T3","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 27T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_CBASS_27T3_vs_27T1_temp_and_gen_model.csv")


#
#### 2  design only RSS spistillata ####
#
# get temp from header
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
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) # for CTL and AMF
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,-c(1:40)] #to remove lib
design2<-design[-c(1:40),]
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
#ddsSKtemprlog<-rlog(ddsSKtemp)

#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")


#to do DEG which one with model and follow up on the simple model deseq2 and then more difficult model limma
#### reorder factor and model
#ddsCTLAMF$condition<-factor(ddsCTLAMF$condition, levels=c("CTL","AMF"))
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
  new[i,1]<-(sum(table[i,] > 1)/8)*100   #### keep only site that have at least 1reads
}
new2<-subset(new,new[,1]>=80) ###### in 100 of treatments
dim(new2)

#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Spist_RSS_table_rowMeans_per_treatments")

resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)

#####results contrast "34T1","27T1" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis1000.t1", intgroup="temp",main="Spist RSS Spis1000.t1 Common gene in recovery for CBASS and RSS for 29,32,34")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_34T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "32T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_32T1_vs_27T1_temp_and_gen_model.csv")


#####results contrast "29T1","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T1","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_29T1_vs_27T1_temp_and_gen_model.csv")





#####results contrast "34T3","27T3" ########
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T3","27T3"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="Spis1000.t1", intgroup="temp",main="Spist RSS Spis1000.t1 Common gene in recovery for CBASS and RSS for 29,32,34")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 34T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_34T3_vs_27T3_temp_and_gen_model.csv")


#####results contrast "32T3","27T3" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T3","27T3"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_32T3_vs_27T3_temp_and_gen_model.csv")


#####results contrast "29T3","27T3" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T3","27T3"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_29T3_vs_27T3_temp_and_gen_model.csv")

#######
####### Stress vs recovery
#######
#####results contrast "34T3","34T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"34T3","34T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_34T3_vs_34T1_temp_and_gen_model.csv")

#####results contrast "32T3","32T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"32T3","32T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 32T1 vs 32T3
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_32T3_vs_32T1_temp_and_gen_model.csv")

#####results contrast "29T3","29T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"29T3","29T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 29T1 vs 29T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_29T3_vs_29T1_temp_and_gen_model.csv")

#####results contrast "27T3","27T1" ####
resutlsSKtemp<-results(resSKtemp2,contrast = c('temp',"27T3","27T1"))
resutlsSKtempOrdered<-resutlsSKtemp[order(resutlsSKtemp$padj),]
resutlsSKtempOrderedp05<-subset(resutlsSKtempOrdered,resutlsSKtempOrdered$padj<0.05) #### 2958 gene change expression
dim(resutlsSKtempOrderedp05)
summary(resutlsSKtemp)
sum(resutlsSKtempOrderedp05$padj < 0.05, na.rm=TRUE)
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange >0))
dim(subset(resutlsSKtempOrderedp05, resutlsSKtempOrderedp05$padj < 0.05 & resutlsSKtempOrderedp05$log2FoldChange <0))
plotCounts(ddsSKtemp, gene="11195", intgroup="condition",main="3-oxoacyl-[acyl-carrier-protein] reductase FabG OS=Vibrio cholerae serotype O1")
which.min(resutlsSKtemp$padj)### save results 0.05 sign from 29T1 vs 27T1
### save results 0.05 sign from 27T1 vs 27T1
dim((annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]))
list1annot<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(resutlsSKtempOrderedp05),][,c(1,4)]
rownames(list1annot)<-list1annot$Query$Name
list1annot<-list1annot[order(rownames(list1annot)),]
resutlsSKtempOrderedp052<-resutlsSKtempOrderedp05[order(rownames(resutlsSKtempOrderedp05)),]
resutlsSKtempOrderedp05<-cbind(list1annot,resutlsSKtempOrderedp052)
resutlsSKtempOrderedp05<-resutlsSKtempOrderedp05[order(resutlsSKtempOrderedp05$padj),]
write.csv(as.data.frame(resutlsSKtempOrderedp05),"/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/Spist_RSS_recovery_27T3_vs_27T1_temp_and_gen_model.csv")


























