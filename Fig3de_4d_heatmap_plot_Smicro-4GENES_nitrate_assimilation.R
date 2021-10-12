##### Fig. 3d and e heatmap of Smic plot four genes Nitrate assimilation as Fig 4d?


setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/")


####################


options(stringsAsFactors = FALSE);
library("WGCNA")
library(DESeq2)
library(impute) 
library(dynamicTreeCut)
library(qvalue) 
library(flashClust)
library(Hmisc)
library("pheatmap")
library("DESeq2")
library("limma")
library("edgeR")
library("ggplot2")
library("tximport")
library("adegenet")
library("pheatmap")
library("reshape2")


########################################
########## get Kallisto abundances ####
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


#
#### 2  design only CBASS Smicro ####
#
# get temp from header

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) #
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[1:49109,] #to remove lib
design2<-design
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
plotPCA(ddsSKtempvst, intgroup="gen")
plotPCA(ddsSKtempvst, intgroup="temp")
a<-plotPCA(ddsSKtempvst, intgroup="gen")
#Make a counts table that is scaled by the size factors

resSKtemp<-DESeq(ddsSKtemp,modelMatrixType="expanded", betaPrior=T )


temp = t(sizeFactors(resSKtemp))
sizematrix<-matrix(data=temp, nrow=nrow(counts(ddsSKtemp)), ncol=ncol(temp), byrow=TRUE)
scaledcounts = counts(ddsSKtemp)/sizematrix
head(scaledcounts)
head(sizematrix)
dim(sizematrix)
table<-cbind(rowMeans(scaledcounts[,1:5]),
             rowMeans(scaledcounts[,6:10]),
             rowMeans(scaledcounts[,11:15]),
             rowMeans(scaledcounts[,16:20]),
             rowMeans(scaledcounts[,21:25]),
             rowMeans(scaledcounts[,26:30]),
             rowMeans(scaledcounts[,31:35]),
             rowMeans(scaledcounts[,36:40]),
             rowMeans(scaledcounts[,41:45]),
             rowMeans(scaledcounts[,46:50]),
             rowMeans(scaledcounts[,51:55]),
             rowMeans(scaledcounts[,56:60]),
             rowMeans(scaledcounts[,61:65]),
             rowMeans(scaledcounts[,66:70]),
             rowMeans(scaledcounts[,71:75]),
             rowMeans(scaledcounts[,76:80]))
colnames(table)<-c("CB-T1-27","CB-T1-29","CB-T1-32","CB-T1-34","CB-T2-27","CB-T2-29","CB-T2-32","CB-T2-34",
                   "RSS-T1-27","RSS-T1-29","RSS-T1-32","RSS-T1-34","RSS-T2-27","RSS-T2-29","RSS-T2-32","RSS-T2-34")

new<-scaledcounts[,1:2]
for(i in 1:dim(table)[1]){
  new[i,1]<-(sum(table[i,] > 5)/16)*100   #### keep only site that have at least 5reads in 88%
}
new2<-subset(new,new[,1]>=80)
dim(new2)

#write.csv(table, "/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/31_Last_analysis_of_DEGs_with_filtering_following_Dan_discussion/02_Save_DEGs_list/05_New_filtering_following_Dan_advices_table/Smicro_CBASS_table_rowMeans_per_treatments")


resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]


####### four genes: 
########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
ddsSKtempvst<-varianceStabilizingTransformation(resSKtemp2)#PCA
plotPCA(ddsSKtempvst, intgroup="gen")
plotPCA(ddsSKtempvst, intgroup="temp")


vstCBASScoralT1<-assay(ddsSKtempvst) # all 80
#### remove gene with zero variance

vstCBASScoral3<-vstCBASScoralT1[apply(vstCBASScoralT1,1,var)>0,] # filtered for low variance

vstCBASScoral4<-t(vstCBASScoral3)


###### removal of lib that are too different (7)
gsg = goodSamplesGenes(vstCBASScoral4, verbose = 3);
sampleTree = hclust(dist(vstCBASScoral4), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 450, col = "red");
clust = cutreeStatic(sampleTree, cutHeight = 450, minSize = 10)
table(clust)
keepSamples = (clust==1)
vstCBASScoral5 = vstCBASScoral4[keepSamples, ]
nGenes = ncol(vstCBASScoral5)
nSamples = nrow(vstCBASScoral5)
vstCBASScoral5

dim(vstCBASScoral5)
rownames(vstCBASScoral5)

vstCBASScoral6<-vstCBASScoral5





########## ALL gene in the nitrate module
data<-t(vstCBASScoral6)

nitmoduley<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/01_Module/Smicro_Allgreenyellow.csv")

data2<-data[rownames(data) %in% nitmoduley$rownames.geneModuleMembershipM.,]

heatmap(data2,scale="column")
pheatmap(data2)


#### Fig 3d ####
pheatmap(data2,scale="row")


#### Fig 3e #####
#####
#### juste the 10genes of nitrate GO

annotyGO<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/Python_first_script/Smicro_genename_Smicro_GO.csv")
annotyGO2<-transform(annotyGO$gene.V2, gene.V2=colsplit(annotyGO$gene.V2, " ", names = c('Gene', 'GO')))
head(annotyGO2)
annotyGO3<-annotyGO2[,c(2:3)]
colnames(annotyGO3)<-c("Gene","GO")
annotyGO4<-annotyGO3[grepl("GO",annotyGO3$GO),]


library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME/02_Symbiodinium/01_CLADE_A_Smicroadriaticum/smic_tabulated_annots-NEW.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query," ", names = c('Name','Namelong')))

nitrate<-annotyGO4[annotyGO4$GO %in% "GO:0042128",]
data3<-data2[rownames(data2) %in% nitrate$Gene,]
nitratename<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(data3),]
nitratename2<-nitratename[order(nitratename$Query$Name),]
data4<-data3[order(rownames(data3)),]
rownames(data4)<-paste(rownames(data4),nitratename2$Hit.description,sep=" ")


par(mar=c(5,3,3,10))
#heatmap(data4,scale="row")
pheatmap(data3, scale="row")



#### Fig 4d? ####

###### Plot four genes always differnetially express between genotypes
data2<-data[rownames(data) %in% c("Smic41892", "Smic33528", "Smic42197", "Smic16438"),]
nitratename<-annotGOSm2[annotGOSm2$Query$Name %in% rownames(data2),]
nitratename2<-nitratename[order(nitratename$Query$Name),]
data4<-data2[order(rownames(data2)),]
rownames(data4)<-paste(rownames(data4),nitratename2$Hit.description,sep=" ")

##### plot expression:
colna<-colnames(data3)
colna2<-gsub("-5",".5",colna)
a<-strsplit(colna2,"-")
df <- data.frame(matrix(NA,nrow=73,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
gen<-df[,1]

#In CBASS G15A_1 and G15A_2 and in RSS G15A_1 and G15A_2
#### 7 => G15_2 following PCA scores
CB-T1-29-5-G15A
CB-T1-32-G15A
CB-T3-27-G15A
CB-T3-29-5-G15A
CB-T3-34-G15A
RSS-T1-29-5-G15A
RSS-T3-27-G15A

colnames(data3)[c(8,13,23,28,38,48,61)]
gen[c(8,13,23,28,38,48,61)]<-"G15A_2"

gen

gen<-gsub("G8A","06G8A",gen)
gen<-gsub("G15A_2","05G15A_2",gen)
gen<-gsub("\\<G15A\\>","04G15A_1",gen)
gen<-gsub("G13A","03G13A",gen)
gen<-gsub("G14A","02G14A",gen)
gen<-gsub("G9A","01G9A",gen)
##### boxplot of genes : Nitrate reductase Smic41892, Smic33528, Nitrate transporter Smic42197, Smic16438
par(mfrow=c(1,4),mar=c(5,5,2,4))
boxplot(data4[1,]~gen,las=2, main=rownames(data4)[1],cex.main=0.5,ylab="vst normalized expression")
boxplot(data4[4,]~gen,las=2, main=rownames(data4)[4],cex.main=0.5,ylab="vst normalized expression")
boxplot(data4[2,]~gen,las=2, main=rownames(data4)[2],cex.main=0.5,ylab="vst normalized expression")
boxplot(data4[3,]~gen,las=2, main=rownames(data4)[3],cex.main=0.5,ylab="vst normalized expression")

