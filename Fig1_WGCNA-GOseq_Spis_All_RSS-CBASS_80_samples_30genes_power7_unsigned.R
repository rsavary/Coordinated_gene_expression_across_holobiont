##### WGCNA for all Spis samples


setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/")


####################
########################################
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

#WGCNA following tuto 2, using only Spis

#######################################
#design only Spis
#######################################
#######################################

########################################
##########Load Kallisto for Stylophora pistillata ####
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
####  design only Spis ####
#
# get temp from header

temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime)
# test inculding genotype
temptime<-paste(temp,time,sep="")
design<-data.frame(temp=temptime, gen=gen)

############# make the DESeq data
ddsTxiSKtemp<-DESeqDataSetFromTximport(txiSK, colData=design, design=formula(~temp+gen)) 
#####get rid of library of the RSS
ddsTxiSKtemp<-counts(ddsTxiSKtemp)[49110:74878,] # keep only Spis abundances
design2<-design
design2<-data.frame(temp=factor(design2$temp),gen=factor(design2$gen))
ddsTxiSKtemp2<-DESeqDataSetFromMatrix(ddsTxiSKtemp, colData = design2, design = formula(~temp+gen))
ddsSKtemp<-ddsTxiSKtemp2[rowSums(counts(ddsTxiSKtemp2)) >1,] #trim no express gene
ddsSKtempvst<-varianceStabilizingTransformation(ddsSKtemp)#PCA
plotPCA(ddsSKtempvst, intgroup="gen")
plotPCA(ddsSKtempvst, intgroup="temp")
#Make a counts table that is scaled by the size factors

#### reorder factor and model
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
dim(new2) # 16756 genes kept


resSKtemp2<-resSKtemp[rownames(counts(resSKtemp)) %in% rownames(new2), ]
########## log fold change shrinkage to remove gene and have a better visualization
resultsNames(ddsSKtemp)
ddsSKtempvst<-varianceStabilizingTransformation(resSKtemp2)#PCA
plotPCA(ddsSKtempvst, intgroup="gen")



#plotPCA(ddsSKtemprlog)
plotPCA(ddsSKtempvst, intgroup="temp")

#ddsSKtemprlog<-rlog(ddsSKtemp)
vstCBASScoralT1<-assay(ddsSKtempvst) # all 80
#### remove gene with zero variance

##### from here we are using only T1
#####
#####building stress modules
#####
#####

vstCBASScoral3<-vstCBASScoralT1[apply(vstCBASScoralT1,1,var)>0,] # filtered for low variance

vstCBASScoral4<-t(vstCBASScoral3)


###### removal of lib
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

Spist_all_vst_vstCBASScoral6<-vstCBASScoral6
save(Spist_all_vst_vstCBASScoral6,file = "/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Spist_all_vst_vstCBASScoral6.RData")

####2a automatic network construction and module detection for Spis
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(vstCBASScoral6, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])

##### We use 7, as it is the first soft threshold that reach 0.8 of scalde free topology.

#WGCNA has a function cor that enter into conflit with cor function of R then you have assign the cor function of 
#WGCNA to the name like this cor <- WGCNA::cor

#### set wd to save session on the right place:
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/")

cor <- WGCNA::cor

###building network with one function
SpistCBASS = blockwiseModules(vstCBASScoral6, power =7,
                              TOMType = "unsigned", minModuleSize =30,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = TRUE,
                              saveTOMFileBase = "Spist ALL",
                              verbose = 3)

cor<-stats::cor # reput function cor 
#####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

table(SpistCBASS$colors)
mean(table(SpistCBASS$colors))
median(table(SpistCBASS$colors))
SpistCBASS$dendrograms[[1]];

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColorsSpistCBASS = labels2colors(SpistCBASS$colors)
table(mergedColorsSpistCBASS)
dim(table(mergedColorsSpistCBASS))
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(SpistCBASS$dendrograms[[1]], mergedColorsSpistCBASS[SpistCBASS$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabelsSpist = SpistCBASS$colors
moduleColorsSpist = labels2colors(SpistCBASS$colors)
MEsSpist = SpistCBASS$MEs;
geneTreeSpist = SpistCBASS$dendrograms[[1]];

#save
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/")
save(MEsSpist, moduleLabelsSpist, moduleColorsSpist, geneTreeSpist,
     file = "Spist_ALL-30genes_p7_unsigned_networkConstruction-auto.RData")
View(table(moduleColorsSpist))
dim(table(moduleColorsSpist))


write.csv(MEsSpist,"MEs_Spist_ALL_30genes_7.csv")



###################
############## Module eigengenes and phenotypic traits
####################

#############LOAD PHENOTYPIC TRAITS
#############
##### load phenotypic data: test one with PAM data make mean data for every samples. after that make a single data frame with the 85 samples and the phenoptypic value associated.
CBT1<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/03_Traits/Physiological_traits/All_physiology_all_samples_CBASS_vs_RSS.csv", header = T)
rownames(CBT1)<-CBT1$Name
CBT1<-CBT1[,-1]
View(CBT1)
CBT1<-CBT1[-c(58,76),] # remove phenotypic trait for the gene expression samples removed.
dim(CBT1)


ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)



modNamesM = substring(names(MEsMT1), 3)
geneModuleMembershipM = cor(CBT1, MEsMT1, use = "p")
MMPvalueM = corPvalueStudent(as.matrix(geneModuleMembershipM), nSamples)

#### make correction for number of test
MMPvalueMadj<-p.adjust(MMPvalueM, method ="bonferroni") # "bonferroni" or "BH"=Benjamini & Hochberg 
#MMPvalueMadj<-MMPvalueM
#Since  we  have  a  moderately  large  number  of  modules  and  traits,  a  suitable  graphical  representation  will  help  inreading the table.  We color code each association by the correlation value:
sizeGrWindow(10,6)# Will display correlations and their p-values

mat<-signif(MMPvalueMadj, 1)
mat[mat>0.05]<-"-"
mat2<-signif(geneModuleMembershipM, 2)
mat2[which(mat =="-")] <- ""


#textMatrix = paste(mat2, "\n(",mat, " adj)", sep = "");
#textMatrix[textMatrix=="\n(- adj)"]<-"-"
mat2[mat2==""]<-"-"
textMatrix = mat2


table(moduleColorsSpist)
namesplot<-as.data.frame(cbind(colnames(MEsMT1),gsub("ME","",colnames(MEsMT1)),c(1:length(colnames(MEsMT1)))))
namesplot<-namesplot[order(namesplot$V2),]
namesplot2<-as.data.frame(cbind(namesplot,as.data.frame(table(moduleColorsSpist))))
namesplot3<-namesplot2[order(as.numeric(namesplot2$V3)),] ## sort by numeric
namesplot4<-paste(namesplot3$V1,namesplot3$Freq,sep=" ")

namescolumn<-gsub("\\(100\\)","",row.names(geneModuleMembershipM))
namescolumn<-gsub("_"," ",namescolumn)

dim(textMatrix) = dim(geneModuleMembershipM)

# Display the correlation values within a heatmap plot
pdf(paste("Module-physio_Spist_All_p7__30genes_38modules_heatmap1",".pdf",sep="")) 
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = t(geneModuleMembershipM),
               xLabels = colnames(CBT1),
               yLabels = names(MEsMT1),
               ySymbols =namesplot4,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships, Spist- All p7,30genes, 38modules vs Phenotypic trait"))

dev.off()

dev.off()
library("pheatmap")
pdf(paste("Module-physio_Spist_ALL_p7_30genes_38modules_heatmap2",".pdf",sep="")) 
#par(mar = c(12, 12, 3, 3))

b<-pheatmap(t(geneModuleMembershipM),
            cluster_cols=T, 
            cluster_row=T,
            ySymbols =namesplot4,
            col= blueWhiteRed(50),
            display_numbers = t(textMatrix),
            #angle_col=45,
            labels_row=namesplot4,
            labels_col=namescolumn,
            main= paste("WGCNA Spist-ALL p7, 30genes, 38modules vs trait")
)
dev.off()

dim(table(moduleColorsSpist))



######
##### annotation and save modules ####
########

library("gtools")
annotGOSm<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/spis_tabulated_annots_csv.csv")
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, " ", names = c('Name','Namelong')))

nSamples=nrow(vstCBASScoral6)
geneModuleMembershipM = as.data.frame(cor(vstCBASScoral6, MEsMT1, use = "p"))
MMPvalueM = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembershipM), nSamples))

###### get module saved with gene and hub genes and annotations


for(i in 1:length(colnames(MEsMT1))){
  a<-as.data.frame(table(moduleColorsSpist))[i,1]
  b<-paste("ME",a,sep="")
  d<-paste("MMPvalueM.ME",a,sep="")
  e<-paste("geneModuleMembershipM.ME",a,sep="")
  probesM = colnames(vstCBASScoral6)[moduleColorsSpist==a]
  modulemodule<-data.frame(rownames(geneModuleMembershipM),geneModuleMembershipM[[b]],MMPvalueM[[b]])
  modulemodule2<-modulemodule[modulemodule$rownames.geneModuleMembershipM. %in% probesM,]
  modulemodule4<-annotGOSm2[annotGOSm2$Query$Name %in% modulemodule2$rownames.geneModuleMembershipM.,][,c(1,4)]
  modulemodule5<-modulemodule4[order(modulemodule4$Query$Name),]
  modulemodule6<-cbind(modulemodule2,modulemodule5)
  modulemodule7<-modulemodule6[order(modulemodule6[,3]),]
  write.csv(modulemodule7, file = paste0("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Spist_All",a,".csv"))
  dim(subset(modulemodule7,modulemodule7[,2]> 0.9 | modulemodule7[,2] < -0.9))
  moduleHUB<-subset(modulemodule7,modulemodule7[,2]> 0.9 | modulemodule7[,2] < -0.9)
  write.csv(moduleHUB, file = paste0("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Spist_ALL_HUB_",a,".csv"))
}



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


###########
########### file 1 length
###########

annotgff3M<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_RNAseq_ANALYSIS/00_GENOME_reefgenomics/01_Coral/Stylophora_pistillata_genome_reefgenomics/Spis.genome.annotation.gff3",header=F)
annotgff3M2<-subset(annotgff3M,annotgff3M$V3=="exon")
head(annotgff3M2)

rm(annotgff3M)
dim(annotgff3M2)
annotgff3M3= transform(annotgff3M2, V9 = colsplit(V9, ";", names = c('ID','Parent')))
annotgff3M4= transform(annotgff3M3$V9, Names = colsplit(Parent, "=", names = c('V9', 'Names')))
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
annotGOSm2<-transform(annotGOSm, Query = colsplit(Query, " ", names = c('Name','Namelong')))
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
annotyGO2<-transform(annotyGO$gene.GO, gene.GO=colsplit(annotyGO$gene.GO, " ", names = c('Gene', 'GO')))
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



########### for list of gene from each Spis module
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/")
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
  DETS2<-as.character(DETS2$rownames.geneModuleMembershipM.)
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
  write.csv2(newfile4,paste("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/08_Last_analysis_modules/03_WGCNA_Spist_all_80/",name,"GO_enrichement_",dimen,"_terms.csv"))
}















































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
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/01_Module/")

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
  DETS2<-as.character(DETS2$rownames.geneModuleMembershipM.)
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
  write.csv2(newfile4,paste("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/02_Go_enrichment_of_module/",name,"GO_enrichement_",dimen,"_terms.csv"))
}


