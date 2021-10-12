#Fig 1a
#Fig 2b
#Fig 2a
#Fig S3 Spis-turquoise module vs Smic-brown module 
#Fig S4b Fv/Fm vs Spis-turquoise
#Fig S5 OTU1 relative abundance vs Smic-brown module eigengenes values

######## Correlation between Module eigengene value of Symbiodinium and Coral
library("WGCNA")
#1) load environnment of Symbiodinium all samples and keep ME value
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_ALL-30genes_p6_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_all_vst_vstCBASScoral6.RData")
moduleColorsSmicro<-moduleColorsSpist
vstCBASScoral6<-Smicro_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSmicro)$eigengenes
MEsMT1 = orderMEs(MEsT1)

Smicro<-MEsMT1

#2) load environnment of Spist all samples and keep ME value

load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Spist_ALL-30genes_p7_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Spist_all_vst_vstCBASScoral6.RData")

vstCBASScoral6<-Spist_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)

Spist<-MEsMT1

Spist2<-Spist[-c(56,75:78),] # remove samples that are not kept in the Smic dataset (7 samples), two were already not used for Spis module construction
######### ajust both matrix
cbind(rownames(Smicro),rownames(Spist2))
nSamples=nrow(Smicro)

#3)

# plot :

modNamesM = substring(names(Smicro), 3)
geneModuleMembershipM = cor(Smicro, Spist2, use = "p")
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

#### Spist module name with number of gene for the heatmap
table(moduleColorsSpist)
namesplot<-as.data.frame(cbind(colnames(Spist2),gsub("ME","",colnames(Spist2)),c(1:length(colnames(Spist2)))))
namesplot<-namesplot[order(namesplot$V2),]
namesplot2<-as.data.frame(cbind(namesplot,as.data.frame(table(moduleColorsSpist))))
namesplot3<-namesplot2[order(as.numeric(namesplot2$V3)),] ## sort by numeric
namesplot4<-paste(namesplot3$V1,namesplot3$Freq,sep=" ")


#### Smic module name with number of gene for the heatmap
table(moduleColorsSmicro)
namesplot<-as.data.frame(cbind(colnames(Smicro),gsub("ME","",colnames(Smicro)),c(1:length(colnames(Smicro)))))
namesplot<-namesplot[order(namesplot$V2),]
namesplot2<-as.data.frame(cbind(namesplot,as.data.frame(table(moduleColorsSmicro))))
namesplot3<-namesplot2[order(as.numeric(namesplot2$V3)),] ## sort by numeric
namesplotSmicro4<-paste(namesplot3$V1,namesplot3$Freq,sep=" ")

dim(textMatrix) = dim(geneModuleMembershipM)



#set
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/")
# Display the correlation values within a heatmap plot
#pdf(paste("Module-physio_Smicro-Spist_All_p6-7__30genes_18vs38modules_heatmap1",".pdf",sep="")) 
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = t(geneModuleMembershipM),
               xLabels = colnames(Smicro),
               yLabels = names(Spist2),
               ySymbols =namesplot4,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships, Smicro-Spist All p6-7,30genes, 18vs38modules"))

#dev.off()

#dev.off()

###### Fig 1a ######

library("pheatmap")
#pdf(paste("Module-physio_Smicro-Spist_ALL_p6_30genes_19vs38modules_heatmap2",".pdf",sep="")) 
#par(mar = c(12, 12, 3, 3))

b<-pheatmap(t(geneModuleMembershipM),
            cluster_cols=T, 
            cluster_row=T,
            ySymbols =namesplot4,
            col= blueWhiteRed(50),
            display_numbers = t(textMatrix),
            #angle_col=45,
            labels_row=namesplot4,
            labels_col=namesplotSmicro4,
            main= paste("WGCNA Smicro-Spist-ALL p6-7, 30genes, 18 vs 38 modules")
)
#dev.off()




#### add colour for the different stress (27 blue, 29.5 yellow, 32 orange, 34.5 red.)
col1<-c("blue","yellow","orange","red")
####
header<-colnames(t(Smicro))
a<-strsplit(header,"-")
df <- data.frame(matrix(NA,nrow=40,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
col<-df[,1]
col<-col[1:73]
col<-gsub("27", col1[1], col)
col<-gsub("29", col1[2], col)
col<-gsub("32", col1[3], col)
col<-gsub("34", col1[4], col)

### symbole for the genotypes
gen<-0
gen[1:20]<-c("22","21","24","23","25")
gen[21:40]<-c("0","1","2","5","6")
gen[41:60]<-c("22","21","24","23","25")
gen[61:80]<-c("0","1","2","5","6")
gen2<-as.numeric(gen)
gen2<-gen2[c(1:55,57,59,60,61:75)]

df<-as.data.frame(cbind(rownames(Spist2),gen2))
####  G15 is a merge colony of two different genotypes (see Savary et al., 2021), we thus assign two different name to the samples that belong to genotype 1 G15A_1 and that belong to genotype 2 G15A_2
#the 7 following samples are from genotype 2 => G15_2 following PCA scores based SNPs calling from RNAseq (Savary et al., 2021, supplementary figure S3A and B)
CB-T1-29-5-G15A
CB-T1-32-G15A
CB-T3-27-G15A
CB-T3-29-5-G15A
CB-T3-32-G15A
RSS-T1-29-5-G15A
RSS-T3-27-G15A

df[c(8,13,23,28,38,48,61),1]

df$V1[c(8,13,23,28,38,48,61)]
levels(df$gen2) <- c(levels(df$gen2), c("3","4")) # give new symboles to colony G15A-2
df$gen2[c(8,13,48)]<-"3"
df$gen2[c(23,28,38,61)]<-"4"

dff<-df

##### Fig S3 ######

par(mfrow=c(1,1))
plot(Smicro$MEbrown,Spist2$MEturquoise,pch=as.numeric(as.character(df$gen2)),bg=col,col=col,xlab="Smic-brown 2299 genes", ylab="Spis-turquoise, 9132 genes")
abline(lm(Spist2$MEturquoise~Smicro$MEbrown))


######heatmap bacteria Smicro modules
#1) reload environnment of Symbiodinium all samples and keep ME value
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_ALL-30genes_p6_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_all_vst_vstCBASScoral6.RData")
moduleColorsSmicro<-moduleColorsSpist
vstCBASScoral6<-Smicro_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)


library("reshape2")



##### load phenotypic data: test one with PAM data make mean data for every samples. after that make a single data frame with the 85 samples and the phenoptypic value associated.
CBT1<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/03_Traits/Physiological_traits/All_physiology_all_samples_CBASS_vs_RSS.csv", header = T)
rownames(CBT1)<-CBT1$Name
CBT1<-CBT1[,-1]
View(CBT1)
CBT1<-CBT1[-c(56,58,76:80),] # remove the physiological data of the 7 RNAseq sample removed
dim(CBT1)



############# Fig 2b ####
par(mfrow=c(1,1))

plot(CBT1$PAM~Smicro$MEbrown,pch=as.numeric(as.character(dff$gen2)),bg=col,col=col,xlab="Smic-brown 2,299 genes", ylab="Fv/Fm")
abline(lm(CBT1$PAM~Smicro$MEbrown),col="darkgrey")

par(mfrow=c(1,1))


### Fig 2a #####
a<-strsplit(rownames(MEsMT1),"-")
df <- data.frame(matrix(NA,nrow=73,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][2]
}
time<-df[,1]

a<-strsplit(rownames(MEsMT1),"-")
df <- data.frame(matrix(NA,nrow=73,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][3]
}
temp<-df[,1]

a<-strsplit(rownames(MEsMT1),"-")
df <- data.frame(matrix(NA,nrow=73,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][1]
}
exp<-df[,1]

a<-strsplit(rownames(MEsMT1),"-")
df <- data.frame(matrix(NA,nrow=44,ncol=1))
for(i in 1:length(a)){
  df[i,]<-sub('.*(?=.{4}$)', '', header[i], perl=T)
}
gen<-df[,1]
gen<-gsub('[-]', '', gen)

a<-c(22,21,24,23,25)
#a<-c(0,1,2,5,6)
b<-rep(a,16)

data<-gsub("T3", "T2", paste(exp,temp,time,sep="_"))
par(mfrow=c(1,1),mar=c(6,4,2,2))
boxplot(MEsMT1$MEbrown~data, las=2,border="white", main="Smic-brown module (2,299 genes), “the heat-stress module”",ylab="Module eigengene expression",ylim=c(-0.4,0.2),outline=FALSE)#col=c("red","green","red","green","red","green","red","green","red","green","red","green","red","green","red") #col=c("red","red","red","red","green","green","green","green","red","red","red","red","green","green","green"),outline=FALSE)
points(MEsMT1$MEbrown~as.factor(data),pch=as.numeric(as.character(dff$gen2)),col=col,bg=col)
