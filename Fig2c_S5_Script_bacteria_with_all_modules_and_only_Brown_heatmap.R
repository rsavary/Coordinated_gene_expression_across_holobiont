######heatmap bacteria Smicro modules

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.12")

#install.packages("BiocManager")
#BiocManager::install("WGCNA") 

library("WGCNA")
#install.packages("compositions")
library(compositions)

#1) load environnment of Symbiodinium all samples and keep ME value
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_ALL-30genes_p6_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_all_vst_vstCBASScoral6.RData")
moduleColorsSmicro<-moduleColorsSpist
vstCBASScoral6<-Smicro_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)


###### related Smicro T1 gene coexpression module to OTUs relative abundance #####
#########

#############LOAD bacterial taxonomic association, and bacterial relative abundance
bact<-read.delim("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/05_Paper_version/12_Version_Nils_and_Chris_comment/05_New_bacteria_ASV_counts_and_clr_transformed_counts/Table_S5_NEW.txt",row.names = 1,header=T)
###### clr transform counts
bact_clr=as.data.frame(t(apply(bact[,1:68],2,clr)))

# number of bact taken for the correlation wiht gene co-expression
bact_nb<-50

fam<-bact$Family
sp<-bact$Species
genus<-bact$Genus
fam[is.na(fam)] <- "Unknown Family"
genus[is.na(genus)] <- "Unknown Genus"
sp[is.na(sp)] <- "sp."

#### bacteria species name
species<-paste(fam,";",genus,sp,sep=" ")[1:bact_nb]
#### Clr transformed counts for 60 species of bacteria
df<-bact_clr[,1:bact_nb]
df2<-df[-c(56,58),] #remove samples of bacterial abundance for which we do not have Smic modules
nSamples=nrow(df2)
#### get on both gene expression module and bacteria transformed counts corresponding samples
MEsMT2<-MEsMT1[-c(65,66,69:73),] # remove samples of Smic modules for which we do not have bacterial abundance 

###### Pearson correlation between Smic eigengenes value of co-expression modules and transformed clr abundance
modNamesM = substring(names(MEsMT2), 3)
#or with all modules
geneModuleMembershipM = cor(df2, MEsMT2, use = "p") 
MMPvalueM = corPvalueStudent(as.matrix(geneModuleMembershipM), nSamples)

#### make correction for number of test
MMPvalueMadj<-p.adjust(MMPvalueM, method ="bonferroni") # "bonferroni" or "BH"=Benjamini & Hochberg 
#MMPvalueMadj<-MMPvalueM
#Since  we  have  a  moderately  large  number  of  modules  and  traits,  a  suitable  graphical  representation  will  help  inreading the table.  We color code each association by the correlation value:
sizeGrWindow(10,6)# Will display correlations and their p-values

mat<-signif(MMPvalueMadj, 1)
mat[mat>0.05]<-""
mat2<-signif(geneModuleMembershipM, 2)
mat2[which(mat =="")] <- ""

mat2[mat2==""]<-""
textMatrix = mat2


table(moduleColorsSpist)
namesplot<-as.data.frame(cbind(colnames(MEsMT2),gsub("ME","",colnames(MEsMT2)),c(1:length(colnames(MEsMT2)))))
namesplot<-namesplot[order(namesplot$V2),]
namesplot2<-as.data.frame(cbind(namesplot,as.data.frame(table(moduleColorsSpist))))
namesplot3<-namesplot2[order(as.numeric(namesplot2$V3)),] ## sort by numeric
namesplot4<-paste(namesplot3$V1,namesplot3$Freq,sep=" ")

namescolumn<-row.names(geneModuleMembershipM)
namescolumn<-paste(namescolumn,species,sep=" ")

dim(textMatrix) = dim(geneModuleMembershipM)

# Display the correlation values within a heatmap plot
#pdf(paste("Module-bactos_Smicro_ALL_p6_30genes_19modules_heatmap1",".pdf",sep="")) 
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = t(geneModuleMembershipM),
               xLabels = namescolumn,
               yLabels = names(MEsMT2),
               ySymbols =namesplot4,
               xSymbols = namescolumn,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships, Smicro ALL p6, 30genes, 19modules vs 60 ASV bacterial clr-transformed relative abund"))
#dev.off()

#dev.off()
#pdf(paste("Module-bactos_Smicro_All_p6_30genes_19modules_heatmap2",".pdf",sep="")) 
library("pheatmap")
par(mar = c(12, 12, 3, 3))

pheatmap(t(geneModuleMembershipM),
         cluster_cols=T, 
         cluster_row=F,
         ySymbols =namesplot4,
         col= blueWhiteRed(50),
         display_numbers = t(textMatrix),
         #angle_col=45,
         labels_row=namesplot4,
         labels_col=namescolumn,
         main= paste("WGCNA Smicro ALL p6, 30genes, 19 modules vs 60 most abundant bacterial ASVs"))
#dev.off()

#### Fig. S5

library("pheatmap")
par(mar = c(12, 12, 3, 3))

pheatmap(geneModuleMembershipM,
         cluster_cols=F, 
         cluster_row=T,
         ySymbols =namescolumn,
         col= blueWhiteRed(50),
         display_numbers = textMatrix,
         #angle_col=45,
         labels_row=namescolumn,
         labels_col=namesplot4,
         main= paste("WGCNA Smicro ALL p6, 30genes, 19 modules vs 50 most abundant bacterial ASVs"))
#dev.off()





#### Fig 2c ###


#with only brown module
geneModuleMembershipM = cor(df2, MEsMT2$MEbrown, use = "p") 

MMPvalueM = corPvalueStudent(as.matrix(geneModuleMembershipM), nSamples)

#### make correction for number of test
MMPvalueMadj<-p.adjust(MMPvalueM, method ="bonferroni") # "bonferroni" or "BH"=Benjamini & Hochberg 
#MMPvalueMadj<-MMPvalueM
#Since  we  have  a  moderately  large  number  of  modules  and  traits,  a  suitable  graphical  representation  will  help  inreading the table.  We color code each association by the correlation value:
sizeGrWindow(10,6)# Will display correlations and their p-values

mat<-signif(MMPvalueMadj, 1)
mat[mat>0.05]<-""
mat2<-signif(geneModuleMembershipM, 2)
mat2[which(mat =="")] <- ""

mat2[mat2==""]<-""
textMatrix = mat2


table(moduleColorsSpist)
namesplot<-as.data.frame(cbind(colnames(MEsMT2),gsub("ME","",colnames(MEsMT2)),c(1:length(colnames(MEsMT2)))))
namesplot<-namesplot[order(namesplot$V2),]
namesplot2<-as.data.frame(cbind(namesplot,as.data.frame(table(moduleColorsSpist))))
namesplot3<-namesplot2[order(as.numeric(namesplot2$V3)),] ## sort by numeric
namesplot4<-paste(namesplot3$V1,namesplot3$Freq,sep=" ")

namescolumn<-row.names(geneModuleMembershipM)
namescolumn<-paste(namescolumn,species,sep=" ")

dim(textMatrix) = dim(geneModuleMembershipM)

# Display the correlation values within a heatmap plot
#pdf(paste("Module-bactos_Smicro_ALL_p6_30genes_19modules_heatmap1",".pdf",sep="")) 
par(mar = c(12, 12, 3, 3))
labeledHeatmap(Matrix = t(geneModuleMembershipM),
               xLabels = namescolumn,
               yLabels = names(MEsMT2)[15],
               ySymbols =namesplot4[15],
               xSymbols = namescolumn,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = t(textMatrix),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships, Smicro ALL p6, 30genes, 19modules vs 60 ASV bacterial clr-transformed relative abund"))
#dev.off()

#dev.off()
#pdf(paste("Module-bactos_Smicro_All_p6_30genes_19modules_heatmap2",".pdf",sep="")) 
library("pheatmap")
par(mar = c(12, 12, 3, 3))

pheatmap(t(geneModuleMembershipM),
         cluster_cols=T, 
         cluster_row=F,
         ySymbols =namesplot4[15],
         col= blueWhiteRed(50),
         display_numbers = t(textMatrix),
         #angle_col=45,
         labels_row=namesplot4[15],
         labels_col=namescolumn,
         main= paste("WGCNA Smicro ALL p6, 30genes, Brown modules vs 60 most abundant bacterial ASVs"))
#dev.off()



library("pheatmap")
par(mar = c(12, 12, 3, 3))

pheatmap((geneModuleMembershipM),
         cluster_cols=F, 
         cluster_row=T,
         ySymbols =namesplot4[15],
         col= blueWhiteRed(50),
         display_numbers = textMatrix,
         #angle_col=45,
         labels_row=namescolumn,
         labels_col=namesplot4[15],
         main= paste("WGCNA Smicro ALL p6, 30genes, Brown modules vs 50 most abundant bacterial ASVs"))
#dev.off()














