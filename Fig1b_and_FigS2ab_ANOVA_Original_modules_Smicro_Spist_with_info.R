######## Correlation between Module eigengene value of Symbiodinium and Coral
library("WGCNA")

#1) load environnment of Symbiodinium all samples and keep ME value
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/Smicro_ALL-30genes_p6_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Smicro_modulesize30/01_Module/data/Smicro_all_vst_vstCBASScoral6.RData")
moduleColorsSmicro<-moduleColorsSpist
vstCBASScoral6<-Smicro_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)

Smicro<-MEsMT1
dim(Smicro)



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

a2<-gsub("-5", ".5", rownames(MEsMT1))
a<-strsplit(a2,"-")
df <- data.frame(matrix(NA,nrow=73,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
gen<-df[,1]

treat1<-paste(time,temp,sep="_")
treat2<-paste(exp,time,temp,sep="_")
######### change gen by spliting G15 in the two genotypes: Indeed according to Savary et al., 2021, SNPs calling and Pictures, genotype G15 is actually the results of two genotypic distinct colonies that have "merged" or that are really close to each other.
rownames(MEsMT1)
#In CBASS G15A_1 and G15A_2 and in RSS G15A_1 and G15A_2
#### 7 => G15_2 following PCA scores
CB-T1-29-5-G15A
CB-T1-32-G15A
CB-T2-27-G15A
CB-T2-29-5-G15A
CB-T2-32-G15A
RSS-T1-29-5-G15A
RSS-T2-27-G15A

rownames(MEsMT1)[c(8,13,23,28,33,48,61)]
gen[c(8,13,23,28,33,48,61)]<-"G15A_2"

gen


##### model selection

#install.packages("AICcmodavg")
library(AICcmodavg)

i=1
model1<-aov(Smicro[,i]~gen+exp+temp+time)
model2<-aov(Smicro[,i]~gen+exp+treat1)
model3<-aov(Smicro[,i]~gen+exp+temp+time+exp/temp/time)
model4<-aov(Smicro[,i]~gen*exp*temp)



colnames(Smicro)[i]
AIC(model1,model2,model3,model3bis,model4)

models<-list(model1,model2,model3,model4)
modelsnames<-c("model1","model2","model3","model4")
aictab(cand.set = models, modnames = modelsnames)



##### save results of anova in files in folder

df <- data.frame(matrix(NA, nrow = dim(Smicro)[2], ncol = 6))
colnames(df)<-c("gene","exp","temp","time","exp-temp","exp-temp-time")

for(i in 1:length(colnames(Smicro))){
  rownames(df)[i]<-colnames(Smicro)[i]
  model<-aov(Smicro[,i]~gen+exp+temp+time+exp/temp/time)
  df[i,1:6]<-summary(model)[[1]][["Pr(>F)"]][1:36]
}

##### plot
pline<- -log10(0.05)
b<- -log10(df)
b<-b[order(-b$gene),]
df<-b
module<-rep(rownames(df),each=6)
module2<-paste(rep(18:1,each=6),module,sep="")
module3<-paste(c(rep("",54),rep("0",54)),module2,sep="")
treat<-rep(paste(c(1:6),colnames(b),sep=""),18)
pval<-c(t(df[1,1:6]),
        t(df[2,1:6]),
        t(df[3,1:6]),
        t(df[4,1:6]),
        t(df[5,1:6]),
        t(df[6,1:6]),
        t(df[7,1:6]),
        t(df[8,1:6]),
        t(df[9,1:6]),
        t(df[10,1:6]),
        t(df[11,1:6]),
        t(df[12,1:6]),
        t(df[13,1:6]),
        t(df[14,1:6]),
        t(df[15,1:6]),
        t(df[16,1:6]),
        t(df[17,1:6]),
        t(df[18,1:6]))


d<-data.frame(a=factor(module3),b=as.numeric(pval),onto=treat)
d<-d[order(d$onto,-d$b),]

p<-ggplot(data=d, aes(x=a, y=b, fill=onto, group = onto)) + geom_bar(stat="identity", fill="darkgrey") + 
  xlab("") + ylab("-log10 (P value)") + ggtitle("ANOVA res module gene expression vs gen+treat+exp") +
  facet_wrap(~onto) 

p + coord_flip() + geom_hline(yintercept = pline, 
                              color = "black", size=0.5)




#2) load environnment of Spist all samples and keep ME value

load("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Spist_modulesize30/01_Module/data/Spist_ALL-30genes_p6_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/02_with_Dan_filtering/02_WGCNA_with_bacteria_and_traits/0X2_test_all80_samples_Spist_modulesize30/01_Module/data/Spist_all_vst_vstCBASScoral6.RData")

vstCBASScoral6<-Spist_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)

Spist<-MEsMT1
dim(Spist)


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

a2<-gsub("-5", ".5", rownames(MEsMT1))
a<-strsplit(a2,"-")
df <- data.frame(matrix(NA,nrow=73,ncol=1))
for(i in 1:length(a)){
  df[i,]<-a[[i]][4]
}
gen<-df[,1]

######### change gen by spliting G15 in the two genotypes:
rownames(MEsMT1)
#In CBASS G15A_1 and G15A_2 and in RSS G15A_1 and G15A_2
#### 7 => G15_2 following PCA scores
CB-T1-29-5-G15A
CB-T1-32-G15A
CB-T3-27-G15A
CB-T3-29-5-G15A
CB-T3-32-G15A
RSS-T1-29-5-G15A
RSS-T3-27-G15A

rownames(MEsMT1)[c(8,13,23,28,33,48,62)]
gen[c(8,13,23,28,33,48,62)]<-"G15A_2"

treat<-paste(time,temp,sep="_")

##### save results of anova in files in folder
setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/05_Paper_version/02_Version/ANOVA_results_genotypes_experiment_model_exp_treat_gen/")

for(i in 1:length(colnames(Spist))){
  name<-colnames(Spist)[i]
  capture.output(summary(aov(Spist[,i]~gen+exp/treat)), file = paste("ANOVA_results_Spist_Module_NESTERD_vs_gene_time-temp_exp_",name,".txt"))
}



####### plot
df <- data.frame(matrix(NA, nrow = dim(Spist)[2], ncol = 6))
colnames(df)<-c("gene","exp","temp","time","exptemp","exptemptime")

for(i in 1:length(colnames(Spist))){
  rownames(df)[i]<-colnames(Spist)[i]
  model<-aov(Spist[,i]~gen+exp+temp+time+exp/temp/time)
  df[i,1:6]<-summary(model)[[1]][["Pr(>F)"]][1:6]
}

##### plot
pline<- -log10(0.05)
b<- -log10(df)
b<-b[order(-b$gene),]
df<-b
module<-rep(rownames(df),each=6)
module2<-paste(rep(38:1,each=6),module,sep="")
module3<-paste(c(rep("",174),rep("0",54)),module2,sep="")
treat<-rep(paste(c(1:6),colnames(b),sep=""),38)
pval<-c(t(df[1,1:6]),
        t(df[2,1:6]),
        t(df[3,1:6]),
        t(df[4,1:6]),
        t(df[5,1:6]),
        t(df[6,1:6]),
        t(df[7,1:6]),
        t(df[8,1:6]),
        t(df[9,1:6]),
        t(df[10,1:6]),
        t(df[11,1:6]),
        t(df[12,1:6]),
        t(df[13,1:6]),
        t(df[14,1:6]),
        t(df[15,1:6]),
        t(df[16,1:6]),
        t(df[17,1:6]),
        t(df[18,1:6]),
        t(df[19,1:6]),
        t(df[20,1:6]),
        t(df[21,1:6]),
        t(df[22,1:6]),
        t(df[23,1:6]),
        t(df[24,1:6]),
        t(df[25,1:6]),
        t(df[26,1:6]),
        t(df[27,1:6]),
        t(df[28,1:6]),
        t(df[29,1:6]),
        t(df[30,1:6]),
        t(df[31,1:6]),
        t(df[32,1:6]),
        t(df[33,1:6]),
        t(df[34,1:6]),
        t(df[35,1:6]),
        t(df[36,1:6]),
        t(df[37,1:6]),
        t(df[38,1:6]))


d<-data.frame(a=factor(module3),b=as.numeric(pval),onto=treat)
d<-d[order(d$onto,-d$b),]

p<-ggplot(data=d, aes(x=a, y=b, fill=onto, group = onto)) + geom_bar(stat="identity", fill="darkgrey") + 
  xlab("") + ylab("-log10 (P value)") + ggtitle("ANOVA res module gene expression vs gen+treat+exp") +
  facet_wrap(~onto) 

p + coord_flip() + geom_hline(yintercept = pline, 
                              color = "black", size=0.5)
















