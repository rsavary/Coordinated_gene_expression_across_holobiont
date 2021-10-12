###### Expression of Nitrate module with physiological traits.
library("WGCNA")
#1) load environnment of Symbiodinium all samples and keep ME value
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_ALL-30genes_p6_unsigned_networkConstruction-auto.RData")
load("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/11_WGCNA_last_analysis_with_script_forGithub/Smicro_all_vst_vstCBASScoral6.RData")
moduleColorsSmicro<-moduleColorsSpist
vstCBASScoral6<-Smicro_all_vst_vstCBASScoral6

ncol(vstCBASScoral6)
nSamples=nrow(vstCBASScoral6)

MEsT1 = moduleEigengenes(vstCBASScoral6, moduleColorsSpist)$eigengenes
MEsMT1 = orderMEs(MEsT1)

Smicro<-MEsMT1




##### physiological trait
##### load phenotypic data: test one with PAM data make mean data for every samples. after that make a single data frame with the 85 samples and the phenoptypic value associated.
CBT1<-read.csv("/Users/romainsavary/Desktop/EPFL-Post-Doc/09_Wally_RNAseq_analysis/37_WGCNA_For_paper_final/03_Traits/Physiological_traits/All_physiology_all_samples_CBASS_vs_RSS.csv", header = T)
rownames(CBT1)<-CBT1$Name
CBT1<-CBT1[,-1]
View(CBT1)
CBT1<-CBT1[-c(56,58,76:80),]
dim(CBT1)



######Newdataframe
CBT2<-as.data.frame(cbind(Smicro$MEgreenyellow,CBT1,rownames(CBT1)))
colnames(CBT2)[10]<-"name"
CBT3<-transform(CBT2,Name=colsplit(name,"_",names=c("exp","time","temp","gen")))
dim(CBT3)


#In CBASS G15A_1 and G15A_2 and in RSS G15A_1 and G15A_2
#### 7 => G15_2 following PCA scores
CB-T1-29-5-G15A
CB-T1-32-G15A
CB-T3-27-G15A
CB-T3-29-5-G15A
CB-T3-34-G15A     ###### PROBLEME TO FIX HERE SHOULD BE CB-T2-32-G15A !!!!!!!!!!!!!!!
RSS-T1-29-5-G15A
RSS-T3-27-G15A

CBT3$Name.gen[c(8,13,23,28,38,48,61)]
CBT3$Name.gen[c(8,13,23,28,38,48,61)]<-"15_2"
symbol<-CBT3$Name.gen
symbol<-gsub("9","25",symbol)
symbol<-gsub("8","23",symbol)
symbol<-gsub("13","22",symbol)
symbol<-gsub("14","21",symbol)
symbol<-gsub("\\<15\\>","24",symbol)
symbol<-gsub("\\<15_2\\>","4",symbol)
CBT4<-as.data.frame(cbind(CBT3,symbol))

col<-rep(rep(c("blue","yellow","orange","red"),each=5),4)
col1<-col[-c(56,57,76:80)]

##### Fig 3a ######

par(mfrow=c(1,1),mar=c(6,4,2,2))
boxplot(MEsMT1$MEgreenyellow~paste(CBT4$Name.exp,CBT4$Name.temp,CBT4$Name.time,sep="_"),border = "white", las=2, main="Smic-greenyellow module (348 genes), “the nitrate assimilation module”",ylim=c(-0.2,0.2))
rect(16,-0.215,0,0.215,col="grey90",lty=0)    
points(MEsMT1$MEgreenyellow~as.factor(paste(CBT4$Name.exp,CBT4$Name.temp,CBT4$Name.time,sep="_")),pch=as.numeric(symbol), col=col1)


##### Fig 3b #####
symbol2<-symbol

symbol4<-gsub("23","5_G8",symbol2)
symbol5<-gsub("25","1_G9A",symbol4)
symbol6<-gsub("24","3_G15A-1",symbol5)
symbol7<-gsub("4","6_G15A-2",symbol6)
symbol8<-gsub("22","4_G13A",symbol7)
symbol9<-gsub("21","2_G14A",symbol8)


boxplot(MEsMT1$MEgreenyellow~symbol9, las=2,outline=FALSE, ylim=c(-0.10,0.22))
stripchart(MEsMT1$MEgreenyellow~symbol9, vertical = TRUE, 
          method = "jitter", add = TRUE, pch = 20, col = 'black')

summary(aov(MEsMT1$MEgreenyellow~symbol9))
res<-TukeyHSD(aov(MEsMT1$MEgreenyellow~symbol9))
res
options(scipen = 999)
res$symbol9[,4]

text(1,-0.045,"a")
text(2,-0.045,"a")
text(3,-0.045,"ab")
text(4,-0.04,"b")
text(5,0.21,"c")
text(6,0.21,"c")


