


###### list all the GO term in all genotypes comparisons of genotypes and do a barplot

setwd("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Spistillata/02_Go_enrichment/Go_enrichment_RSS_witout_34-5C_samples/")

filesT1 <- list.files(pattern =  "\\Spist")

#a<-read.csv2("/Users/romainsavary/Desktop/EPFL-Post-Doc/00_0_Paper_writting/03_Paper_WGCNA_Traits_Smicro/03_Deseq2_genotypes/Smicro/02_Go_enrichment/ Smicro_CBASS_G13A_vs_14A_temp_and_gen_model GO_enrichement_ 4 _terms.csv",header=T)
b<-0
for (i in 1:length(filesT1))
{
  a<-read.csv2(filesT1[i],header=T)
  a<-head(a,-10)
  b<-c(b,as.character(a$term))
}
View(sort(table(b)))



#####Only biological process

b<-0
for (i in 1:length(filesT1))
{
  a<-read.csv2(filesT1[i],header=T)
  a<-head(a,-10)
  a<-subset(a,a$ontology=="BP")
  b<-c(b,a$term)
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



#### data frame with all compariosn with name
b<-data.frame()
for (i in 1:length(filesT1))
{
  a<-read.csv2(filesT1[i],header=T)
  a<-head(a,-10)
  #a[,10]<-filesT1[i]
  d<-head(subset(a,a$term=="regulation of I-kappaB kinase/NF-kappaB signaling"),n=1)
  if(dim(d)[1]>0){
  d[10]<-filesT1[i]
  b<-rbind(b,d)
  }
}
View(b)


