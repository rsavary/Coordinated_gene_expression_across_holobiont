##### barplot for each module how many gene are connected to smicro
library("ggplot2")
####
a<-(c("01Greenyellow","02Cyan","03Brown","04Magenta","05Purple","06Lightcyan","07Tan","08Grey60","09Red","10Midnightblue","11Turquoise","12Blue","13Grey","14Black","15Pink","16Salmon","17Green","18Yellow"))
b<-(c(2465,1791,12943,2133,56,0,134,0,191,0,134,0,0,0,0,133,0,616))
b2<-(b/16756)*100 # 25768 new number, only the gene that passed the filtering step
d<-data.frame(a=factor(a),b=as.numeric(b2))


p<-ggplot(data=d, aes(x=a, y=b)) + geom_bar(stat="identity", fill="darkgrey") + 
  xlab("") + ylab("Spist Gene nb")
p 


##### barplot for each module how many gene are connected to spist
a<-(c("01Paleturquoise","02Saddlebrown","03Greenyellow","04Turquoise","05Black","06Blue",
      "07Violet","08Green","09Darkorange","10Red","11Tan","12Pink","13Darkgreen",
      "14Midnightblue","15darkgrey","16Brown","17Yellowgreen","18Darkmagenta","19Yellow","20Sienna",
      "21skyblue3","22skyblue","23grey60","24royalblue","25Cyan","26Magenta","27Steelblue",
      "28Darkturquoise","29Salmon","30orange","31grey","32white","33lightcyan","34lightgreen",
      "35lightyellow","36darkred","37darkolivegreen","38purple"))
b<-(c(3180,0,2299,2299,2299,2299,3433,5552,348,348,1431,
      0,1134,786,1134,1431,2647,0,2299,0,0,2299,2299,2299,2299,4121,
      1054,0,16069,2604,0,0,0,0,0,2299,0,0))
b2<-(b/29557)*100 #new number, only the gene that passed the filtering step

d<-data.frame(a=factor(a),b=as.numeric(b2))
d<-d[order(-d$a),]

p<-ggplot(data=d, aes(x=a, y=b)) + geom_bar(stat="identity", fill="darkgrey") + 
  xlab("") + ylab("Smicro Gene nb")
p # + coord_flip()
