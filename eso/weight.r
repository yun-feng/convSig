setwd("/data/ted/multi/eso")
cos_bg<-read.table("bg_cos400.txt",sep="\t",header=F,row.names=NULL)
cos<-read.table("seq_open_cos_approx400.txt",sep="\t",header=F,row.names=NULL)
wkdata<-cos_bg/cos
wkdata<-wkdata/1000.0
write.table(wkdata,"weight_open400.txt",sep="\t",quote=F,row.names = F,col.names =F)
