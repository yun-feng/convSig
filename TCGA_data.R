library(maftools)
library(TCGAmutations)

tcga_load(study = "COAD")
getSampleSummary(tcga_coad)

mut<-subsetMaf(maf = tcga_coad)[,c(9,2,3,5,6)]


mut_sorted<-mut[order(mut[,c(3)]),]
mut_sorted<-mut_sorted[order(mut_sorted[,c(2)]),]
mut_sorted<-mut_sorted[order(mut_sorted[,c(1)]),]

setwd("C:\\MutationSignature\\convolve\\deep\\multisample")
write.table(mut_sorted,"TCGA_COAD.txt",sep="\t",quote=F,row.names = F,col.names =F)

names<-mut_sorted[,1]
name<-getSampleSummary(tcga_coad)[,1]

sample=matrix(rep(0,nrow(name)*nrow(names)),nrow=nrow(names),ncol=nrow(name))

for(i in 1:nrow(names)){
  sample[i,as.numeric(names[i])]=1
}
write.table(sample,"TCGA_COAD_sample.txt",sep="\t",quote=F,row.names = F,col.names =F)


res<-c()
for(i in 2:13757){
  if(mut_sorted[i,2]==mut_sorted[i-1,2]){
  res<-c(res,as.numeric(mut_sorted[i,3]-mut_sorted[i-1,3]))
  }
}
plot((log10(res)))
plot(density(log10(res),n=1000))
