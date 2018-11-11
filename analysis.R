setwd("C:\\MutationSignature\\TCGA\\eso\\deep")
pre_v<-read.table("res_mean_ac.txt",header=F,sep="\t")
pre_l<-read.table("label_mean_ac.txt",header=F,sep="\t")
library(ggplot2)
pre_v<-exp(pre_v)
lab<-rep("other",1540)
tar=1
lab[which(pre_l[,1]==tar)]<-"target"
prob<-data.frame(type=lab,probability=pre_v[,tar])
p_a<-ggplot(data=prob)+geom_density(aes(x=(probability),color=type,fill=type),alpha=0.5)+
  theme(text = element_text(size=28))+
  theme(axis.text = element_text(size = 16))
p_a

for(i in 1:10){
  grad<-read.table(paste("bg_grad_rc2_sm_",i,".txt",sep=""),header=F,sep="\t")
  grad_pos<-data.frame(gradient=unlist(c(grad)),type=rep(c("A","C","G","T"),each=401),pos=rep(c(-200:200),4))
  p<-ggplot(data=grad_pos)+
    geom_point(aes(x=pos,y=gradient,color=type),alpha=0.5,size=1)+
    geom_smooth(aes(x=pos,y=gradient,color=type),alpha=0.5,size=1,method="loess",span=0.3)#+
    #scale_y_continuous(limits=c(-0.25,0.25))
  
  ggsave(paste("bg_grad_rc2_sm_",i,".png",sep=""),p,width=13,height=7)
}

wkdata<-read.table("loss_mean_ac.txt",header=F,sep="\t")
plot(wkdata[,1])
