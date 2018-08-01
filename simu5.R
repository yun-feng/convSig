#theta<-read.table("C:\\MutationSignature\\comp\\XRXS_base3\\XRXSsig.txt",sep="\t",header=F)
setwd("C:\\MutationSignature\\comp\\base5_simu_sup_weak")

temp<-read.table("signature.txt",sep="\t",header=T,row.names=1)

temp<-as.matrix(temp[,c(1,5,6,18,20)])

mutation<-read.table("mutation.txt",sep="\t",header=F)[,1]
temp<-temp[rank(mutation),]

Nsig=dim(temp)[2]

theta<-matrix(0,ncol=Nsig,nrow=96*16)

for(i in 1:Nsig){
  temp_m=matrix(temp[,i],nrow=3,ncol=96/3)
  theta[,i]=rep(rbind(temp_m,temp_m,temp_m,temp_m),4)/16
  
}

theta[seq(915,960,3),1]=theta[seq(915,960,3),1]*2
theta[c(1057,1058,1059,1069,1070,1071,1081,1082,1083,1093,1094,1095),5]=theta[c(1057,1058,1059,1069,1070,1071,1081,1082,1083,1093,1094,1095),5]*2
theta=theta/apply(theta,2,sum)


iter=20;
K=10;


#W=3*10^9/96

#W=rep(W,32)

#write.table(t(W),"Background.txt",sep="\t",quote=F,row.names = F,col.names =F)

W=read.table("Background.txt",sep="\t",header=F)
#W=rep(W,each=3)

#EMu=rep(W,each=K)
#EMu=matrix(EMu,ncol=96,nrow=K)
#write.table(EMu,"Background_EMu.txt",sep="\t",quote=F,row.names = F,col.names =F)





for(s in 1:iter){
  
  Activity=-5#-runif(Nsig)*3-3
  
  X=c()
  p=c()
  for(i in 1:K){
    temp_p=c(runif(Nsig)*2*10^Activity)
    
    temp_X<-rep(0,96*16)
    for(j in 1:Nsig){
      temp_X=temp_X+temp_p[j]*theta[,j]
    }
    
    temp_X<-(W*temp_X)
    for(j in 1:(96*16)){
      temp_X[j]<-rpois(1,temp_X[1,j])
      #temp_X[j]<-rnorm(1,temp_X[j],2)
      #if(temp_X[j]<0){temp_X[j]=0}
    }
    X<-rbind(X,as.integer(temp_X+0.5))
    p<-rbind(p,temp_p)
  }
  
  write.table(X,paste("Simu",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
  write.table(p,paste("p",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
  
}

##exp
err<-c()
for(s in 1:iter){
  conv<-as.matrix(read.table(paste("exp\\conv",s,".txt",sep=""),sep="\t",header=F))
  temp_err<-c()
  for(i in 1:(Nsig)){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),2,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*conv,2,sum)/sqrt(apply(conv^2,2,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}


##ReLU
#err<-c()
for(s in 1:iter){
  conv<-as.matrix(read.table(paste("ReLU\\conv_",s,".txt",sep=""),sep="\t",header=F))
  #conv=conv/apply(conv,2,sum)
  temp_err<-c()
  for(i in 1:(Nsig)){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),2,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*conv,2,sum)/sqrt(apply(conv^2,2,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}


##pmsignature
for(s in 1:iter){
  X=read.table(paste("Simu",s,".txt",sep=""),sep="\t",header=F)
  
  SamName=c()
  data=c()
  for(i in 1:ncol(X)){
    type=c((i-1)%%3)
    temp=(i-1)%/%3
    for(j in 1:5){
      if(j!=3){
        type=c(type,temp%%4);
        temp=temp%/%4
      }
      else{
        type=c(type,temp%%2);
        temp=temp%/%2
      }
    }
    for(j in 1:nrow(X)){
      temp_SamName=rep(paste("sample_",j,sep=""),X[j,i])
      temp_data=rep(c(type[4]*3+type[1]+1,type[6]+1,type[5]+1,type[3]+1,type[2]+1,1),X[j,i])
      temp_data=t(matrix(temp_data,nrow=6,ncol=X[j,i]))
      SamName=c(SamName,temp_SamName)
      data=rbind(data,temp_data)
    }
  }
  wkdata<-data.frame(SamName,data)
  write.table(wkdata,paste("C:\\MutationSignature\\comp\\base5_simu_sup\\pmsignature\\simu",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
}

library(pmsignature)

for(s in 1:iter){
  G <- readMFVFile(paste("C:\\MutationSignature\\comp\\base5_simu_sup\\pmsignature\\simu",s,".txt",sep=""), numBases = 5, type="independent", trDir=TRUE)
  
  Param <- getPMSignature(G, K = Nsig,numInit=1)
  
  conv=c()
  for(i in 1:ncol(X)){
    type=c((i-1)%%3)
    temp=(i-1)%/%3
    for(j in 1:5){
      if(j!=3){
        type=c(type,temp%%4);
        temp=temp%/%4
      }
      else{
        type=c(type,temp%%2);
        temp=temp%/%2
      }
    }
    temp_conv=c()
    for (k in 1:Nsig){
      feature=getSignatureValue(Param, k)
      temp_conv=c(temp_conv,feature[1,type[4]*3+type[1]+1]*feature[2,type[6]+1]*feature[3,type[5]+1]*feature[4,type[3]+1]*feature[5,type[2]+1])
    }
    conv=cbind(conv,temp_conv)
  }
  write.table(conv,paste("C:\\MutationSignature\\comp\\base5_simu_sup\\pmsignature\\feature",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
}


#err<-c()
for(s in 1:iter){
  conv<-read.table(paste("pmsignature\\feature",s,".txt",sep=""),sep="\t",header=F)
  temp_err<-c()
  for(i in 1:(Nsig)){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),1,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*t(conv),2,sum)/sqrt(apply(conv^2,1,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}


#summary
err=c(err)
sig_id=rep(paste("Sig",1:Nsig),each=3*iter)
Method=rep(rep(c("convSig (Exp)","convSig (ReLU)","pmsignature"),each=iter),Nsig)
Iter=rep(1:iter,3*Nsig)
wkdata<-data.frame(Cos=err,Sig_ID=sig_id,Method=Method,Iter=Iter)
#write.table(wkdata,"wkdata.txt",sep="\t",quote=F,row.names = F,col.names =F)
setwd("C:\\MutationSignature\\comp\\base5_simu_sup_weak")
#wkdata<-read.table("wkdata.txt",sep="\t",header=F)
colnames(wkdata)<-c("Cos","Sig_ID","Method","Iter")
wkdata$Method<-factor(wkdata$Method,levels=c("pmsignature","convSig (Exp)","convSig (ReLU)"))

library(ggplot2)
p<-ggplot(wkdata)+geom_boxplot(aes(x=Method,y=Cos,fill=factor(Method)))+
  facet_wrap(~Sig_ID,ncol=5)+
  theme(text = element_text(size=25))+
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank())+
  labs(x=NULL,y="Cosine Similarity",fill="Method")


p

ggsave("Base5_simu_sup_weak_f.png",p,width=13,height=7)


x<-err
y<-err
p_value=c()
for (i in 1:Nsig){
  res<-t.test(x[,i],y[,i],alternative ="greater")
  p_value<-c(p_value,res$p.value)
}


