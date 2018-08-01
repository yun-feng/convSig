#simu
#theta<-read.table("C:\\MutationSignature\\comp\\XRXS_base3\\XRXSsig.txt",sep="\t",header=F)
setwd("C:\\MutationSignature\\comp\\base3_simu_sup_CRC")

theta<-read.table("signature.txt",sep="\t",header=T,row.names=1)

theta<-as.matrix(theta[,c(1,5,6,18,20)])

Nsig=dim(theta)[2]

iter=20;
K=10;


#W=3*10^9/96

#W=rep(W,32)

#write.table(t(W),"Background.txt",sep="\t",quote=F,row.names = F,col.names =F)

W=read.table("Background.txt",sep="\t",header=F)[[1]]
W=rep(W,each=3)

EMu=rep(W,each=K)
EMu=matrix(EMu,ncol=96,nrow=K)
write.table(EMu,"Background_EMu.txt",sep="\t",quote=F,row.names = F,col.names =F)







mutation<-read.table("mutation.txt",sep="\t",header=F)[,1]
theta<-theta[rank(mutation),]

for(s in 1:iter){
  
  X=c()
  p=c()
  for(i in 1:K){
    temp_p=c(runif(Nsig)*2*10^-6)
    
    temp_X<-rep(0,96)
    for(j in 1:Nsig){
      temp_X=temp_X+temp_p[j]*theta[,j]
    }
    
    temp_X<-W*temp_X
    for(j in 1:96){
      temp_X[j]<-rpois(1,temp_X[j])
    }
    X<-rbind(X,temp_X)
    p<-rbind(p,temp_p)
  }
  
  write.table(X,paste("Simu",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
  write.table(p,paste("p",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
  
}

##exp
err<-c()
for(s in 1:iter){
  conv<-as.matrix(read.table(paste("exp\\conv",s,".txt",sep=""),sep="\t",header=F))
  conv<-conv/apply(conv,2,sum)[1]
  temp_err<-c()
  for(i in 1:(Nsig)){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),2,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*conv,2,sum)/sqrt(apply(conv^2,2,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}
write.table(err,"exp\\sigerr.txt",sep="\t",quote=F,row.names = F,col.names =F)


##ReLU
err<-c()
for(s in 1:iter){
  conv<-as.matrix(read.table(paste("ReLU\\conv_",s,".txt",sep=""),sep="\t",header=F))
  #conv=conv/apply(conv,2,sum)
  temp_err<-c()
  for(i in 1:Nsig){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),2,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*conv,2,sum)/sqrt(apply(conv^2,2,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}
write.table(err,"ReLU\\sigerr.txt",sep="\t",quote=F,row.names = F,col.names =F)

##BayesNMF
BayesNMF.L1KL <- function(V0,n.iter,a0,tol,K,K0,phi) {
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[,active_nodes]
  V <- V0-min(V0) + eps
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  W <- matrix(runif(N * K)*sqrt(Vmax),ncol=K)
  H <- matrix(runif(M * K)*sqrt(Vmax),ncol=M)
  V.ap <- W %*% H + eps
  I <- array(1,dim=c(N,M))
  
  C <- N + M + a0 + 1
  b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
  lambda.bound <- b0/C
  lambda <- (colSums(W)+rowSums(H)+b0)/C
  lambda.cut <- 1.5*lambda.bound
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  
  iter <- 2
  count <- 1
  while ((del >= tol) & (iter < n.iter)) {
    H <- H * (t(W) %*% (V/V.ap))/(matrix(rep(colSums(W)+phi/lambda,M),ncol=M) + eps)
    V.ap <- W %*% H + eps
    W <- W * ((V/V.ap) %*% t(H))/t(matrix(rep(rowSums(H)+phi/lambda,N),ncol=N) + eps)
    V.ap <- W %*% H + eps
    lambda <- (colSums(W) + rowSums(H) + b0) / C
    del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
    like <- sum(V*log(V/V.ap)+V.ap-V)
    n.like[[iter]] <- like
    n.evid[[iter]] <- like+phi*sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V-V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
    }
    iter <- iter+1
  }
  return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

n.iter <- 10 
Kcol <- 95   
tol <- 1.e-07
a0 <- 10

for(s in 1:iter){
  lego96<-as.matrix(read.table(paste("C:\\MutationSignature\\comp\\base3_simu_sup\\Simu",s,".txt",sep=""),sep="\t",header=F))
  num=0;
  while(num<(Nsig-3)){
    res <- BayesNMF.L1KL(as.matrix(lego96),100000,a0,tol,Kcol,Kcol,1.0)
    num=length(which(apply(res[[1]],2,sum)>0))
  }
  feature=res[[2]][which(apply(res[[1]],2,sum)>0),]
  
  temp=matrix(0*96*num,nrow=num,ncol=96)
  temp[paste("V",1:96,sep="") %in% colnames(feature)]=feature
  
  p=res[[1]][,which(apply(res[[1]],2,sum)>0)]
  write.table(p,paste("C:\\MutationSignature\\comp\\base3_simu_sup\\BayesNMF\\p",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
  write.table(temp,paste("C:\\MutationSignature\\comp\\base3_simu_sup\\BayesNMF\\theta",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
}

err<-c()
nmf_num<-c()
for(s in 1:iter){
  conv<-as.matrix(read.table(paste("BayesNMF\\theta",s,".txt",sep=""),sep="\t",header=F))
  conv=conv/apply(conv,1,sum)
  temp_err<-c()
  nmf_num<-c(nmf_num,nrow(conv))
  for(i in 1:Nsig){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),1,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*t(conv),2,sum)/sqrt(apply(conv^2,1,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}
write.table(err,"BayesNMF\\sigerr.txt",sep="\t",quote=F,row.names = F,col.names =F)

##EMu
#err<-c()
#for(s in 1:iter){
#  conv<-as.matrix(read.table(paste("C:\\MutationSignature\\comp\\base3_simu_sup\\EMu\\Simu",s,"_ml_spectra.txt",sep=""),sep=" ",header=F))[,1:96]
#  conv=conv/apply(conv,1,sum)
#  temp_err<-c()
#  for(i in 1:Nsig){
#    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),1,sum))/96)
#    temp_err<-c(temp_err,max(apply(theta[,i]*t(conv),2,sum)/sqrt(apply(conv^2,1,sum)*sum(theta[,i]^2))))
#  }
#  err<-rbind(err,temp_err)
#}


##pmsignature
for(s in 1:iter){
  X=read.table(paste("Simu",s,".txt",sep=""),sep="\t",header=F)
  
  SamName=c()
  data=c()
  for(i in 1:ncol(X)){
    type=c((i-1)%%3)
    temp=(i-1)%/%3
    for(j in 1:3){
      if(j!=2){
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
      temp_data=rep(c(type[3]*3+type[1]+1,type[4]+1,type[2]+1,1),X[j,i])
      temp_data=t(matrix(temp_data,nrow=4,ncol=X[j,i]))
      SamName=c(SamName,temp_SamName)
      data=rbind(data,temp_data)
    }
  }
  wkdata<-data.frame(SamName,data)
  write.table(wkdata,paste("C:\\MutationSignature\\comp\\base3_simu_sup\\pmsignature\\simu",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
}

library(pmsignature)

for(s in 1:iter){
  inputFile <-paste("C:\\MutationSignature\\comp\\base3_simu_sup\\pmsignature\\simu",s,".txt",sep="")
  
  G <- readMFVFile(paste("C:\\MutationSignature\\comp\\base3_simu_sup\\pmsignature\\simu",s,".txt",sep=""), numBases = 3, type="independent", trDir=TRUE)
  
  
  BG_prob <- readBGFile(G)
  Param <- getPMSignature(G, K =Nsig,numInit=1)
  
  conv=c()
  for(i in 1:96){
    type=c((i-1)%%3)
    temp=(i-1)%/%3
    for(j in 1:3){
      if(j!=2){
        type=c(type,temp%%4);
        temp=temp%/%4
      }
      else{
        type=c(type,temp%%2);
        temp=temp%/%2
      }
    }
    temp_conv=c()
    for (k in 1:(Nsig)){
      feature=getSignatureValue(Param, k)
      temp_conv=c(temp_conv,feature[1,type[3]*3+type[1]+1]*feature[2,type[4]+1]*feature[3,type[2]+1])
    }
    conv=cbind(conv,temp_conv)
  }
  write.table(conv,paste("C:\\MutationSignature\\comp\\base3_simu_sup\\pmsignature2\\feature",s,".txt",sep=""),sep="\t",quote=F,row.names = F,col.names =F)
}


err<-c()
for(s in 1:iter){
  conv<-read.table(paste("pmsignature\\feature",s,".txt",sep=""),sep="\t",header=F)
  temp_err<-c()
  for(i in 1:(Nsig)){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),1,sum))/96)
    temp_err<-c(temp_err,max(apply(theta[,i]*t(conv),2,sum)/sqrt(apply(conv^2,1,sum)*sum(theta[,i]^2))))
  }
  err<-rbind(err,temp_err)
}
write.table(err,"pmsignature\\sigerr.txt",sep="\t",quote=F,row.names = F,col.names =F)

#summary
err=c(err)
sig_id=rep(paste("Sig",1:Nsig),each=4*iter)
Method=rep(rep(c("convSig (Exp)","convSig (ReLU)","BayesNMF","pmsignature"),each=iter),Nsig)
Iter=rep(1:iter,4*Nsig)
wkdata<-data.frame(Cos=err,Sig_ID=sig_id,Method=Method,Iter=Iter)
#write.table(wkdata,"wkdata.txt",sep="\t",quote=F,row.names = F,col.names =F)
setwd("C:\\MutationSignature\\comp\\base3_simu_sup_CRC")
#wkdata<-read.table("wkdata.txt",sep="\t",header=F)
colnames(wkdata)<-c("Cos","Sig_ID","Method","Iter")
wkdata$Method<-factor(wkdata$Method,levels=c("BayesNMF","pmsignature","convSig (Exp)","convSig (ReLU)"))

library(ggplot2)
p<-ggplot(wkdata)+geom_boxplot(aes(x=Method,y=Cos,fill=factor(Method)))+
  facet_wrap(~Sig_ID,ncol=5)+
  theme(text = element_text(size=25))+
  theme(axis.text = element_text(size = 16),axis.text.x = element_blank())+
  labs(x=NULL,y="Cosine Similarity",fill="Method")

p

ggsave("Base3_simu_sup_CRC_f.png",p,width=13,height=7)


x<-err
y<-err
p_value=c()
for (i in 1:Nsig){
  res<-t.test(x[,i],y[,i],alternative ="less")
  p_value<-c(p_value,res$p.value)
}

test<-read.table("exp\\test_loss.txt",sep="\t")
test<-test[,1:6]+test[,7:12]
test<-test-apply(test,1,min)
loss<-unlist(c(test))
simu<-rep(1:20,6)
num<-rep(2:7,each=20)
wkdata<-data.frame(loss=loss,simu=simu,num=num)
library(ggplot2)
p<-ggplot(wkdata)+geom_point(aes(x=num,y=loss,color=factor(simu)))+
  geom_line(aes(x=num,y=loss,color=factor(simu),group=factor(simu)))+
  theme(text = element_text(size=25))+
  theme(axis.text = element_text(size = 16))+
  labs(x="Number",y="Relatice Loss",color="Simulations")

p

num<-rep(0,6)
for(i in 1:nrow(test)){
  num[which(test[i,]<0.00001)]=num[which(test[i,]<0.00001)]+1
}
num2<-rep(0,6)
for(i in 1:nrow(test)){
  num2[which(test[i,]<0.00001)]=num2[which(test[i,]<0.00001)]+1
}
num<-c(num,num2)
pos<-rep(2:7,2)
method=rep(c("ConvSig(ReLU)","ConvSig(exp)"),each=6)
wkdata<-data.frame(method=method,num=num,pos)
p<-ggplot(wkdata)+geom_bar(aes(x=pos,y=num,fill=factor(method)),stat="identity",position = "dodge",show.legend=F)+
  #geom_rect(data=NULL,aes(xmin=4.5,xmax=5.5,ymin=-Inf,ymax=Inf),fill="gold",alpha=0.01)+
  facet_wrap(~method)+
  theme(text = element_text(size=25))+
  theme(axis.text = element_text(size = 16))+
  labs(x="Number of signatures",y="Simulations")
  

p
ggsave("Base3_number.png",p,width=13,height=7)
