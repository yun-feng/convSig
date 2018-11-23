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

for(s in 1:10){
step=0.2
lambda_1=3e-2
lambda_2=9e-3
niter=6000

grad<-read.table(paste("bg_grad_rc_sm_si3_",s,".txt",sep=""),header=F,sep="\t")
new_grad=grad
for(c in 1:4){
  target=grad[,c]
  L=length(target)
  u_1=rep(0,L)
  u_2=u_1[2:L]
  
   for(i in 1:niter){
     x=u_1+target+c(u_2,0)-c(0,u_2)
     u_1=u_1-step*x
     u_2=u_2-step*(x[1:(L-1)]-x[2:L])
     u_1[which(u_1>lambda_1)]=lambda_1
     u_1[which(u_1<(-lambda_1))]=-lambda_1
     u_2[which(u_2>lambda_2)]=lambda_2
     u_2[which(u_2<(-lambda_2))]=lambda_2
     
   }
  new_grad[,c]=x
}
grad=new_grad
grad_pos<-data.frame(gradient=unlist(c(grad)),type=rep(c("A","C","G","T"),each=1001),pos=rep(c(-500:500),4))
p<-ggplot(data=grad_pos)+
  geom_point(aes(x=pos,y=gradient,color=type),alpha=0.5,size=1)+
  geom_smooth(aes(x=pos,y=gradient,color=type),alpha=0.5,size=1,method="loess",span=0.3)+
scale_y_continuous(limits=c(-0.25,0.25))
ggsave(paste("bg_grad_rc_sm_si3_post_",s,".png",sep=""),p,width=13,height=7)
}


grad<-read.table(paste("bg_grad_rc_sm_si3_",s,".txt",sep=""),header=F,sep="\t")
mu=grad
L=nrow(grad)
for(c in 1:4){
  target=grad[,c]
  mu[,c]=runif(L)
  sigma=runif(1)
  lambda=runif(1)
  max_range=max(abs(grad))
  for(i in 1:niter){
    mu[,c]=dnorm(target,mean=0,sd=sigma)
    mu[,c]=lambda*mu[,c]/(lambda*mu[,c]+(1-lambda)/(2*max_range))
    
    lambda=sum(mu[,c])/L
    sigma=sqrt((sum(mu[,c]*(target^2)))/sum(mu[,c]))
    
  }
  new_grad[,c]=1-mu[,c]
}
grad=new_grad
grad_pos<-data.frame(gradient=unlist(c(grad)),type=rep(c("A","C","G","T"),each=1001),pos=rep(c(-500:500),4))
p<-ggplot(data=grad_pos)+
  geom_point(aes(x=pos,y=gradient,color=type),alpha=0.5,size=1)+
  geom_smooth(aes(x=pos,y=gradient,color=type),alpha=0.5,size=1,method="loess",span=0.3)
p



#failed normal mixture, for the second normal will have small mean as well
for(c in 1:4){
  target=grad[,c]
  mu[,c]=runif(L)
  sigma=runif(1)
  lambda=runif(1)
  mu_2=runif(1)
  sigma2=runif(1)
  for(i in 1:niter){
    mu[,c]=dnorm(target,mean=0,sd=sigma)
    temp=dnorm(target,mean=mu_2,sd=sigma2)
    mu[,c]=lambda*mu[,c]/(lambda*mu[,c]+(1-lambda)*temp)
    
    lambda=sum(mu[,c])/L
    sigma=sqrt((sum(mu[,c]*(target^2)))/sum(mu[,c]))
    mu_2=sum((1-mu[,c])*target)/sum(1-mu[,c])
    sigma2=sqrt((sum((1-mu[,c])*((target-mu_2)^2)))/sum(1-mu[,c]))
  }
  new_grad[,c]=1-mu[,c]
}

#multivariat guassiam mixture
niter=1000
target=as.matrix(grad)
mu=runif(L)
sigma=matrix(runif(16),ncol=4,nrow=4)
sigma=sigma%*%t(sigma)
lambda=runif(1)
max_range=apply(grad,2,max)-apply(grad,2,min)
for(i in 1:niter){
    temp=sqrt(det(sigma))
    for(j in 1:L){
      mu[j]=(temp/(2*pi)^2)*exp(-0.5*t(target[j,])%*%sigma%*%(target[j,]))
    }
    mu=lambda*mu/(lambda*mu+(1-lambda)/prod(max_range))
    
    lambda=sum(mu)/L
    sigma=solve((t(target)%*%diag(mu)%*%target)/sum(mu))
  }
new_grad=1-mu


