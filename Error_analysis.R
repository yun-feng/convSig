frag_length=1001
Num_frag=100

base=c(1,2,3,4)
err_array=c()
err_array_2=c()
cos_array=c()
for (iter in 1:20){
  T=matrix(rep(0,4*frag_length*Num_frag),ncol=4*frag_length,nrow=Num_frag)
  for(i in 1:Num_frag){
    for(j in 1:frag_length){
      T[i,(j-1)*4+base[floor(runif(1,max=4))+1]]=1
    }
  }

  target=rep(0,4*frag_length)

  for(j in 1:frag_length){
    target[(j-1)*4+base[floor(runif(1,max=6))+1]]=1
  }

  T_inv=solve(T%*%t(T)/frag_length)
  S=1/((T%*%t(T)/frag_length)^2%*%rep(1,Num_frag))
  S=diag(S[,1])

  err_array=c(err_array,sum((T_inv-S)^2)/sum(T_inv^2))

  err_array_2=c(err_array_2,sum((diag(1,Num_frag)-S%*%T%*%t(T)/frag_length)^2)/Num_frag)

  weight_exact=T_inv%*%(T%*%target/frag_length)
  weight_exact[which(weight_exact<0)]=0
  weight_cal=S%*%(T%*%target/frag_length)
  cos_sim=sum(weight_exact*weight_cal)/(sum(weight_cal^2)*sum(weight_exact^2))^(0.5)
  cos_array=c(cos_array,cos_sim)
}

