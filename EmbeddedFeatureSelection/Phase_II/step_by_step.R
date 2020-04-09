maxiter=300; max_feature_num=8; threshold1=4e-3; threshold2=0.01 ; option=2

xx = data_info$data
c = data_info$continous
b = data_info$binary
m = data_info$multi
s = data_info$s
K = data_info$numcluster

D=c+b+m+s;N=nrow(xx);delete=NULL;keep=1:D;xx_new=xx; delta_matrix=matrix(0,maxiter,D)
c1=c;b1=b;m1=m;s1=s; cat1=cat
i=1
x=xx[,1:(c+b)]

if (option==1){
  obj1 = data_info$obj1
  gammaKn = e.step(x,obj1)
  df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
  gammaKn = t(df)
} else{
  obj2 = data_info$obj2
  gammaKn = t(obj2$z)
}

v <- Mstep(c,b,m,s,cat,xx,gammaKn);
xx_new=xx
cur.loglik <- loglike(c,b,m,s,cat,K,xx_new,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]],v[[6]],v[[7]])
loglik.vector <- cur.loglik
for(h in 1:5){
  gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]);
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
}

i = 2
  # Repeat E and M steps till convergence
  gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  if(c>=1){GammaKn1=lapply(1:c,function(j){Estep(c-1,b,m,s,cat,xx_new[,-j],K,v[[1]],v[[2]][-j,],v[[3]][-j,-j,],v[[4]],v[[5]],v[[6]],v[[7]])})}
  if(c==0) {GammaKn1=NULL}
  if(b>=1){GammaKn2=lapply(1:b,function(j){Estep(c,b-1,m,s,cat,xx_new[,-(j+c)],K,v[[1]],v[[2]],v[[3]],v[[4]][-j],v[[5]],v[[6]],v[[7]])})}
  if(b==0) {GammaKn2=NULL}
  if(m>=1){GammaKn3=lapply(1:m,function(j){
    if(j==1){cut=(1:cat[1])+c+b}
    if(j>1){cut=c((sum(cat[1:(j-1)])+1): sum(cat[1:j]))+c+b}
    Estep(c,b,m-1,s,cat[-j],xx_new[,-cut],K,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]][-j],v[[6]],v[[7]])})}
  if(m==0) {GammaKn3=NULL}
  if(s>=1){GammaKn4=lapply(1:s,function(j){Estep(c,b,m,s-1,cat,xx_new[,-(j+c+b+sum(cat[1:m]))],K,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]],v[[6]][-j],v[[7]][-j])})}
  if(s==0){GammaKn4=NULL}
  GammaKn=c(GammaKn1,GammaKn2,GammaKn3,GammaKn4)
  d=c+b+m
  delta_new=sapply(1:d,function(j){Mdiff1(gammaKn,GammaKn[[j]])})
  delta_matrix[i,keep]=delta_new
  index=which.min(delta_new)
  if((abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))<threshold1)&(delta_new[index]<threshold2)&(i>3))
  {
    delete=c(delete,keep[index])
    c=c-1*(keep[index]<=c1); b=b-1*(keep[index]<=b1+c1 & keep[index]>c1);
    m=m-1*(keep[index]>c1+b1 & keep[index]<=c1+b1+m1);s=s-1*(keep[index]>c1+b1+m1)
    if(keep[index]>c1+b1 & keep[index]<=c1+b1+m1){
      l=index-c-b
      if(l==1){cut1=(1:cat[1])+c+b}
      if(l>1){cut1=c((sum(cat[1:(l-1)])+1): sum(cat[1:l]))+c+b}
      xx_new=xx_new[,-cut1];cat=cat[-l]}
    if(keep[index]<=c1+b1){xx_new=xx_new[,-index]}
    if(keep[index]>c1+b1+m1){xx_new=xx_new[,-(sum(cat[1:m])+c+b)]}
    keep=keep[-index]
    for(h in 1:10){
      gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]);
      v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
    }
  }
  i = i+1
  
  keep
  abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))
  index
  delta_new[index]
  
  
  
 # v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  a <-loglike(c,b,m,s,cat,K,xx_new,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]);
  loglik.vector <- c(loglik.vector, a)
  loglik.diff <- abs((cur.loglik - a)) 
  if((loglik.diff < 1e-15 & length(keep)<(max_feature_num+1))|| ncol(xx_new)<=2) {
    for(h in 1:40){
      v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
      gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
    } 
    break } 
  cur.loglik <- a
}