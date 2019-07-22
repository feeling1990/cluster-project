truelabel=GOSE90
K=5; maxiter=400; 
featurenum=5; 

threshold1=1e-5; threshold2=0.01

D=c+b+m+s;N=nrow(xx);delete=NULL;keep=1:D;xx_new=xx; delta_matrix=matrix(0,maxiter,D)
c1=c;b1=b;m1=m;s1=s; cat1=cat
# Initialization
x=xx[,1:(c+b)]
obj=init.EM(x,nclass=K)
#load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart_obj.RData")

i=1
gammaKn = e.step(x,obj)
df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
gammaKn=t(df)
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
arandi(mylabel, truelabel, adjust = TRUE)



for(h in 1:20){
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
}
v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
v[[1]]
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
arandi(mylabel, truelabel, adjust = FALSE)


  # Repeat E and M steps till convergence
  i=i+1
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
  delta_new
  delta_matrix[i,keep]=delta_new
  index=which.min(delta_new)
  if((abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))<threshold1)&(delta_new[index]<threshold2)&(i>2))
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
    for(h in 1:20){
      v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
      gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
    }
  }
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  mylabel=classassign(t(gammaKn))
  confusion=table(truelabel,mylabel)
  confusion
  myaccuracy(confusion)
  arandi(mylabel, truelabel, adjust = TRUE)
  index
  keep[index]
  keep
  delta_new
  abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))

