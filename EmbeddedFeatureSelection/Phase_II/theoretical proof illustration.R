delta_new=sapply(1:d,function(j){Mdiff1(gammaKn,GammaKn[[j]])})
f=function(p,q){
  if(p==0){f=0}
  if(q==0){f=1}
  if(p!=0 & q!=0){f=p*log(p/q)}
  return(f)
}

Mdiff1=function(x,y){
  h=x
  for (i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      h[i,j]=f(x[i,j],y[i,j])
    }
  }
  return(sum(h)/(dim(x)[1]*dim(x)[2]))
}

gik <-function(c,b,m,s,cat,xn,i,k,mu,sigma,pi,pb,pm,q,lambda)
gk <-function(c,b,m,s,cat,xn,k,mu,sigma,pi,pb,pm,q,lambda)

part=function(data,c,b,m,s,gammaKn,GammaKn,j,v){
  pi=v[[1]];mu=v[[2]];sigma=v[[3]];pb=v[[4]];pm=v[[5]];q=v[[6]];lambda=v[[7]]
  x=gammaKn;y=GammaKn[[j]]
  p1=x;h=x
  for(k in 1:dim(x)[1]){
    for(n in 1:dim(x)[2]){
      xn=data[n,]
      h[k,n]=f(x[k,n],y[k,n])
      if(j>c){p1=x[k,n]*log(gik(c,b,m,s,cat,xn,j,k,mu,sigma,pi,pb,pm,q,lambda))}
      if(j<=c){p1=x[k,n]*log(dmvnorm(xn[(1:c)], mu[,k], sigma[,,k])/dmvnorm(xn[(1:c)][-j], mu[-j,k], sigma[-j,-j,k]))}
    }
  }
  part1=sum(p1)/(dim(x)[1]*dim(x)[2])
  part2=sum(h)/(dim(x)[1]*dim(x)[2])-part1
  return(c(part1,part2,part1+part2))
}


##simulation data
sim1=simulation1_1(n=600)
X=sim1[[1]]; truelabel=sim1[[2]]
c=10;b=1;m=4;s=0;K=2
Y=X[,-c(1:(c+b))]; 
cat=c(3,3,4,10)
Y1=model.matrix(~ . + 0, data=Y, contrasts.arg = lapply(Y, contrasts, contrasts=FALSE))
X1=data.frame(X[,1:(c+b)],Y1)
X1=as.matrix(X1)
xx=X1
x=xx[,1:(c+b)]
obj=init.EM(x,nclass=K)
gammaKn = e.step(x,obj)
df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
gammaKn=t(df)
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
myaccuracy(confusion)
xx_new=xx
for(h in 1:40){
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
}
v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
myaccuracy(confusion)
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


FI=matrix(,15,3)
for(i in 1:15){
  FI[i,]=part(data=xx,c=10,b=1,m=4,s=0,gammaKn,GammaKn,i,v)
}
FI
FI1=FI


FI=as.vector(FI1[,-3])
type=as.factor(c(rep("part1",15),rep("part2",15)))
feature=rep(c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15"),2)
feature=factor(feature,levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15"))
mydata=data.frame(FI,type,feature)
dat1 <- subset(mydata,FI >= 0)
dat2 <- subset(mydata,FI < 0)
figure0_1=ggplot() + 
  geom_bar(data = dat1, aes(x=feature, y=FI),stat = "identity") +
  theme_bw()+ggtitle("values for second term")
figure0_1
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure0_1.tiff', units="in", width=6, height=4, res=300)
figure0_1
dev.off()
figure0_2=ggplot() + 
  geom_bar(data = dat2, aes(x=feature, y=FI),stat = "identity") +
  theme_bw()+ggtitle("values for the first term")
figure0_2
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure0_2.tiff', units="in", width=6, height=4, res=300)
figure0_2
dev.off()

  