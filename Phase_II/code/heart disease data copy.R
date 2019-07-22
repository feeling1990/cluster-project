library(ggplot2)
####Heart data
x=read.table("/Users/yinlinfu/Desktop/benchmark data/heart.dat",head=FALSE)
x_label=x[,14]
A1=as.numeric(as.character(x[,1]))
A4=as.numeric(as.character(x[,4]))
A5=as.numeric(as.character(x[,5]))
A8=as.numeric(as.character(x[,8]))
A10=as.numeric(as.character(x[,10]))
x_continuous=data.frame(A1,A4,A5,A8,A10)
x_binary=x[,c(2,6,9)]
A7=as.factor(as.character(x[,7]))
A3=as.factor(as.character(x[,3]))
A11=as.factor(as.character(x[,11]))
A12=as.factor(as.character(x[,12]))
A13=as.factor(as.character(x[,13]))
cate=data.frame(A7,A3,A11,A12,A13)           
x_cate=model.matrix(~ . + 0, data=cate, contrasts.arg = lapply(cate, contrasts, contrasts=FALSE))
heart=data.frame(x_continuous, x_binary,x_cate)
cat=c(3,4,3,3,3)
c=5;b=3;m=5;s=0
xx=heart
clusternum=2; maxiter=300; truelabel=x_label; featurenum=4; threshold1=4e-5;threshold2=0.01
K=2


######## Exp 1. test classical LCA on the simulated data
#########E step and M step basic function###########
gik <-function(c,b,m,cat,xn,i,k,mu,sigma,pi,pb,pm){
  if(i>c & i<=c+b){f=max(pb[[i-c]][k],0.001)^xn[i]*(1-min(pb[[i-c]][k],0.999))^(1-xn[i])}
  if(i<=c){
    if(!is.vector(mu)){f=dnorm(xn[i],mu[i,k],sigma[i,k])}
    else{f=dnorm(xn[i],mu[k],sigma[k])}
  }
  if(i>c+b){
    h=i-c-b;
    if(h==1){cut=(1:cat[1])+c+b}
    if(h>1){cut=c((sum(cat[1:(h-1)])+1): sum(cat[1:h]))+c+b;}
    f=prod(pm[[h]][,k]^xn[cut])}
  return(f)
}

gk <-function(c,b,m,cat,xn,k,mu,sigma,pi,pb,pm){
  g=1;
  for(i in 1:(c+b+m)){
    g<-g*gik(c,b,m,cat,xn,i,k,mu,sigma,pi,pb,pm)
  }
  return(g)
}

responsibility <- function(c,b,m,cat,xn,k,K,mu,sigma,pi,pb,pm){
  a1=pi[k]*gk(c,b,m,cat,xn,k,mu,sigma,pi,pb,pm)
  b1=sum(sapply(1:K,function(j){pi[j]*gk(c,b,m,cat,xn,j,mu,sigma,pi,pb,pm)}))
  a1/b1
}

Estep <- function(c,b,m,cat,xx,K,pi,mu,sigma,pb,pm) {
  apply(xx, 1, function(x) {
    sapply(1:K, function(k) {responsibility(c,b,m,cat,x,k,K,mu,sigma,pi,pb,pm)})
  })
}

nK=function(k,gammaKn){
  sum(gammaKn[k,])
}

muNew=function(c,xx,gammaKn){
  K=nrow(gammaKn)
  nKList <- sapply(1:K, function(k) { nK(k, gammaKn); });
  if(c>1){
    mu=matrix(,c,K)
    for(k in 1:K){
      xxx=data.frame(gammaKn[k,],xx[,1:c])
      mu[,k]=rowSums(apply(xxx,1,function(x){x[1]*x[-1]}))/nKList[k]
    }
  }
  else{
    mu=rep(0,K)
    for(k in 1:K){
      xxx=data.frame(gammaKn[k,],xx[,1:c])
      mu[k]=t(gammaKn[k,])%*%xx[,1]/nKList[k]
    }
  }
  return(mu)
}

pbNew=function(c,b,m,xx,gammaKn){
  K <- nrow(gammaKn);
  if(b!=0){
    pb=list()
    nKList <- sapply(1:K, function(k) { nK(k, gammaKn)});
    for(i in 1:b){
      pb[[i]]=rep(0,K)
      for(k in 1:K){
        xxx=data.frame(gammaKn[k,],xx[,(c+1):(c+b)])
        pb[[i]][k]=t(gammaKn[k,])%*%xx[,(c+i)]/nKList[k]
      }} 
  }
  else{pb=NULL}
  return(pb)
}

pmNew=function(c,b,m,cat,xx,gammaKn){
  K <- nrow(gammaKn);
  if(m!=0){
    pm=list()
    nKList <- sapply(1:K, function(k) { nK(k, gammaKn) });
    for (i in 1:m){
      pm[[i]]=matrix(,cat[i],K)
      if(i==1){cut=(1:cat[1])+c+b}
      if(i>1){cut=c((sum(cat[1:(i-1)])+1): sum(cat[1:i]))+c+b}
      for (k in 1:K){
        pm[[i]][,k]=t(gammaKn[k,])%*%as.matrix(xx[,cut])/nKList[k]
      }
    }
  }
  else{pm=NULL}
  return(pm)}

sigmaNew=function(c,xx,gammaKn,muKNext){
  K <- nrow(gammaKn);
  nKList <- sapply(1:K, function(k) { nK(k, gammaKn) });
  if(c>1){
    muKNext=as.matrix(muKNext)
    N=nrow(xx);sigma=matrix(,c,K);
    for(i in 1:c){
      for(k in 1:K){
        sigma[i,k]=(t(gammaKn[k,])%*%((xx[,i]-rep(muKNext[i,k],N))^2)/nKList[k])^(1/2)
      }
    }}
  else{
    N=nrow(xx);sigma=rep(1,K);
    for(k in 1:K){
      sigma[k]=(t(gammaKn[k,])%*%((xx[,1]-rep(muKNext[k],N))^2)/nKList[k])^(1/2)
    }
  }
  return(sigma)
}

piNew <- function(xx, k, gammaKn, nK) {
  if(is.vector(xx)){nK/length(xx)}
  else{nK / nrow(xx);}
}

Mstep <- function(c,b,m,cat,xx,gammaKn) {
  K <- nrow(gammaKn);
  nKList <- sapply(1:K, function(k) { nK(k, gammaKn); });
  piNext <- sapply(1:K, function(k) { piNew(xx, k, gammaKn, nKList[k]); });
  muNext <- muNew(c,xx,gammaKn);
  pbNext <- pbNew(c,b,m,xx,gammaKn);
  pmNext <- pmNew(c,b,m,cat,xx,gammaKn);
  sigmaNext <- sigmaNew(c,xx,gammaKn,muNext);
  list(piNext, muNext, sigmaNext,pbNext,pmNext); 
}

loglike <- function(c,b,m,cat,K,xx,pi,mu,sigma,pb,pm) {
  sum(apply(xx, 1, function(xn) {
    log( sum(sapply(1:K, function(j) {pi[j]*gk(c,b,m,cat,xn,j,mu,sigma,pi,pb,pm)})),base=exp(1)); }))
}

#########experiments
xx=as.matrix(heart);xx_new=xx
x=xx[,1:(c+b)]
obj=init.EM(x,nclass=K)
gammaKn = e.step(x,obj)
df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
gammaKn=t(df)
mylabel=classassign(t(gammaKn));confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
for(h in 1:20){
  v <- Mstep(c,b,m,cat,xx_new, gammaKn);
  gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
}
v <- Mstep(c,b,m,cat,xx_new, gammaKn);
gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]);
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
save(obj,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart_obj.RData")
load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart_obj.RData")


###2. test on improved LCA using multivariate normal No Feature Selection

########run experiments
xx=heart;xx_new=xx
cat=c(3,4,3,3,3)
c=5;b=3;m=5;s=0; K=2
xx=heart
x=xx[,1:(c+b)]

load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart_obj.RData")
obj=init.EM(x,nclass=K)
gammaKn = e.step(x,obj)
df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
gammaKn=t(df)
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
for(h in 1:20){
  v <- Mstep(c,b,m,cat,xx_new, gammaKn);
  gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
}
v <- Mstep(c,b,m,cat,xx_new, gammaKn);
gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)

##3. test on improved LCA using multivariate normal After Feature Selection
cat=c(3,4,3,3,3)
c=5;b=3;m=5;s=0
xx=heart
clusternum=2; maxiter=300; truelabel=x_label; featurenum=4; threshold1=4e-5;threshold2=0.01
K=2
load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart_obj.RData")
heart1=ECSM(xx,obj,c=5,b=3,m=5,s=0,cat=c(3,4,3,3,3),K=2,maxiter=400,truelabel,featurenum=4,threshold1=8e-5, threshold2=0.04)
heart1
save(heart1,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart1.RData")


##RI plot
iter=17; 
#delta_matrix=heart1[[1]]
load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/delta_matrix for heart")
iterations=1:(iter)
delta=delta_matrix[2:(iter+1),]
delta_v=matrix(delta,ncol=1)
iter_v=rep(iterations, 13)
variable=c(rep("f1",iter),rep("f2",iter),rep("f3",iter),rep("f4",iter),rep("f5",iter),rep("f6",iter),
           rep("f7",iter),rep("f8",iter),rep("f9",iter),rep("f10",iter),rep("f11",iter),rep("f12",iter),
           rep("f13",iter))
feature=factor(variable, levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13"))
Delta=data.frame(iter_v, delta_v, feature)
colnames(Delta)=c("iterations","FI","feature")

figure2=ggplot(Delta,aes(x=iterations, y=FI,shape=feature))+geom_point(size=1)+
  geom_line(size=0.3)+
  scale_x_discrete(name="iterations", breaks=seq(from=0,to=iter,by=5))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  scale_shape_manual(values=1:nlevels(Delta$feature),labels=c(expression(f[1]),expression(f[2]),expression(f[3]),expression(f[4]),expression(f[5]),
                                                              expression(f[6]),expression(f[7]),expression(f[8]),expression(f[9]),expression(f[10]),
                                                              expression(f[11]),expression(f[12]),expression(f[13])))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text",x=17.5,y=0.160,label="f[4]",parse=TRUE,size=4)+
  annotate("text",x=17.5,y=0.116,label="f[5]",parse=TRUE,size=4)+
  annotate("text",x=17.6,y=0.025,label="f[11]",parse=TRUE,size=4)+
  annotate("text",x=14,y=-0.005,label="others(f[1]-f[3],f[6]-f[10],f[12]-f[13])",parse=TRUE,size=4)
figure2
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure2_2.tiff', units="in", width=6, height=4, res=300)
figure2
dev.off()


