library(mvtnorm)
library(EMCluster)
library(mcclust)
library(mvnfast)
#functions needed
myaccuracy=function(confusion){
  if(dim(confusion)[1]==3){
    type1=(confusion[1,1]+confusion[2,2]+confusion[3,3])/sum(confusion) #1-1;2-2;3-3
    type2=(confusion[1,1]+confusion[2,3]+confusion[3,2])/sum(confusion) #1-1;2-3;3-2
    type3=(confusion[1,2]+confusion[2,1]+confusion[3,3])/sum(confusion) #1-2;2-1;3-3
    type4=(confusion[1,2]+confusion[2,3]+confusion[3,1])/sum(confusion) #1-2;2-3;3-1
    type5=(confusion[1,3]+confusion[2,2]+confusion[3,1])/sum(confusion) #1-3;2-2;3-1
    type6=(confusion[1,3]+confusion[2,1]+confusion[3,2])/sum(confusion) #1-3;2-1;3-2
    acc=max(type1,type2,type3,type4,type5,type6) 
  }
  if(dim(confusion)[1]==2){
    type1=(confusion[1,1]+confusion[2,2])/sum(confusion)
    type2=(confusion[1,2]+confusion[2,1])/sum(confusion)
    acc=max(type1,type2)
  }
  return(acc)
}

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
classassign=function(gammaKn){
  apply(gammaKn,1,which.max)
}

#parallel computing
library(parallel)
no_cores=detectCores()-1
cl=makeCluster(no_cores)

#########E step and M step basic function###########
gik <-function(c,b,m,s,cat,xn,i,k,mu,sigma,pi,pb,pm,q,lambda){
  mm=sum(cat[1:m])
  if(i<=c){
    if(!is.vector(mu)){
      f=dmvn(xn[(1:c)], mu[,k], sigma[,,k])}
    else{f=dnorm(xn[(1:c)],mu[k],sigma[k])}
  }
  if(i>c & i<=c+b){f=max(pb[[i-c]][k],0.001)^xn[i]*(1-min(pb[[i-c]][k],0.999))^(1-xn[i])}
  if(i>c+b & i<=c+b+m){
    h=i-c-b;
    if(h==1){cut=(1:cat[1])+c+b}
    if(h>1){cut=c((sum(cat[1:(h-1)])+1): sum(cat[1:h]))+c+b;}
    f=prod(pm[[h]][,k]^xn[cut])}
  if(i>c+b+m){
    l=i-m-c-b;h=i-m+mm
    if(xn[h]==0){f=1-q[[l]][k]+q[[l]][k]*exp(-lambda[[l]][k])}
    if(xn[h]!=0){f=q[[l]][k]*exp(-lambda[[l]][k])*lambda[[l]][k]^xn[h]/factorial(xn[h])}
  }
  return(f)
}

gk <-function(c,b,m,s,cat,xn,k,mu,sigma,pi,pb,pm,q,lambda){
  g=gik(c,b,m,s,cat,xn,1,k,mu,sigma,pi,pb,pm,q,lambda)#all continuous
  if(b+m+s>0){
    for(i in (c+1):(c+b+m+s)){
      g<-g*gik(c,b,m,s,cat,xn,i,k,mu,sigma,pi,pb,pm,q,lambda)
    }
  }
  return(g)
}

responsibility <- function(c,b,m,s,cat,xn,k,K,mu,sigma,pi,pb,pm,q,lambda){
  a1=pi[k]*gk(c,b,m,s,cat,xn,k,mu,sigma,pi,pb,pm,q,lambda)
  b1=sum(sapply(1:K,function(j){pi[j]*gk(c,b,m,s,cat,xn,j,mu,sigma,pi,pb,pm,q,lambda)}))
  a1/b1
}

Estep <- function(c,b,m,s,cat,xx,K,pi,mu,sigma,pb,pm,q,lambda) {
  apply(xx, 1, function(x) {
    sapply(1:K, function(k) {responsibility(c,b,m,s,cat,x,k,K,mu,sigma,pi,pb,pm,q,lambda)});
  }); 
}

nK=function(k,gammaKn){
  sum(gammaKn[k,])
}

muNew=function(xx,k,gammaKn,nK){
  if(is.vector(xx)){
    gammaKn[k,]%*%xx/nK
  }
  else{
    xxx=data.frame(gammaKn[k,],xx)
    rowSums(apply(xxx,1,function(x){x[1]*x[-1]}))/nK
  }
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

sigmaNew=function(xx,k,gammaKn,muKNew,nK){
  if(is.vector(xx)){
  matrix(gammaKn[k,]%*%((xx-rep(muKNew[k],length(xx)))^2)/nK,ncol=1)
  }
  else{
    xxx=data.frame(gammaKn[k,],xx)
    d=ncol(xxx)-1
    matrix(rowSums(apply(xxx,1,function(x){x[1]*((x[-1]-muKNew[,k])%*%t(x[-1]-muKNew[,k]))}))/nK,ncol=d)
  }
}

rK=function(k,gammaKn,yn){index1=which(yn==0);sum(gammaKn[k,index1])}

lambdaNew<-function(c,b,m,s,xx,gammaKn){
  if(s!=0){lambda=list();K=nrow(gammaKn);
  for(i in 1:s){
    lambda[[i]]=rep(0,K);
    begin=sum(cat[1:m])+c+b
    yi=xx[,(begin+i)]
    for(k in 1:K){
      Nk=nK(k,gammaKn)
      rk=rK(k,gammaKn,yi)
      f=function(x){x*(Nk-rk)-(1-exp(-x))*sum(yi)}
      sol=uniroot.all(f,c(0.001, 1000),tol = .Machine$double.eps^0.2, maxiter = 1000, n = 100)
      if(length(sol)==0){lambda[[i]][k]=gammaKn[k,]%*%yi/Nk}
      if(length(sol)!=0){lambda[[i]][k]=sol}
    }
  }}
  if(s==0){lambda=NULL}
  return(lambda)
}

qNew<- function(c,b,m,s,xx,gammaKn,lambda){
  if(s!=0){q=list(); K <- nrow(gammaKn);
  nKList <- sapply(1:K, function(k) { nK(k, gammaKn) });
  for(i in 1:s){
    q[[i]]=rep(0,K);
    begin=sum(cat[1:m])+c+b
    yi=xx[,(begin+i)]
    for(k in 1:K){
      rk=rK(k,gammaKn,yi)
      q[[i]][k]=(nKList[k]-rk)/nKList[k]/(1-exp(-lambda[[i]][k]))
    }
  }}
  else{q=NULL}
  return(q)
}

piNew <- function(xx, k, gammaKn, nK) {
  if(is.vector(xx)){nK/length(xx)}
  else{nK / nrow(xx);}
}


Mstep <- function(c,b,m,s,cat,xx,gammaKn) {
  K <- nrow(gammaKn);
  nKList <- sapply(1:K, function(k) { nK(k, gammaKn); });
  piNext <- sapply(1:K, function(k) { piNew(xx, k, gammaKn, nKList[k]); });
  muNext <- sapply(1:K, function(k) { muNew(xx[,1:c], k, gammaKn, nKList[k]); });
  pbNext <- pbNew(c,b,m,xx,gammaKn);
  pmNext <- pmNew(c,b,m,cat,xx,gammaKn);
  sigmaNext <- lapply(1:K, function(k) { sigmaNew(xx[,1:c], k, gammaKn, muNext, nKList[k]); });
  sigmaNext=array(unlist(sigmaNext), dim = c(nrow(sigmaNext[[1]]), ncol(sigmaNext[[1]]), length(sigmaNext)))
  lambdaNext<-lambdaNew(c,b,m,s,xx,gammaKn)
  qNext <-qNew(c,b,m,s,xx,gammaKn,lambdaNext);
  list(piNext, muNext, sigmaNext,pbNext,pmNext,qNext,lambdaNext); 
}

loglike <- function(c,b,m,s,cat,K,xx,pi,mu,sigma,pb,pm,q,lambda) {
  sum(apply(xx, 1, function(xn) {
    log( sum(sapply(1:K, function(j) {pi[j]*gk(c,b,m,s,cat,xn,j,mu,sigma,pi,pb,pm,q,lambda)})),base=exp(1)); }))
}

# Initialization

#x=xx[,1:(c+b)]
#obj=init.EM(x,nclass=K)
####ECSM algorithm##################################
ECSM=function(xx,obj,c,b,m,s,cat,K,maxiter,truelabel,featurenum,threshold1, threshold2){
  D=c+b+m+s;N=nrow(xx);delete=NULL;keep=1:D;xx_new=xx; delta_matrix=matrix(0,maxiter,D)
  c1=c;b1=b;m1=m;s1=s; cat1=cat
  i=1
  gammaKn = e.step(x,obj)
  df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
  v <- Mstep(c,b,m,s,cat,xx,t(df));
  xx_new=xx
  cur.loglik <- loglike(c,b,m,s,cat,K,xx_new,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]],v[[6]],v[[7]])
  loglik.vector <- cur.loglik
  for(h in 1:20){
    gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
    v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  }
  for (i in 2:maxiter) {
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
      for(h in 1:20){
        v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
        gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
      }
    }
    v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
    a <-loglike(c,b,m,s,cat,K,xx_new,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]);
    loglik.vector <- c(loglik.vector, a)
    loglik.diff <- abs((cur.loglik - a)) 
    if((loglik.diff < 1e-15 & length(keep)<(featurenum+1))|| length(keep)<=featurenum) {
      for(h in 1:40){
      v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
      gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
      } 
      break } 
    cur.loglik <- a
  }
  mylabel=classassign(t(gammaKn))
  confusion=table(truelabel,mylabel)
  acc=myaccuracy(confusion)
  ARI=arandi(mylabel, truelabel, adjust = TRUE)
  results=list(delta_matrix,keep,acc,confusion,ARI,delete,i,v)
  return(results)
}

