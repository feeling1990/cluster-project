library(mvtnorm)
library(EMCluster)
library(mcclust)

#functions needed
classassign=function(gammaKn){
  apply(gammaKn,1,which.max)
}

myaccuracy = function(truelabel,mylabel){
  return(1-classError(mylabel,truelabel)$errorRate)
}

KL_dist=function(p,q){
  if(p==0){
    f=0} else if(q==0){
      f=1} else if(p!=0 & q!=0){
        f=p*(log(p)-log(q))} else{
          f=NULL
        }
  return(f)
}

Mdiff1=function(x,y){
  h=x
  for (i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      h[i,j]=KL_dist(x[i,j],y[i,j])
    }
  }
  return(sum(h)/(dim(x)[1]*dim(x)[2]))
}

#parallel computing
# library(parallel)
# no_cores=detectCores()-1
# cl=makeCluster(no_cores)

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

####ECSM algorithm##################################
ECSM_basic=function(xx,obj, c,b,m,cat,K,maxiter,truelabel,featurenum,threshold1, threshold2){
  D=c+b+m;N=nrow(xx);delete=NULL;keep=1:D;xx_new=xx; delta_matrix=matrix(0,maxiter,D)
  c1=c;b1=b;m1=m; cat1=cat
  # Initialization
  i=1
  x=xx[,1:(c+b)]
  gammaKn = e.step(x,obj)
  df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
  v <- Mstep(c,b,m,cat,xx,t(df));
  xx_new=xx
  cur.loglik <- loglike(c,b,m,cat,K,xx_new,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]])
  loglik.vector <- cur.loglik
  for(h in 1:40){
    v <- Mstep(c,b,m,cat,xx_new, gammaKn);
    gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
  }
  for (i in 2:maxiter) {
    # Repeat E and M steps till convergence
    gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
    v <- Mstep(c,b,m,cat,xx_new, gammaKn);
    if(c>=1){GammaKn1=lapply(1:c,function(j){Estep(c-1,b,m,cat,xx_new[,-j],K,v[[1]],v[[2]][-j,],v[[3]][-j,],v[[4]],v[[5]])})}
    if(c==0) {GammaKn1=NULL}
    if(b>=1){GammaKn2=lapply(1:b,function(j){Estep(c,b-1,m,cat,xx_new[,-(j+c)],K,v[[1]],v[[2]],v[[3]],v[[4]][-j],v[[5]])})}
    if(b==0) {GammaKn2=NULL}
    if(m>=1){GammaKn3=lapply(1:m,function(j){
      if(j==1){cut=(1:cat[1])+c+b}
      if(j>1){cut=c((sum(cat[1:(j-1)])+1): sum(cat[1:j]))+c+b}
      Estep(c,b,m-1,cat[-j],xx_new[,-cut],K,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]][-j])})}
    if(m==0) {GammaKn3=NULL}
    GammaKn=c(GammaKn1,GammaKn2,GammaKn3)
    d=c+b+m
    delta_new=sapply(1:d,function(j){Mdiff1(gammaKn,GammaKn[[j]])})
    delta_matrix[i,keep]=delta_new
    index=which.min(delta_new)
    if((abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))<threshold1)&(delta_new[index]<threshold2)&(i>3))
    {
      delete=c(delete,keep[index])
      c=c-1*(keep[index]<=c1); b=b-1*(keep[index]<=b1+c1 & keep[index]>c1);m=m-1*(keep[index]>c1+b1)
      if(keep[index]>c1+b1){
        l=index-c-b
        if(l==1){cut1=(1:cat[1])+c+b}
        if(l>1){cut1=c((sum(cat[1:(l-1)])+1): sum(cat[1:l]))+c+b}
        xx_new=xx_new[,-cut1];cat=cat[-l]}
      if(keep[index]<=c1+b1){xx_new=xx_new[,-index]}
      keep=keep[-index]
      for(h in 1:20){
        v <- Mstep(c,b,m,cat,xx_new, gammaKn);
        gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
      }
    }
    v <- Mstep(c,b,m,cat,xx_new, gammaKn);
    a <-loglike(c,b,m,cat,K,xx_new,v[[1]],v[[2]],v[[3]],v[[4]],v[[5]]);
    loglik.vector <- c(loglik.vector, a)
    loglik.diff <- abs((cur.loglik - a)) 
    if((loglik.diff < 1e-15)|| length(keep)<=featurenum) {
      for(h in 1:40){
        v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
        gammaKn = Estep(c,b,m,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]]); 
      } 
      break } 
    cur.loglik <- a
  }
  mylabel=classassign(t(gammaKn))
  confusion=table(truelabel,mylabel)
  acc=myaccuracy(confusion)
  results=list(delta_matrix,keep,acc,confusion,delete,i,v)
  return(results)
}

