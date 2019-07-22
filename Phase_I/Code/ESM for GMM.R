##EM feature selection algorithm
install.packages("EMCluster","mvtnorm","mvnfast","parallel")
library(EMCluster)
library(mvtnorm)
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
  if(p!=0 & q!=0){f=p*(log(p)-log(q))}
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

Mdiff=function(x,y){Reduce("+",abs(y-x))/(dim(x)[1]*dim(x)[2])}

#parallel computing
library(parallel)
no_cores=detectCores()-1
cl=makeCluster(no_cores)

responsibility <- function(xn, k, K, pi, mu, sigma) {
  if(is.vector(mu)){
    a <- pi[k] * dnorm(xn, mu[k], sigma[k]);
    b <- sum(sapply(1:K, function(j) { pi[j] * dnorm(xn, mu[j], sigma[j]) }));
    a/b;
  }
  else{
    a <- pi[k] * dmvnorm(xn, mu[,k], sigma[,,k]);
    b <- sum(sapply(1:K, function(j) { pi[j] * dmvn(xn, mu[,j], sigma[,,j]) }));
    a / b;
  } 
}

Estep <- function(xx, pi, mu, sigma) {
  xx=as.matrix(xx)
  if(is.vector(mu)){ K=length(mu) }
  else{K <-ncol(mu);}
  apply(xx, 1, function(x) {
    sapply(1:K, function(k) {
      responsibility(x, k, K, pi, mu, sigma);
    });
  }); 
}

nK=function(xx,k,gammaKn){
  sum(gammaKn[k,])
}

muNew=function(xx,k,gammaKn,nK){
  xxx=data.frame(gammaKn[k,],xx)
  rowSums(apply(xxx,1,function(x){x[1]*x[-1]}))/nK
}

sigmaNew=function(xx,k,gammaKn,muKNew,nK){
  xxx=data.frame(gammaKn[k,],xx)
  d=ncol(xxx)-1
  matrix(rowSums(apply(xxx,1,function(x){x[1]*((x[-1]-muKNew)%*%t(x[-1]-muKNew))}))/nK,ncol=d)
}

piNew <- function(xx, k, gammaKn, nK) {
  if(is.vector(xx)){nK/length(xx)}
  else{nK / nrow(xx);}
}

Mstep <- function(xx, gammaKn) {
  K <- nrow(gammaKn);
  nKList <- sapply(1:K, function(k) { nK(xx, k, gammaKn); });
  piNext <- sapply(1:K, function(k) { piNew(xx, k, gammaKn, nKList[k]); });
  muNext <- sapply(1:K, function(k) { muNew(xx, k, gammaKn, nKList[k]); });
  sigmaNext <- lapply(1:K, function(k) { sigmaNew(xx, k, gammaKn, muNext[,k], nKList[k]); });
  sigmaNext=array(unlist(sigmaNext), dim = c(nrow(sigmaNext[[1]]), ncol(sigmaNext[[1]]), length(sigmaNext)))
  list(piNext, muNext, sigmaNext);
}

loglike <- function(xx, pi, mu, sigma) {
  if(is.vector(mu)){
    sum(apply(xx, 1, function(xn) {
      log( sum(sapply(1:K, function(j) {pi[j] * dnorm(xn, mu[j], sigma[j])})),base=exp(1)); }))
  }
  else{
    K <- ncol(mu);
    sum(apply(xx, 1, function(xn) {
      log( sum(sapply(1:K, function(j) {pi[j] * dmvn(xn, mu[,j], sigma[,,j])})),base=exp(1));
    }));
  }
}

pi=v[[1]];mu=v[[2]];sigma=v[[3]]
Estep1=function(xx,pi, mu,sigma){
  N=nrow(xx);p=ncol(xx);K=length(pi)
  d=matrix(,N,K);dx=array(,dim=c(N,K,p));RI=rep(0,p)
  gammaKn=matrix(,K,N);GammaKn=array(,dim=c(K,N,p))
  for(k in 1:K){
    d[,k]=apply(xx,1,function(xn){dmvn(xn,mu[,k],sigma[,,k])})
    for(j in 1:p){
      dx[,k,j]=apply(xx[,-j],1,function(xn){dmvn(xn,mu[-j,k],sigma[-j,-j,k])})
    }
  }
  for(k in 1:K){
    gammaKn[k,]=t(pi[k]*d[,k]/(d%*%pi))
    for(j in 1:p){
      GammaKn[k,,j]=t(pi[k]*dx[,k,j]/(dx[,,j]%*%pi))
    }
  }
  for(j in 1:p){
    RI[j]=Mdiff1(gammaKn,GammaKn[,,j])
  }
  #for(j in 1:p){
    #part2==1/(N*K)*sum(gammaKn%*%(log(dx[,,j]%*%pi)-log(d%*%pi))) ##for one j
    #part1=1/(N*K)*sum(t(gammaKn)*(log(d)-log(dx[,,j])))}
  return(list(gammaKn,RI))
}


#obj=init.EM(xx,nclass=clusternum)
####ESM for GMM#####
ESM=function(xx,obj, clusternum,maxiter,truelabel,featurenum,threshold1,threshold2){
  p=ncol(xx);N=nrow(xx);delete=NULL;keep=1:p;xx_new=xx; delta_matrix=matrix(0,maxiter,p)
  cl=makeCluster(no_cores);
  cur.loglik=0
  # Initialization
  gammaKn = e.step(xx,obj)
  df <- data.frame(matrix(unlist(gammaKn), ncol=clusternum, byrow=F),stringsAsFactors=FALSE)
  v <- Mstep(xx, t(df));
  cur.loglik <- loglike(xx,v[[1]],v[[2]],v[[3]])
  loglikvector <- cur.loglik
  delta_cur=numeric(p)
  for (i in 2:maxiter) {
    # Repeat E and M steps till convergence
    temper=Estep1(xx_new,v[[1]],v[[2]],v[[3]])
    gammaKn=temper[[1]];delta_new=temper[[2]]
    delta_matrix[i,keep]=temper[[2]]
    index=which.min(delta_new)
    if((abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))<threshold1)&(delta_new[index]<threshold2)&(i>5))
    {delete=c(delete,keep[index])
    xx_new=xx_new[,-index]
    keep=keep[-index]}
    v <- Mstep(xx_new, gammaKn);
    a <-loglike(xx_new,v[[1]],v[[2]],v[[3]])
    loglikvector <- c(loglikvector, a)
    loglikdiff <- abs((cur.loglik - a))
    #loglik.vector
    if((loglik.diff < 1e-15 & length(keep)<(featurenum+1))|| length(keep)<=featurenum) {
      for(h in 1:10){
        v <- Mstep(xx_new, gammaKn);
        gammaKn = Estep(xx_new,v[[1]],v[[2]],v[[3]]) 
      } 
      break } 
    else {
      cur.loglik <- a
    } 
  }
  mylabel=classassign(t(gammaKn))
  confusion=table(truelabel,mylabel)
  acc=myaccuracy(confusion)
  results=list(delta_matrix,keep,acc,confusion,delete,i)
  return(results)
}



##test on simulation data
sim3=simulation1(n=400)
start_time <- Sys.time()
obj=obj=init.EM(sim3[[1]],nclass=2)
myresults2=ESM(sim3[[1]],obj,clusternum=2,maxiter=100,truelabel=sim3[[2]],featurenum=2,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time

myresults2


