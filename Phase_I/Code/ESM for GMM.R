##EM feature selection algorithm 

## pre-install library packages
#install.packages("EMCluster","mvnfast","parallel","ggplot2")
#install.packages("ggplot2","mclust","VarSelLCM","vscc","SelvarMix","clustvarsel")
library(EMCluster)
library(mvnfast)
library(ggplot2)
library(mclust)
library(VarSelLCM)
library(vscc)
library(SelvarMix)
library(clustvarsel)

## functions needed
myaccuracy = function(truelabel,mylabel){
  return(1-classError(mylabel,truelabel)$errorRate)
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
    a <- pi[k] * dmvn(xn, mu[,k], sigma[,,k]);
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
      if (p>2){ dx[,k,j]=apply(xx[,-j],1,function(xn){dmvn(xn,mu[-j,k],sigma[-j,-j,k])})}
      else{
        dx[,k,j] = sapply(xx[,-j],function(xn){dnorm(xn,mu[-j,k],sigma[-j,-j,k])})
      }
      
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
  return(list(gammaKn,RI))
}

##plot function##
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

####ESM for GMM#####
ESM=function(xx,obj, clusternum,maxiter,truelabel,max_feature_num,threshold1,threshold2){
  p=ncol(xx);N=nrow(xx);delete=NULL;keep=1:p;xx_new=xx; delta_matrix=matrix(0,maxiter,p)
  cl=makeCluster(no_cores);
  cur.loglik=0
  # Initialization way1
  #gammaKn = e.step(xx,obj)
  gammaKn = obj$z
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
    if((loglikdiff < 1e-15 & length(keep)<(max_feature_num+1)) || ncol(xx_new)<max_feature_num) {
      print('yes')
      for(h in 1:10){
        v <- Mstep(xx_new, gammaKn);
        gammaKn = Estep(xx_new,v[[1]],v[[2]],v[[3]]) 
      } 
      break } 
    else {
      cur.loglik <- a
    } 
  }
  i = i + 1
  mylabel=classassign(t(gammaKn))
  confusion=table(truelabel,mylabel)
  acc=myaccuracy(truelabel,mylabel)
  results=list(mylabel,delta_matrix,keep,acc,confusion,delete)
  return(results)
}

##Example 
##test on simulation data
## prepare simulation data
simulation1=function(n){
  S1 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
  mu1 <- c(-1,1)
  S2 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
  mu2 <- c(2,-1)
  n1 <- n/2 ;n2 <- n/2
  label=c(rep(1,n1),rep(2,n2))
  val1 <- mvrnorm(n1,mu=mu1,Sigma=S1)
  val2 <- mvrnorm(n2,mu=mu2,Sigma=S2)
  allval <- rbind(val1,val2)      ## combine
  x3=rnorm(n, mean=1.5, sd=1)
  x4=rnorm(n, mean=3, sd=0.5)
  x5=rnorm(n, mean=1.8, sd=0.9)
  x6=rnorm(n, mean=2.7, sd=1.5)
  x7=rnorm(n, mean=0.3, sd=0.5)
  x8=rnorm(n, mean=0.8, sd=0.9)
  x9=rnorm(n, mean=-2, sd=0.5)
  x10=rnorm(n, mean=-3, sd=0.9)
  sim=data.frame(label,allval,x3,x4,x5,x6,x7,x8,x9,x10)
  order=sample(n,n)
  newlabel=sim[order,1]
  sim1=sim[order,-1]
  results=list(sim1,newlabel)
  return(results)
}



sim3=simulation1(400)
start_time <- Sys.time()
obj = Mclust(sim3[[1]],G=2) # a more robust way of initialization
##another way of initialization is #obj=init.EM(sim3[[1]],nclass=2) and change ESM function initialization way1
myresults2=ESM(sim3[[1]],obj,clusternum=2,maxiter=100,truelabel=sim3[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time

myresults2


