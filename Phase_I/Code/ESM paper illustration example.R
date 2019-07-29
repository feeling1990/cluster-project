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
  #for(j in 1:p){
  #part2==1/(N*K)*sum(gammaKn%*%(log(dx[,,j]%*%pi)-log(d%*%pi))) ##for one j
  #part1=1/(N*K)*sum(t(gammaKn)*(log(d)-log(dx[,,j])))}
  return(list(gammaKn,RI))
}



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
    if((loglikdiff < 1e-15 & length(keep)<(featurenum+1)) || ncol(xx_new)<featurenum) {
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
  acc=myaccuracy(confusion)
  results=list(delta_matrix,keep,acc,confusion,delete,i)
  return(results)
}
##illustration
set.seed(999)
n=1000
n1 <- n/2 ;n2 <- n/2

x1=rnorm(n,mean=4,sd=1)
x11=rnorm(n1, mean=1.5, sd=1)
x12=rnorm(n2,mean=1.5,sd=5)
x2=c(x11,x12)
print(paste('the mean of the first feature:', mean(x1)))
print(paste('the mean of the second feature:', mean(x2)))
label=c(rep(1,n1),rep(2,n2))
data=data.frame(label,x1,x2)

## reverse the order 
order=sample(n,n)
newlabel=data[order,1]
data1=data[order,-1]

obj=obj=init.EM(data1,nclass=2)
myresults2=ESM(data1,obj,clusternum=2,maxiter=100,truelabel=newlabel,featurenum=1,threshold1=1e-5,threshold2=0.05)
#myresults2=EM_Minus(data1,clusternum=2,maxiter=100,truelabe=newlabel,featurenum=1,threshold=1e-5)

myresults2
v[[1]]=c(0.5,0.5)
v[[2]]=matrix(c(1.5,-1,1.5,-1),2,2)
v[[3]][,,1]=matrix(c(1,0,0,1))
v[[3]][,,2]=matrix(c(1,0,0,10))