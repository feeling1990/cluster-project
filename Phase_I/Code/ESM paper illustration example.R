responsibility <- function(xn, k, K, pi, mu, sigma) {
  if(is.vector(mu)){
    a <- pi[k] * dnorm(xn, mu[k], sigma[k]);
    b <- sum(sapply(1:K, function(j) { pi[j] * dnorm(xn, mu[j], sigma[j]) }));
    a/b;
  }
  else{
    a <- pi[k] * dmvnorm(xn, mu[,k], sigma[,,k]);
    b <- sum(sapply(1:K, function(j) { pi[j] * dmvnorm(xn, mu[,j], sigma[,,j]) }));
    a / b;
  } 
}


##illustration
n=100
n1 <- n/2 ;n2 <- n/2
x11=rnorm(n1, mean=1.5, sd=1);x12=rnorm(n2,mean=1.5,sd=10)
x2=c(x11,x12)
x1=rnorm(n,mean=-1,sd=1)
label=c(rep(1,n1),rep(2,n2))
data=data.frame(label,x1,x2)
order=sample(n,n)
newlabel=data[order,1]
data1=data[order,-1]

myresults2=EM_Minus(data1,clusternum=2,maxiter=100,truelabe=newlabel,featurenum=1,threshold=1e-5)

v[[1]]=c(0.5,0.5)
v[[2]]=matrix(c(1.5,-1,1.5,-1),2,2)
v[[3]][,,1]=matrix(c(1,0,0,1))
v[[3]][,,2]=matrix(c(1,0,0,10))