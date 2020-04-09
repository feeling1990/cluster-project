
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

myresults2
