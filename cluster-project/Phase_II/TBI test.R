x=read.csv("/Users/yinlinfu/Desktop/benchmark data/TBI/track_X_CompX.csv")[,-1]
y=read.csv("/Users/yinlinfu/Desktop/benchmark data/TBI/track_Y_CompY.csv")
#####           
x_cont=data.frame(x[,1])
x_binary=data.frame(x[,2],x[,3],x[,5],x[,6],x[,c(8:15)],x[,17],x[,c(19:28)],x[,34],x[,42])
cate1=data.frame(x[,4],x[,7],x[,16],x[,18],x[,29],x[,30],x[,31],x[,32],x[,33],x[,35],
                 x[,36],x[,37],x[,38],x[,39],x[,40],x[,41])
cate=as.data.frame(apply(cate1,2,as.factor))
x_cat=model.matrix(~ . + 0, data=cate, contrasts.arg = lapply(cate, contrasts, contrasts=FALSE))
cat=c(4,2,8,3,3,3,4,6,3,2,2,3,3,2,3,3)
GOSE90=x[,43] ;GOSE180=x[,44]                              
c=1;b=ncol(x_binary);m=length(cat);s=0
xx=as.matrix(data.frame(x_cont,x_binary,x_cat))
K=2

clusternum=2; maxiter=300; truelabel=GOSE90; featurenum=5; threshold1=1e-5;threshold2=0.005
TBIresults1=ECSM(xx,c=1,b=ncol(x_binary),m=length(cat),s=0,cat=c(4,2,8,3,3,3,4,6,3,2,2,3,3,2,3,3),
                 K=2,maxiter=400,truelabel=GOSE90,featurenum=6,threshold1=1e-5, threshold2=0.05)


###Bing's paper data
x=read.csv("/Users/yinlinfu/Dropbox/benchmark data/TBI/track_X_CompX.csv")[,-1]
GOSE90=x[,43] ;GOSE180=x[,44] 
x=read.csv("/Users/yinlinfu/Dropbox/benchmark data/TBI/2017_5_7_track_complete_3.csv")[,-1]
x1=as.data.frame(apply(x,2,as.factor))
x_cont=data.frame(x[,1],x[,9])
A7=as.factor(x[,7]-1)
A25=x[,25]-2
A26=x[,26]-2
A28=x[,28]-1 ##very unbalanced can be removed only 2 "0"
A30=x[,30]-2
A31=x[,31]-2
A33=x[,33]-1 ##1/484
x_binary=data.frame(x[,2:3],x[,5:6],x[,8],x[,c(10:18)],x[,21],x[,23:24],x[,35:46],x[,52:53],
                    A7,A25,A26,A28,A30,A31,A33)
cate1=data.frame(x[,4],x[,19],x[,20],x[,22],x[,27],x[,29],x[,32],x[,34],x[,47:51])
cate=as.data.frame(apply(cate1,2,as.factor))
x_cat=model.matrix(~ . + 0, data=cate, contrasts.arg = lapply(cate, contrasts, contrasts=FALSE))
cat=c(4,5,7,3,3,3,3,3,3,3,4,6,3)
c=2;b=ncol(x_binary);m=length(cat);s=0
xx=as.matrix(data.frame(x_cont,x_binary,x_cat))
K=5

x=xx[,1:(c+b)]
obj=init.EM(x,nclass=K)
clusternum=5; maxiter=300; truelabel=GOSE90; featurenum=5; threshold1=1e-3;threshold2=0.005

TBIresults5=ECSM(as.matrix(xx),obj,c=2,b=ncol(x_binary),m=length(cat),s=0,cat=c(4,2,5,7,3,2,2,3,2,3,2,2,3,2,3,3,3,4,6,3),
                 K=5,maxiter=400,truelabel=GOSE90,featurenum=12,threshold1=1e-3, threshold2=0.05)


###cluster on 12 critical features
x=read.csv("/Users/yinlinfu/Desktop/benchmark data/TBI/12 critical features.csv")
xx=data.frame(x_cont,x)
c=2;b=12;m=0;s=0
K=5
clusternum=5; maxiter=300; truelabel=GOSE90; featurenum=5; threshold1=1e-3;threshold2=0.005

