#install.packages('tsne')
library(tsne)

## simulation data2
simulation2=function(n){
  S1 <- matrix(c(1,0.1,0.1,1),nrow=2,byrow=TRUE)
  mu1 <- c(-1, 2)
  S2 <- matrix(c(1,0.1,0.1,1),nrow=2,byrow=TRUE)
  mu2 <- c(2,-1)
  S=array(0, dim=c(2,2,2))
  S[,,1]=S1; S[,,2]=S2
  #overlap(c(0.5,0.5), Mu=cbind(mu1,mu2), S, eps = 1e-06, lim = 1e06)
  n1 <- n/2 ;n2 <- n/2
  label=c(rep(1,n1),rep(2,n2))
  val1 <- mvrnorm(n1,mu=mu1,Sigma=S1)
  val2 <- mvrnorm(n2,mu=mu2,Sigma=S2)
  allval <- rbind(val1,val2)      ## combine
  x1=allval[,1];
  x2=allval[,2]
  x3=rnorm(n, mean=1.5, sd=1)
  x4=rnorm(n, mean=3, sd=0.5)
  x5=rnorm(n, mean=1.8, sd=0.9)
  x6=rnorm(n, mean=0.3, sd=0.5)
  x7=rnorm(n, mean=2, sd=1)
  x8=rnorm(n, mean=-2, sd=2)
  x9=rnorm(n, mean=4, sd=3)
  x10=rnorm(n, mean=1.5, sd=0.1)
  x11=rnorm(n, mean=-4,sd=2)
  x12_13=mvrnorm(n, mu=c(2,2),Sigma=matrix(c(1,0.3,0.3,1),nrow=2,byrow=TRUE))
  x14_15=mvrnorm(n, mu=c(-3,2),Sigma=matrix(c(1,0.1,0.1,1),nrow=2,byrow=TRUE))
  sim=data.frame(label,allval,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12_13,x14_15)
  order=sample(n,n)
  newlabel=sim[order,1]
  sim2=sim[order,-1]
  #add an outlier
  # sim2 = rbind(sim2,c(-1,2,rnorm(13, 1.1,1)))
  # sim2 = rbind(sim2,c(2,2,rnorm(13, 0.4,0.5)))
  # newlabel = append(newlabel,c(2,2))
  results=list(sim2,newlabel)
  return(results)
}

sim2 = simulation2(n=600)
obj = Mclust(sim2[[1]],G=2) # a more robust way of initialization
start_time <- Sys.time()
myresults2=ESM_fast(sim2[[1]],obj,clusternum=2,maxiter=100,truelabel=sim2[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time
myresults2


misclassified = myresults2[['misclassified']]
deltabig = myresults2[['outliers']]
truelabel = as.factor(sim2[[2]])
data = sim2[[1]][,myresults2[['keep']]] #use selected features
n = length(truelabel)

type = rep('normal', n)
type[misclassified] = 'misclassified' 
type[deltabig] = 'outliers_delta'

data_plot = as.data.frame(cbind(data,truelabel,as.factor(type)))
colnames(data_plot) = c('x1','x2', 'truelabel','type')
ggplot(data_plot, aes(x=x1, y=x2,color=type,  shape=truelabel)) +
  geom_point()+
  ggtitle('Detect Outliers by Feature Selection')+
  theme(plot.title = element_text(hjust = 0.5))


# misclassified_vector = rep(0, n)
# misclassified_vector[misclassified] = 1
# deltabig_vector = rep(0,n)
# deltabig_vector[deltabig] = 1
# 
# data_plot = as.data.frame(cbind(sim2[[1]][,1:2],truelabel,as.factor(deltabig_vector), 
#                                 as.factor(misclassified_vector)))
# colnames(data_plot) = c('x1','x2', 'truelabel','deltabig','misclassified')
# 
# g = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel, size=misclassified, shape=deltabig)) +
#   geom_point()

library(M3C)
tsne(pollen$data,labels=as.factor(pollen$celltypes))
