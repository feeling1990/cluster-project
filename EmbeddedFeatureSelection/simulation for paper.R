##Sim Exp1
library(ggplot2)
# simulate data first example, 2 relevant and 5 irrelvent variables (independent)
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


RFS10=matrix(,10,2);
acc10=NULL
for (i in 1:10){
sim10=simulation1(n=1000)
obj = Mclust(sim10[[1]],G=2) # a more robust way of initialization
R1000=ESM(sim10[[1]],obj,clusternum=2,maxiter=100,truelabel=sim10[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
acc10=c(acc10,R1000[[4]]);RFS10[i,]=R1000[[3]]}
mean10=mean(acc10); sd10=sd(acc10)

RFS9 =matrix(,10,2);acc9 =NULL
for (i in 1:10){
  sim9=simulation1(n=900)
  obj = Mclust(sim9[[1]],G=2)
  R900=ESM(sim9[[1]],obj,clusternum=2,maxiter=100,truelabel=sim9[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc9=c(acc9,R900[[4]]);RFS9[i,]=R900[[3]]}
mean9=mean(acc9); sd9=sd(acc9)

RFS8 =matrix(,10,2);acc8 =NULL
for (i in 1:10){
  sim8=simulation1(n=800)
  obj = Mclust(sim8[[1]],G=2)
  R800=ESM(sim8[[1]],obj,clusternum=2,maxiter=100,truelabel=sim8[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc8=c(acc8,R800[[4]]);RFS8[i,]=R800[[3]]}
mean8=mean(acc8); sd8=sd(acc8)

RFS7 =matrix(,10,2);acc7 =NULL
for (i in 1:10){
  sim7=simulation1(n=700)
  obj = Mclust(sim7[[1]],G=2)
  R700=ESM(sim7[[1]],obj,clusternum=2,maxiter=100,truelabel=sim7[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc7=c(acc7,R700[[4]]);RFS7[i,]=R700[[3]]}
mean7=mean(acc7); sd7=sd(acc7)

RFS6 =matrix(,10,2);acc6 =NULL
for (i in 1:10){
  sim6=simulation1(n=600)
  obj = Mclust(sim6[[1]],G=2)
  R600=ESM(sim6[[1]],obj,clusternum=2,maxiter=100,truelabel=sim6[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc6=c(acc6,R600[[4]]);RFS6[i,]=R600[[3]]}
mean6=mean(acc6); sd6=sd(acc6)

RFS5 =matrix(,10,2);acc5 =NULL
for (i in 1:10){
  sim5=simulation1(n=500)
  obj = Mclust(sim5[[1]],G=2)
  R500=ESM(sim5[[1]],obj,clusternum=2,maxiter=100,truelabel=sim5[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc5=c(acc5,R500[[4]]);RFS5[i,]=R500[[3]]}
mean5=mean(acc5); sd5=sd(acc5)

RFS4 =matrix(,10,2);acc4 =NULL
for (i in 1:10){
  sim4=simulation1(n=400)
  obj = Mclust(sim4[[1]],G=2)
  R400=ESM(sim4[[1]],obj,clusternum=2,maxiter=100,truelabel=sim4[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc4=c(acc4,R400[[4]]);RFS4[i,]=R400[[3]]}
mean4=mean(acc4); sd4=sd(acc4)

RFS3 =matrix(,10,2);acc3 =NULL
for (i in 1:10){
  sim3=simulation1(n=300)
  obj = Mclust(sim3[[1]],G=2)
  R300=ESM(sim3[[1]],obj,clusternum=2,maxiter=100,truelabel=sim3[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc3=c(acc3,R300[[4]]);RFS3[i,]=R300[[3]]}
mean3=mean(acc3); sd3=sd(acc3)

RFS2 =matrix(,10,2);acc2 =NULL
for (i in 1:10){
  sim2=simulation1(n=200)
  obj = Mclust(sim2[[1]],G=2)
  R200=ESM(sim2[[1]],obj,clusternum=2,maxiter=100,truelabel=sim2[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc2=c(acc2,R200[[4]]);RFS2[i,c(1,2)]=R200[[3]][c(1,2)]}
mean2=mean(acc2); sd2=sd(acc2)

RFS1 =matrix(,10,2);acc1 =NULL
for (i in 1:10){
  sim1=simulation1(n=100)
  obj = Mclust(sim1[[1]],G=2)
  R100=ESM(sim1[[1]],obj,clusternum=2,maxiter=100,truelabel=sim1[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
  acc1=c(acc1,R100[[4]]);RFS1[i,c(1,2)]=R100[[3]][c(1,2)]}
mean1=mean(acc1); sd1=sd(acc1)

##keep a record in excel
RFS=cbind(RFS1,RFS2,RFS3,RFS4,RFS5,RFS6,RFS7,RFS8,RFS9,RFS10)
#write.csv(RFS,"C:/Users/yinlinfu/Desktop/EM implementation/RFS in the first experiment for diff data sizes.csv")

##data size plot
ACC=c(mean1,mean2,mean3,mean4,mean5,mean6,mean7,mean8,mean9,mean10)
SD=c(sd1,sd2,sd3,sd4,sd5,sd6,sd7,sd8,sd9,sd10)
minimum=c(min(acc1),min(acc2),min(acc3),min(acc4),min(acc5),min(acc6),min(acc7),min(acc8),min(acc9),min(acc10))
maximum=c(max(acc1),max(acc2),max(acc3),max(acc4),max(acc5),max(acc6),max(acc7),max(acc8),max(acc9),max(acc10))
size=seq(from=100,to=1000,by=100)
sizeplot=data.frame(size,ACC,SD,minimum,maximum)
#write.csv(sizeplot,"c:/Users/yinlinfu/Desktop/EM implementation/sizeplotdata.csv")
#sizeplot=read.csv("/Users/yinlinfu/Desktop/EM implementation/sizeplotdata.csv")

figure1=ggplot(sizeplot,aes(x=size,y=ACC))+geom_point()+scale_x_discrete(name="Data Size", breaks=seq(from=0,to=1000,by=100))+
  geom_errorbar(aes(ymin=ACC-SD,ymax=maximum),width=15)+
  ggtitle("Clustering Performance vs. Data Size")+scale_y_continuous(name="Clustering Accuracy", limits=c(0.3,1))+
  geom_line()+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(size=12))

tiff('/Users/yinlinfu/Desktop/figures/figure1.tiff', units="in", width=6, height=4, res=1000)
figure1
dev.off()

png('/Users/yinlinfu/Desktop/figures/figure1.png', units="in", width=6, height=4, res=1000)
figure1
dev.off()

pdf('/Users/yinlinfu/Desktop/figures/figure1.pdf', width=6, height=4)
figure1
dev.off()

##for n=300 the accuracy with full feature set
  xx=sim3[[1]]; label.true=sim3[[2]]
  obj=init.EM(xx,nclass=2)
  model1=emcluster(xx,obj)
  summary(model1)
  ret.new <- assign.class(xx, model1, return.all = FALSE)
  mylabel1=ret.new$class
  confusion=table(label.true,ret.new$class)
  acc1=(confusion[1,1]+confusion[2,2])/nrow(xx); acc=max(acc1,1-acc1)
  confusion
  acc
  
##RI plot
sim3=simulation1(n=300)
obj = Mclust(sim3[[1]],G=2)
R300=ESM(sim3[[1]],obj,clusternum=2,maxiter=100,truelabel=sim3[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)

#save(R300,file="/Users/yinlinfu/Desktop/figures/R300.RData")
#load("/Users/yinlinfu/Desktop/figures/R300.RData")

iter=12
iterations=1:12
delta=R300[[2]][2:13,]
delta_v=matrix(delta,ncol=1)
iter_v=rep(iterations, 10)
variable=c(rep("f1",iter),rep("f2",iter),rep("f3",iter),rep("f4",iter),rep("f5",iter),rep("f6",iter),
rep("f7",iter),rep("f8",iter),rep("f9",iter),rep("f10",iter))
feature=factor(variable, levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10"))
Delta=data.frame(iter_v, delta_v, feature)
colnames(Delta)=c("iterations","RI","feature")

figure2=ggplot(Delta,aes(x=iterations, y=RI,color=feature))+geom_line()+
  scale_x_discrete(name="", breaks=seq(from=0,to=iter,by=1))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  scale_color_discrete(breaks=levels(feature),labels=c(expression(f[1]),expression(f[2]),expression(f[3]),expression(f[4]),expression(f[5]),
                                                       expression(f[6]),expression(f[7]),expression(f[8]),expression(f[9]),expression(f[10])))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text",x=10,y=0.19,label="f[1]",parse=TRUE,size=4)+
  annotate("text",x=10,y=0.07,label="f[2]",parse=TRUE,size=4)+
  annotate("text",x=9,y=0.02,label="f[3]-f[10]",parse=TRUE,size=4)
figure2

tiff('/Users/yinlinfu/Desktop/figures/figure2.tiff', units="in", width=6, height=4, res=300)
figure2
dev.off()

png('/Users/yinlinfu/Desktop/figures/figure2.png', units="in", width=6, height=4, res=300)
figure2
dev.off()

pdf('/Users/yinlinfu/Desktop/figures/figure2.pdf', width=6, height=4)
figure2
dev.off()
##############################Second simulation example#####################
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
  results=list(sim2,newlabel)
  return(results)
}

sim=simulation2(n=300)
obj = Mclust(sim[[1]], G=2)
R=ESM(sim[[1]],obj,clusternum=2,maxiter=100,truelabel=sim[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)

#save(R,file="/Users/yinlinfu/Desktop/figures/R.RData")
#load("/Users/yinlinfu/Desktop/figures/R.RData")

##delta plot
iter=17
iterations=1:17
delta=R[[2]][2:18,]
delta_v=matrix(delta,ncol=1)
iter_v=rep(iterations, 15)
variable=c(rep("f1",iter),rep("f2",iter),rep("f3",iter),rep("f4",iter),rep("f5",iter),rep("f6",iter),
           rep("f7",iter),rep("f8",iter),rep("f9",iter),rep("f10",iter),rep("f11",iter),rep("f12",iter),rep("f13",iter),rep("f14",iter),rep("f15",iter))
feature=factor(variable, levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15"))
Delta=data.frame(iter_v, delta_v, feature)
colnames(Delta)=c("iterations","RI","feature")
figure3=ggplot(Delta,aes(x=iterations, y=RI,color=feature))+geom_line()+
  scale_x_discrete(name="", breaks=seq(from=0,to=17,by=1))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  scale_color_discrete(breaks=levels(feature),labels=c(expression(f[1]),expression(f[2]),expression(f[3]),expression(f[4]),expression(f[5]),
                                                       expression(f[6]),expression(f[7]),expression(f[8]),expression(f[9]),expression(f[10]),
                                                       expression(f[11]),expression(f[12]),expression(f[13]),expression(f[14]),expression(f[15])))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text",x=17,y=0.097,label="f[1]",parse=TRUE,size=4)+
  annotate("text",x=17,y=0.078,label="f[2]",parse=TRUE,size=4)+
  annotate("text",x=16,y=0.008,label="f[3]-f[15]",parse=TRUE,size=4)
figure3
tiff('/Users/yinlinfu/Desktop/figures/figure3.tiff', units="in", width=6, height=4.2, res=300)
figure3
dev.off()

png('/Users/yinlinfu/Desktop/figures/figure3.png', units="in", width=6, height=4.2, res=300)
figure3
dev.off()

pdf('/Users/yinlinfu/Desktop/figures/figure3.pdf', width=6, height=4)
figure3
dev.off()

###AD data
x=read.csv("C:/Users/yinlinfu/Desktop/AD/AD data for clustering.csv")
labelv1=x[,1];labelv2=x[,2]
xx=x[,-c(1,2,3)]
AD_v=data.frame(labelv1,x[,3],xx)
colnames(AD_v)=c("label","Convert","Age","MMSE","CDR","Hippocampus","Ventricles","WholeBrain","Entorhinal","ICV","HCI","sROI","mcSUVRcere","mcSUVRwm")
ggplot(AD_v,aes(x=mcSUVRcere,y=Hippocampus,color=label,shape=Convert,size=Convert))+geom_point()
