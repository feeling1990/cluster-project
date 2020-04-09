#install.packages("glmpath")
library(glmpath)
library(ggplot2)
####Heart dat
x <- read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/heart-disease/processed.cleveland.data",header=FALSE,sep=",",na.strings = '?')
names(x) <- c( "age", "sex", "cp", "trestbps", "chol","fbs", "restecg",
                        "thalach","exang", "oldpeak","slope", "ca", "thal", "num")

x=read.table("/Volumes/primary/Research/cluster-project/Phase_II/code/heart.dat",head=FALSE)
x_label=x[,14]
A1=as.numeric(as.character(x[,1]))
A4=as.numeric(as.character(x[,4]))
A5=as.numeric(as.character(x[,5]))
A8=as.numeric(as.character(x[,8]))
A10=as.numeric(as.character(x[,10]))
x_continuous=data.frame(A1,A4,A5,A8,A10)
x_binary=x[,c(2,6,9)]
A7=as.factor(as.character(x[,7]))
A3=as.factor(as.character(x[,3]))
A11=as.factor(as.character(x[,11]))
A12=as.factor(as.character(x[,12]))
A13=as.factor(as.character(x[,13]))
cate=data.frame(A7,A3,A11,A12,A13)           
x_cate=model.matrix(~ . + 0, data=cate, contrasts.arg = lapply(cate, contrasts, contrasts=FALSE))
heart=data.frame(x_continuous, x_binary,x_cate)
cat=c(3,4,3,3,3)
c=5;b=3;m=5;s=0
xx=heart
clusternum=2; maxiter=300; truelabel=x_label; featurenum=4; threshold1=4e-5;threshold2=0.01
K=2


###2. test on improved LCA using multivariate normal No Feature Selection

########run experiments
xx=heart;xx_new=xx
cat=c(3,4,3,3,3)
c=5;b=3;m=5;s=0; K=2
xx=heart
x=xx[,1:(c+b)]

obj=init.EM(x,nclass=K)
gammaKn = e.step(x,obj)
df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
gammaKn=t(df)

for(h in 1:20){
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
}
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(truelabel,mylabel)

##3. test on improved LCA using multivariate normal After Feature Selection
cat=c(3,4,3,3,3)
c=5;b=3;m=5;s=0
xx=heart
clusternum=2; maxiter=300; truelabel=x_label; featurenum=4; threshold1=4e-5;threshold2=0.01
K=2
#load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart_obj.RData")
heart1=ECSM_v0(xx,obj,c=5,b=3,m=5,s=0,cat=c(3,4,3,3,3),K=2,maxiter=400,truelabel,featurenum=4,threshold1=4e-5, threshold2=0.04)
heart1
#save(heart1,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/heart1.RData")


##RI plot
iter=17; 
#delta_matrix=heart1[[1]]
load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/delta_matrix for heart")
iterations=1:(iter)
delta=delta_matrix[2:(iter+1),]
delta_v=matrix(delta,ncol=1)
iter_v=rep(iterations, 13)
variable=c(rep("f1",iter),rep("f2",iter),rep("f3",iter),rep("f4",iter),rep("f5",iter),rep("f6",iter),
           rep("f7",iter),rep("f8",iter),rep("f9",iter),rep("f10",iter),rep("f11",iter),rep("f12",iter),
           rep("f13",iter))
feature=factor(variable, levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13"))
Delta=data.frame(iter_v, delta_v, feature)
colnames(Delta)=c("iterations","FI","feature")

figure2=ggplot(Delta,aes(x=iterations, y=FI,shape=feature))+geom_point(size=1)+
  geom_line(size=0.3)+
  scale_x_discrete(name="iterations", breaks=seq(from=0,to=iter,by=5))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  scale_shape_manual(values=1:nlevels(Delta$feature),labels=c(expression(f[1]),expression(f[2]),expression(f[3]),expression(f[4]),expression(f[5]),
                                                              expression(f[6]),expression(f[7]),expression(f[8]),expression(f[9]),expression(f[10]),
                                                              expression(f[11]),expression(f[12]),expression(f[13])))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text",x=17.5,y=0.160,label="f[4]",parse=TRUE,size=4)+
  annotate("text",x=17.5,y=0.116,label="f[5]",parse=TRUE,size=4)+
  annotate("text",x=17.6,y=0.025,label="f[11]",parse=TRUE,size=4)+
  annotate("text",x=14,y=-0.005,label="others(f[1]-f[3],f[6]-f[10],f[12]-f[13])",parse=TRUE,size=4)
figure2
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure2_2.tiff', units="in", width=6, height=4, res=300)
figure2
dev.off()


