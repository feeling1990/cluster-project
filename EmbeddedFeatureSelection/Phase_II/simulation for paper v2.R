library(ggplot2)
library(pkgmaker)
#### simulation Data
simulation1_1=function(n){
  S1 <- matrix(c(1,0.2,0.2,1),nrow=2,byrow=TRUE)
  mu1 <- c(-1,1)
  S2 <- matrix(c(1,0.2,0.2,1),nrow=2,byrow=TRUE)
  mu2 <- c(1,0)
  n1 <- n*2/3 ;n2 <- n/3
  label=c(rep(1,n1),rep(2,n2))
  val1 <- mvrnorm(n1,mu=mu1,Sigma=S1)
  val2 <- mvrnorm(n2,mu=mu2,Sigma=S2)
  allval <- rbind(val1,val2)      ## combine
  z11=rbinom(n1,1,p=0.7);z12=rbinom(n2,1,p=0.2)
  z1=c(z11,z12)
  z21 = sample(letters[1:3], n1,replace=TRUE, prob=c(0.7,0.2,0.2))
  z22 = sample(letters[1:3], n2,replace=TRUE, prob=c(0.2,0.4,0.5))
  z2 = c(z21, z22)
  z3 = sample(letters[1:3], n,replace=TRUE, prob=c(1/3,1/3,1/3))
  z4 = sample(letters[1:4], n,replace=TRUE, prob=c(1/4,1/4,1/4,1/4))
  z5 = sample(letters[1:10], n,replace=TRUE, prob=c(1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10))
  # z21=t(rmultinom(n1,size=1,prob=c(0.7,0.2,0.2))); z21=charmap(z21,maps=c("a","b","c"))
  # z22=t(rmultinom(n2,size=1,prob=c(0.2,0.4,0.5))); z22=charmap(z22,maps=c("a","b","c"))
  # z2=c(z21,z22)
  # z31=t(rmultinom(n,size=1,prob=c(1/3,1/3,1/3))); z3=charmap(z31,maps=c("a","b","c"))
  # z41=t( rmultinom(n,size=1,prob=c(1/4,1/4,1/4,1/4)));z4=charmap(z41,maps=c("a","b","c","d"))
  # z51=t( rmultinom(n,size=1,prob=c(1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10)));z5=charmap(z51,maps=c("a","b","c","d","e","f","g","h","i","j"))
  x3=rnorm(n, mean=1.5, sd=1)
  x4=rnorm(n, mean=3, sd=0.5)
  x5=rnorm(n, mean=1.8, sd=0.9)
  x6=rnorm(n, mean=2.7, sd=1.5)
  x7_8=mvrnorm(n, mu=c(2,2),Sigma=matrix(c(1,0.3,0.3,1),nrow=2,byrow=TRUE))
  x9_10=mvrnorm(n, mu=c(-3,2),Sigma=matrix(c(1,0.1,0.1,1),nrow=2,byrow=TRUE))
  sim=data.frame(label,allval,x3,x4,x5,x6,x7_8,x9_10,z1,z2,z3,z4,z5)
  order=sample(n,n)
  newlabel=sim[order,1]
  sim1=sim[order,-1]
  results=list(sim1,newlabel)
  return(results)
}

## prepare simulation
sim1=simulation1_1(n=600)
X=sim1[[1]]; truelabel=sim1[[2]]
c=10;b=1;m=4;s=0
Y=X[,-c(1:(c+b))]; 
cat=c(3,3,4,10)
Y1=model.matrix(~ . + 0, data=Y, contrasts.arg = lapply(Y, contrasts, contrasts=FALSE))
X1=data.frame(X[,1:(c+b)],Y1)
X1=as.matrix(X1)
xx=X1


###2. test on improved LCA using multivariate normal No Feature Selection
xx=X1;xx_new=xx;K=2
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

##3. Proposed alogrithm ECSM: test on improved LCA using multivariate normal After Feature Selection 
xx=X1
test3_1=ECSM(xx,obj,c=10,b=1,m=4,s=0,cat=c(3,3,4,10),K=2,maxiter=400,truelabel,featurenum=4,threshold1=8e-5, threshold2=0.04)
test3_1


#save(test3_1,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/test3_1.RData")

# xx=X1
# test3_2=ECSM(xx,c=10,b=1,m=4,s=0,cat=c(3,3,4,10),K=2,maxiter=400,truelabel,featurenum=4,threshold1=8e-5, threshold2=0.04)
# test3_2
#save(test3_2,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/test3_2.RData")

##RI plot
# load("/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/test3_2.RData")
# iter=25; delta_matrix=test3_2[[1]]
# iterations=1:(iter)
# delta=delta_matrix[2:(iter+1),]
# delta_v=matrix(delta,ncol=1)
# iter_v=rep(iterations, 15)
# variable=c(rep("f1",iter),rep("f2",iter),rep("f3",iter),rep("f4",iter),rep("f5",iter),rep("f6",iter),
#            rep("f7",iter),rep("f8",iter),rep("f9",iter),rep("f10",iter),rep("f11",iter),rep("f12",iter),
#            rep("f13",iter),rep("f14",iter),rep("f15",iter))
# feature=factor(variable, levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15"))
# Delta=data.frame(iter_v, delta_v, feature)
# colnames(Delta)=c("iterations","FI","feature")
# figure1=ggplot(Delta,aes(x=iterations, y=FI, shape=feature))+geom_point(size=1)+
#   geom_line(size=0.3)+
#   scale_x_discrete(name="iterations", breaks=seq(from=0,to=iter,by=5))+
#   theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
#   scale_shape_manual(values=1:nlevels(Delta$feature),labels=c(expression(f[1]),expression(f[2]),expression(f[3]),expression(f[4]),expression(f[5]),
#                                                                expression(f[6]),expression(f[7]),expression(f[8]),expression(f[9]),expression(f[10]),
#                                                                expression(f[11]),expression(f[12]),expression(f[13]),expression(f[14]),expression(f[15])))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   annotate("text",x=26,y=0.078,label="f[1]",parse=TRUE,size=4)+
#   annotate("text",x=26,y=0.030,label="f[2]",parse=TRUE,size=4)+
#   annotate("text",x=26,y=0.022,label="f[11]",parse=TRUE,size=4)+
#   annotate("text",x=26,y=0.015,label="f[12]",parse=TRUE,size=4)+
#   annotate("text",x=21.5,y=-0.004,label="others(f[3]-f[10],f[13]-f[15])",parse=TRUE,size=4)
# figure1
# tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure1_1.tiff', units="in", width=6, height=4, res=300)
# figure1
# dev.off()
