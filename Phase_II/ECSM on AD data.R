x=read.csv("/Users/yinlinfu/Desktop/benchmark data/AD_removemissing.csv")
labelv1=x[,1];labelv2=x[,2]
APOE=x$APOE.Status
Gender=as.factor(x$PTGENDER); levels(Gender)=c("0","1")
AV45=as.factor(x$AV45.v1.Pos118);levels(AV45)=c("0","1")
Convert=x[,3]
Convert=factor(Convert,levels=c("Non-converter","Conv: NC-MCI","Conv: NC-AD","Conv: MCI-AD"))
xx=as.matrix(x[,c(4:15)])
AD_v=cbind(labelv1,Convert,xx,Gender,AV45,APOE)
colnames(AD_v)=c("label","Convert","Age","MMSE","CDR","Hippocampus","Ventricles","WholeBrain","Entorhinal",
                 "ICV","HCI","sROI","mcSUVRcere","mcSUVRwm","Gender","AV45","APOE")

label_num=numeric()
for(i in 1:length(labelv1)){
  if(labelv1[i]=="AD"){label_num[i]=3}
  else if(labelv1[i]=="MCI"){label_num[i]=2}
  else label_num[i]=1
}

cate=data.frame(Gender,AV45,APOE)           
x_cate=model.matrix(~ . + 0, data=cate, contrasts.arg = lapply(cate, contrasts, contrasts=FALSE))
AD_v=cbind(labelv1,Convert,xx,x_cate)
clusternum=3; maxiter=300; truelabel=labelv1; featurenum=5; threshold1=1e-5;threshold2=0.005
K=3

##test on continuous+gender+AV45+APOE
c=12;b=0;m=3;s=0;cat=c(2,2,3);K=3;threshold1=1e-5;threshold2=0.01
xx=AD_v[,-c(1,2)];xx_new=xx
x=xx[,1:(c+b)]
obj=init.EM(x,nclass=K)
gammaKn = e.step(x,obj)
df <- data.frame(matrix(unlist(gammaKn), ncol=K, byrow=F),stringsAsFactors=FALSE)
gammaKn=t(df)
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
for(h in 1:20){
  v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
  gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]); 
}
v <- Mstep(c,b,m,s,cat,xx_new, gammaKn);
gammaKn = Estep(c,b,m,s,cat,xx_new,K, v[[1]], v[[2]], v[[3]],v[[4]],v[[5]],v[[6]],v[[7]]);  
mylabel=classassign(t(gammaKn))
confusion=table(truelabel,mylabel)
confusion
myaccuracy(confusion)
save(obj,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/AD_obj.RData")
load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/AD_obj.RData")


####
results2=ECSM(xx=AD_v[,-c(1,2)],obj,c=12,b=0,m=3,s=0,cat=c(2,2,3),K=3,maxiter=400,truelabel=labelv1,featurenum=3,threshold1=1e-5, threshold2=0.1)
save(results2,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/results2 for AD.RData")



results3=ECSM(xx=AD_v[,-c(1,2)],obj,c=12,b=0,m=3,s=0,cat=c(2,2,3),K=3,maxiter=400,truelabel=labelv1,featurenum=3,threshold1=1e-5, threshold2=0.05)
save(results3,file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/results3 for AD.RData")


##RI plot
load(file="/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/results2 for AD.RData")

iter=35; delta_matrix=results2[[1]]
#i=35; v=results3[[8]];x=AD_v[,-c(1,2)];xx=x[,c(3,9,11)];xx_new=xx;c=3;b=0;m=0
iterations=1:(iter)
delta=delta_matrix[2:(iter+1),]
delta_v=matrix(delta,ncol=1)
iter_v=rep(iterations, 15)
variable=c(rep("f1",iter),rep("f2",iter),rep("f3",iter),rep("f4",iter),rep("f5",iter),rep("f6",iter),
           rep("f7",iter),rep("f8",iter),rep("f9",iter),rep("f10",iter),rep("f11",iter),rep("f12",iter),
           rep("f13",iter),rep("f14",iter),rep("f15",iter))
feature=factor(variable, levels=c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13","f14","f15"))
Delta=data.frame(iter_v, delta_v, feature)
colnames(Delta)=c("iterations","FI","feature")

figure3=ggplot(Delta,aes(x=iterations, y=FI,shape=feature))+geom_point(size=1)+
  geom_line(size=0.3)+
  scale_x_discrete(name="iterations", breaks=seq(from=0,to=iter+5,by=5))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  scale_shape_manual(values=1:nlevels(Delta$feature),labels=c(expression(f[1]),expression(f[2]),expression(f[3]),expression(f[4]),expression(f[5]),
                                                              expression(f[6]),expression(f[7]),expression(f[8]),expression(f[9]),expression(f[10]),
                                                              expression(f[11]),expression(f[12]),expression(f[13]),expression(f[14]),expression(f[15])))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  annotate("text",x=36,y=0.095,label="f[3]",parse=TRUE,size=4)+
  annotate("text",x=36.5,y=0.075,label="f[11]",parse=TRUE,size=4)+
  annotate("text",x=36,y=0.045,label="f[9]",parse=TRUE,size=4)+
  annotate("text",x=29,y=-0.005,label="others(f[1],f[2],f[4]-f[8],f[10],f[12]-f[15])",parse=TRUE,size=3)
figure3
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure3_3.tiff', units="in", width=6, height=4, res=300)
figure3
dev.off()

png('/Users/yinlinfu/Desktop/figures/figure4.png', units="in", width=5, height=4, res=300)
figure4
dev.off()

pdf('/Users/yinlinfu/Desktop/figures/figure4.pdf', width=6, height=4)
figure4
dev.off()

### AD visulization
override.shape <- c(16, 17, 16, 17, 16)
override.linetype <- c(1, 1, 3, 3, 4)

g <- ggplot(df.merged, aes(x, y, colour = interaction(Type, Method), linetype = Method, shape = Type)) + geom_line() + geom_point()
g <- g + guides(colour = guide_legend(override.aes = list(shape = override.shape, linetype = override.linetype)))
g <- g + scale_shape(guide = FALSE)
g <- g + scale_linetype(guide = FALSE)
print(g)

override.color <- c("LINE1"="#f04546","LINE2"="#3591d1","BAR"="#62c76b")
cols <- c("LINE1"="#f04546","LINE2"="#3591d1","BAR"="#62c76b")

figure5=ggplot(AD_v,aes(x=mcSUVRcere,y=CDR,color=label, shape=Convert,size=Convert))+geom_point()+
  scale_alpha_discrete(range=c(0.6,1,1,1))+scale_size_discrete(range=c(2,5,5,4))+
  xlab(expression(mcSUVRcere(f[11])))+ylab(expression(CDR(f[3])))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  guides(shape=guide_legend(override.aes=list(color=override.color)))
figure5
tiff('/Users/yinlinfu/Dropbox/ESM Phase II/paper/figure/figure5.tiff', units="in", width=6, height=4, res=300)
figure5
dev.off()

pdf('/Users/yinlinfu/Desktop/figures/figure5.pdf',width=7,height=4)
figure5
dev.off()

summary(Convert)
##statistical summary of AD data
summary(xx)
mean=apply(xx,2,mean)
std=apply(xx,2,sd)
