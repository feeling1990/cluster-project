###Benchmark data 'Crab'
benchmark_compare=function(data_X,Class,num_class){
  #initial model with full features 
  mod1=Mclust(X,G=1:5)
  ARI_1=ARI(mod1$classification,Class)
  Error_1=classError(Class,mod1$classification)$errorRate
  # Model 2
  start_time <- Sys.time()
  out=clustvarsel(X,G=1:5)
  mod2=Mclust(X[,out$subset],G=1:5)
  end_time <- Sys.time()
  time2=end_time - start_time
  ARI_2=ARI(Class,mod2$classification)
  Error_2=classError(Class,mod2$classification)$errorRate
  # Model 3
  start_time <- Sys.time()
  mod3=VarSelCluster(X,gvals=num_class,nbcores=2,initModel=1000,crit.varsel="BIC")
  end_time <- Sys.time()
  time3=end_time - start_time
  ARI_3=ARI(Class, fitted(mod3))
  Error_3=classError(Class,fitted(mod3))$errorRate
  # Model 4
  start_time <- Sys.time()
  mod4 = vscc(X,G=1:5)
  end_time <- Sys.time()
  time4=end_time - start_time
  ARI_4=ARI(Class,mod4$bestmodel$classification)
  Error_4=classError(Class,mod4$bestmodel$classification)$errorRate
  # Model 5
  start_time <- Sys.time()
  mod5 = SelvarClustLasso(X,nbcluster=1:4,nbcores=4)
  end_time <- Sys.time()
  time5=end_time - start_time
  ARI_5=ARI(Class, mod5$partition)
  Error_5=classError(Class,mod5$partition)$errorRate
  # ESM
  start_time <- Sys.time()
  mod6 = ESM_v2(X,mod1,clusternum=num_class,maxiter=100,truelabel=Class,featurenum=4,threshold1=0.01,threshold2=0.05)
  end_time <- Sys.time()
  time6 = end_time - start_time
  ARI_6=ARI(mod6[[1]], Class)
  Error_6=classError(Class,mod6[[1]])$errorRate
  algorithm=c("clustvarsel","VarSelCluster","vscc","selvarclustLasso","Proposed")
  algorithm=factor(algorithm,levels=c("clustvarsel","VarSelCluster","vscc","selvarclustLasso","Proposed"))
  ARI=c(ARI_2,ARI_3, ARI_4,ARI_5, ARI_6)
  ErrorRate=c(Error_2,Error_3,Error_4,Error_5,Error_6)
  Time = c(time2,time3,time4,time5,time6)
  benchmark_results=data.frame(algorithm,ARI,ErrorRate,Time)
  return (benchmark_results)
}


data("crabs",package="MASS")
X=crabs[,4:8]
Class=with(crabs,paste(sp,sex,sep="|"))
table(Class)
result1 = benchmark_compare(X,Class,num_class = 4)
result1$ACC = round(1- result1$ErrorRate,4)*100
figure6=ggplot(data=result1,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on Crab Data")

result1$Time = round(result1$Time,3)
figure7=ggplot(data=result1,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,6.5))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on Crab Data")

tiff('/Users/yinlin/Desktop/ESM phase I/ESM paper/figures/crab.tiff', units="in", width=6, height=8, res=300)
multiplot(figure6, figure7, cols=1)
dev.off()

##################################################################################

## Data 2 'Wine'
data("wine",package='SelvarMix')
summary(wine)
X=wine[,1:27]
Class = wine[,28]
result2 = benchmark_compare(X,Class,num_class = 3)
result2$ACC = round(1- result2$ErrorRate,4)*100
result2$ARI = round(result2$ARI,3)
figure8=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on Wine Data")

result2$Time = round(result2$Time,3)
figure9=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on Wine Data")


tiff('/Users/yinlin/Desktop/ESM phase I/ESM paper/figures/wine.tiff', units="in", width=6, height=8, res=300)
multiplot(figure8, figure9, cols=1)
dev.off()