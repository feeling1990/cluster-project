
###Benchmark data 'Crab'
data("crabs",package="MASS")
X=crabs[,4:8]
Class=with(crabs,paste(sp,sex,sep="|"))
max_feature_num = 4
num_class = 4
truelabel = Class
table(Class)
results_crab = benchmark_compare(X,Class,num_class = 4,max_feature_num = 3)
result1 = results_crab[[1]]
result1$ACC = round(1- result1$ErrorRate,4)*100
figure6=ggplot(data=result1,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on Crab Data")
figure6

result1$Time = round(result1$Time,3)
figure7=ggplot(data=result1,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,6.5))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on Crab Data")
figure7

##save the plot
tiff('/Users/yinlin/Github/cluster-project/Phase_I/Figures/crab.tiff', units="in", width=6, height=8, res=300)
multiplot(figure6, figure7, cols=1)
dev.off()

##################################################################################

## Data 2 'Wine'
data("wine",package='SelvarMix')
summary(wine)
X=wine[,1:27]
Class = wine[,28]
results_wine = benchmark_compare(X,Class,num_class = 3,max_feature_num = 8)
result2 = results_wine[[1]]
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
figure8

result2$Time = round(result2$Time,3)
figure9=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on Wine Data")
figure9

## save the plots
# tiff('/Users/yinlin/Github/cluster-project/Phase_I/Figures/wine.tiff', units="in", width=6, height=8, res=300)
# multiplot(figure8, figure9, cols=1)
# dev.off()