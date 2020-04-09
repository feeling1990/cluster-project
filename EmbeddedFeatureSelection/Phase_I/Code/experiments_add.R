library(foreign)
###Benchmark data 'Crab'
set.seed(123)
data("crabs",package="MASS")
X=crabs[,4:8]
Class=with(crabs,paste(sp,sex,sep="|"))
max_feature_num = 4
num_class = 4
truelabel = Class
table(Class)
results_crab = benchmark_compare(X,Class,num_class = 4,max_feature_num = 3)
result1 = results_crab[[1]]
result1$ACC = round(1- result1$ErrorRate,3)*100
result1$ARI = round(result1$ARI,3)
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
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/crab.tiff', units="in", width=6, height=8, res=300)
multiplot(figure6, figure7, cols=1)
dev.off()

##################################################################################

## Data 2 'Wine'
data("wine",package='SelvarMix')
summary(wine)
X=wine[,1:27]
Class = wine[,28]
results_wine = benchmark_compare(X,Class,num_class = 3,max_feature_num = 20)
result2 = results_wine[[1]]
result2$ACC = round(1- result2$ErrorRate,3)
result2$ARI = round(result2$ARI,3)
figure8=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,1))+
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
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/wine.tiff', units="in", width=6, height=8, res=300)
multiplot(figure8, figure9, cols=1)
dev.off()


## Data 3 "tetragonula"
# data = read.csv("/Volumes/primary/Research/cluster-project/Tetragonula.csv")
# summary(data)
# X = data[,1:ncol(data)-1]
# Class = data[,ncol(data)]
# table(Class)
# dim(X)
# G = 9
# results_tet = benchmark_compare(X,Class,num_class = 9,max_feature_num = 10)
# obj = Mclust(X,G=10) # a more robust way of initialization
# start_time <- Sys.time()
# myresults2=ESM_fast(X,obj,clusternum=10,maxiter=100,truelabel=Class,max_feature_num=10,threshold1=0.001,threshold2=0.05)
# end_time <- Sys.time()
# end_time - start_time
# myresults2[[3]]

## Data 4 "G2"
x = read.table("/Volumes/primary/Research/cluster-project/g2-txt/g2-128-40.txt")
y = read.table("/Volumes/primary/Research/cluster-project/g2-gt-pa/g2-128-40-gt.pa", head=FALSE, skip =4)
X=x
Class = y[,1]
table(Class)
dim(X)
max_feature_num = 64
num_class = 2
results_g2_128_40 = benchmark_compare(X,Class,num_class =2,max_feature_num = 40, include_model2=FALSE)# good keep
#results_g2_128_40_full = benchmark_compare(X,Class,num_class =2,max_feature_num = 40, include_model2=TRUE)
result2 = results_g2_128_40[[1]]
result2$ACC = round(1- result2$ErrorRate,4)*100
result2$ARI = round(result2$ARI,3)
figure10=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on g2 Data")

result2$Time = round(result2$Time,3)
figure11=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on G2 Data")

## save the plots
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/g2.tiff', units="in", width=6, height=8, res=300)
multiplot(figure10, figure11, cols=1)
dev.off()

## iris
data = read.arff("/Volumes/primary/Research/cluster-project/clustering-benchmark/src/main/resources/datasets/real-world/iris.arff")
summary(data)
X = data[,1:ncol(data)-1]
Class = data[, ncol(data)]
table(Class)
max_feature_num = 4
G = 3
results_iris = benchmark_compare(X,Class=Class,num_class =G,
                                      max_feature_num = max_feature_num, include_model2=TRUE)# good keep

result2 = results_iris[[1]]
result2$ACC = round(1- result2$ErrorRate,4)*100
result2$ARI = round(result2$ARI,3)
figure10=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on iris Data")

result2$Time = round(result2$Time,3)
figure11=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on iris Data")


## save the plots
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/iris.tiff', units="in", width=6, height=8, res=300)
multiplot(figure10, figure11, cols=1)
dev.off()


##thy
data = read.arff("/Volumes/primary/Research/cluster-project/clustering-benchmark/src/main/resources/datasets/real-world/thy.arff")
summary(data)
X = data[,1:ncol(data)-1]
Class = data[, ncol(data)]
table(Class)
ncol(X)
max_feature_num = 4
G = 3
results_thy = benchmark_compare(X,Class=Class,num_class =G,
                             max_feature_num = max_feature_num, include_model2=TRUE)# good keep
result2 = results_thy[[1]]
result2$ACC = round(1- result2$ErrorRate,4)*100
result2$ARI = round(result2$ARI,3)

figure10=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on thy Data")

result2$Time = round(result2$Time,3)
figure11=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on thy Data")
## save the plots
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/thy.tiff', units="in", width=6, height=8, res=300)
multiplot(figure10, figure11, cols=1)
dev.off()

##sonar results not very good
data = read.arff("/Volumes/primary/Research/cluster-project/clustering-benchmark/src/main/resources/datasets/real-world/sonar.arff")
summary(data)
X = data[,1:ncol(data)-1]
Class = data[, ncol(data)]
table(Class)
dim(X)
max_feature_num = 45
G = 2
results_soar = benchmark_compare(X,Class=Class,num_class =G,
                            max_feature_num = max_feature_num, include_model2=TRUE)
result2 = results_soar[[1]]
result2$ACC = round(1- result2$ErrorRate,4)*100
result2$ARI = round(result2$ARI,3)
figure10=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on Sonar Data")

result2$Time = round(result2$Time,3)
figure11=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on Sonar Data")
## save the plots
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/soar.tiff', units="in", width=6, height=8, res=300)
multiplot(figure10, figure11, cols=1)
dev.off()

##vowel results ok
data = read.arff("/Volumes/primary/Research/cluster-project/clustering-benchmark/src/main/resources/datasets/real-world/vowel.arff")
summary(data)
X = data[,5:ncol(data)-1]
Class = data[, ncol(data)]
table(Class)
dim(X)
max_feature_num = 9
G = 11
results_vowel = benchmark_compare(X,Class=Class,num_class =G,
                            max_feature_num = max_feature_num, include_model2=TRUE)
result2 = results_vowel[[1]]
result2$ACC = round(1- result2$ErrorRate,3)
result2$ARI = round(result2$ARI,3)
figure10=ggplot(data=result2,aes(x=algorithm,y=ACC))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,100))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(ACC,'%')), vjust=-0.3, size=3.5)+
  ggtitle("Accuracy on Vowel Data")

result2$Time = round(result2$Time,3)
figure11=ggplot(data=result2,aes(x=algorithm,y=Time))+geom_bar(stat="identity",width=0.5)+
  scale_y_continuous()+coord_cartesian(ylim=c(0,71))+
  theme(axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  #guides(shape=guide_legend(override.aes=list(color=override.color)))+
  geom_text(aes(label=paste0(Time)), vjust=-0.3, size=3.5)+
  ggtitle("Running Time on Vowel Data")
## save the plots
tiff('/Volumes/primary/Research/cluster-project/Phase_I/Figures/vowel.tiff', units="in", width=6, height=8, res=300)
multiplot(figure10, figure11, cols=1)
dev.off()