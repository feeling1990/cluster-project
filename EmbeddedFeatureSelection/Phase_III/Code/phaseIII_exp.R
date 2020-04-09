#install.packages('tsne')
library(tsne)

## simulation data 1
set.seed(123)
sim1=simulation1(300)
X = sim1[[1]]
truelabel = as.factor(sim1[[2]])

obj = Mclust(X,G=2) # a more robust way of initialization
start_time <- Sys.time()
myresults1=ESM_fast_confidence(X,obj,clusternum=2,maxiter=100,truelabel=truelabel,max_feature_num=3,threshold1=0.001,threshold2=0.06)
end_time <- Sys.time()
end_time - start_time
myresults1

# generate plot
misclassified = myresults1[['misclassified']]
outlier_candidates = myresults1[['outliers']]
frequency = table(outlier_candidates)
outliers = as.numeric(names(which(frequency>2)))
frequency
misclassified
intersect(outliers, misclassified)
intersect(outlier_candidates, misclassified)


data = X[,myresults1[['keep']]] #use selected features
n = length(truelabel)
type = rep('normal', n)
type[misclassified] = 'misclassified' 
type[outliers] = 'outliers_delta'
outlier_name = rep(0,n)
outlier_name[outliers] = outliers
data_plot = as.data.frame(cbind(data,truelabel,as.factor(type),outlier_name))
colnames(data_plot) = c('x1','x2', 'truelabel','type')
ggplot(data_plot, aes(x=x1, y=x2,color=type,  shape=truelabel)) +
  geom_point()+
  ggtitle('Detect Outliers by Feature Selection')+
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label=ifelse(outlier_name>0,as.character(outlier_name),'')),hjust=0,vjust=0,size=2.5)


#+===================================simulation data2==================================
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
  #add extra outlier
  # sim2 = rbind(sim2,c(-1,2,rnorm(13, 1.1,1)))
  # sim2 = rbind(sim2,c(2,2,rnorm(13, 0.4,0.5)))
  # newlabel = append(newlabel,c(2,2))
  results=list(sim2,newlabel)
  return(results)
}

sim2 = simulation2(n=600)
X = sim2[[1]]
truelabel = as.factor(sim2[[2]])
n = length(truelabel)

obj = Mclust(X,G=2) # a more robust way of initialization
start_time <- Sys.time()
myresults2=ESM_fast_confidence(X,obj,clusternum=2,maxiter=100,truelabel=truelabel,max_feature_num=3,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time
#myresults2

# generate plot
misclassified = myresults2[['misclassified']]
outlier_candidates = myresults2[['outliers']]
frequency = table(outlier_candidates)
outliers = as.numeric(names(which(frequency>1)))
frequency
misclassified
intersect(outliers, misclassified)
intersect(outlier_candidates, misclassified)


data = X[,myresults2[['keep']]] #use selected features
type = rep('normal', n)
type[misclassified] = 'misclassified' 
type[outliers] = 'outliers_delta'
outlier_name = rep(0,n)
outlier_name[outliers] = outliers
data_plot = as.data.frame(cbind(data,truelabel,as.factor(type),outlier_name))
colnames(data_plot) = c('x1','x2', 'truelabel','type')
ggplot(data_plot, aes(x=x1, y=x2,color=type,  shape=truelabel)) +
  geom_point()+
  ggtitle('Detect Outliers by Feature Selection')+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label=ifelse(outlier_name>0,as.character(outlier_name),'')),hjust=0,vjust=0,size=2.5)



## =======================crab data==============================
data("crabs",package="MASS")
X=crabs[,4:8]
Class=with(crabs,paste(sp,sex,sep="|"))
table(Class)

max_feature_num = 4
num_class = 4
truelabel = Class
obj = Mclust(X,G=num_class) # a more robust way of initialization
start_time <- Sys.time()
myresults3=ESM_fast_confidence(X,obj,clusternum=num_class,maxiter=100,truelabel=truelabel,
                               max_feature_num=max_feature_num,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time
myresults3

#compare plot: use full features vs. use selected features
tsne_crab_full = tsne(X) 
data <- data.frame(apply(tsne_crab_full, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel))
colnames(data_plot) = c('x1','x2', 'truelabel')
p1 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel)) +
  geom_point()+
  ggtitle('t-SNE on full features')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

tsne_crab_selected = tsne(X[,myresults3$keep])
data <- data.frame(apply(tsne_crab_selected, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel))
colnames(data_plot) = c('x1','x2', 'truelabel')
p2 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel)) +
  geom_point()+
  ggtitle('t-SNE on selected features')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

tiff('/Volumes/primary/Research/cluster-project/Phase_III/Figures/t-SNE_crab.tiff', 
     units="in", width=10, height=4.5, res=300)
multiplot(p1, p2, cols=2)
dev.off()

####### compare detected outliers with misclassified ================
truelabel = as.character(truelabel)
n = length(truelabel)
misclassified = myresults3[['misclassified']]
outlier_candidates = myresults3[['outliers']]
frequency = table(outlier_candidates)
outliers = as.numeric(names(which(frequency>1)))
frequency
misclassified
intersect(outliers, misclassified)
intersect(outlier_candidates, misclassified)

group_outlier = rep('normal', n)
group_outlier[outliers] = 'outliers'
outlier_index = rep(0,n)
outlier_index[outliers] = outliers
group_misclassified = rep('correct',n)
group_misclassified[misclassified] = 'misclassified'
misclassified_index = rep(0,n)
misclassified_index[misclassified] = misclassified

p3 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel,size=group_outlier)) +
  geom_point()+
  ggtitle('detected outlier')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(outlier_index>0,as.character(outlier_index),'')),hjust=0.9,vjust=-1,size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3

p4 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel,size=group_misclassified)) +
  geom_point()+
  ggtitle('misclassified')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(misclassified_index>0,as.character(misclassified_index),'')),hjust=0.9,vjust=-1,size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4

tiff('/Volumes/primary/Research/cluster-project/Phase_III/Figures/crab_outlier_misclassified.tiff', 
     units="in", width=10, height=4.9, res=300)
multiplot(p3, p4, cols=2)
dev.off()

#==================data Wine====================================================
data("wine",package='SelvarMix')
summary(wine)
X=wine[,1:27]
Class = wine[,28]
max_feature_num = 20
num_class = 3
truelabel = Class
table(Class)

obj = Mclust(X,G=num_class) # a more robust way of initialization
start_time <- Sys.time()
myresults4=ESM_fast_confidence(X,obj,clusternum=num_class,maxiter=100,truelabel=truelabel,max_feature_num=max_feature_num,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time
myresults4

#compare plot: use full features vs. use selected features
truelabel = as.character(truelabel)
n = length(truelabel)
tsne_full = tsne(X) 
data <- data.frame(apply(tsne_full, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel))
colnames(data_plot) = c('x1','x2', 'truelabel')
p1 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel)) +
  geom_point()+
  ggtitle('t-SNE on full features')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1
tsne_selected = tsne(X[,myresults4$keep])
data <- data.frame(apply(tsne_selected, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel))
colnames(data_plot) = c('x1','x2', 'truelabel')
p2 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel)) +
  geom_point()+
  ggtitle('t-SNE on selected features')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2
tiff('/Volumes/primary/Research/cluster-project/Phase_III/Figures/t-SNE_wine_v3.tiff', 
     units="in", width=10, height=4.5, res=300)
multiplot(p1, p2, cols=2)
dev.off()

# outliers 
truelabel = as.character(truelabel)
misclassified = myresults4[['misclassified']]
outlier_candidates = myresults4[['outliers']]
frequency = table(outlier_candidates)
outliers = as.numeric(names(which(frequency>50)))
frequency
misclassified
intersect(outliers, misclassified)
intersect(outlier_candidates, misclassified)

group_outlier = rep('normal', n)
group_outlier[outliers] = 'outliers'
outlier_index = rep(0,n)
outlier_index[outliers] = outliers
group_misclassified = rep('correct',n)
group_misclassified[misclassified] = 'misclassified'
misclassified_index = rep(0,n)
misclassified_index[misclassified] = misclassified

data <- data.frame(apply(tsne_full, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel,as.factor(group_outlier),outlier_index,group_misclassified, misclassified_index))
colnames(data_plot) = c('x1','x2', 'truelabel','group_outlier', 'outlier_index', 'group_misclassified', 'misclassified_index')
p3 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel,size=group_outlier)) +
  geom_point()+
  ggtitle('detected outlier')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(outlier_index>0,as.character(outlier_index),'')),hjust=0.9,vjust=-1,size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3

p4 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel,size=group_misclassified)) +
  geom_point()+
  ggtitle('misclassified')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(misclassified_index>0,as.character(misclassified_index),'')),hjust=0.9,vjust=-1,size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4

tiff('/Volumes/primary/Research/cluster-project/Phase_III/Figures/wine_outlier_misclassified_v3.tiff', 
     units="in", width=10, height=4.9, res=300)
multiplot(p3, p4, cols=2)
dev.off()

##=================Vowel================================
data = read.arff("/Volumes/primary/Research/cluster-project/datasets/vowel.arff")
summary(data)
X = data[,5:ncol(data)-1]
Class = data[, ncol(data)]
table(Class)
dim(X)
max_feature_num = 9
num_class = 11
truelabel = as.character(Class)

obj = Mclust(X,G=num_class) # a more robust way of initialization
start_time <- Sys.time()
myresults5=ESM_fast_confidence(X,obj,clusternum=num_class,maxiter=100,truelabel=truelabel,max_feature_num=max_feature_num,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time
myresults5

#compare plot: use full features vs. use selected features
truelabel = as.character(truelabel)
n = length(truelabel)
tsne_full = tsne(X) 
data <- data.frame(apply(tsne_full, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel))
colnames(data_plot) = c('x1','x2', 'truelabel')
p1 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel)) +
  geom_point()+
  ggtitle('t-SNE on full features')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p1
tsne_selected = tsne(X[,myresults5$keep])
data <- data.frame(apply(tsne_selected, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel))
colnames(data_plot) = c('x1','x2', 'truelabel')
p2 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel)) +
  geom_point()+
  ggtitle('t-SNE on selected features')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p2
tiff('/Volumes/primary/Research/cluster-project/Phase_III/Figures/t-SNE_vowel.tiff', 
     units="in", width=10, height=4.5, res=300)
multiplot(p1, p2, cols=2)
dev.off()

# outliers 
truelabel = as.character(truelabel)
misclassified = myresults5[['misclassified']]
outlier_candidates = myresults5[['outliers']]
frequency = table(outlier_candidates)
outliers = as.numeric(names(which(frequency>1)))
frequency
misclassified
intersect(outliers, misclassified)
intersect(outlier_candidates, misclassified)

group_outlier = rep('normal', n)
group_outlier[outliers] = 'outliers'
outlier_index = rep(0,n)
outlier_index[outliers] = outliers
group_misclassified = rep('correct',n)
group_misclassified[misclassified] = 'misclassified'
misclassified_index = rep(0,n)
misclassified_index[misclassified] = misclassified

data <- data.frame(apply(tsne_selected, 2, as.numeric))
data_plot = as.data.frame(cbind(data,truelabel,as.factor(group_outlier),outlier_index,group_misclassified, misclassified_index))
colnames(data_plot) = c('x1','x2', 'truelabel','group_outlier', 'outlier_index', 'group_misclassified', 'misclassified_index')
p3 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel,size=group_outlier)) +
  geom_point()+
  ggtitle('detected outlier')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(outlier_index>0,as.character(outlier_index),'')),hjust=0.9,vjust=-1,size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p3

p4 = ggplot(data_plot, aes(x=x1, y=x2,color=truelabel,size=group_misclassified)) +
  geom_point()+
  ggtitle('misclassified')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none") +
  geom_text(aes(label=ifelse(misclassified_index>0,as.character(misclassified_index),'')),hjust=0.9,vjust=-1,size=4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
p4

tiff('/Volumes/primary/Research/cluster-project/Phase_III/Figures/vowel_outlier_misclassified.tiff', 
     units="in", width=10, height=4.9, res=300)
multiplot(p3, p4, cols=2)
dev.off()

