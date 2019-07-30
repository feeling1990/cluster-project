load_mnist <- function() {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <<- load_image_file('/Volumes/primary/Research/Clustering_v1/VAE-Tensorflow/MNIST_data/train-images-idx3-ubyte')
  test <<- load_image_file('/Volumes/primary/Research/Clustering_v1/VAE-Tensorflow/MNIST_data/t10k-images-idx3-ubyte')
  
  train$y <<- load_label_file('/Volumes/primary/Research/Clustering_v1/VAE-Tensorflow/MNIST_data/train-labels-idx1-ubyte')
  test$y <<- load_label_file('/Volumes/primary/Research/Clustering_v1/VAE-Tensorflow/MNIST_data/t10k-labels-idx1-ubyte')  
  list(train=train, test=test) 
}

data= load_mnist()

## Run ESM on MNIST data
sample_set = sample(60000,1000)
X = data$train$x[sample_set,]
y = data$train$y[sample_set]

pca = prcomp(X,retx=TRUE)
X_pca = pca$x[,1:50]
mod1=Mclust(X_pca,G=10)
Class = y
max_feature_num = 30
num_class = 10
X = X_pca

print(Error_4)
benchmark_compare=function(X,Class,num_class,max_feature_num){
  #initial model with full features 
  mod1=Mclust(X,G=num_class)
  ARI_1=ARI(mod1$classification,Class)
  Error_1=classError(Class,mod1$classification)$errorRate
  print("Finished model1")
  # Model 2
  start_time <- Sys.time()
  out=clustvarsel(X,G=num_class)
  mod2=Mclust(X[,out$subset],G=num_class)
  end_time <- Sys.time()
  time2=end_time - start_time
  ARI_2=ARI(Class,mod2$classification)
  Error_2=classError(Class,mod2$classification)$errorRate
  print("Finished model2")
  # Model 3
  start_time <- Sys.time()
  mod3=VarSelCluster(X,gvals=num_class,nbcores=2,initModel=100,crit.varsel="BIC")
  end_time <- Sys.time()
  time3=end_time - start_time
  ARI_3=ARI(Class, fitted(mod3))
  Error_3=classError(Class,fitted(mod3))$errorRate
  print("Finished model3")
  # Model 4
  start_time <- Sys.time()
  mod4 = vscc(X,G=num_class)
  end_time <- Sys.time()
  time4=end_time - start_time
  ARI_4=ARI(Class,mod4$bestmodel$classification)
  Error_4=classError(Class,mod4$bestmodel$classification)$errorRate
  print("Finished model4")
  # Model 5
  start_time <- Sys.time()
  mod5 = SelvarClustLasso(X,nbcluster=num_class,nbcores=4)
  end_time <- Sys.time()
  time5=end_time - start_time
  ARI_5=ARI(Class, mod5$partition)
  Error_5=classError(Class,mod5$partition)$errorRate
  print("Finished model5")
  # ESM
  start_time <- Sys.time()
  mod6 = ESM(X,mod1,clusternum=num_class,maxiter=100,truelabel=Class,featurenum=feature_num,threshold1=0.01,threshold2=0.05)
  end_time <- Sys.time()
  time6 = end_time - start_time
  ARI_6=ARI(mod6[[1]], Class)
  Error_6=classError(Class,mod6[[1]])$errorRate
  print("Finished model6")
  
  algorithm=c("clustvarsel","VarSelCluster","vscc","selvarclustLasso","Proposed")
  algorithm=factor(algorithm,levels=c("clustvarsel","VarSelCluster","vscc","selvarclustLasso","Proposed"))
  ARI=c(ARI_2,ARI_3, ARI_4,ARI_5, ARI_6)
  ErrorRate=c(Error_2,Error_3,Error_4,Error_5,Error_6)
  Time = c(time2,time3,time4,time5,time6)
  benchmark_results=data.frame(algorithm,ARI,ErrorRate,Time)
  return (benchmark_results)
}

result_mnist = benchmark_compare(X_pca,y,num_class = 10,feature_num = 78)


## other possible datasets
data(scenarioCor)












show_digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}