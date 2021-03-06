##EM feature selection algorithm 

## pre-install library packages
#install.packages("EMCluster","mvnfast","parallel","ggplot2","progress","matrixcalc","robustHD")
#install.packages("ggplot2","mclust","VarSelLCM","vscc","SelvarMix","clustvarsel")
library(matrixcalc)
library(EMCluster)
library(mvnfast)
library(ggplot2)
library(mclust)
library(VarSelLCM)
library(vscc)
library(SelvarMix)
library(clustvarsel)
library(progress)
library(robustHD)
## functions needed
myaccuracy = function(truelabel,mylabel){
  return(1-classError(mylabel,truelabel)$errorRate)
}

KL_dist=function(p,q){
  if(p==0){
    f=0} else if(q==0){
      f=1} else if(p!=0 & q!=0){
        f=p*(log(p)-log(q))} else{
          f=NULL
        }
  return(f)
}

Mdiff1=function(x,y){
  h=x
  for (i in 1:dim(x)[1]){
    for(j in 1:dim(x)[2]){
      h[i,j]=KL_dist(x[i,j],y[i,j])
    }
  }
  return(sum(h)/(dim(x)[1]*dim(x)[2]))
}

classassign=function(gammaKn){
  apply(gammaKn,1,which.max)
}


#parallel computing
#library(parallel)
#no_cores=detectCores()-1
#cl=makeCluster(no_cores)

##plot function##
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

adjust_gammaKn_function2 = function(xx_new,gammaKn){
  G = ncol(gammaKn)
  for (k in 1:G){
    if (sum(gammaKn[,k])==0){
      print("restart")
      obj=init.EM(xx_new,nclass=G)
      gammaKn_new = e.step(xx_new,obj)$Gamma
      return(gammaKn_new)
    }
  }
  return(gammaKn)
}

CalculateRI = function(data, parameters, esEst){
  p = ncol(data)
  N = nrow(data)
  pro = parameters$pro
  mu = parameters$mean
  sigma = parameters$variance$sigma
  G = length(pro)
  GammaKn = array(,dim=c(G,N,p))
  parameters_reduced = parameters
  for (j in 1:p){
    mu_reduced = mu[-j,]
    sigma_reduced = sigma[-j,-j,]
    parameters_reduced$mean = mu_reduced
    parameters_reduced$variance$sigma = sigma_reduced
    parameters_reduced$variance$cholsigma = parameters_reduced$variance$sigma
    G = dim(parameters_reduced$variance$sigma)[3]
    for(k in 1:G) {
      sigma_cur = parameters_reduced$variance$sigma[,,k]
      if(is.positive.definite(sigma_cur, tol=1e-8)){ sigma_chol = chol(sigma_cur)}
      else{
        sigma_nearPD = as.matrix(nearPD(sigma_cur)$mat)
        sigma_chol = chol(sigma_nearPD)
      }
      parameters_reduced$variance$cholsigma[,,k] = sigma_chol
    }
    data_reduced = data[,-j]
    temp1 = estep(modelName = "VVV", data = data_reduced, parameters = parameters_reduced)
    GammaKn[,,j] = t(temp1$z)
  }
  RI=rep(0,p)
  gammaKn = t(esEst$z)
  for(j in 1:p){
    RI[j]=Mdiff1(gammaKn,GammaKn[,,j])
  }
  return (RI)}

### ESM_fast
# obj = init.EM(xx,nclass=clusternum)
# gammaKn = e.step(xx,obj)$Gamma

ESM_fast=function(xx,obj,clusternum,maxiter,truelabel,max_feature_num,threshold1,threshold2){
  p=ncol(xx);
  N=nrow(xx);
  delete=NULL;
  keep=1:p;
  RI_matrix=matrix(0,maxiter,p)
 # cl=makeCluster(no_cores);
  loglik_cur=0
  gammaKn = obj$z
  parameters = obj$parameters
  loglik_cur = obj$loglik
  loglikvector <- loglik_cur
  xx_new=xx; 
  msEst = obj
  esEst = obj
  ###begin loop for E and M step
  pb <- progress_bar$new(total = 100)
  for (i in 2:maxiter){
    pb$tick()
    #gammaKn = t(esEst$z)
    RI = CalculateRI(xx_new, msEst$parameters, esEst)
    RI_matrix[i, keep] = RI
    index = which.min(RI)
    #print (abs(mean(RI_matrix[i,]-RI_matrix[i-1,])))
    #print (RI[index])
    if((abs(mean(RI_matrix[i,]-RI_matrix[i-1,]))<threshold1)&(RI[index]<threshold2)&(i>5)){
      delete = c(delete, keep[index])
      xx_new = xx_new[,-index]
      keep = keep[-index]
    }
    # msEst = mstep(modelName="VVV", data=xx_new, z=esEst$z)
    # esEst = estep(modelName="VVV", data=xx_new, parameters=msEst$parameters)
    # 
    msEst = mstep(modelName="VVV", data=xx_new, z=esEst$z,warnings=TRUE)
    for(k in 1:clusternum){
      sigma_cur = msEst$parameters$variance$sigma[,,k]
      if(is.positive.definite(sigma_cur, tol=1e-8)==FALSE){
        sigma_adjust = sigma_cur + diag(nrow(sigma_cur))*1e-4
        msEst$parameters$variance$sigma[,,k] = sigma_adjust
      }
    }
    esEst = estep(modelName="VVV", data=xx_new, parameters=msEst$parameters,warnings=TRUE)
    esEst$z =  adjust_gammaKn_function2(xx_new, gammaKn = esEst$z)
    
    a = esEst$loglik
    loglikvector = c(loglikvector, a)
    loglikdiff = abs((loglik_cur - a))
    
    if((loglikdiff < 1e-15 & length(keep)<(max_feature_num+1)) || ncol(xx_new)<=2) {
      print(paste('converged at iteration:', i))
      for(h in 1:10){
        msEst = mstep(modelName="VVV", data=xx_new, z=esEst$z)
        esEst = estep(modelName="VVV", data=xx_new, parameters=msEst$parameters)
      }
      break}
    else {
      loglik_cur = a
    }
  }
  mylabel = classassign(esEst$z)
  confusion = table(truelabel, mylabel)
  acc = myaccuracy(truelabel, mylabel)
  results=list(mylabel,RI_matrix,keep,acc,confusion,delete)
  return(results)
}

ESM=function(xx,obj, clusternum,maxiter,truelabel,max_feature_num,threshold1,threshold2){
  p=ncol(xx);N=nrow(xx);delete=NULL;keep=1:p;xx_new=xx; delta_matrix=matrix(0,maxiter,p)
#  cl=makeCluster(no_cores);
  cur.loglik=0
  # Initialization way1
  #gammaKn = e.step(xx,obj)
  gammaKn = obj$z
  df <- data.frame(matrix(unlist(gammaKn), ncol=clusternum, byrow=F),stringsAsFactors=FALSE)
  v <- Mstep(xx, t(df));
  cur.loglik <- loglike(xx,v[[1]],v[[2]],v[[3]])
  loglikvector <- cur.loglik
  delta_cur=numeric(p)
  pb <- progress_bar$new(total = 100)
  for (i in 2:maxiter){
    pb$tick()
    # Repeat E and M steps till convergence
    temper=Estep1(xx_new,v[[1]],v[[2]],v[[3]])
    gammaKn=temper[[1]];delta_new=temper[[2]]
    delta_matrix[i,keep]=temper[[2]]
    index=which.min(delta_new)
    if((abs(mean(delta_matrix[i,]-delta_matrix[i-1,]))<threshold1)&(delta_new[index]<threshold2)&(i>5))
    {delete=c(delete,keep[index])
    xx_new=xx_new[,-index]
    keep=keep[-index]}
    v <- Mstep(xx_new, gammaKn);
    a <-loglike(xx_new,v[[1]],v[[2]],v[[3]])
    loglikvector <- c(loglikvector, a)
    loglikdiff <- abs((cur.loglik - a))
    #loglik.vector
    if((loglikdiff < 1e-15 & length(keep)<(max_feature_num+1)) || ncol(xx_new)<=2) {
      print(paste('converged at iteration:', i))
      for(h in 1:10){
        v <- Mstep(xx_new, gammaKn);
        gammaKn = Estep(xx_new,v[[1]],v[[2]],v[[3]]) 
      } 
      break } 
    else {
      cur.loglik <- a
    } 
  }
  i = i + 1
  mylabel=classassign(t(gammaKn))
  confusion=table(truelabel,mylabel)
  acc=myaccuracy(truelabel,mylabel)
  results=list(mylabel,delta_matrix,keep,acc,confusion,delete)
  return(results)
}


benchmark_compare=function(X,Class,num_class,max_feature_num, include_model2=TRUE){
  #initial model with full features 
  mod0=Mclust(X,G=num_class)
  ARI_0=ARI(mod0$classification,Class)
  Error_0=classError(Class,mod0$classification)$errorRate
  print(paste("Finished model0",Error_0))
  
  # ESM Model1
  start_time <- Sys.time()
  mod1 = ESM_fast(X,mod0,clusternum=num_class,maxiter=100,truelabel=Class,
                  max_feature_num=max_feature_num,threshold1=0.001,threshold2=0.05)
  end_time <- Sys.time()
  time1 = end_time - start_time
  ARI_1=ARI(mod1[[1]], Class)
  Error_1=classError(Class,mod1[[1]])$errorRate
  print(paste("Finished model1",Error_1))
  
  # Model 2
  if(include_model2==TRUE){
    start_time <- Sys.time()
    out=clustvarsel(X,G=num_class)
    mod2=Mclust(X[,out$subset],G=num_class)
    end_time <- Sys.time()
  } else{mod2 = mod0}
  time2=end_time - start_time
  ARI_2=ARI(Class,mod2$classification)
  Error_2=classError(Class,mod2$classification)$errorRate
  print(paste("Finished model2",Error_2))
  # Model 3
  # mod3 = mod0
  start_time <- Sys.time()
  mod3=VarSelCluster(X,gvals=num_class,nbcores=4,initModel=1000,crit.varsel="BIC")
  end_time <- Sys.time()
  time3=end_time - start_time
  ARI_3=ARI(Class, fitted(mod3))
  Error_3=classError(Class,fitted(mod3))$errorRate
  print(paste("Finished model3",Error_3))
  # Model 4
  start_time <- Sys.time()
  mod4 = vscc(X,G=num_class)
  end_time <- Sys.time()
  time4=end_time - start_time
  ARI_4=ARI(Class,mod4$bestmodel$classification)
  Error_4=classError(Class,mod4$bestmodel$classification)$errorRate
  print(paste("Finished model4",Error_4))
  # Model 5
  start_time <- Sys.time()
  mod5 = SelvarClustLasso(X,nbcluster=num_class,nbcores=4)
  end_time <- Sys.time()
  time5=end_time - start_time
  ARI_5=ARI(Class, mod5$partition)
  Error_5=classError(Class,mod5$partition)$errorRate
  print(paste("Finished model5",Error_5))
  
  # Kmeans Model 6
  start_time <- Sys.time()
  model_kmeans = kmeans(X, centers=num_class, iter.max=100)
  end_time <- Sys.time()
  time6 = end_time - start_time
  ARI_6=ARI(Class, model_kmeans$cluster)
  Error_6=classError(Class,model_kmeans$cluster)$errorRate
  print(paste("Finished model6 Kmeans",Error_6))
  
  algorithm=c("Proposed","clustvarsel","VarSelCluster","vscc","selvarclustLasso","KMeans")
  algorithm=factor(algorithm,levels=c("Proposed","clustvarsel","VarSelCluster","vscc","selvarclustLasso","KMeans"))
  ARI=c(ARI_1, ARI_2,ARI_3, ARI_4,ARI_5, ARI_6)
  ErrorRate=c(Error_1, Error_2,Error_3,Error_4,Error_5,Error_6)
  Time = c(time1,time2,time3,time4,time5,time6)
  metrics=data.frame(algorithm,ARI,ErrorRate,Time)
  models = list(mod1,mod2, mod3, mod4, mod5, model_kmeans)
  benchmark_results = list()
  benchmark_results[[1]] = metrics
  benchmark_results[[2]] = models
  return (benchmark_results)
}


##Example 
##test on simulation data
## prepare simulation data
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


set.seed(999)
sim3=simulation1(400)

obj = Mclust(sim3[[1]],G=2) # a more robust way of initialization
start_time <- Sys.time()
myresults2=ESM_fast(sim3[[1]],obj,clusternum=2,maxiter=100,truelabel=sim3[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
end_time <- Sys.time()
end_time - start_time
myresults2[[3]]

result2 = benchmark_compare(sim3[[1]],sim3[[2]],num_class =2,max_feature_num = 3)

