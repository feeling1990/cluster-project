##EM feature selection algorithm 

## pre-install library packages
#install.packages("EMCluster","mvnfast","parallel","ggplot2","progress","matrixcalc","robustHD")
#install.packages("ggplot2","mclust","VarSelLCM","vscc","SelvarMix","clustvarsel")
library(matrixcalc)
library(EMCluster)
library(mclust)
library(mvnfast)
library(ggplot2)
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

CalculateDeltaVector=function(x, y){
  h = x
  for (i in 1:dim(x)[1]){
    for(j in 1:dim(y)[2]){
      h[i, j] = abs(x[i,j]-y[i,j])
    }
  }
  return(h)
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

GetMaxIndex = function(data_matrix){
  m = nrow(data_matrix)
  n = ncol(data_matrix)
  max_value = max(data_matrix)
  index_value = which(data_matrix==max_value, arr.ind = TRUE)
  if (dim(index_value)[1]>1){
    index_value = index_value[1,]}
#  index_j = colnames(data_matrix[index_value[2]])
#  if (is.null(index_j)){return (index_value[[2]])}
  index_j = index_value[[2]]
  return(index_j)
}



shuffle = function(x){
  y = NULL
  for (i in 1:length(x)){
    if (x[i]==1){y = c(y,2)}
    else{y=c(y,1)}
  }
  return (y)
}
data_matrix = matrix(9:1, nrow = 3, ncol = 3)
index = GetMaxIndex(data_matrix)
print(index)


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
  
  gammaKn = t(esEst$z)
  # for(j in 1:p){
  #   RI[j]=Mdiff1(gammaKn,GammaKn[,,j])
  # }
  # exp on variance of deltas
  RI_mean=rep(0,p)
  RI_sd=rep(0,p)
  
  Outliers = NULL
  DeltaVector = list()
  for(j in 1:p){
    DeltaVector[[j]] = CalculateDeltaVector(gammaKn, GammaKn[,,j])
    delta_j_mean = mean(DeltaVector[[j]])
    delta_j_sd = sd(DeltaVector[[j]])
    delta_j_max = max(DeltaVector[[j]])
    delta_j_min = min(DeltaVector[[j]])
    delta_j_median = median(DeltaVector[[j]])
    RI_mean[j] = delta_j_mean
    RI_sd[j] = delta_j_sd
    # print("============================")
    # print(paste("feature:", j))
    # print(paste("delta_j_mean", delta_j_mean))
    # print(paste("delta_j_standard_deviation", delta_j_sd))
    # print(paste("delta_j_max", delta_j_max))
    # print(paste("delta_j_min", delta_j_min))
    # print(paste("delta_j_median", delta_j_median))
    #print(DeltaVector[[j]])
    #hist(as.vector(DeltaVector[[j]]))
  }
 RI_upper = RI_mean + RI_sd
  for(j in 1:p){
    if(RI_upper[j] < max(RI_upper)){
      max_index = GetMaxIndex(DeltaVector[[j]])
      Outliers = c(Outliers, max_index)
    }
  }
  return (list(RI_mean,RI_sd, Outliers,DeltaVector))
  }
  
SummarizeDelta = function(DeltaVector,keep){
  p = length(DeltaVector)
  summary_matrix = matrix(0,p,5)
  for (j in 1:p){
    delta_j_mean = mean(DeltaVector[[j]])
    delta_j_sd = sd(DeltaVector[[j]])
    delta_j_max = max(DeltaVector[[j]])
    delta_j_min = min(DeltaVector[[j]])
    summary_matrix[j,] = c(keep[j], delta_j_mean, delta_j_sd, delta_j_min, delta_j_max)
  }
  colnames(summary_matrix) = c('feature', 'mean', 'sd', 'min', 'max')
  return (summary_matrix)
}

match_label = function(list1, list2){
  names = sort(unique(list2))
  confusion = table(list1,factor(list2, levels=names))
  column_order = apply(confusion, 1, which.max)
  list1_reorder = rep('a', length(list1)) 
  for (j in 1:length(list1)){
    list1_reorder[j] = names[column_order[list1[j]]]
  }
  #confusion_reorder = table(list1_reorder, list2)
  return(list1_reorder)
}


### ESM_fast_confidence
ESM_fast_confidence=function(xx,obj,clusternum,maxiter,truelabel,max_feature_num,threshold1,threshold2){
  p=ncol(xx);
  N=nrow(xx);
  delete=NULL;
  keep=1:p;
  RI_matrix=matrix(0,maxiter,p)
  Outlier_Candidates = NULL
  loglik_cur=0
  gammaKn = obj$z
  parameters = obj$parameters
  loglik_cur = obj$loglik
  loglikvector <- loglik_cur
  xx_new=xx; 
  msEst = obj
  esEst = obj
  delta_summary_matrix = list()
  ###begin loop for E and M step
  pb <- progress_bar$new(total = 100)
  for (i in 2:maxiter){
    pb$tick()
    RI_info = CalculateRI(xx_new, msEst$parameters, esEst)
    delta_summary_matrix[[i]] = SummarizeDelta(RI_info[[4]],keep)
    print ("=========RI_info==========")
    print(RI_info)
    RI_mean = RI_info[[1]]
    RI_sd = RI_info[[2]]
    RI_matrix[i, keep] = RI_mean
    RI_upper = RI_mean + RI_sd
    index = which.min(RI_upper)
    if(abs(mean(RI_matrix[i,]-RI_matrix[i-1,])) < threshold1){
      Outlier_Candidates = c(Outlier_Candidates,RI_info[[3]])
    }
    if((abs(mean(RI_matrix[i,]-RI_matrix[i-1,]))<threshold1)&(RI_upper[index]<threshold2)&(i>5)){
      delete = c(delete, keep[index])
      xx_new = xx_new[,-index]
      keep = keep[-index]
    }
    
    
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
  mylabel_raw = classassign(esEst$z)
  print(mylabel_raw)
  print(truelabel)
  mylabel = match_label(mylabel_raw, truelabel)
  print(mylabel)
  #mylabel = classassign(esEst$z)
  confusion = table(mylabel,truelabel)
  acc = myaccuracy(mylabel, truelabel)
  misclassified = which(truelabel != mylabel) 
  ari = ARI(mylabel, truelabel)
  results = list()
  results[['mylabel']]= mylabel
  results[['RI_matrix']] = RI_matrix
  results[['delta_summary_matrix']] = delta_summary_matrix
  results[['keep']] = keep
  results[['acc']] = acc
  results[['ari']] = ari
  results[['confusion']]= confusion
  results[['misclassified']] = misclassified
  results[['delete']] = delete
  results[['outliers']] = Outlier_Candidates
  #results=list(mylabel,RI_matrix,keep,acc,confusion,misclassified, delete,Outliers)
  return(results)
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


# set.seed(999)
# sim3=simulation1(300)
# 
# obj = Mclust(sim3[[1]],G=2) # a more robust way of initialization
# start_time <- Sys.time()
# myresults2=ESM_fast_confidence(sim3[[1]],obj,clusternum=2,maxiter=100,truelabel=sim3[[2]],max_feature_num=3,threshold1=0.001,threshold2=0.05)
# end_time <- Sys.time()
# end_time - start_time
# myresults2

#result2 = benchmark_compare(sim3[[1]],sim3[[2]],num_class =2,max_feature_num = 3)

