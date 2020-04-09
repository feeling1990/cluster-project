
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

adjust_gammaKn_function2 = function(xx_new,gammaKn){
  for (k in 1:num_class){
    if (sum(gammaKn[,k])==0){
      print("restart")
      obj=init.EM(xx_new,nclass=num_class)
      gammaKn_new = e.step(xx_new,obj)$Gamma
      return(gammaKn_new)
    }
  }
  return(gammaKn)
}

### ESM_fast
xx = X
clusternum = 9
mod0 = Mclust(xx,G=clusternum)
obj = mod0

truelabel = Class
print(dim(X))
print(length(truelabel))
maxiter = 100
max_feature_num = 10
threshold1 = 0.001
threshold2 = 0.05

#ESM_fast=function(xx,obj,clusternum,maxiter,truelabel,max_feature_num,threshold1,threshold2){
  p=ncol(xx);
  N=nrow(xx);
  delete=NULL;
  keep=1:p;
  RI_matrix=matrix(0,maxiter,p)
  loglik_cur=0
  loglik_cur = obj$loglik
  loglikvector <- loglik_cur
  xx_new=xx; 
  msEst = obj
  esEst = obj
  ###begin loop for E and M step
  #pb <- progress_bar$new(total = 100)
  #for (i in 2:maxiter){
    i = 2
   # pb$tick()
    #gammaKn = t(esEst$z)
    RI = CalculateRI(xx_new, msEst$parameters, esEst)
    RI_matrix[i, keep] = RI
    index = which.min(RI)
    #print (abs(mean(RI_matrix[i,]-RI_matrix[i-1,])))
    #print (RI[index])
    if((abs(mean(RI_matrix[i,]-RI_matrix[i-1,]))<threshold1)&(RI[index]<threshold2)&(i>1)){
      delete = c(delete, keep[index])
      xx_new = xx_new[,-index]
      keep = keep[-index]
    }
   
    msEst = mstep(modelName="VVV", data=xx_new, z=esEst$z,warnings=TRUE)
    for(k in 1:num_class){
      sigma_cur = msEst$parameters$variance$sigma[,,k]
      if(is.positive.definite(sigma_cur, tol=1e-8)==FALSE){
        sigma_adjust = sigma_cur + diag(nrow(sigma_cur))*1e-2
        msEst$parameters$variance$sigma[,,k] = sigma_adjust
      }
    }
    is.positive.definite(msEst$parameters$variance$sigma[,,k])
    esEst = estep(modelName="VVV", data=xx_new, parameters=msEst$parameters,warnings=TRUE)
    esEst$z =  adjust_gammaKn_function2(xx_new, gammaKn = esEst$z)
    #esEst$z
    a = esEst$loglik
    loglikvector = c(loglikvector, a)
    loglikdiff = abs((loglik_cur - a))
    i = i+1
    delete

    # if((loglikdiff < 1e-15 & length(keep)<(max_feature_num+1)) || ncol(xx_new)<=2) {
    #   print(paste('converged at iteration:', i))
    #   for(h in 1:10){
    #     msEst = mstep(modelName="VVV", data=xx_new, z=esEst$z)
    #     esEst = estep(modelName="VVV", data=xx_new, parameters=msEst$parameters)
    #   }}
    # else {loglik_cur = a}
   
 # }
  mylabel = classassign(esEst$z)
  confusion = table(truelabel, mylabel)
  acc = myaccuracy(truelabel, mylabel)
  results=list(mylabel,RI_matrix,keep,acc,confusion,delete)
#   return(results)
# }

parameters = msEst$parameters
mu=parameters$mean
pro = parameters$pro
G = 10
K = G
data = X
n = nrow(data)
p = ncol(data)
Vinv = NULL
temp <- .Fortran("esvvv", as.logical(1), as.double(data), 
                 as.double(mu), as.double(parameters$variance$cholsigma), 
                 as.double(pro), as.integer(n), as.integer(p), as.integer(G), 
                 as.double(if (is.null(Vinv)) -1 else Vinv), double(p), 
                 double(1), double(n * K), PACKAGE = "mclust")[10:12]


temp <- .Fortran("msvvvp", as.double(data), as.double(z), 
                 as.integer(n), as.integer(p), as.integer(G), as.double(priorParams$shrinkage), 
                 as.double(priorParams$mean), as.double(if (any(priorParams$scale != 
                                                                0)) chol(priorParams$scale) else priorParams$scale), 
                 as.double(priorParams$dof), double(p), double(p * 
                                                                 G), double(p * p * G), double(G), double(p * 
                                                                                                            p), PACKAGE = "mclust")[11:13]




