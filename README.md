# cluster-project
Experiments for clustering
This is my Phase I reasearch study. The study focuses on doing feature selection for Gaussian Mixture Model (GMM).

## Pre-required Packages in R
install.packages("EMCluster","mvtnorm","mvnfast","parallel")

## How to run the code
The core and basic functions are in file ESM for GMM.R. You need to run all these functions before doing other experiments.
The core function is ESM=function(xx,obj, clusternum,maxiter,truelabel,featurenum,threshold1,threshold2). You can call this function for any new dataset. 

## Example to run
sim3=simulation1(400) # the function is in [simulation for paper.R]
data = sim3[[1]]
truelabel = sim3[[2]]

obj=obj=init.EM(data,nclass=2)
myresults2=ESM(data,obj,clusternum=2,maxiter=100,truelabel=truelabel,featurenum=3,threshold1=0.001,threshold2=0.05)

myresults2 returns a list of results:
--delta_matrix, the difference between keeping and removing each feature at each iteration.
--keep, final selected features
--acc, accuracy of the final clustering results
--confusion, confusion matrix of clustering label and true label
--delete, removed set of features
--i, number of iterations when finished.

## Purpose of each file
1. ESM for GMM.R provides core functions for the ESM algorithm.
2. simulation for paper provides the code for the simulation data experiment in paper.
3. experiments_add.R provides additional expeirments on benchmarkdata like Hearts/Crab/Wine.
4. ESM paper illustration example.R provides code for the illustration example in paper.
