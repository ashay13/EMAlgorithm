#Libraries
#install.packages('readr')
#install.packages('mvtnorm')
#install.packages('ggplot2')
#install.packages('reshape2')
library(reshape2)
library(ggplot2)
library(mvtnorm)
library(readr)

#Import Data
getDataSet = function(){
  ionosphere = read_csv("https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data", 
                        col_names = FALSE, col_types = cols(X2 = col_skip(),X35 = col_character()))
  return(ionosphere)
}


multivariateNormalDistribution = function(x,mean,co,prior){
  dmv = dmvnorm(x = x,mean = mean,sigma = co)
  return(dmv*prior)
}

findWeightedProbability = function(j,i,ionosphere_matrix,clusterMean,covarianceMatrix,priors){
  
  mvnd_currentCluster = 0
  mvnd_allCluster = 0
  for(a in 1:nrow(clusterMean)){
    mvnd = multivariateNormalDistribution(ionosphere_matrix[j,],clusterMean[a,],covarianceMatrix[[a]],priors[a])
    if(a == i){
      mvnd_currentCluster = mvnd
    }
    mvnd_allCluster = mvnd + mvnd_allCluster
  }
  
  f = mvnd_currentCluster/mvnd_allCluster
  
  return(f)
}

updateClusterMean = function(i,weightedProbability,ionosphere_matrix){
  numerator = 0
  denominator = 0
  
  numerator = sum(weightedProbability[,i]*ionosphere_matrix[i,])
  denominator = sum(weightedProbability[,i])
  
  return(numerator/denominator)
}

updatePriors = function(i,weightedProbability,n){
  s = sum(weightedProbability[,i])
  p = s/n
  return(p)
}

findMeanDifference = function(clusterMean_old,clusterMean_new){
  return(sum((clusterMean_new - clusterMean_old)^2))
}


EvaluateEMAlgorithm = function(ionosphere_matrix, k,threshold,initialClusterMean){
  dimension = ncol(ionosphere_matrix)
  
  #Initialize parameters
  
  #Cluster Mean
  clusterMean_new = matrix(nrow = k,ncol = dimension)
  clusterMean_old = matrix(nrow = k,ncol = dimension)
  clusterMean_new = initialClusterMean
  
  #Covariance Matrix for each cluster
  covarianceMatrix = rep(list(data.matrix(diag(dimension))), k)
  
  #Uniformly distributed prior probabilities
  priors = rep((1/k),k)
  
  #Weighted Probability (Wij)
  weightedProbability = matrix(nrow = nrow(ionosphere_matrix),ncol = k)
  
  #threshold
  t = 9999999
  iter = 1
  
  while(t > threshold){
    
    clusterMean_old = clusterMean_new
    
    #Begin Expectation Step
    for(i in 1:k){
      for(j in 1:nrow(ionosphere_matrix)){
        weightedProbability[j,i] = findWeightedProbability(j,i,ionosphere_matrix,clusterMean_old,covarianceMatrix,priors) + 0.00000005
      }
    }
    
    #Begin Maximization Step
    
    #Re-estimate cluster means
    for(i in 1:k){
      clusterMean_new[i,] = updateClusterMean(i,weightedProbability,ionosphere_matrix)
    }
    
    #Re-estimate Covariance Matrix
    for(i in 1:k){
      temp = cov.wt(x = ionosphere_matrix,wt = weightedProbability[,i],center = clusterMean_new[i,])
      covarianceMatrix[[i]] = temp$cov
      
    }
    
    #Re-estimate priors
    for(i in 1:k){
      priors[i] = updatePriors(i,weightedProbability,nrow(ionosphere_matrix))    
    }
    
    #Compute cluster mean difference
    t = findMeanDifference(clusterMean_old,clusterMean_new)
    #cat("Mean difference = ",t,"\n")
    iter = iter + 1
  }
  
  ## EM Alogorithm is successfully. Now, let's find labels for each data
  cluster_label = rep(0,nrow(ionosphere_matrix))
  for(i in 1:nrow(weightedProbability)){
    cluster_label[i] = which.max(weightedProbability[i,])
  }
  
  returnList = list("clusterLabel" = cluster_label, "iter" = iter-1, "mean" = clusterMean_new)
  return(returnList)
}

#Calculate Euclidean Distance function
calculate_Euclidean_Dist = function(TrueMean,PredictedMean){
  sum = 0
  for(i in 1:(length(PredictedMean))){
    sum = (TrueMean[i]-PredictedMean[i])^2 + sum
  }
  
  return (sqrt(sum))
}

#Calculate Error when K=2
findErrorForK2 = function(trueLabels,clusterLabel){
  
  #cluster_error = rep(0,2)
  Total_error = 0
  for(i in 1:2){
    b=0
    g=0
    for(j in 1:nrow(trueLabels)){
      if(clusterLabel[j] == i){
        if(isTRUE(all.equal.character(as.character(trueLabels[j,]),'b'))){
          b = b + 1
        }else if(isTRUE(all.equal.character(as.character(trueLabels[j,]),'g'))){
          g = g + 1
        }
      }
    }
    if(b != 0 & g != 0){
      #Find Cluster-wise error
      if(g > b){
        error = b/(g+b)
        #cluster_error[i] = error
        Total_error = Total_error+error
      }else{
        error = g/(g+b)
        #cluster_error[i] = error
        Total_error = Total_error+error
      }  
    }
    
  }
  return(Total_error)  
}

#This Function Calculates Number of Good and Bad in clusters 
calculateNumberOfGoodAndBad = function(k,clusterLabel,trueLabels){
  g = 0
  b = 0
  for(i in 1:length(clusterLabel)){
    if(clusterLabel[i] == k){
      if(isTRUE(all.equal.character(as.character(trueLabels[i,]),'b'))){
        b = b + 1
      }else if(isTRUE(all.equal.character(as.character(trueLabels[i,]),'g'))){
        g = g + 1
      }
    }
  }
  returnList = list("good" = g, "bad" = b)
  return(returnList) 
}


#Calculate Error for K = 3,4,5
findErrorForOtherK = function(trueGood,trueBad,predictMean,clusterLabel,trueLabels){
  
  #TotalError = 0
  newClusterLabel = rep(NA,length(clusterLabel))
  
  for(k in 1:nrow(predictMean)){
    good_dist = calculate_Euclidean_Dist(trueGood,predictMean[k,])
    bad_dist = calculate_Euclidean_Dist(trueBad,predictMean[k,])
    
    if(good_dist > bad_dist){
      for(i in 1:length(clusterLabel)){
        if(clusterLabel[i] == k){
          newClusterLabel[i] = 1
        }
      }
    }
    else if(good_dist < bad_dist){
      for(i in 1:length(clusterLabel)){
        if(clusterLabel[i] == k){
          newClusterLabel[i] = 2
        }
      }
    }
  }
  TotalError = findErrorForK2(trueLabels = trueLabels, clusterLabel = newClusterLabel)
  
  return(TotalError)
}

#=======================================================================================================================================
#  K-Means Plus Plus Initialization 
#

matequal <- function(x, y){
  x = unname(x,force = TRUE)
  if(length(x) == length(y)){
    flag = TRUE
    for(i in 1:length(x)){
      if(isTRUE(all.equal(x[i],y[i]))){
        
      }else{
        return(FALSE)
      }
    }
  }else{
    return(FALSE)
  }
  return(flag)
}

#This function will return next Centroid Vector
find_next_centroid = function(A,centroid,n){
  min_point_centroid_dist = rep(0,nrow(A))
  possible_centroid = matrix(nrow = nrow(A) , ncol = (ncol(A)))
  k = 1
  for(i in 1:nrow(A)){
    repeat_flag=0
    dist_nearest_cluster = rep(0,n)
    for(j in 1:n){
      if(matequal(A[i,],centroid[j,])){
        repeat_flag = 1
        break
      }
      dist_nearest_cluster[j]=calculate_Euclidean_Dist(A[i,],centroid[j,])
    }
    if(repeat_flag == 1){
      next
    }
    min_point_centroid_dist[k]=dist_nearest_cluster[which.min(dist_nearest_cluster)]
    possible_centroid[k,] = A[i,]
    k = k + 1
  }
  
  sum_of_points = sum(min_point_centroid_dist)
  probability = (sum_of_points)^-1 * min_point_centroid_dist
  return(possible_centroid[which.max(probability),])
}

initializeKMPPClusters = function(A,k){
  centroid = matrix(nrow = k , ncol = (ncol(A)))
  centroid[1,] = A[sample(nrow(A), 1),]
  #centroid[1,] = A[1,]
  for(i in 2:k){
    centroid[i,] = find_next_centroid(A,centroid,i-1)
  }
  return(centroid)
}
#########################################################################################################

# Lets start with comparing EM Algorithm with random initialization and KM++ initialization
#
# First Download the Dataset and preprocess it.
ionosphere = getDataSet()
ionosphere_matrix = data.matrix(subset(ionosphere, select = -c(34)))
ionosphere_label_matrix = as.vector(subset(ionosphere, select = c(34)))

# Here I'm calculating true centroids based on labels
ionosphere = as.data.frame(ionosphere)

good = ionosphere[ionosphere$X35 == 'g', ]
bad = ionosphere[ionosphere$X35 == 'b', ]

good = data.matrix(subset(good, select = -c(34)))
true_good = colMeans(good)


#Reporting Variables initialization
All.iter = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Initialization", "K", "Iteration"))
All.error = setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("Initialization", "K", "Error"))



bad = data.matrix(subset(bad, select = -c(34)))
true_bad = colMeans(bad)

####################################################################################################  

# 1. First evaluate EM Algorithm with Random Initialization of Gaussians for k=2,3..5 20 Runs Each. Let's set Threshold value to 0.0005
threshold = 0.001

EM_iter = matrix(nrow = 4 , ncol = 20)
EM_error = matrix(nrow = 4 , ncol = 20)
for(k in 2:5){
  cat("EM Random k = ",k,"\n")
  for(run in 1:20){
    cat("Run = ",run,", ")
    initialClusterMean = ionosphere_matrix[sample(nrow(ionosphere_matrix),size=k,replace=TRUE),]
    returnList = EvaluateEMAlgorithm(ionosphere_matrix,k,threshold,initialClusterMean)
    
    cluster_label = returnList$clusterLabel
    iter = returnList$iter
    clusterMean = returnList$mean
    
    EM_iter[k-1,run] = iter
    
    if(k == 2){
      EM_error[k-1,run] = findErrorForK2(ionosphere_label_matrix,cluster_label)
    }
    else{
      EM_error[k-1,run] = findErrorForOtherK(true_good,true_bad,clusterMean,cluster_label,ionosphere_label_matrix)
    }
  }
}

cat("******************* Number of Iterations in EM - ",EM_iter)
cat("******************* Error in EM - ",EM_error)


####################################################################################################  

# 2. Now evaluate EM Algorithm with K-means plus plus initialization for k=2,3..5 20 Runs Each.


EMPP_iter = matrix(nrow = 4 , ncol = 20)
EMPP_error = matrix(nrow = 4 , ncol = 20)
for(k in 2:5){
  cat("EM KM++ k = ",k,"\n")
  for(run in 1:20){
    cat("Run = ",run,", ")
    initialClusterMean = initializeKMPPClusters(ionosphere_matrix,k)
    returnList = EvaluateEMAlgorithm(ionosphere_matrix,k,threshold,initialClusterMean)
    
    cluster_label = returnList$clusterLabel
    iter = returnList$iter
    clusterMean = returnList$mean
    
    EMPP_iter[k-1,run] = iter
    
    if(k == 2){
      EMPP_error[k-1,run] = findErrorForK2(ionosphere_label_matrix,cluster_label)
    }
    else{
      EMPP_error[k-1,run] = findErrorForOtherK(true_good,true_bad,clusterMean,cluster_label,ionosphere_label_matrix)
    }
  }
}

cat("******************* Number of Iterations in EM with KM++  - ",EMPP_iter)
cat("******************* Error in EM with KM++ - ",EMPP_error)



#______________________________________________________________________________________________________________  

##Whisker plot for Iteration

trans1 = t(EM_iter)
colnames(trans1) <- c("K2","K3","K4","K5")
temp1 = melt(trans1)
colnames(temp1) <- c("Run","Cluster","Iteration")
temp1 = data.frame(temp1)
temp1$Algorithm <- "EM"

trans2 = t(EMPP_iter)
colnames(trans2) <- c("K2","K3","K4","K5")
temp2 <- melt(trans2)
colnames(temp2) <- c("Run","Cluster","Iteration")
temp2 <- data.frame(temp2)
temp2$Algorithm <- "EM with K-Means++"
temp3 = rbind(temp1,temp2)

ggplot(data = temp3, aes(x=Cluster, y=Iteration)) + geom_boxplot(aes(fill=Algorithm))


#-----------------------------------------------------------------------------------------

##Whisker plot for Error

trans3 = t(EM_error)
colnames(trans3) <- c("K2","K3","K4","K5")
temp3 = melt(trans3)
colnames(temp3) <- c("Run","Cluster","Error")
temp3 = data.frame(temp3)
temp3$Algorithm <- "EM"

trans4 = t(EMPP_error)
colnames(trans4) <- c("K2","K3","K4","K5")
temp4 <- melt(trans4)
colnames(temp4) <- c("Run","Cluster","Error")
temp4 <- data.frame(temp4)
temp4$Algorithm <- "EM with K-Means++"
temp5 = rbind(temp3,temp4)

ggplot(data = temp5, aes(x=Cluster, y=Error)) + geom_boxplot(aes(fill=Algorithm))

