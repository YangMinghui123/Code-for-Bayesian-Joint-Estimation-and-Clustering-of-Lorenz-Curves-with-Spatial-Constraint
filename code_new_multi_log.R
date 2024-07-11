#install.packages("gtools")
require(gtools)
# p.assumed<-rep(0.3,10)
# 
# #y = as.vector(sapply(1:length(p.assumed), function(x){sapply(1:30, function(y) {rbinom(1,1,p.assumed[x])})}))
# y = sapply(1:length(p.assumed), function(x){rbinom(1,20,p.assumed[x])})

#fit2 = CDMFM_new1(datasimple, niterations = 500, alpha = 1, beta = 1, GAMMA = 1, LAMBDA = 1, initNClusters = 5,VN = VN)
mbeta = function(Alpha)
{
  exp(sum(lgamma(Alpha))-lgamma(sum(Alpha)))
}

lmbeta = function(Alpha)
{
 sum(lgamma(Alpha))-lgamma(sum(Alpha))
}

cut = function(d,c)
{as.numeric(d<=c)}

## Dahl's method to summarize the samples from the MCMC
getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from CDMFM_new ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

###### data is an array, with the third dim denoting the year, the first dim representing the region index i, the second dim representing the component q 

## function for Collapsed sampler
CDMFM_new1 <- function(data,niterations, lambda1, neighbour,distance, alpha, GAMMA, LAMBDA, initNClusters,VN) # lambda1 is spatial parameter, neighbour=1, distance is graph distance calculated by adjmatrix.
{ 
  n = dim(data)[1]
  #precomputation for prespecified coefficient VN
  lambda <- LAMBDA
  gamma <- GAMMA
  N=n ## n is the number of oberservations
  
  #initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  ## a matrix with each row denoting the probability for a cluster
  phi<-rdirichlet(initNClusters, alpha = alpha)
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          cost = exp(lambda1*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+c.counts.noi[x])*cost*prod(sapply(1:T, function(t) {dmultinom(x = data[i,,t], size = NULL, prob = phi[x,] , log = FALSE)}))### data is an array
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*exp(sum(sapply(1:T, function(t) {lfactorial(sum(data[i,,t]))-sum(lfactorial(data[i,,t]))}))+lmbeta(alpha+rowSums(data[i,,]))-lmbeta(alpha) + VN[nClusters+1]-VN[nClusters])#*prod(sapply(1:T, function(t) {exp(lfactorial(sum(data[i,,t]))-sum(lfactorial(data[i,,t])))}))*mbeta(alpha+rowSums(data[i,,]))/mbeta(alpha)*exp(VN[nClusters+1]-VN[nClusters])
                              
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          phinew = matrix(0,nClusters+1,dim(phi)[2])
          phinew[1:nClusters,] = phi
          phinew[nClusters+1,] = rdirichlet(1, alpha = alpha)
          phi = phinew
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {phi = phi
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
        
      } else {
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          cost = exp(lambda1*sum(cut(distance[i,clusterAssign_temp==x],neighbour)))
          (GAMMA+c.counts.noi[x])*cost*prod(sapply(1:T, function(t) {dmultinom(x = data[i,,t], size = NULL, prob = phi[x,] , log = FALSE)}))
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*exp(sum(sapply(1:T, function(t) {lfactorial(sum(data[i,,t]))-sum(lfactorial(data[i,,t]))}))+lmbeta(alpha+rowSums(data[i,,]))-lmbeta(alpha) + VN[nClusters+1]-VN[nClusters])#*prod(sapply(1:T, function(t) {exp(lfactorial(sum(data[i,,t]))-sum(lfactorial(data[i,,t])))}))*mbeta(alpha+rowSums(data[i,,]))/mbeta(alpha)*exp(VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else
        {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          phi = phi[-cur.cluster.i,]}
      }
    }
    # end for loop over subjects i
    ## update phi ##
    phi = matrix(0, nClusters,length(alpha))
    for (r in 1:nClusters){
      if (sum(clusterAssign == r) == 1) {
        phi[r,] = rdirichlet(1,alpha + rowSums(data[clusterAssign == r,,]))
      } else {
        phi[r,] = rdirichlet(1,alpha + colSums(rowSums(data[clusterAssign == r,,],dims = 2))) 
      }
    }
    History[[iter]] <- list(zout = clusterAssign,phiout = phi)
    #cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}


gamma = 1; lambda = 1; N = dim(multinomial_data)[1]
VN<-0
tmax = N+10
for (t in 1:tmax)
{
  r = log(0)
  for (k in t:500)
  {
    b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
    m = max(b,r)
    r = log(exp(r-m) + exp(b-m)) + m
  }
  VN[t] = r
}

#CDMFM_new1(data,100, rep(1,10), GAMMA = 1, LAMBDA = 1, initNClusters = 2,VN)
