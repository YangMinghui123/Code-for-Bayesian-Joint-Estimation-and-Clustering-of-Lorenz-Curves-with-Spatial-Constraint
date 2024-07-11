#---Design 1 聚类方法比较---#

library(DirichletReg)
library(usmap)
library(ggplot2)
library(DescTools)
library(magrittr)
library(dplyr)
library(ggpubr)
library(future)
library(parallel)
library(purrr)
library(furrr)
library(fdasrvf)
library(ineq)
library(pracma)
library(factoextra)
library(cluster)
library(flexclust)
library(mclust)
library(dbscan)
library(apcluster)
library(MASS)
library(fossil)

n=51#洲的个数
Q=10# 收入等级的个数
T=5#5年
prob=matrix(0,Q,3)#将美国各洲根据收入平等程度划分为3类

disjoint_simulation_data1=array(0,c(n,Q,T))
disjoint_simulation_data2<-array(0,c(n,Q,T))
loglikelihood1=array(0,c(n,T,300))
loglikelihood2=array(0,c(n,T,300))
loglikelihood_sum1=c()
loglikelihood_sum2=c()
loglikelihood_dahlans_1=matrix(0,n,T)
loglikelihood_dahlans_2=matrix(0,n,T)
deviance1=c()
deviance2=c()
deviance_bar1=c()
deviance_bar2=c()
deviance_bar_Dahl1=c()
deviance_bar_Dahl2=c()
lambda1_disjoint_simulation_data1=c()
lambda1_disjoint_simulation_data2=c()

x=c()
k=c()
k_best=c()
cluster_CDMFM1_design1=c()
ARI_CDMFM1_design1=c()
cluster_CDMFM2_design1=c()
ARI_CDMFM2_design1=c()
cluster_dbscan1_design1=c()
ARI_dbscan1_design1=c()
cluster_dbscan2_design1=c()
ARI_dbscan2_design1=c()
cluster_Mclust1_design1=c()
ARI_Mclust1_design1=c()
cluster_Mclust2_design1=c()
ARI_Mclust2_design1=c()
cluster_kmeans1_design1=c()
ARI_kmeans1_design1=c()
cluster_kmeans2_design1=c()
ARI_kmeans2_design1=c()
cluster_optics1_design1=c()
ARI_optics1_design1=c()
cluster_optics2_design1=c()
ARI_optics2_design1=c()
cluster_ap1_design1=c()
ARI_ap1_design1=c()
labels_pred_ap1=c()
cluster_ap2_design1=c()
ARI_ap2_design1=c()
labels_pred_ap2=c()

ARI_CDMFM1_design1_RI=c()
ARI_CDMFM2_design1_RI=c()
ARI_CDMFM1_design1_ARI=c()
ARI_CDMFM2_design1_ARI=c()

ARI_CDMFM1_non_design1_RI=c()
ARI_CDMFM2_non_design1_RI=c()
ARI_CDMFM1_non_design1_ARI=c()
ARI_CDMFM2_non_design1_ARI=c()

ARI_dbscan1_design1_RI=c()
ARI_dbscan2_design1_RI=c()
ARI_dbscan1_design1_ARI=c()
ARI_dbscan2_design1_ARI=c()

ARI_Mclust1_design1_RI=c()
ARI_Mclust2_design1_RI=c()
ARI_Mclust1_design1_ARI=c()
ARI_Mclust2_design1_ARI=c()

ARI_kmeans1_design1_RI=c()
ARI_kmeans2_design1_RI=c()
ARI_kmeans1_design1_ARI=c()
ARI_kmeans2_design1_ARI=c()

ARI_optics1_design1_RI=c()
ARI_optics2_design1_RI=c()
ARI_optics1_design1_ARI=c()
ARI_optics2_design1_ARI=c()

ARI_ap1_design1_RI=c()
ARI_ap2_design1_RI=c()
ARI_ap1_design1_ARI=c()
ARI_ap2_design1_ARI=c()

#---without spatial---#
cluster_CDMFM1_non_design1=c()
cluster_CDMFM2_non_design1=c()
ARI_CDMFM1_non_design1=c()
ARI_CDMFM2_non_design1=c()


#--- design 1 ---#

stateinfo_disjoint_init <- readRDS("D:/Edge Download/state_fips_name.rds")
c1 <- c("AK", "HI", "CA", "OR", "WA", "ID", "NV", "AZ", "UT", "MT",
        "WY", "ME", "NH", "VT", "NY", "CT", "MA", "RI", "PA", "DE",
        "DC", "MD", "VA", "NC", "WV", "NJ", "SC", "OH")
c2 <- c("ND", "SD", "MN", "IA", "WI", "MI", "IL", "IN", "MO", "KS",
        "NE", "CO")
c3 <-  setdiff(stateinfo_disjoint_init$abbr, c(c1, c2))
stateinfo_disjoint_init$cluster <- as.factor(1 * (stateinfo_disjoint_init$abbr %in% c1) + 2 * (stateinfo_disjoint_init$abbr %in% c2)+3 * (stateinfo_disjoint_init$abbr %in% c3))
z=as.numeric(stateinfo_disjoint_init$cluster)#3类Cluster
labels_true=z


seed=20240219:20240318

for (time in 1:100) {#生成50次模拟数据
  set.seed(seed[time])
  
  disjoint_simulation_data1=array(0,c(n,Q,T))
  disjoint_simulation_data2=array(0,c(n,Q,T))
  
  #disjoint_simulation_data1
  probability=c()
  prob[,1]=c(0.23,0.14,0.13,0.12,0.10,0.09,0.07,0.05,0.04,0.03)
  prob[,2]=c(0.24,0.14,0.13,0.12,0.10,0.09,0.07,0.05,0.03,0.03)
  prob[,3]=c(0.25,0.15,0.13,0.12,0.10,0.09,0.06,0.05,0.03,0.02)
  mean=rep(0,Q-1)
  sigma=diag(8e-8,Q-1)#rand index和正确的类别数没有MFM高7e-8
                      #模拟10次，9e-8暂时没问题
  for (i in 1:n){
    probability[1:Q-1]=prob[1:Q-1,z[i]]+mvrnorm(1,mean,sigma)
    probability[Q]=1-sum(probability[1:Q-1])#probability需要加和为1
    for (j in 1:T){
      disjoint_simulation_data1[i,,j]<-rmultinom(1,sample(5000:15000,1),probability)#随机设定人数
    }
  }
  #disjoint_simulation_data2
  stateinfo <- readRDS("D:/Edge Download/state_fips_name.rds")#导入洲名数据
  p1=0.2#rand index比MFM高，但正确率太低6e2
        #2e2没问题
  for (s in 1:T){
  income_data_simu1<-matrix(0,10000,51)
  income_data_simu1[, stateinfo$abbr %in% c1] <-
    rgamma(length(c1) * 10000, shape = 1.15, scale = 50000) + 
    rbinom(length(c1) * 10000, 1, p1) *
    rgamma(length(c1) * 10000, shape = 0.3, scale = 50000)
  income_data_simu1[, stateinfo$abbr %in% c2] <-
    rgamma(length(c2) * 10000, shape = 1.2, scale = 50000) + 
    rbinom(length(c2) * 10000, 1, p1) *
    rgamma(length(c2) * 10000, shape = 0.3, scale = 50000)
  income_data_simu1[, stateinfo$abbr %in% c3] <-
    rgamma(length(c3) * 10000, shape = 1.25, scale = 50000) + 
    rbinom(length(c3) * 10000, 1, p1) *
    rgamma(length(c3) * 10000, shape = 0.3, scale = 50000)
  state_fips_name <- readRDS("D:/Edge Download/state_fips_name.rds")#载入洲名称数据
  state_fips<-as.numeric(state_fips_name$fips)#将洲数据转化成数值型
  
    for(i in 1:n){
      income_data_simu1_i<-as.data.frame(income_data_simu1[,i])
      income_data_simu1_i$income<-income_data_simu1[,i]+rnorm(10000,0,2e2)
      income_data_simu1_i<-subset(income_data_simu1_i,select = income)%>%#选第i个洲的收入数据
        arrange(income)%>%#给收入数据从小到大排序
        mutate(accumulate_population = rank(income)/n())%>%#计算人口累计百分比
        mutate(accumulate_Income = cumsum(as.numeric(income))/sum(income)) #计算收入累计百分比
      for (j in 1:Q) {#求各洲的Q个收入等级的人口数
        disjoint_simulation_data2[i,j,s]<-sum(income_data_simu1_i$accumulate_Income<j/Q&income_data_simu1_i$accumulate_Income>=(j-1)/Q)
      }
    }
  }
  
  print(disjoint_simulation_data1)
  print(disjoint_simulation_data2)
  
  #---DIC---#
  
  loglikelihood1=array(0,c(n,T,300))
  loglikelihood2=array(0,c(n,T,300))
  loglikelihood_sum1=c()
  loglikelihood_sum2=c()
  loglikelihood_dahlans_1=matrix(0,n,T)
  loglikelihood_dahlans_2=matrix(0,n,T)
  deviance1=c()
  deviance2=c()
  deviance_bar1=c()
  deviance_bar2=c()
  deviance_bar_Dahl1=c()
  deviance_bar_Dahl2=c()
  
  for (lambda_test in seq(0,3,0.1)) {#31个模型参数
    result_disjoint1=
      CDMFM_new1(data= disjoint_simulation_data1 , niterations=500, lambda1=lambda_test , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
    result_disjoint2=
      CDMFM_new1(data= disjoint_simulation_data2 , niterations=500, lambda1=lambda_test , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=10 , VN=VN)
    #写loglikelihood function，将每次迭代的参数估计值代入
    for (q in 201:500) { #从第201次迭代开始计算,舍去前200次
      for (j in 1:T) {
        for (i in 1:n) {
          loglikelihood1[i,j,q-200]=dmultinom(x=disjoint_simulation_data1[i,,j],
                                              size = NULL, 
                                              prob=result_disjoint1[["Iterates"]][[q]][["phiout"]][result_disjoint1[["Iterates"]][[q]][["zout"]][i],],
                                              log = TRUE)#计算每次迭代时，洲数据对应概率估计的对数似然值
          loglikelihood2[i,j,q-200]=dmultinom(x=disjoint_simulation_data2[i,,j],
                                              size = NULL, 
                                              prob=result_disjoint2[["Iterates"]][[q]][["phiout"]][result_disjoint2[["Iterates"]][[q]][["zout"]][i],],
                                              log = TRUE)
        }
      }
    }
    for (q in 201:500) {
      loglikelihood_sum1[q-200]=sum(loglikelihood1[,,q-200])
      loglikelihood_sum2[q-200]=sum(loglikelihood2[,,q-200])
    }
    #Deviance=-2*loglikelihood
    deviance1=-2*loglikelihood_sum1#迭代300次的deviance被记录下来
    deviance2=-2*loglikelihood_sum2
    #D(theta)_bar
    deviance_bar1[10*lambda_test+1]=sum(deviance1)/300#求迭代的deviance均值
    deviance_bar2[10*lambda_test+1]=sum(deviance2)/300
    #D(theta_bar)#求theta_bar的deviance
    DahlAns_disjoint1<-getDahl(result_disjoint1,burn=200)
    DahlAns_disjoint2<-getDahl(result_disjoint2,burn=200)
    for (j in 1:T) {
      for (i in 1:n) {
        loglikelihood_dahlans_1[i,j]=dmultinom(x=disjoint_simulation_data1[i,,j],
                                               size = NULL, 
                                               prob= DahlAns_disjoint1[["phiout"]][DahlAns_disjoint1[["zout"]][i],],
                                               log = TRUE)#计算在最优迭代下，每个洲数据对应概率估计的对数似然值
        loglikelihood_dahlans_2[i,j]=dmultinom(x=disjoint_simulation_data2[i,,j],
                                               size = NULL, 
                                               prob= DahlAns_disjoint2[["phiout"]][DahlAns_disjoint2[["zout"]][i],],
                                               log = TRUE)
      }
    }
    for (j in 1:T) {
      for (i in 1:n) {
        loglikelihood_sum_dahlans_1=sum(loglikelihood_dahlans_1)
        loglikelihood_sum_dahlans_2=sum(loglikelihood_dahlans_2)
      }
    }
    deviance_bar_Dahl1[10*lambda_test+1]=-2*loglikelihood_sum_dahlans_1
    deviance_bar_Dahl2[10*lambda_test+1]=-2*loglikelihood_sum_dahlans_2
  }
  
  #pD——参数复杂度
  pD1=deviance_bar1-deviance_bar_Dahl1
  pD2=deviance_bar2-deviance_bar_Dahl2
  #DIC
  DIC1=deviance_bar1+pD1
  DIC2=deviance_bar2+pD2
  #根据最小的DIC选出最优的lambda1
  lambda1_disjoint_simulation_data1[time]=seq(0,3,0.1)[which.min(DIC1)]
  lambda1_disjoint_simulation_data2[time]=seq(0,3,0.1)[which.min(DIC2)]
  
  #---调用本文算法(CDMFM_new1)---#
  CDMFM_result_disjoint1 = CDMFM_new1(data=disjoint_simulation_data1, niterations=500, lambda1=lambda1_disjoint_simulation_data1[time] , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
  DahlAns_disjoint1<-getDahl(CDMFM_result_disjoint1,burn=200)
  CDMFM_result_disjoint2 = CDMFM_new1(data=disjoint_simulation_data2, niterations=500, lambda1=lambda1_disjoint_simulation_data2[time] , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
  DahlAns_disjoint2<-getDahl(CDMFM_result_disjoint2,burn=200)
  cluster_CDMFM1_design1[time]=length(unique(DahlAns_disjoint1[["zout"]])) #记录聚类类别数
  cluster_CDMFM2_design1[time]=length(unique(DahlAns_disjoint2[["zout"]]))
  labels_pred_CDMFM1 = DahlAns_disjoint1[["zout"]] #记录聚类结果
  labels_pred_CDMFM2 = DahlAns_disjoint2[["zout"]]
  ARI_CDMFM1_design1[time]=adjustedRandIndex(labels_true,labels_pred_CDMFM1) #计算调整兰德指数
  ARI_CDMFM2_design1[time]=adjustedRandIndex(labels_true,labels_pred_CDMFM2)
  
  ARI_CDMFM1_design1_RI[time]=rand.index(labels_true,labels_pred_CDMFM1)
  ARI_CDMFM2_design1_RI[time]=rand.index(labels_true,labels_pred_CDMFM2)
  ARI_CDMFM1_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_CDMFM1)
  ARI_CDMFM2_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_CDMFM2)
  
  #---调用本文算法(CDMFM_new1)(without spatial)---#
  CDMFM_result_disjoint1_non = CDMFM_new1(data=disjoint_simulation_data1, niterations=500, lambda1=0 , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
  DahlAns_disjoint1_non<-getDahl(CDMFM_result_disjoint1_non,burn=200)
  CDMFM_result_disjoint2_non = CDMFM_new1(data=disjoint_simulation_data2, niterations=500, lambda1=0 , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
  DahlAns_disjoint2_non<-getDahl(CDMFM_result_disjoint2_non,burn=200)
  cluster_CDMFM1_non_design1[time]=length(unique(DahlAns_disjoint1_non[["zout"]])) #记录聚类类别数
  cluster_CDMFM2_non_design1[time]=length(unique(DahlAns_disjoint2_non[["zout"]]))
  labels_pred_CDMFM1_non = DahlAns_disjoint1_non[["zout"]] #记录聚类结果
  labels_pred_CDMFM2_non = DahlAns_disjoint2_non[["zout"]]
  ARI_CDMFM1_non_design1[time]=adjustedRandIndex(labels_true,labels_pred_CDMFM1_non) #计算调整兰德指数
  ARI_CDMFM2_non_design1[time]=adjustedRandIndex(labels_true,labels_pred_CDMFM2_non)
  
  ARI_CDMFM1_non_design1_RI[time]=rand.index(labels_true,labels_pred_CDMFM1_non)
  ARI_CDMFM2_non_design1_RI[time]=rand.index(labels_true,labels_pred_CDMFM2_non)
  ARI_CDMFM1_non_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_CDMFM1_non)
  ARI_CDMFM2_non_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_CDMFM2_non)
 
  #处理模拟数据
  disjoint_simulation_data1_processing=as.data.frame(disjoint_simulation_data1)
  for (i in 1:n) {
    disjoint_simulation_data1_processing[i,]=disjoint_simulation_data1_processing[i,]/sum(disjoint_simulation_data1_processing[i,])
  }
  disjoint_simulation_data2_processing=as.data.frame(disjoint_simulation_data2)
  for (i in 1:n) {
    disjoint_simulation_data2_processing[i,]=disjoint_simulation_data2_processing[i,]/sum(disjoint_simulation_data2_processing[i,])
  }

  
  
  #disjoint_simulation_data1_processing<-as.data.frame.array(disjoint_simulation_data1)#转化模拟数据
  #disjoint_simulation_data1_processing<-scale(disjoint_simulation_data1_processing)#标准化模拟数据
  #disjoint_simulation_data1_processing<-disjoint_simulation_data1_processing[,is.finite(colSums(disjoint_simulation_data1_processing))]#去掉模拟数据中的NAN非数值数据
  
  #disjoint_simulation_data2_processing<-as.data.frame.array(disjoint_simulation_data2)#转化模拟数据
  #disjoint_simulation_data2_processing<-scale(disjoint_simulation_data2_processing)#标准化模拟数据
  #disjoint_simulation_data2_processing<-disjoint_simulation_data2_processing[,is.finite(colSums(disjoint_simulation_data2_processing))]#去掉模拟数据中的NAN非数值数据
  
  #---调用DBSCAN---#
  dbscan_result1 <- dbscan(disjoint_simulation_data1_processing, eps = 4 ,MinPts = 5)
  cluster_dbscan1_design1[time] = length(unique(dbscan_result1$cluster))#记录聚类类别数
  labels_pred_dbscan1 <- dbscan_result1$cluster#记录聚类结果
  ARI_dbscan1_design1[time]=adjustedRandIndex(labels_true,labels_pred_dbscan1) #计算调整兰德指数
  
  dbscan_result2 <- dbscan(disjoint_simulation_data2_processing, eps = 4 ,MinPts = 5)
  cluster_dbscan2_design1[time] = length(unique(dbscan_result2$cluster))#记录聚类类别数
  labels_pred_dbscan2 <- dbscan_result2$cluster#记录聚类结果
  ARI_dbscan2_design1[time]=adjustedRandIndex(labels_true,labels_pred_dbscan2) #计算调整兰德指数
 
  ARI_dbscan1_design1_RI[time]=rand.index(labels_true,labels_pred_dbscan1) 
  ARI_dbscan2_design1_RI[time]=rand.index(labels_true,labels_pred_dbscan2)
  ARI_dbscan1_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_dbscan1)
  ARI_dbscan2_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_dbscan2)

  
  #---调用基于模型的聚类算法---#
  Mclust_result1 <- Mclust(disjoint_simulation_data1_processing)
  cluster_Mclust1_design1[time]=length(unique(Mclust_result1[["classification"]]))#记录类别数
  labels_pred_Mclust1 = Mclust_result1[["classification"]]#记录聚类结果
  ARI_Mclust1_design1[time]=adjustedRandIndex(labels_true,labels_pred_Mclust1)#计算调整兰德指数
  
  Mclust_result2 <- Mclust(disjoint_simulation_data2_processing)
  cluster_Mclust2_design1[time]=length(unique(Mclust_result2[["classification"]]))#记录类别数
  labels_pred_Mclust2 = Mclust_result2[["classification"]]#记录聚类结果
  ARI_Mclust2_design1[time]=adjustedRandIndex(labels_true,labels_pred_Mclust2)#计算调整兰德指数
  
  ARI_Mclust1_design1_RI[time]=rand.index(labels_true,labels_pred_Mclust1) 
  ARI_Mclust2_design1_RI[time]=rand.index(labels_true,labels_pred_Mclust2) 
  ARI_Mclust1_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_Mclust1)
  ARI_Mclust2_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_Mclust2)
  
  #---调用kmeans---#
  fviz_nbclust1<-fviz_nbclust(disjoint_simulation_data1_processing, kmeans, method = "wss")#画手肘图
  for (a in 1:10) {
    x[a]=fviz_nbclust1[["data"]][a,2]/fviz_nbclust1[["data"]][a+1,2]
  }
  for (b in 1:10) {
    k[1]=0
    k[b+1]=abs((x[b]-x[b+1])/x[b+1])
  }
  k_best1=which.max(k)#设置初始类别数
  kmeans_result1 <- kmeans(disjoint_simulation_data1_processing, centers = k_best1, nstart = 1)#调用k-means聚类算法
  cluster_kmeans1_design1[time]=length(unique(kmeans_result1[["cluster"]]))#记录聚类类别数
  labels_pred_kmeans1 = kmeans_result1[["cluster"]]#记录聚类结果
  ARI_kmeans1_design1[time]=adjustedRandIndex(labels_true,labels_pred_kmeans1)#计算调整兰德指数
  
  fviz_nbclust2<-fviz_nbclust(disjoint_simulation_data2_processing, kmeans, method = "wss")#画手肘图
  for (a in 1:10) {
    x[a]=fviz_nbclust2[["data"]][a,2]/fviz_nbclust2[["data"]][a+1,2]
  }
  for (b in 1:10) {
    k[1]=0
    k[b+1]=abs((x[b]-x[b+1])/x[b+1])
  }
  k_best2=which.max(k)#设置初始类别数
  kmeans_result2 <- kmeans(disjoint_simulation_data2_processing, centers = k_best2, nstart = 1)#调用k-means聚类算法
  cluster_kmeans2_design1[time]=length(unique(kmeans_result2[["cluster"]]))#记录聚类类别数
  labels_pred_kmeans2 = kmeans_result2[["cluster"]]#记录聚类结果
  ARI_kmeans2_design1[time]=adjustedRandIndex(labels_true,labels_pred_kmeans2)#计算调整兰德指数
  
  ARI_kmeans1_design1_RI[time]=rand.index(labels_true,labels_pred_kmeans1) 
  ARI_kmeans2_design1_RI[time]=rand.index(labels_true,labels_pred_kmeans2) 
  ARI_kmeans1_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_kmeans1) 
  ARI_kmeans2_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_kmeans2) 
  
  #---调用OPTICS聚类算法---#
  optics_result1<-optics(disjoint_simulation_data1_processing, minPts = 10)
  optics_result1<-extractDBSCAN(optics_result1, eps_cl = 0.07)
  optics_result1<-extractXi(optics_result1, xi = 0.01)
  cluster_optics1_design1[time]=length(unique(optics_result1[["cluster"]]))
  labels_pred_optics1 = optics_result1[["cluster"]]
  ARI_optics1_design1[time]=adjustedRandIndex(labels_true,labels_pred_optics1)
  
  optics_result2<-optics(disjoint_simulation_data2_processing, minPts = 10)
  optics_result2<-extractDBSCAN(optics_result2, eps_cl = 0.07)
  optics_result2<-extractXi(optics_result2, xi = 0.01)
  cluster_optics2_design1[time]=length(unique(optics_result2[["cluster"]]))
  labels_pred_optics2 = optics_result2[["cluster"]]
  ARI_optics2_design1[time]=adjustedRandIndex(labels_true,labels_pred_optics2)
  
  ARI_optics1_design1_RI[time]=rand.index(labels_true,labels_pred_optics1)
  ARI_optics2_design1_RI[time]=rand.index(labels_true,labels_pred_optics2)
  ARI_optics1_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_optics1)
  ARI_optics2_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_optics2)
  
  #---调用Affinity Propagation(AP)聚类算法---#
  
  labels_pred_ap1=c()
  labels_pred_ap2=c()
  
  ap_result1 <- apcluster(negDistMat(r=2),disjoint_simulation_data1_processing)#计算相似度矩阵，调用AP函数
  cluster_ap1_design1[time]=length(ap_result1@clusters)#记录聚类类别数
  #转换ap_result的聚类结果，与labels_true格式相对应
  for (ap_i in 1:51) {
    for (ap_j in 1:length(ap_result1@clusters)) {
      if( ap_i %in% ap_result1@clusters[[ap_j]]){labels_pred_ap1[ap_i]=ap_j}
    }
  } 
  labels_pred_ap1=as.numeric(labels_pred_ap1)#记录聚类结果
  ARI_ap1_design1[time]=adjustedRandIndex(labels_true,labels_pred_ap1)
  
  ap_result2 <- apcluster(negDistMat(r=2),disjoint_simulation_data2_processing)#计算相似度矩阵，调用AP函数
  cluster_ap2_design1[time]=length(ap_result2@clusters)#记录聚类类别数
  #转换ap_result的聚类结果，与labels_true格式相对应
  for (ap_i in 1:51) {
    for (ap_j in 1:length(ap_result2@clusters)) {
      if( ap_i %in% ap_result2@clusters[[ap_j]]){labels_pred_ap2[ap_i]=ap_j}
    }
  } 
  labels_pred_ap2=as.numeric(labels_pred_ap2)#记录聚类结果
  ARI_ap2_design1[time]=adjustedRandIndex(labels_true,labels_pred_ap2)
  
  ARI_ap1_design1_RI[time]=rand.index(labels_true,labels_pred_ap1)
  ARI_ap2_design1_RI[time]=rand.index(labels_true,labels_pred_ap2)
  ARI_ap1_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_ap1)
  ARI_ap2_design1_ARI[time]=adj.rand.index(labels_true,labels_pred_ap2)
  
  print(time)
 
}


ARI_avar_CDMFM1_design1=sum(ARI_CDMFM1_design1)/100 #看100次模拟的调整兰德系数的均值（越接近1说明聚类准确度越好）                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       cc
ARI_avar_CDMFM2_design1=sum(ARI_CDMFM2_design1)/100
ARI_avar_CDMFM1_non_design1=sum(ARI_CDMFM1_non_design1)/100
ARI_avar_CDMFM2_non_design1=sum(ARI_CDMFM2_non_design1)/100
ARI_avar_dbscan1_design1=sum(ARI_dbscan1_design1)/100
ARI_avar_dbscan2_design1=sum(ARI_dbscan2_design1)/100
ARI_avar_Mclust1_design1=sum(ARI_Mclust1_design1)/100
ARI_avar_Mclust2_design1=sum(ARI_Mclust2_design1)/100
ARI_avar_kmeans1_design1=sum(ARI_kmeans1_design1)/100
ARI_avar_kmeans2_design1=sum(ARI_kmeans2_design1)/100
ARI_avar_optics1_design1=sum(ARI_optics1_design1)/100
ARI_avar_optics2_design1=sum(ARI_optics2_design1)/100
ARI_avar_ap1_design1=sum(ARI_ap1_design1)/100
ARI_avar_ap2_design1=sum(ARI_ap2_design1)/100


ARI_avar_CDMFM1_design1_RI=sum(ARI_CDMFM1_design1_RI)/100 #看100次模拟的调整兰德系数的均值（越接近1说明聚类准确度越好）                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       cc
ARI_avar_CDMFM2_design1_RI=sum(ARI_CDMFM2_design1_RI)/100
ARI_avar_CDMFM1_non_design1_RI=sum(ARI_CDMFM1_non_design1_RI)/100
ARI_avar_CDMFM2_non_design1_RI=sum(ARI_CDMFM2_non_design1_RI)/100
ARI_avar_dbscan1_design1_RI=sum(ARI_dbscan1_design1_RI)/100
ARI_avar_dbscan2_design1_RI=sum(ARI_dbscan2_design1_RI)/100
ARI_avar_Mclust1_design1_RI=sum(ARI_Mclust1_design1_RI)/100
ARI_avar_Mclust2_design1_RI=sum(ARI_Mclust2_design1_RI)/100
ARI_avar_kmeans1_design1_RI=sum(ARI_kmeans1_design1_RI)/100
ARI_avar_kmeans2_design1_RI=sum(ARI_kmeans2_design1_RI)/100
ARI_avar_optics1_design1_RI=sum(ARI_optics1_design1_RI)/100
ARI_avar_optics2_design1_RI=sum(ARI_optics2_design1_RI)/100
ARI_avar_ap1_design1_RI=sum(ARI_ap1_design1_RI)/100
ARI_avar_ap2_design1_RI=sum(ARI_ap2_design1_RI)/100


ARI_avar_CDMFM1_design1_ARI=sum(ARI_CDMFM1_design1_ARI)/100 #看100次模拟的调整兰德系数的均值（越接近1说明聚类准确度越好）                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       cc
ARI_avar_CDMFM2_design1_ARI=sum(ARI_CDMFM2_design1_ARI)/100
ARI_avar_CDMFM1_non_design1_ARI=sum(ARI_CDMFM1_non_design1_ARI)/100
ARI_avar_CDMFM2_non_design1_ARI=sum(ARI_CDMFM2_non_design1_ARI)/100
ARI_avar_dbscan1_design1_ARI=sum(ARI_dbscan1_design1_ARI)/100
ARI_avar_dbscan2_design1_ARI=sum(ARI_dbscan2_design1_ARI)/100
ARI_avar_Mclust1_design1_ARI=sum(ARI_Mclust1_design1_ARI)/100
ARI_avar_Mclust2_design1_ARI=sum(ARI_Mclust2_design1_ARI)/100
ARI_avar_kmeans1_design1_ARI=sum(ARI_kmeans1_design1_ARI)/100
ARI_avar_kmeans2_design1_ARI=sum(ARI_kmeans2_design1_ARI)/100
ARI_avar_optics1_design1_ARI=sum(ARI_optics1_design1_ARI)/100
ARI_avar_optics2_design1_ARI=sum(ARI_optics2_design1_ARI)/100
ARI_avar_ap1_design1_ARI=sum(ARI_ap1_design1_ARI)/100
ARI_avar_ap2_design1_ARI=sum(ARI_ap2_design1_ARI)/100
