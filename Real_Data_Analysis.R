#---处理income data(转化成矩阵，51行20列)---#
#51代表美国洲的数量，20代表将收入等级划分为10类#
library(DirichletReg)
library(usmap)
library(ggplot2)
library(DescTools)
library(magrittr)
library(dplyr)
library(ggpubr)

#----将真实收入数据转换成收入等级矩阵---#

load("D:/Edge Download/income_data.Rdata")
income_data<-data
state_fips_name <- readRDS("D:/Edge Download/state_fips_name.rds")#载入洲名称数据
state_fips<-as.numeric(state_fips_name$fips)#将洲数据转化成数值型
Q<-25#收入等级超过30则会报错：too few positive probabilities
multinomial_data<-matrix(0,51,Q)#将收入等级设置为Q
for(i in 1:51){
  income_data_new<-subset(income_data,data.ST==state_fips[i],select = data.HINCP)%>%#选第i个洲的收入数据
    arrange(data.HINCP)%>%#给收入数据从小到大排序
    mutate(accumulate_population = rank(data.HINCP)/n())%>%#计算人口累计百分比
    mutate(accumulate_Income = cumsum(as.numeric(data.HINCP))/sum(data.HINCP)) #计算收入累计百分比
  for (j in 1:Q) {#求各洲的Q个收入等级的人口数
    multinomial_data[i,j]<-sum(income_data_new$accumulate_Income<j/Q&income_data_new$accumulate_Income>=(j-1)/Q)
  }
}






#---给real data 选最优lambda1---#

Q=25#real data的等级设置大一点，这样画lorenz curve会更平滑
n=51
T=5
loglikelihood=array(0,c(n,T,300))
loglikelihood_sum=c()
deviance_bar=c()
loglikelihood_dahlans=matrix(0,n,T)
deviance_bar_Dahl=c()
#---DIC---#
for (lambda1 in seq(0,3,0.1)) {#31个模型参数
  result_real=
    CDMFM_new1(data= multinomial_data_1 , niterations=500, lambda1=lambda1 , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
  #写loglikelihood function，将每次迭代的参数估计值代入
  for (q in 201:500) { #从第201次迭代开始计算,舍去前200次
    for (j in 1:T) {
      for (i in 1:n) {
        loglikelihood[i,j,q-200]=dmultinom(x=multinomial_data_1[i,,j],
                                            size = NULL, 
                                            prob=result_real[["Iterates"]][[q]][["phiout"]][result_real[["Iterates"]][[q]][["zout"]][i],],
                                            log = TRUE)#计算每次迭代时，洲数据对应概率估计的对数似然值
      }
    }
  }
  for (q in 201:500) {
    loglikelihood_sum[q-200]=sum(loglikelihood[,,q-200])
  }
  #Deviance=-2*loglikelihood
  deviance=-2*loglikelihood_sum#迭代300次的deviance被记录下来
  #D(theta)_bar
  deviance_bar[10*lambda1+1]=sum(deviance)/300#求迭代的deviance均值
  #D(theta_bar)#求theta_bar的deviance
  DahlAns<-getDahl(result_real,burn=200)
  for (j in 1:T) {
    for (i in 1:n) {
      loglikelihood_dahlans[i,j]=dmultinom(x=multinomial_data_1[i,,j],
                                             size = NULL, 
                                             prob= DahlAns[["phiout"]][DahlAns[["zout"]][i],],
                                             log = TRUE)#计算在最优迭代下，每个洲数据对应概率估计的对数似然值
     }
  }
  for (j in 1:T) {
    for (i in 1:n) {
      loglikelihood_sum_dahlans=sum(loglikelihood_dahlans)
    }
  }
  deviance_bar_Dahl[10*lambda1+1]=-2*loglikelihood_sum_dahlans
}

#pD——参数复杂度
pD=deviance_bar-deviance_bar_Dahl
#DIC
DIC=deviance_bar+pD
#根据最小的DIC选出最优的lambda1
lambda1_real_data=seq(0,3,0.1)[which.min(DIC)]

#---调用本文算法(CDMFM_new1)---#

result=
  CDMFM_new1(data=multinomial_data_1 , niterations=500, lambda1=lambda1_real_data , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
DahlAns<-getDahl(result,burn=200)

result_0=
  CDMFM_new1(data=multinomial_data_1 , niterations=500, lambda1=0 , neighbour=1 , distance=distance , alpha=rep(1:Q) , GAMMA = 1 , LAMBDA = 1 , initNClusters=5 , VN=VN)
DahlAns_0<-getDahl(result_0,burn=200)






#---画各种图---#

#install.packages("REAT")
#install.packages("patchwork")
library(REAT)
library(reshape2)
library(patchwork)


#--- 对real data计算Gini index、画lorenz curve ---#

#计算各洲的基尼系数
Gini=c()
for (i in 1:51) {
  income_data_new<-subset(income_data,data.ST==state_fips[i],select = data.HINCP)%>%#选第i个洲的收入数据
    arrange(data.HINCP)%>%#给收入数据从小到大排序
    mutate(accumulate_population = rank(data.HINCP)/n())%>%#计算人口累计百分比
    mutate(accumulate_Income = cumsum(as.numeric(data.HINCP))/sum(data.HINCP))#计算收入累计百分比
    Gini[i] = sum(2*(income_data_new$accumulate_population - income_data_new$accumulate_Income)/nrow(income_data_new))
}
#计算每个类别的平均基尼系数
gini=rep(0,length(unique(DahlAns$zout)))
t1=t2=t3=t4=t5=0
for (j in 1:51) {
  if (DahlAns$zout[j]==1) {
    gini[1]=gini[1]+Gini[j]
    t1=t1+1
  }
  if (DahlAns$zout[j]==2) {
    gini[2]=gini[2]+Gini[j]
    t2=t2+1
  }
  if (DahlAns$zout[j]==3) {
    gini[3]=gini[3]+Gini[j]
    t3=t3+1
  }
  if (DahlAns$zout[j]==4) {
    gini[4]=gini[4]+Gini[j]
    t4=t4+1
  }
  if (DahlAns$zout[j]==5) {
    gini[5]=gini[5]+Gini[j]
    t5=t5+1
  }
}
gini=c(gini[1]/t1,gini[2]/t2,gini[3]/t3,gini[4]/t4,gini[5]/t5)


#计算各洲的基尼系数
Gini_0=c()
for (i in 1:51) {
  income_data_new<-subset(income_data,data.ST==state_fips[i],select = data.HINCP)%>%#选第i个洲的收入数据
    arrange(data.HINCP)%>%#给收入数据从小到大排序
    mutate(accumulate_population = rank(data.HINCP)/n())%>%#计算人口累计百分比
    mutate(accumulate_Income = cumsum(as.numeric(data.HINCP))/sum(data.HINCP))#计算收入累计百分比
  Gini_0[i] = sum(2*(income_data_new$accumulate_population - income_data_new$accumulate_Income)/nrow(income_data_new))
}
#计算每个类别的平均基尼系数
gini_0=rep(0,length(unique(DahlAns_0$zout)))
t1=t2=t3=t4=t5=0
for (j in 1:51) {
  if (DahlAns_0$zout[j]==1) {
    gini_0[1]=gini_0[1]+Gini_0[j]
    t1=t1+1
  }
  if (DahlAns_0$zout[j]==2) {
    gini_0[2]=gini_0[2]+Gini_0[j]
    t2=t2+1
  }
  if (DahlAns_0$zout[j]==3) {
    gini_0[3]=gini_0[3]+Gini_0[j]
    t3=t3+1
  }
  if (DahlAns_0$zout[j]==4) {
    gini_0[4]=gini_0[4]+Gini_0[j]
    t4=t4+1
  }
  if (DahlAns_0$zout[j]==5) {
    gini_0[5]=gini_0[5]+Gini_0[j]
    t5=t5+1
  }
}
gini_0=c(gini_0[1]/t1,gini_0[2]/t2,gini_0[3]/t3,gini_0[4]/t4,gini_0[5]/t5)

#---根据类别画lorenz curve---#

lorenz_data=matrix(0,Q,length(unique(DahlAns$zout)))
for (i in 1:length(unique(DahlAns$zout))) {
  lorenz_data[,i]=cumsum(DahlAns$phiout[i,])
}#求累计概率
#为画图调整数据
lorenz_data=rbind(rep(0,length(unique(DahlAns$zout))),lorenz_data)
lorenz_data=cbind(lorenz_data,seq(0,1,0.04))
colnames(lorenz_data)=c("Cluster 1 (0.43)","Cluster 2 (0.46)","Cluster 3 (0.48)","Absolute average line")#调整将要画lorenz曲线的数据
lorenz_data=melt(lorenz_data)
colnames(lorenz_data)=c("y","Cluster","value")
lorenz_data$y=seq(0,1,0.04)
color=c("#f0e54b","#2b9f78", "#5cb6ea","#000000")
lorenz_curve=ggplot(lorenz_data)+
  geom_line(aes(x=value,y=y,color=Cluster,linetype=Cluster),size=0.6)+
  scale_color_manual(values = color)+
  scale_linetype_manual(values = c("dashed","dotdash","twodash","solid"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x = "% of Population", y = "% of Income",title = "(b)")


lorenz_data_0=matrix(0,Q,length(unique(DahlAns_0$zout)))
for (i in 1:length(unique(DahlAns_0$zout))) {
  lorenz_data_0[,i]=cumsum(DahlAns_0$phiout[i,])
}#求累计概率
#为画图调整数据
lorenz_data_0=rbind(rep(0,length(unique(DahlAns_0$zout))),lorenz_data_0)
lorenz_data_0=cbind(lorenz_data_0,seq(0,1,0.04))
colnames(lorenz_data_0)=c("Cluster 1 (0.43)","Cluster 2 (0.50)","Cluster 3 (0.46)","Cluster 4 (0.49)","Cluster 5 (0.48)","Absolute average line")#调整将要画lorenz曲线的数据
lorenz_data_0=melt(lorenz_data_0)
colnames(lorenz_data_0)=c("y","Cluster","value")
lorenz_data_0$y=seq(0,1,0.04)
color_0=c("#f0e54b","#e5a11c","#2b9f78", "#cc7daa","#5cb6ea","#000000")
lorenz_curve_0=ggplot(lorenz_data_0)+
  geom_line(aes(x=value,y=y,color=Cluster,linetype=Cluster),size=0.6)+
  scale_color_manual(values = color_0)+
  scale_linetype_manual(values = c("dashed","dotdash","twodash","dotted","longdash","solid"))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  labs(x = "% of Population", y = "% of Income",title = "(d)") 



#---根据收入等级矩阵的聚类结果画图---#

stateinfo <- readRDS("D:/Edge Download/state_fips_name.rds")#导入洲名数据
stateinfo$Cluster <- as.factor(DahlAns[["zout"]])#根据聚类结果给各洲赋上类别
usmap=plot_usmap(data = stateinfo, values = "Cluster", labels = TRUE) +
  scale_fill_manual(values = color)+#根据聚类结果画美国地图
  labs(title="(a)")

color_0=c("#f0e54b","#e5a11c","#2b9f78","#cc7daa","#5cb6ea")
stateinfo_0 <- readRDS("D:/Edge Download/state_fips_name.rds")#导入洲名数据
stateinfo_0$Cluster <- as.factor(DahlAns_0[["zout"]])#根据聚类结果给各洲赋上类别
usmap0=plot_usmap(data = stateinfo_0, values = "Cluster", labels = TRUE) +
  scale_fill_manual(values = color_0)+#根据聚类结果画美国地图
  labs(title="(c)")


#把四张图画在一起（两张地图+两张洛伦兹曲线）
usmap+lorenz_curve+usmap0+lorenz_curve_0+plot_layout(heights = c(1, 1),widths=c(2,1))
