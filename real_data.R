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
