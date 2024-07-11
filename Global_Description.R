#---洛伦兹曲线的比较---#
library(REAT)
par(mfrow = c(1,2))
lorenz(data$data.HINCP,  bg.col = "white", lctitle = "(a)")
## Minnesota
lorenz(data$data.HINCP[data$data.ST == 27], lctitle = "", bg.col = "white",
       lc.col =  "#E69F00", add.lc = TRUE, ltype = "dashed")
## New York
lorenz(data$data.HINCP[data$data.ST == 36], lctitle = "", bg.col = "white",
       lc.col = "#56B4E9", add.lc = TRUE, ltype = "dotdash")
legend("topleft", legend = c("National", "Minnesota", "New York"),
       col = c("black", "#E69F00", "#56B4E9"),
       lty = c("solid", "dashed", "longdash"))

lorenz(data$data.HINCP, lctitle = "(b)", bg.col = "white")
for (i in unique(data$data.ST)) {
  lorenz(data$data.HINCP[data$data.ST == i], lctitle = "", bg.col = "white",
         lc.col =  "grey", add.lc = TRUE, ltype = "dashed")
}
lorenz(data$data.HINCP, lctitle = "", bg.col = "white", add.lc = TRUE)

#---画gini系数图和平均收入图，用地图表示---#
library(purrr)
library(DescTools)
index <- unique(data$data.ST)
medincome <- map_dbl(index, ~median(data$data.HINCP[data$data.ST == .x]))
Gini_coefficient <- rep(0, length(index))
for (i in 1:51) {
  income_data_new<-subset(income_data,data.ST==state_fips[i],select = data.HINCP)%>%#选第i个洲的收入数据
    arrange(data.HINCP)%>%#给收入数据从小到大排序
    mutate(accumulate_population = rank(data.HINCP)/n())%>%#计算人口累计百分比
    mutate(accumulate_Income = cumsum(as.numeric(data.HINCP))/sum(data.HINCP))#计算收入累计百分比
  Gini_coefficient[i] = sum(2*(income_data_new$accumulate_population - income_data_new$accumulate_Income)/nrow(income_data_new))
}
stateinfo <- stateinfo %>% mutate(fips = as.numeric(fips)) %>%
  mutate(Gini = Gini_coefficient) %>%
  mutate(medincome = medincome)

gini <- plot_usmap(data = stateinfo, value = "Gini", labels = TRUE) +
  scale_fill_continuous(low = "white", high = "#004D40", name = "Gini coefficient") +
  ggtitle("(a)")

income <- plot_usmap(data = stateinfo, value = "medincome", labels = TRUE) +
  scale_fill_continuous(low = "white", high = "#004D40", name = "Median income") +
  ggtitle("(b)")
gridExtra::grid.arrange(gini, income, nrow = 1)


#---两种模拟设置---#
color=c("#2b9f78", "#5cb6ea", "#e5a11c")
plot_disjoint_map=plot_usmap(data = stateinfo_disjoint_init, values = "cluster", labels = TRUE) +
  scale_fill_manual(values = color)+ggtitle("(a)")
plot_joint_map=plot_usmap(data = stateinfo_joint_init, values = "cluster", labels = TRUE) +
  scale_fill_manual(values = color)+ggtitle("(b)")
gridExtra::grid.arrange(plot_disjoint_map, plot_joint_map, nrow = 1)

