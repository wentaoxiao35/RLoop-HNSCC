rm(list = ls())
library(dplyr)
library(survival)
rt=read.table(file = "output//riskoutside1.txt", header=T, sep="\t", check.names=F)
# rownames(rt)=rt[,1]
# rt=rt[,-1]
rt <- rt %>% 
  dplyr::select(-c(futime,risk,riskScore))

#结果变量变为因子型
rt$fustat <- factor(rt$fustat)
dim(rt)

###使用经典的randomForest建立随机森林模型
library(randomForest)

# set.seed(1)
fit <- randomForest(fustat~., data = rt)

fit
###结果给出了树的数量：500颗；OOB错误率：23.91%；还给出了混淆矩阵。

# 下面是可视化整体错误率和树的数量的关系，可以看到随着树的数量增加，
# 错误率逐渐降低并渐趋平稳，中间的黑色线条是整体的错误率，
# 上下两条是结果变量中两个类别的错误率。

plot(fit)

# 可以看到结果有一个类别的错误率竟然是逐渐增加的，
# 因为我们这个数据的存在严重的类不平衡问题，
# 也就是结果变量中的两种类别差异很大：

table(rt$fustat)

# 类别0有226个，类别1只有71个，
# 模型为了提高整体准确率，就会牺牲掉类别为1的准确性~

# 查看整体错误率最小时有几棵树：

which.min(fit$err.rate[,1])

# 查看各个变量的重要性，这里给出了mean decrease gini，数值越大说明变量越重要：
importance(fit)

# 可视化变量重要性，可以使用默认的方法，如下所示，
# 也可以使用其他专门的R包

varImpPlot(fit)


library(export)
graph2ppt(file= paste0("output/随机森林（1）.ppt"),width=5,height=5)
# # 通过变量重要性，大家就可以选择比较重要的变量了。
# # 你可以选择前5个，前10个，
# # 或者大于所有变量性平均值(中位数，百分位数等)的变量等等。
# 
# # randomForest还提供了使用交叉验证法进行递归特征消除，
# # 筛选变量的方法：rfcv，下面是使用5折交叉验证进行递归特征消除：
# 
# set.seed(647)
# res <- rfcv(trainx = rt[,-1],trainy = rt[,1],
#             cv.fold = 5,
#             recursive = T
# )
# res$n.var #变量个数
# res$error.cv #错误率
# 
# # 可以看到在变量个数为7的时候，错误率是最小的(和59一样，但是肯定选简单的)。
# # 
# # 可视化这个结果，很明显变量个数为7(和59)的时候错误率最小：
# 
# with(res, plot(n.var, error.cv, type="o", lwd=2))
# # 结合上面的变量重要性，你可以选择前7个最重要的变量。
# 
# library(export)
# # graph2ppt(file= paste0("output/随机森林（2）.ppt"),width=5,height=5)















#https://www.jianshu.com/p/a369437cc334
#https://blog.csdn.net/weixin_52486108/article/details/136773664
#https://blog.csdn.net/dege857/article/details/135494384
rm(list = ls())
# if (!require(randomForestSRC)) install.packages("randomForestSRC")
# # 设置包的URL，这里假设你需要安装版本2.9.3
# package_url <- "https://cran.r-project.org/src/contrib/Archive/randomForestSRC/randomForestSRC_2.9.3.tar.gz"
# # 下载并安装旧版本的randomForestSRC包
# install.packages(package_url, repos = NULL, type = "source")
# if (!require(survival)) install.packages("survival")
# install.packages("ggRandomForests")
library(dplyr)
library(randomForestSRC)
library(survival)

rt=read.table(file = "output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)
rt <- rt %>%
  dplyr::select(-c(risk,riskScore))

#建立随机生存森林模型
###旧的安装包装好了，但是环境不兼容一跑就报错，去win上跑吧！
###旧版 randomForestSRC 和当前 R 环境不兼容
rfsrc_pbcmy <- rfsrc(Surv(futime, fustat) ~ .,
                     data = rt,
                     nsplit = 10,
                     na.action = "na.impute",
                     tree.err = TRUE,
                     splitrule='logrank',
                     proximity = T,
                     forest = T,
                     #ntime = 60,#时间节点
                     #mtry = 10,#变量数
                     ntree=1000,#树的数量
                     importance = TRUE#变量重要性VIMP
                     #block.size=100,#第几颗树的错误率
)

library(ggRandomForests)
#首先查看VIMP法，排名前15的变量重要性
gg_dta <- gg_vimp(rfsrc_pbcmy, nvar=15)
plot(gg_dta)
# library(export)
# graph2ppt(file = "output/Figure/randomForest1",width=5, height=5)
#最小深度法查看变量重要性
gg_dta <- gg_minimal_depth(rfsrc_pbcmy,lbls = st.labs)
plot(gg_dta)
# library(export)
# graph2ppt(file = "output/Figure/randomForest2",width=5, height=5)
#两种方法的结合  VIMP+min_depth
gg_dta <- gg_minimal_vimp(rfsrc_pbcmy)
plot(gg_dta)
# library(export)
# graph2ppt(file = "output/Figure/randomForest3",width=5, height=5)





library(randomForestSRC)
library(survival)

fit_test <- rfsrc(
  Surv(futime, fustat) ~ .,
  data = rt,
  ntree = 50
)

