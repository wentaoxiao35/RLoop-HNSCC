##清空数据
rm(list = ls())
library(survminer)
library(survival)
outside <- read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)     
data <- outside %>% 
  dplyr::select(futime,fustat,SERPINE1)

res.cut <- surv_cutpoint(data, time = "futime",
                         event = "fustat",
                         variables = names(data)[3:ncol(data)],
                         minprop = 0.1) #minprop = 0.3更好一些
res.cat <- surv_categorize(res.cut)
summary(res.cut)
str(res.cat)
###测试看看
genes <- "SERPINE1"
your.surv <- Surv(res.cat$futime, res.cat$fustat) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- res.cat[,genes] #分组
  survival_dat <- data.frame(group = group)
  group <- factor(group, levels = c("low", "high")) 
  fit <- survfit(your.surv ~ group)
  sdf <- survdiff(your.surv ~ group,rho=0)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
  p.val
  photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
                        legend.title = genes,#定义图例的名称
                        legend.labs = c("low","high"), #所以上面要因子化分组顺序
                        #legend = "top",#图例位置
                        pval = T, #在图上添加log rank检验的p值
                        #pval.method = TRUE,#添加p值的检验方法
                        conf.int = TRUE,#添加置信区间
                        risk.table = F, #在图下方添加风险表
                        #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                        risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                        #linetype = "strata", #改变不同组别的生存曲线的线型
                        #surv.median.line = "hv", #标注出中位生存时间
                        xlab = "Time in years", #x轴标题
                        xlim = c(0,max(res.cat$futime)+1), #展示x轴的范围
                        break.time.by = 2, #x轴间隔
                        size = 1.5, #线条大小
                        #ggtheme = theme_bw(), #为图形添加网格
                        palette = c("#2878b5", "#c82423")#图形颜色风格
  ) 
  
  photo2  #看一下图
  # 修改图例
  # 修改风险表的图例名称 
  photo2$table <- photo2$table 
  photo2  #再看一下图
}
##2.测试函数功能
your.km.plot("SERPINE1")
library(export)
graph2ppt(file = "output/Figure/outside_SERPINE1_KM.ppt",width=6.5, height=5.5)

##清空数据
rm(list = ls())
library(survminer)
library(survival)
outside <- read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1)     
data <- outside %>% 
  dplyr::select(futime,fustat,SERPINE1)

res.cut <- surv_cutpoint(data, time = "futime",
                         event = "fustat",
                         variables = names(data)[3:ncol(data)],
                         minprop = 0.1) #minprop = 0.3更好一些
res.cat <- surv_categorize(res.cut)
summary(res.cut)
str(res.cat)
###测试看看
genes <- "SERPINE1"
your.surv <- Surv(res.cat$futime, res.cat$fustat) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- res.cat[,genes] #分组
  survival_dat <- data.frame(group = group)
  group <- factor(group, levels = c("low", "high")) 
  fit <- survfit(your.surv ~ group)
  sdf <- survdiff(your.surv ~ group,rho=0)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
  p.val
  photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
                        legend.title = genes,#定义图例的名称
                        legend.labs = c("low","high"), #所以上面要因子化分组顺序
                        #legend = "top",#图例位置
                        pval = T, #在图上添加log rank检验的p值
                        #pval.method = TRUE,#添加p值的检验方法
                        conf.int = TRUE,#添加置信区间
                        risk.table = F, #在图下方添加风险表
                        #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                        risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                        #linetype = "strata", #改变不同组别的生存曲线的线型
                        #surv.median.line = "hv", #标注出中位生存时间
                        xlab = "Time in years", #x轴标题
                        xlim = c(0,max(res.cat$futime)+1), #展示x轴的范围
                        break.time.by = 2, #x轴间隔
                        size = 1.5, #线条大小
                        #ggtheme = theme_bw(), #为图形添加网格
                        palette = c("#2878b5", "#c82423")#图形颜色风格
  ) 
  
  photo2  #看一下图
  # 修改图例
  # 修改风险表的图例名称 
  photo2$table <- photo2$table 
  photo2  #再看一下图
}
##2.测试函数功能
your.km.plot("SERPINE1")
library(export)
graph2ppt(file = "output/Figure/outside2_SERPINE1_KM.ppt",width=6.5, height=5.5)

##清空数据
rm(list = ls())
library(survminer)
library(survival)
outside <- read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1)     
data <- outside %>% 
  dplyr::select(futime,fustat,SERPINE1)

res.cut <- surv_cutpoint(data, time = "futime",
                         event = "fustat",
                         variables = names(data)[3:ncol(data)],
                         minprop = 0.1) #minprop = 0.3更好一些
res.cat <- surv_categorize(res.cut)
summary(res.cut)
str(res.cat)
###测试看看
genes <- "SERPINE1"
your.surv <- Surv(res.cat$futime, res.cat$fustat) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- res.cat[,genes] #分组
  survival_dat <- data.frame(group = group)
  group <- factor(group, levels = c("low", "high")) 
  fit <- survfit(your.surv ~ group)
  sdf <- survdiff(your.surv ~ group,rho=0)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
  p.val
  photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
                        legend.title = genes,#定义图例的名称
                        legend.labs = c("low","high"), #所以上面要因子化分组顺序
                        #legend = "top",#图例位置
                        pval = T, #在图上添加log rank检验的p值
                        #pval.method = TRUE,#添加p值的检验方法
                        conf.int = TRUE,#添加置信区间
                        risk.table = F, #在图下方添加风险表
                        #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                        risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                        #linetype = "strata", #改变不同组别的生存曲线的线型
                        #surv.median.line = "hv", #标注出中位生存时间
                        xlab = "Time in years", #x轴标题
                        xlim = c(0,max(res.cat$futime)+1), #展示x轴的范围
                        break.time.by = 2, #x轴间隔
                        size = 1.5, #线条大小
                        #ggtheme = theme_bw(), #为图形添加网格
                        palette = c("#2878b5", "#c82423")#图形颜色风格
  ) 
  
  photo2  #看一下图
  # 修改图例
  # 修改风险表的图例名称 
  photo2$table <- photo2$table 
  photo2  #再看一下图
}
##2.测试函数功能
your.km.plot("SERPINE1")
library(export)
graph2ppt(file = "output/Figure/outside3_SERPINE1_KM.ppt",width=6.3, height=5.5)


