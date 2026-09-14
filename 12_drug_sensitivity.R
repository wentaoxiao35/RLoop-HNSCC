##清空数据
rm(list = ls())

##加载R包
library(tidyr)
library(dplyr)
library(tibble)
library(survival)
library(caret)
library("glmnet")
library(survminer)
library(export)
options(stringsAsFactors = F)


##读取药物相关数据
##首先，使用readxl包中的read_excel()函数，读取药物相关的数据。
##由于前7行为注释信息，因此使用参数skip进行跳过前7行。
library(readxl)
rt1 <- read_excel(path = "resource/DTP_NCI60_ZSCORE.xlsx", skip = 7)

##同时，将第一行作为列名，并去除末尾两列其余信息。
colnames(rt1) <- rt1[1,]
rt1 <- rt1[-1,-c(67,68)]

##筛选药物标准
##使用table()函数药物标准，结果显示，其中574种经过了临床试验，218种经过FDA批准。
table(rt1$`FDA status`)
##为了保证分析结果的可靠性，选取经过临床试验（Clinical trial）和FDA批准（FDA approved）的药物结果。
##当然，你也可以选择将所有的药物进行保留。最终，一共得到792个药物结果，并将其保存为txt文件用于后续分析。
rt1 <- rt1[rt1$`FDA status` %in% c("FDA approved", "Clinical trial"),]
rt1 <- rt1[,-c(1, 3:6)]
write.table(rt1, file = "output/drug.txt",sep = "\t",row.names = F,quote = F)


##同样的，读取表达数据，并对其进行整理和保存，用于后续的分析。
rt2 <- read_excel(path = "resource/RNA__RNA_seq_composite_expression.xls", skip = 9)
colnames(rt2) <- rt2[1,]
rt2 <- rt2[-1,-c(2:6)]
write.table(rt2, file = "output/geneExp.txt",sep = "\t",row.names = F,quote = F)

##药物敏感性分析
rm(list = ls())
library(impute)
library(limma)

##读取药物输入文件
##首先，读取前面保存的药物敏感性结果，设置相应的行名，并将其转换为矩阵形式
rt <- read.table("output/drug.txt",sep="\t",header=T,check.names=F)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
drug <- rt[,2:ncol(rt)]
dimnames <- list(rownames(drug),colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)
##考虑到药物敏感性数据中存在部分NA缺失值，通过impute.knn()函数来评估并补齐药物数据。
##其中，impute.knn()函数是一个使用最近邻平均来估算缺少的表达式数据的函数。
# 你好，正如报错提示，因为缺失值太多，knn法不能准确地补充缺失值，这个时候你需要先把缺失值过多的基因删除，
# 比如我会把缺失值多于60%的基因删除
# 可以用Na_ratio <- apply(data,2,function(x) sum(is.na(x))/length(x))>0.6计算得到各列确实值比例
# 然后data[-Na_ratio]即可

mat <- impute.knn(data)
drug <- mat$data
drug <- avereps(drug)


##读取表达输入文件
##同时，读取整理完成的NCI-60细胞系中基因表达情况。
##结果显示：其中包含了60种不同肿瘤细胞系，23808个基因的表达情况。
exp <- read.table("output/geneExp.txt", sep="\t", header=T, row.names = 1, check.names=F)
dim(exp)
exp[1:4, 1:4]
#4.提取特定基因表达
gene <- read.table("output/riskTrain.txt", header=T, row.names = 1, check.names=F)
gene <- gene %>% 
  dplyr::select(-c(futime,fustat,risk))
genelist <- as.vector(colnames(gene))
genelist
#将提前准备的目标基因列表进行读取。

genelist <- gsub(" ","",genelist)
genelist <- intersect(genelist,row.names(exp))
exp1 <- exp[genelist,]

exp1 <- as.data.frame(t(exp1))


#计算riskscore
rt=read.table("output/lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)  

#使用train组构建cox模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型相关信息
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
rs=cbind(id=row.names(outTab),outTab)

rs <- as.data.frame(rs)

cor <- rs %>% 
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
# 将 cor 数据框转换为列表，方便后续查找 cor 值
cor_list <- setNames(cor$coef, cor$id)

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))

# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}

# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore)

exp <- as.data.frame(t(exp1))



# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ##药物敏感性计算
# ##首先，新建一个空的数据框，用于保存后续分析结果。
# outTab <- data.frame()
# ##随后，使用for循环，分别计算每个基因表达与不同药物之间的Pearson相关系数。
# ##根据P值<0.05为界，对分析结果进行筛选，并将结果保存给变量outTab。
# ##结果显示：最终得到了63个相关性分析结果。
# for(Gene in row.names(exp)){
#   x <- as.numeric(exp[Gene,])
#   #对药物循环
#   for(Drug in row.names(drug)){
#     y <- as.numeric(drug[Drug,])
#     if (sd(y)!=0){
#       corT <- cor.test(x,y,method="pearson")
#       cor <- corT$estimate
#       pvalue <- corT$p.value
#       if(pvalue < 0.05){
#         outVector <- cbind(Gene,Drug,cor,pvalue)
#         outTab <- rbind(outTab,outVector)
#       }
#     }
#   }
# }
# 
# ##最后，将相关性分析结果进行输出保存。
# outTab <- outTab[order(as.numeric(as.vector(outTab$pvalue))),]
# write.table(outTab, file="output/drugCor.txt", sep="\t", row.names=F, quote=F)
# 
# 
# ##可视化
# ##在可视化分析之前，首先加载绘图使用的ggplot2包和ggpubr包。
# ##下面，我们提供两种可视化的方法。
# library(ggplot2)
# library(ggpubr)
# 
# ## 方法一：散点图 
# ##定义一个空的列表plotList_1用于保存输出结果。提取分析结果中最显著的前16个结果；
# ##当然，如果结果outTab的行数少于corPlotNum的话，则将其行数赋值给corPlotNum。
# plotList_1 <- list()
# outTab$cor <- as.numeric(outTab$cor)
# outTab$pvalue <- as.numeric(outTab$pvalue)
# 
# outTab <- outTab %>% 
#   arrange(cor)
# 
# 
# corPlotNum <- 20
# 
# 
# ##使用ggplot()函数结合for循环，逐个绘制散点图，进行可视化展示。
# for(i in 1:corPlotNum){
#   Gene <- outTab[i,1]
#   Drug <- outTab[i,2]
#   x <- as.numeric(exp[Gene,])
#   y <- as.numeric(drug[Drug,])
#   cor <- sprintf("%.03f",as.numeric(outTab[i,3]))
#   pvalue=0
#   if(as.numeric(outTab[i,4])<0.001){
#     pvalue="p<0.001"
#   }else{
#     pvalue=paste0("p=",sprintf("%.03f",as.numeric(outTab[i,4])))
#   }
#   df1 <- as.data.frame(cbind(x,y))
#   p1=ggplot(data = df1, aes(x = x, y = y))+
#     geom_point(size=1)+
#     stat_smooth(method="lm",se=FALSE, formula=y~x)+
#     labs(x="Expression",y="IC50",title = paste0(Gene,", ",Drug),subtitle = paste0("Cor=",cor,", ",pvalue))+
#     theme(axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank())+
#     theme_bw()
#   plotList_1[[i]]=p1
# }
# 
# 
# ##将绘制得到的结果按4x4的分布进行排列组合，并输出结果。
# ggarrange(plotlist=plotList_1,nrow=5,ncol=4)
# 
# ## 自由选择要展示的图的索引
# selected_indices <- c(1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)  ## 自由选择的图编号
# selected_plots <- plotList_1[selected_indices]
# 
# ## 将绘制得到的结果按3x3的分布进行排列组合
# ggarrange(plotlist = selected_plots, nrow = 6, ncol = 3)
# # graph2ppt(file = "output/Figure/drug",width=8.8, height=7.2)
# 
# library(export)
# # graph2ppt(file = "output/Figure/drug",width=11.5, height=2.5)
# 
# 
# 
# 
# ## 方法二：箱线图
# plotList_2<- list()
# corPlotNum<- 20
# if(nrow(outTab)<corPlotNum){corPlotNum=nrow(outTab)}
# 
# 
# ##同样，定义一个空的列表plotList_2用于保存输出结果。
# for(i in 1:corPlotNum){
#   Gene <- outTab[i,1]
#   Drug <- outTab[i,2]
#   x <- as.numeric(exp[Gene,])
#   y <- as.numeric(drug[Drug,])
#   df1 <- as.data.frame(cbind(x,y))
#   colnames(df1)[2] <- "IC50"
#   df1$group <- ifelse(df1$x > median(df1$x), "high", "low")
#   compaired <- list(c("low", "high"))
#   p1 <- ggboxplot(df1,
#                   x = "group", y = "IC50",
#                   fill = "group", palette = c("#00AFBB", "#E7B800"),
#                   add = "jitter", size = 0.5,
#                   xlab = paste0("The_expression_of_", Gene),
#                   ylab = paste0("IC50_of_", Drug)) +
#     stat_compare_means(comparisons = compaired,
#                        method = "wilcox.test", #设置统计方法
#                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#                                         symbols = c("***", "**", "*", "ns")))
#   plotList_2[[i]]=p1
# }
# 
# ##通过for循环，结合使用ggpubr包中的ggboxplot()函数，逐个绘制箱线图，
# ##并使用stat_compare_means()函数进行两组间统计分析，从而进行可视化展示。
# nrow<- ceiling(sqrt(corPlotNum))
# ncol<- ceiling(corPlotNum/nrow)
# ggarrange(plotlist=plotList_2,nrow=nrow,ncol=ncol)
# 
# 
# ## 自由选择要展示的图的索引
# selected_indices <- c(3,4,6,9,10)  ## 自由选择的图编号
# selected_plots <- plotList_2[selected_indices]
# 
# ## 将绘制得到的结果按3x3的分布进行排列组合
# ggarrange(plotlist = selected_plots, nrow = 1, ncol = 5)
# 
# library(export)
# # graph2ppt(file = "output/Figure/drug2",width=11.5, height=4)
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
# ################################################################################
index1 <- colnames(exp)
index2 <- colnames(drug)

drugdata <- rbind(exp,drug)

##首先，新建一个空的数据框，用于保存后续分析结果。
outTab <- data.frame()
##随后，使用for循环，分别计算每个基因表达与不同药物之间的Pearson相关系数。
##根据P值<0.05为界，对分析结果进行筛选，并将结果保存给变量outTab。
##结果显示：最终得到了63个相关性分析结果。
for(Gene in row.names(exp)){
  x <- as.numeric(exp[Gene,])
  #对药物循环
  for(Drug in row.names(drug)){
    y <- as.numeric(drug[Drug,])
    if (sd(y)!=0){
      corT <- cor.test(x,y,method="pearson")
      cor <- corT$estimate
      pvalue <- corT$p.value
      if(pvalue < 0.05){
        outVector <- cbind(Gene,Drug,cor,pvalue)
        outTab <- rbind(outTab,outVector)
      }
    }
  }
}

drugdata <- drugdata %>% 
  t() %>% 
  as.data.frame() 
###第353列有问题，太长了
drugdata <- drugdata[,-353]

median_value <- median(drugdata$riskscore, na.rm = TRUE)
drugdata$riskscore <- ifelse(drugdata$riskscore <= median_value, "low", "high")


library(ggpubr)
ggboxplot(
  drugdata, x = "riskscore", y = "Pluripotin",
  color = "riskscore", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(method = "wilcox.test")


my.wilcox = function(x){
  dd <- wilcox.test(drugdata[,x] ~ riskscore, data = drugdata)
  data.frame(gene=x,p.value=dd$p.value)
}
### 测试函数功能
my.wilcox("Pluripotin")

lapplylist = lapply(colnames(drugdata)[-c(1:1)],my.wilcox)
wilcox_data <- do.call(rbind,lapplylist)

outTab$pvalue <- as.numeric(outTab$pvalue)
outTab$cor <- as.numeric(outTab$cor)
wilcox_data$p.value <- as.numeric(wilcox_data$p.value)

wilcox_data_p <- wilcox_data %>% 
  filter(p.value<0.01)
index1 <- wilcox_data_p$gene
outTab_p <- outTab %>% 
  filter(pvalue<0.01) %>% 
  filter(cor<0)
index2 <- outTab_p$Drug
index <- intersect(index1,index2)
index <- c("riskscore",index)

################################################箱图相关性
data <- drugdata[,index]

Tamoxifen <- data[,c(2,1)]
colnames(Tamoxifen) <- c("riskScore","Group")
##因子化，以便排序
Tamoxifen$Group <- factor(Tamoxifen$Group, levels = c("low", "high"))
##比较P值，选取比较对象
levels(Tamoxifen$Group) <- c("Low-risk", "High-risk")

comparisons <- list(c("Low-risk", "High-risk"))
e <- ggplot(Tamoxifen, aes(x = Group, y = riskScore))

e + 
  geom_boxplot(
    width = 0.5,
    fill = c("#00AFBB", "#E7B800"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "black",
    shape = 16,
    position = position_jitter(0.15),
    alpha = 1,
    size = 3
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    label.y = max(Tamoxifen$riskScore, na.rm = TRUE) + 0.4
  ) +
  labs(
    x = "Group",
    y = "IC50 of Tamoxifen"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")
  )


graph2ppt(file = "output/Figure/9999_durg1_Tamoxifen_box.ppt",width=6, height=5.4)




####Drug2
HYPOTHEMYCIN <- data[,c(3,1)]
colnames(HYPOTHEMYCIN) <- c("riskScore","Group")
##因子化，以便排序
HYPOTHEMYCIN$Group <- factor(HYPOTHEMYCIN$Group, levels = c("low", "high"))
levels(HYPOTHEMYCIN$Group) <- c("Low-risk", "High-risk")

comparisons <- list(c("Low-risk", "High-risk"))


e <- ggplot(HYPOTHEMYCIN, aes(x = Group, y = riskScore))

e + 
  geom_boxplot(
    width = 0.5,
    fill = c("#00AFBB", "#E7B800"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "black",
    shape = 16,
    position = position_jitter(0.15),
    alpha = 1,
    size = 3
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    label.y = max(HYPOTHEMYCIN$riskScore, na.rm = TRUE) + 0.4
  ) +
  labs(
    x = "Group",
    y = "IC50 of HYPOTHEMYCIN"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")
  )
graph2ppt(file = "output/Figure/9999_durg2_HYPOTHEMYCIN_box.ppt",width=6, height=5.4)







###Drug3
DOLASTATIN_10 <- data[,c(4,1)]
colnames(DOLASTATIN_10) <- c("riskScore","Group")
##因子化，以便排序
DOLASTATIN_10$Group <- factor(DOLASTATIN_10$Group, levels = c("low", "high"))
##比较P值，选取比较对象
levels(DOLASTATIN_10$Group) <- c("Low-risk", "High-risk")

comparisons <- list(c("Low-risk", "High-risk"))

e <- ggplot(DOLASTATIN_10, aes(x = Group, y = riskScore))

e + 
  geom_boxplot(
    width = 0.5,
    fill = c("#00AFBB", "#E7B800"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "black",
    shape = 16,
    position = position_jitter(0.15),
    alpha = 1,
    size = 3
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    label.y = max(DOLASTATIN_10$riskScore, na.rm = TRUE) + 0.4
  ) +
  labs(
    x = "Group",
    y = "IC50 of DOLASTATIN_10"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")
  )
graph2ppt(file = "output/Figure/9999_durg3_DOLASTATIN_10_box.ppt",width=6, height=5.4)



####Drug4
Vorinostat <- data[,c(5,1)]
colnames(Vorinostat) <- c("riskScore","Group")
##因子化，以便排序
Vorinostat$Group <- factor(Vorinostat$Group, levels = c("low", "high"))
##比较P值，选取比较对象
levels(Vorinostat$Group) <- c("Low-risk", "High-risk")
comparisons <- list(c("Low-risk", "High-risk"))

e <- ggplot(Vorinostat, aes(x = Group, y = riskScore))
e + geom_boxplot(width = 0.5, fill = c("#00AFBB", "#E7B800"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(color = "black", shape = 16, position = position_jitter(0.2), alpha = 1, size = 3)  +
  labs(
    x = "Group",
    y = "IC50 of Vorinostat"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")  # 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test")
graph2ppt(file = "output/Figure/9999_durg4_Vorinostat_box.ppt",width=6, height=5.4)





####Drug5
Nilotinib <- data[,c(6,1)]
colnames(Nilotinib) <- c("riskScore","Group")
##因子化，以便排序
Nilotinib$Group <- factor(Nilotinib$Group, levels = c("low", "high"))
##比较P值，选取比较对象
levels(Nilotinib$Group) <- c("Low-risk", "High-risk")
comparisons <- list(c("Low-risk", "High-risk"))

e <- ggplot(Nilotinib, aes(x = Group, y = riskScore))
e + 
  geom_boxplot(
    width = 0.5,
    fill = c("#00AFBB", "#E7B800"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "black",
    shape = 16,
    position = position_jitter(0.15),
    alpha = 1,
    size = 3
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    label.y = max(Nilotinib$riskScore, na.rm = TRUE) + 0.4
  ) +
  labs(
    x = "Group",
    y = "IC50 of Nilotinib"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")
  )


graph2ppt(file = "output/Figure/9999_durg5_Nilotinib_box.ppt",width=6, height=5.4)





####Drug6
Domatinostat <- data[,c(7,1)]
colnames(Domatinostat) <- c("riskScore","Group")
##因子化，以便排序
Domatinostat$Group <- factor(Domatinostat$Group, levels = c("low", "high"))
##比较P值，选取比较对象
levels(Domatinostat$Group) <- c("Low-risk", "High-risk")
comparisons <- list(c("Low-risk", "High-risk"))

e <- ggplot(Domatinostat, aes(x = Group, y = riskScore))

e + 
  geom_boxplot(
    width = 0.5,
    fill = c("#00AFBB", "#E7B800"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "black",
    shape = 16,
    position = position_jitter(0.15),
    alpha = 1,
    size = 3
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    label.y = max(Domatinostat$riskScore, na.rm = TRUE) + 0.4
  ) +
  labs(
    x = "Group",
    y = "IC50 of Domatinostat"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")
  )
graph2ppt(file = "output/Figure/9999_durg6_Domatinostat_box.ppt",width=6, height=5.4)



####Drug7

AMG_900 <- data[,c(8,1)]
colnames(AMG_900) <- c("riskScore","Group")
##因子化，以便排序
AMG_900$Group <- factor(AMG_900$Group, levels = c("low", "high"))
##比较P值，选取比较对象
levels(AMG_900$Group) <- c("Low-risk", "High-risk")
comparisons <- list(c("Low-risk", "High-risk"))

e <- ggplot(AMG_900, aes(x = Group, y = riskScore))

e + 
  geom_boxplot(
    width = 0.5,
    fill = c("#00AFBB", "#E7B800"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "black",
    shape = 16,
    position = position_jitter(0.15),
    alpha = 1,
    size = 3
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format",
    label.y = max(AMG_900$riskScore, na.rm = TRUE) + 0.4
  ) +
  labs(
    x = "Group",
    y = "IC50 of AMG_900"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")
  )
graph2ppt(file = "output/Figure/9999_durg7_AMG_900_box.ppt",width=6, height=5.4)

# 
# Carmustine <- data[,c(13,1)]
# colnames(Carmustine) <- c("riskScore","clinical")
# ##因子化，以便排序
# Carmustine$clinical <- factor(Carmustine$clinical, levels = c("low", "high"))
# ##比较P值，选取比较对象
# comparisons <- list(c("Low-risk", "High-risk"))
# 
# e <- ggplot(Carmustine, aes(x = clinical, y = riskScore))
# e + geom_boxplot(width = 0.5, fill = c("#00AFBB", "#E7B800"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(color = "black", shape = 16, position = position_jitter(0.2), alpha = 1, size = 3) +
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),  # 将背景设为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     axis.line = element_line(colour = "black"),  # 添加坐标轴线条
#     axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
#     axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
#   )+
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test")
# # graph2ppt(file = "output/Figure/9999_durg8_box.ppt",width=6, height=5.4)
# 
# 
# Barasertib <- data[,c(24,1)]
# colnames(Barasertib) <- c("riskScore","clinical")
# ##因子化，以便排序
# Barasertib$clinical <- factor(Barasertib$clinical, levels = c("low", "high"))
# ##比较P值，选取比较对象
# comparisons <- list(c("Low-risk", "High-risk"))
# 
# e <- ggplot(Barasertib, aes(x = clinical, y = riskScore))
# e + geom_boxplot(width = 0.5, fill = c("#00AFBB", "#E7B800"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(color = "black", shape = 16, position = position_jitter(0.2), alpha = 1, size = 3) +
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),  # 将背景设为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     axis.line = element_line(colour = "black"),  # 添加坐标轴线条
#     axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
#     axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
#   )+
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test")
# # graph2ppt(file = "output/Figure/9999_durg9_box.ppt",width=6, height=5.4)
# 
# 
# Actinomycin_D <- data[,c(2,1)]
# colnames(Actinomycin_D) <- c("riskScore","clinical")
# ##因子化，以便排序
# Actinomycin_D$clinical <- factor(Actinomycin_D$clinical, levels = c("low", "high"))
# ##比较P值，选取比较对象
# comparisons <- list(c("Low-risk", "High-risk"))
# 
# e <- ggplot(Actinomycin_D, aes(x = clinical, y = riskScore))
# e + geom_boxplot(width = 0.5, fill = c("#00AFBB", "#E7B800"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(color = "black", shape = 16, position = position_jitter(0.2), alpha = 1, size = 3) +
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),  # 将背景设为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     axis.line = element_line(colour = "black"),  # 添加坐标轴线条
#     axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
#     axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
#   )+
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test")
# # graph2ppt(file = "output/Figure/9999_durg10_box.ppt",width=6, height=5.4)
# 
# 
# DOLASTATIN_10 <- data[,c(12,1)]
# colnames(DOLASTATIN_10) <- c("riskScore","clinical")
# ##因子化，以便排序
# DOLASTATIN_10$clinical <- factor(DOLASTATIN_10$clinical, levels = c("low", "high"))
# ##比较P值，选取比较对象
# comparisons <- list(c("Low-risk", "High-risk"))
# 
# e <- ggplot(DOLASTATIN_10, aes(x = clinical, y = riskScore))
# e + geom_boxplot(width = 0.5, fill = c("#00AFBB", "#E7B800"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(color = "black", shape = 16, position = position_jitter(0.2), alpha = 1, size = 3) +
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),  # 将背景设为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     axis.line = element_line(colour = "black"),  # 添加坐标轴线条
#     axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
#     axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
#   )+
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test")
# # graph2ppt(file = "output/Figure/9999_durg11_box.ppt",width=6, height=5.4)
# 
# AMG_900 <- data[,c(29,1)]
# colnames(AMG_900) <- c("riskScore","clinical")
# ##因子化，以便排序
# AMG_900$clinical <- factor(AMG_900$clinical, levels = c("low", "high"))
# ##比较P值，选取比较对象
# comparisons <- list(c("Low-risk", "High-risk"))
# 
# e <- ggplot(AMG_900, aes(x = clinical, y = riskScore))
# e + geom_boxplot(width = 0.5, fill = c("#00AFBB", "#E7B800"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(color = "black", shape = 16, position = position_jitter(0.2), alpha = 1, size = 3) +
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),  # 将背景设为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     axis.line = element_line(colour = "black"),  # 添加坐标轴线条
#     axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
#     axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
#   )+
#   stat_compare_means(comparisons = comparisons, method = "wilcox.test")
# # graph2ppt(file = "output/Figure/9999_durg12_box.ppt",width=6, height=5.4)

############################################绘制risk immune 相关性散点图
############################################绘制risk immune 相关性散点图
############################################绘制risk immune 相关性散点图
### 测试数据和相关性分析的方法,药物自己选，一个一个改吧！
drugdata <- rbind(exp,drug)
drugdata <- drugdata %>% 
  t() %>% 
  as.data.frame() 
drugdata <- drugdata[,-353]
data <- drugdata[,index]
colnames(data) <- gsub("-","_",colnames(data))
colnames(data) <- gsub(" ","_",colnames(data))
# ##R-loop里面改一下第七种药的名字
# colnames(data)[7] <- "Domatinostat"
ggplot(data,aes(riskscore,AMG_900))+
  geom_point(col="#45b2cb",size=3,alpha=0.7,stroke=1)+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=1.5,col="#45b2cb")+
  stat_cor(method = "pearson", digits = 3, size=5)+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    axis.title = element_text(size = 15),  # 调整坐标轴标题的字体大小
    axis.text = element_text(size = 15) ,
    axis.title.x = element_text(colour = "black", size = 16, face = "bold"),
    axis.title.y = element_text(colour = "black", size = 16, face = "bold")   # 调整坐标轴刻度标签的字体大小
  )
library(export)
graph2ppt(file = "output/Figure/99999_drug_AMG_900_jitter",width=6,height=5.4)
