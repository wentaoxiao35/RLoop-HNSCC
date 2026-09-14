##清空数据
rm(list = ls())

##加载R包
library(tidyr)
library(dplyr)
library(tibble)
library(survminer)
library(survival)
library(ggplot2)
library(ggprism)
library(tidyverse)
library(timeROC)
library(data.table)
library(pheatmap)
library(ggpubr)


############################################绘制SNV TMB survival 39
############################################绘制SNV TMB survival 39
############################################绘制SNV TMB survival 39
rm(list = ls())
##数据下载
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install(version='devel')
# 
# BiocManager::install("TCGAbiolinks")
# library(TCGAbiolinks)
# #
# query <- GDCquery(
#   project = "TCGA-HNSC",
#   data.category = "Simple Nucleotide Variation",
#   data.type = "Masked Somatic Mutation",
#   access = "open"
# )
# GDCdownload(query)
# GDCprepare(query, save = T,save.filename = "output/TCGA-HNSC_SNP.Rdata")


rm(list = ls())
#BiocManager::install("maftools")
library(maftools)
library(export)

load(file = "output/TCGA-HNSC_SNP.Rdata")

##此处有bug，同一患者可能有多个标本，取突变最多的标本
# data$ID <- substr(data$Tumor_Sample_Barcode,1,12)

new_data <- data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Count = n()) %>% 
  mutate(ID=substr(Tumor_Sample_Barcode,1,12)) %>% 
  arrange(desc(Count)) %>% 
  distinct(ID,.keep_all = T)

index_Tumor_Sample_Barcode <- new_data$Tumor_Sample_Barcode
data <- subset(data,Tumor_Sample_Barcode %in% index_Tumor_Sample_Barcode)
data$ID <- substr(data$Tumor_Sample_Barcode,1,12)


risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt <- risk %>% 
  dplyr::select(riskScore)


##此处也有bug有可能之前的高低分组中，有一些并没有突变数据，所以先把有突变数据的提出来，再分高低两组
index_all <- rownames(rt)
dat <- subset(data,ID %in% index_all)
index_res <- substr(dat$ID,1,12)
index_res <- unique(unlist(index_res))
rt <- rt %>% 
  mutate(TCGA_ID=rownames(rt)) %>% 
  dplyr::select(TCGA_ID,riskScore)

rt <- rt[index_res,]
rownames(rt) <- NULL
rt <- rt %>% 
  column_to_rownames("TCGA_ID")
median_value <- median(rt$riskScore, na.rm = TRUE)
rt$risk <- ifelse(rt$riskScore <= median_value, "low", "high")


high <- rt[rt$risk == "high",]
low <- rt[rt$risk == "low",]

index_high <- rownames(high)
dat <- subset(data,ID %in% index_high)
dat <- dplyr::select(dat,-ID)
dat <- dplyr::select(dat,-X1)
maf.coad <- dat

maf <- read.maf(maf.coad)
oncoplot(maf = maf, top = 15,fontSize = 0.8)
graph2ppt(file = "output/Figure/39_SNV_high.ppt",width=10, height=5.5)

tmb_high = tmb(maf = maf)
graph2ppt(file = "output/Figure/40_TMB_high.ppt",width=4.5, height=8)


index_low <- rownames(low)
dat <- subset(data,ID %in% index_low)
dat <- dplyr::select(dat,-ID)
maf.coad <- dat

maf <- read.maf(maf.coad)
oncoplot(maf = maf, top = 15,fontSize = 0.8)
graph2ppt(file = "output/Figure/41_SNV_low.ppt",width=10, height=5.5)

tmb_low = tmb(maf = maf)
graph2ppt(file = "output/Figure/42_TMB_low.ppt",width=4.5, height=8)


# 为两个数据框添加来源标识
tmb_high$Source <- 'high'
tmb_low$Source <- 'low'
# 合并数据框
combined_df <- rbind(tmb_low, tmb_high)
# 进行Wilcoxon检验
test_result <- wilcox.test(total_perMB_log ~ Source, data = combined_df, exact = FALSE)
# 输出检验结果
print(test_result)

# ##测试画图
# library(ggpubr)
# library(ggplot2)
# ggboxplot(
#   combined_df, x = "Source", y = "total_perMB_log",
#   color = "Source", palette = c("#2878b5", "#c82423"), 
#   add = "jitter"
# )+
#   stat_compare_means(method = "wilcox.test")
# 
# boxplot=ggviolin(combined_df, x="Source", y="total_perMB_log", fill="Source",
#                  xlab="",
#                  ylab="Tumor tmbation burden (log10)",
#                  legend.title="",
#                  palette = c("#2878b5","#c82423"),
#                  add = "boxplot", add.params = list(fill="white"))+
#   stat_compare_means(method = "wilcox.test",label.x = 1.35,label.y = 3)
# print(boxplot)
# library(export)
# # graph2ppt(file = "output/Figure/43_TMB_compare.ppt",width=5.4, height=5.4)

combined_df$Source <- factor(combined_df$Source, levels = c("low", "high"))
comparisons <- list(c("low", "high"))
e <- ggplot(combined_df, aes(x = Source, y = total_perMB_log))
e + geom_violin(aes(fill = Source), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = Source), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 4) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

graph2ppt(file = "output/Figure/43_TMB_compare.ppt",width=4.85, height=5.8)





#################TMB_survival
risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1) 
risk <- risk %>% 
  rownames_to_column("ID") %>% 
  dplyr::select(ID,riskScore,futime,fustat)

combined_df$ID <- substr(combined_df$Tumor_Sample_Barcode,1,12)
combined_df <- combined_df %>% 
  arrange(desc(total_perMB)) %>%
  distinct(ID,.keep_all = T) %>% 
  dplyr::select(ID,total_perMB)

data <- merge(risk,combined_df,by="ID")
colnames(data)[5] <- "TMB"

data <- data %>%
  mutate(risk = if_else(riskScore <= median(riskScore, na.rm = TRUE), "low", "high")) %>% 
  column_to_rownames("ID")

res.cut=surv_cutpoint(data, time = "futime", event = "fustat", variables =c("TMB"))
cutoff=as.numeric(res.cut$cutpoint[1])
tmbType=ifelse(data[,"TMB"]<=cutoff, "L-TMB", "H-TMB")

scoreType=ifelse(data$risk=="low", "L-Risk", "H-Risk")
mergeType=paste0(tmbType, " + ", scoreType)

bioSurvival=function(surData=null, outFile=null){
  diff=survdiff(Surv(futime, fustat) ~ group, data=surData)
  length=length(levels(factor(surData[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  bioCol=c("#223D6C","#D20A13","#6E568C","#088247","#7CC767","#c82423","#2878b5","#FFD121","#11AA4D")
  bioCol=bioCol[1:length]
  surPlot=ggsurvplot(fit, 
                     data=surData,
                     conf.int=F,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Group",
                     legend.labs=levels(factor(surData[,"group"])),
                     font.legend=10,
                     legend = "top",
                     xlab="Time(years)",
                     break.time.by = 2,
                     palette = bioCol,
                     #surv.median.line = "hv",
                     risk.table=F,
                     cumevents=F,
                     risk.table.height=.25)
  print(surPlot)
}

data$group <- tmbType
bioSurvival(surData = data)
graph2ppt(file = "output/Figure/44_TMB_KM.ppt",width=6.3, height=5.5)

data$group <- mergeType
bioSurvival(surData = data)
graph2ppt(file = "output/Figure/45_TMB+risk_KM.ppt",width=6.3, height=5.5)







############################高低突变森林图
rm(list = ls())
#BiocManager::install("maftools")
library(maftools)
library(export)

load(file = "output/TCGA-HNSC_SNP.Rdata")

##此处有bug，同一患者可能有多个标本，取突变最多的标本
# data$ID <- substr(data$Tumor_Sample_Barcode,1,12)

new_data <- data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Count = n()) %>% 
  mutate(ID=substr(Tumor_Sample_Barcode,1,12)) %>% 
  arrange(desc(Count)) %>% 
  distinct(ID,.keep_all = T)

index_Tumor_Sample_Barcode <- new_data$Tumor_Sample_Barcode
data <- subset(data,Tumor_Sample_Barcode %in% index_Tumor_Sample_Barcode)
data$ID <- substr(data$Tumor_Sample_Barcode,1,12)

risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt <- risk %>% 
  dplyr::select(riskScore)


##此处也有bug有可能之前的高低分组中，有一些并没有突变数据，所以先把有突变数据的提出来，再分高低两组
index_all <- rownames(rt)
dat <- subset(data,ID %in% index_all)
index_res <- substr(dat$ID,1,12)
index_res <- unique(unlist(index_res))
rt <- rt %>% 
  mutate(TCGA_ID=rownames(rt)) %>% 
  dplyr::select(TCGA_ID,riskScore)

rt <- rt[index_res,]
rownames(rt) <- NULL
rt <- rt %>% 
  column_to_rownames("TCGA_ID")
median_value <- median(rt$riskScore, na.rm = TRUE)
rt$risk <- ifelse(rt$riskScore <= median_value, "low", "high")

rt <- rt %>% 
  rownames_to_column("id")

high <- rt[rt$risk == "high",]
low <- rt[rt$risk == "low",]

index_high <- high[,"id"]
dat <- subset(data,ID %in% index_high)
dat <- dplyr::select(dat,-ID)
maf.coad <- dat

high_maf <- read.maf(maf.coad)


index_low <- low[,"id"]
dat <- subset(data,ID %in% index_low)
dat <- dplyr::select(dat,-ID)
maf.coad <- dat

low_maf <- read.maf(maf.coad)

#两个队列比较(MAFs)

##输入另一个 MAF 文件

#比较最少Mut个数为5的基因
pt.vs.rt <- mafCompare(m1 = high_maf, m2 = low_maf, m1Name = 'High', m2Name = 'Low', minMut = 10)
print(pt.vs.rt)

#1） Forest plots

forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('#2878b5', '#c82423'), geneFontSize = 0.8)

library(export)
graph2ppt(file = "output/Figure/46_SNP_compare.ppt",width=6.5, height=9)

