##清空数据
rm(list = ls())
library(tidyr)
library(dplyr)
library(tibble)
library(export)
library(DESeq2)
library(ggpubr)
library(ggplot2)
library(survminer)
library(survival)
library(ggprism)
library(tidyverse)
library(timeROC)
library(data.table)
library(pheatmap)
library(rms)
library(foreign)
library(pec)
library(export)

cli=read.table("resource/TCGA.HNSC.CLI.txt",header=T,sep="\t",check.names=F,row.names=1)     
outside=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1) 


outside <- outside %>%  
  rownames_to_column("ID") %>% 
  dplyr::select(ID,fustat,riskScore)
cli <- cli %>% 
  rownames_to_column("ID")

clinical_final <- merge(outside,cli,by="ID")

median_value <- median(clinical_final$age)

clinical_final$age <- ifelse(clinical_final$age > 60, ">60", "<=60")

rt <- clinical_final %>% 
  column_to_rownames("ID")

age <- rt[,c(2,4)]
age$age <- gsub("<=60", "age<=60", age$age)
age$age <- gsub(">60", "age>60", age$age)
age <- na.omit(age)
age_less_60 <- age[age$age == "age<=60",]
age_more_60 <- age[age$age == "age>60",]
age <- rbind(age_less_60,age_more_60)
colnames(age) <- c("riskScore","clinical")
##因子化，以便排序
age$clinical <- factor(age$clinical, levels = c("age<=60", "age>60"))
##比较P值，选取比较对象
comparisons <- list(c("age<=60", "age>60"))


# e <- ggplot(age, aes(x = clinical, y = riskScore))
# 
# e + geom_violin(aes(fill = clinical), trim = FALSE) + 
#   geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
#   scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
#   scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
#   theme(
#     legend.position = "none",
#     panel.background = element_blank(),  # 将背景设为透明
#     panel.grid.major = element_blank(),  # 去掉主要网格线
#     panel.grid.minor = element_blank(),  # 去掉次要网格线
#     axis.line = element_line(colour = "black"),  # 添加坐标轴线条
#     axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
#     axis.title = element_text(colour = "black")  # 设置坐标轴标题颜色
#   )+
#   stat_compare_means(comparisons = comparisons,method = "wilcox.test")
###GPT
e + 
  geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(
    width = 0.5,
    fill = c("#4DBBD5", "#E64B35"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = clinical),
    shape = 16,
    position = position_jitter(0.2),
    alpha = 0.65,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498")) +
  scale_color_manual(values = c("#4DBBD5", "#E64B35")) +
  labs(
    x = "Age (years)",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    size = 6
  )

graph2ppt(file = "output/Figure/99_age_compare.ppt",width=4.85, height=5.8)



#####

gender <- rt[,c(2,3)]
gender <- na.omit(gender)
colnames(gender) <- c("riskScore","clinical")
##因子化，以便排序
gender$clinical <- factor(gender$clinical, levels = c("female", "male"))
##比较P值，选取比较对象
comparisons <- list(c("female", "male"))

e <- ggplot(gender, aes(x = clinical, y = riskScore))
e +  geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(
    width = 0.5,
    fill = c("#4DBBD5", "#E64B35"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = clinical),
    shape = 16,
    position = position_jitter(0.2),
    alpha = 0.65,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498")) +
  scale_color_manual(values = c("#4DBBD5", "#E64B35")) +
  labs(
    x = "Gender",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    size = 6
  )
 graph2ppt(file = "output/Figure/99_gender_compare.ppt",width=4.85, height=5.8)

lymphovascular_invasion <- rt[,c(2,13)]
lymphovascular_invasion <- lymphovascular_invasion %>% 
  filter(lymphovascular_invasion != "unknown")
lymphovascular_invasion <- na.omit(lymphovascular_invasion)
colnames(lymphovascular_invasion) <- c("riskScore","clinical")
##因子化，以便排序
lymphovascular_invasion$clinical <- factor(lymphovascular_invasion$clinical, levels = c("Negative", "Positive"))
##比较P值，选取比较对象
comparisons <- list(c("Negative", "Positive"))

e <- ggplot(lymphovascular_invasion, aes(x = clinical, y = riskScore))
e +  geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(
    width = 0.5,
    fill = c("#4DBBD5", "#E64B35"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = clinical),
    shape = 16,
    position = position_jitter(0.2),
    alpha = 0.65,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498")) +
  scale_color_manual(values = c("#4DBBD5", "#E64B35")) +
  labs(
    x = "lymphovascular invasion",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    size = 6
  )

graph2ppt(file = "output/Figure/99_lymphovascular_invasion_compare.ppt",width=4.85, height=5.8)


perineural_invasion <- rt[,c(2,14)]
perineural_invasion <- perineural_invasion %>% 
  filter(perineural_invasion != "unknown")
perineural_invasion <- na.omit(perineural_invasion)
colnames(perineural_invasion) <- c("riskScore","clinical")
##因子化，以便排序
perineural_invasion$clinical <- factor(perineural_invasion$clinical, levels = c("Negative", "Positive"))
##比较P值，选取比较对象
comparisons <- list(c("Negative", "Positive"))

e <- ggplot(perineural_invasion, aes(x = clinical, y = riskScore))
e +  geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(
    width = 0.5,
    fill = c("#4DBBD5", "#E64B35"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = clinical),
    shape = 16,
    position = position_jitter(0.2),
    alpha = 0.65,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498")) +
  scale_color_manual(values = c("#4DBBD5", "#E64B35")) +
  labs(
    x = "perineural invasion",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    size = 6
  )


graph2ppt(file = "output/Figure/99_perineural_invasion_compare.ppt",width=4.85, height=5.8)


fustat <- rt[,c(1,2)]
fustat$fustat <- gsub("0", "Alive", fustat$fustat)
fustat$fustat <- gsub("1", "Dead", fustat$fustat)
fustat <- na.omit(fustat)
Alive <- fustat[fustat$fustat == "Alive",]
Dead <- fustat[fustat$fustat == "Dead",]
fustat <- rbind(Alive,Dead)
colnames(fustat) <- c("clinical","riskScore")
##因子化，以便排序
fustat$clinical <- factor(fustat$clinical, levels = c("Alive", "Dead"))
##比较P值，选取比较对象
comparisons <- list(c("Alive", "Dead"))

e <- ggplot(fustat, aes(x = clinical, y = riskScore))

e +  geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(
    width = 0.5,
    fill = c("#4DBBD5", "#E64B35"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = clinical),
    shape = 16,
    position = position_jitter(0.2),
    alpha = 0.65,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498")) +
  scale_color_manual(values = c("#4DBBD5", "#E64B35")) +
  labs(
    x = "fustat",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    size = 6
  )


graph2ppt(file = "output/Figure/99_fustat_compare.ppt",width=4.85, height=5.8)


status <- rt[,c(2,16)]
status <- status %>% 
  filter(tumor_status != "unknown")
status$tumor_status <- gsub(" ","_",status$tumor_status)
WITH_TUMOR <- status[status$tumor_status == "With_tumor",]
TUMOR_FREE <- status[status$tumor_status == "Tumor_free",]
status <- rbind(WITH_TUMOR,TUMOR_FREE)
colnames(status) <- c("riskScore","clinical")
##因子化，以便排序
status$clinical <- factor(status$clinical, levels = c("Tumor_free", "With_tumor"))
##比较P值，选取比较对象
comparisons <- list(c("Tumor_free", "With_tumor"))

e <- ggplot(status, aes(x = clinical, y = riskScore))
e +  geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(
    width = 0.5,
    fill = c("#4DBBD5", "#E64B35"),
    color = "black",
    alpha = 1,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(color = clinical),
    shape = 16,
    position = position_jitter(0.2),
    alpha = 0.65,
    size = 4.2
  ) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498")) +
  scale_color_manual(values = c("#4DBBD5", "#E64B35")) +
  labs(
    x = "Tumor status",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.text.y = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    size = 6
  )


graph2ppt(file = "output/Figure/99_status_compare.ppt",width=4.85, height=5.8)

pathologic_t <- rt[,c(2,9)]
pathologic_t  <- pathologic_t  %>% 
  filter(pathologic_T != "unknown")

pathologic_t$pathologic_T <- gsub("T1", "I", pathologic_t$pathologic_T)
pathologic_t$pathologic_T <- gsub("T2", "II", pathologic_t$pathologic_T)
pathologic_t$pathologic_T <- gsub("T3", "III", pathologic_t$pathologic_T)
pathologic_t$pathologic_T <- gsub("T4", "IV", pathologic_t$pathologic_T)
pathologic_t <- na.omit(pathologic_t)

I <- pathologic_t[pathologic_t$pathologic_T == "I",]
II <- pathologic_t[pathologic_t$pathologic_T == "II",]
III <- pathologic_t[pathologic_t$pathologic_T == "III",]
IV <- pathologic_t[pathologic_t$pathologic_T == "IV",]
pathologic_t <- rbind(I,II,III,IV)
colnames(pathologic_t) <- c("riskScore","clinical")

##因子化，以便排序
pathologic_t$clinical <- factor(pathologic_t$clinical, levels = c("I","II","III","IV"))
##比较P值，选取比较对象
comparisons <- list(c("I","II"), c("I","III"), c("I","IV"), c("II","III"), c("II","IV"), c("III", "IV"))

e <- ggplot(pathologic_t, aes(x = clinical, y = riskScore))


e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35", "#6E7C7D", "#9D84BE"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498", "#BABEBF", "#D1C7E2"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35", "#6E7C7D", "#9D84BE"))+
  labs(
    x = "pathologic T",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")









graph2ppt(file = "output/Figure/99_pathologic_t_compare.ppt",width=7.2, height=6.8)



pathologic_n <- rt[,c(2,10)]
pathologic_n  <- pathologic_n  %>% 
  filter(pathologic_N != "unknown")
pathologic_n$pathologic_N <- gsub("N0", "N-", pathologic_n$pathologic_N)
pathologic_n$pathologic_N <- gsub("N1", "N+", pathologic_n$pathologic_N)
pathologic_n$pathologic_N <- gsub("N2", "N+", pathologic_n$pathologic_N)
pathologic_n$pathologic_N <- gsub("N3", "N+", pathologic_n$pathologic_N)


N_0 <- pathologic_n[pathologic_n$pathologic_N == "N-",]
N_1 <- pathologic_n[pathologic_n$pathologic_N == "N+",]
pathologic_n <- rbind(N_0,N_1)
colnames(pathologic_n) <- c("riskScore","clinical")

##因子化，以便排序
pathologic_n$clinical <- factor(pathologic_n$clinical, levels = c("N-","N+"))
##比较P值，选取比较对象
comparisons <- list(c("N-","N+"))

e <- ggplot(pathologic_n, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  labs(
    x = "pathologic N",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

 graph2ppt(file = "output/Figure/99_pathologic_n_compare.ppt",width=4.85, height=5.8)



pathologic_m <- rt[,c(2,11)]
pathologic_m  <- pathologic_m  %>% 
  filter(pathologic_M != "unknown")
pathologic_m$pathologic_M <- gsub("M0", "M-", pathologic_m$pathologic_M)
pathologic_m$pathologic_M <- gsub("M1", "M+", pathologic_m$pathologic_M)

M_0 <- pathologic_m[pathologic_m$pathologic_M == "M-",]
M_1 <- pathologic_m[pathologic_m$pathologic_M == "M+",]
pathologic_m <- rbind(M_0,M_1)
colnames(pathologic_m) <- c("riskScore","clinical")

##因子化，以便排序
pathologic_m$clinical <- factor(pathologic_m$clinical, levels = c("M-","M+"))
##比较P值，选取比较对象
comparisons <- list(c("M-","M+"))

e <- ggplot(pathologic_m, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  labs(
    x = "pathologic M",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

 graph2ppt(file = "output/Figure/99_pathologic_m_compare.ppt",width=4.85, height=5.8)

pathologic_stage <- rt[,c(2,12)]
pathologic_stage  <- pathologic_stage  %>% 
  filter(pathologic_stage != "unknown")
pathologic_stage$pathologic_stage <- gsub(" ", "", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("StageIV", "4", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("StageIII", "3", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("StageII", "2", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("StageI", "1", pathologic_stage$pathologic_stage)

pathologic_stage$pathologic_stage <- gsub("1", "I+II", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("2", "I+II", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("3", "III+IV", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("4", "III+IV", pathologic_stage$pathologic_stage)

I_II <- pathologic_stage[pathologic_stage$pathologic_stage == "I+II",]
III_IV <- pathologic_stage[pathologic_stage$pathologic_stage == "III+IV",]
pathologic_stage <- rbind(I_II,III_IV)
colnames(pathologic_stage) <- c("riskScore","clinical")

##因子化，以便排序
pathologic_stage$clinical <- factor(pathologic_stage$clinical, levels = c("I+II","III+IV"))
##比较P值，选取比较对象
comparisons <- list(c("I+II", "III+IV"))

e <- ggplot(pathologic_stage, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  labs(
    x = "pathologic Stage",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

 graph2ppt(file = "output/Figure/99_pathologic_stage_compare.ppt",width=4.85, height=5.8)




clinical_t <- rt[,c(2,5)]
clinical_t  <- clinical_t  %>% 
  filter(clinical_T != "unknown")

clinical_t$clinical_T <- gsub("T1", "I", clinical_t$clinical_T)
clinical_t$clinical_T <- gsub("T2", "II", clinical_t$clinical_T)
clinical_t$clinical_T <- gsub("T3", "III", clinical_t$clinical_T)
clinical_t$clinical_T <- gsub("T4", "IV", clinical_t$clinical_T)
clinical_t <- na.omit(clinical_t)

I <- clinical_t[clinical_t$clinical_T == "I",]
II <- clinical_t[clinical_t$clinical_T == "II",]
III <- clinical_t[clinical_t$clinical_T == "III",]
IV <- clinical_t[clinical_t$clinical_T == "IV",]
clinical_t <- rbind(I,II,III,IV)
colnames(clinical_t) <- c("riskScore","clinical")

##因子化，以便排序
clinical_t$clinical <- factor(clinical_t$clinical, levels = c("I","II","III","IV"))
##比较P值，选取比较对象
comparisons <- list(c("I","II"), c("I","III"), c("I","IV"), c("II","III"), c("II","IV"), c("III", "IV"))

e <- ggplot(clinical_t, aes(x = clinical, y = riskScore))


e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35", "#6E7C7D", "#9D84BE"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498", "#BABEBF", "#D1C7E2"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35", "#6E7C7D", "#9D84BE"))+
  labs(
    x = "clinical T",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

graph2ppt(file = "output/Figure/99_clinical_t_compare.ppt",width=7.2, height=6.8)



clinical_n <- rt[,c(2,6)]
clinical_n  <- clinical_n  %>% 
  filter(clinical_N != "unknown")
clinical_n$clinical_N <- gsub("N0", "N-", clinical_n$clinical_N)
clinical_n$clinical_N <- gsub("N1", "N+", clinical_n$clinical_N)
clinical_n$clinical_N <- gsub("N2", "N+", clinical_n$clinical_N)
clinical_n$clinical_N <- gsub("N3", "N+", clinical_n$clinical_N)


N_0 <- clinical_n[clinical_n$clinical_N == "N-",]
N_1 <- clinical_n[clinical_n$clinical_N == "N+",]
clinical_n <- rbind(N_0,N_1)
colnames(clinical_n) <- c("riskScore","clinical")

##因子化，以便排序
clinical_n$clinical <- factor(clinical_n$clinical, levels = c("N-","N+"))
##比较P值，选取比较对象
comparisons <- list(c("N-","N+"))

e <- ggplot(clinical_n, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  labs(
    x = "clinical N",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

 graph2ppt(file = "output/Figure/99_clinical_n_compare.ppt",width=4.85, height=5.8)



clinical_m <- rt[,c(2,7)]
clinical_m  <- clinical_m  %>% 
  filter(clinical_M != "unknown")
clinical_m$clinical_M <- gsub("M0", "M-", clinical_m$clinical_M)
clinical_m$clinical_M <- gsub("M1", "M+", clinical_m$clinical_M)

M_0 <- clinical_m[clinical_m$clinical_M == "M-",]
M_1 <- clinical_m[clinical_m$clinical_M == "M+",]
clinical_m <- rbind(M_0,M_1)
colnames(clinical_m) <- c("riskScore","clinical")

##因子化，以便排序
clinical_m$clinical <- factor(clinical_m$clinical, levels = c("M-","M+"))
##比较P值，选取比较对象
comparisons <- list(c("M-","M+"))

e <- ggplot(clinical_m, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  labs(
    x = "clinical M",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

 graph2ppt(file = "output/Figure/99_clinical_m_compare.ppt",width=4.85, height=5.8)

clinical_stage <- rt[,c(2,8)]
clinical_stage  <- clinical_stage  %>% 
  filter(clinical_stage != "unknown")
clinical_stage$clinical_stage <- gsub(" ", "", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("StageIV", "4", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("StageIII", "3", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("StageII", "2", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("StageI", "1", clinical_stage$clinical_stage)

clinical_stage$clinical_stage <- gsub("1", "I+II", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("2", "I+II", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("3", "III+IV", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("4", "III+IV", clinical_stage$clinical_stage)

I_II <- clinical_stage[clinical_stage$clinical_stage == "I+II",]
III_IV <- clinical_stage[clinical_stage$clinical_stage == "III+IV",]
clinical_stage <- rbind(I_II,III_IV)
colnames(clinical_stage) <- c("riskScore","clinical")

##因子化，以便排序
clinical_stage$clinical <- factor(clinical_stage$clinical, levels = c("I+II","III+IV"))
##比较P值，选取比较对象
comparisons <- list(c("I+II", "III+IV"))

e <- ggplot(clinical_stage, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 5) +
  scale_fill_manual(values = c("#A5DDEA", "#F1A498"))+
  scale_color_manual(values = c("#4DBBD5", "#E64B35"))+
  labs(
    x = "clinical Stage",
    y = "riskScore"
  ) +
  theme(
    legend.position = "none",
    panel.background = element_blank(),  # 将背景设为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    axis.line = element_line(colour = "black"),  # 添加坐标轴线条
    axis.text = element_text(colour = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色
    axis.title.x = element_text(colour = "black", face = "bold", size = 16),
    axis.title.y = element_text(colour = "black", face = "bold", size = 16)# 设置坐标轴标题颜色
  )+
  stat_compare_means(comparisons = comparisons,method = "wilcox.test")

 graph2ppt(file = "output/Figure/99_clinical_stage_compare.ppt",width=4.85, height=5.8)




#####################################################生存亚组分析####
##清空数据
rm(list = ls())
risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)     
clinical=read.table("resource/TCGA.HNSC.CLI.txt",header=T,sep="\t",check.names=F,row.names=1) 
index <- rownames(risk) 
clinical_new <- clinical[index,]
clinical_new <- clinical_new %>% 
  rownames_to_column("TCGA_ID")
risk <- risk %>% 
  rownames_to_column("TCGA_ID")
clinical_final <- merge(risk,clinical_new,by="TCGA_ID")
clinical_final$age <- ifelse(clinical_final$age > 60, ">60", "<=60")
clinical <- clinical_final %>% 
  dplyr::select(TCGA_ID,futime,fustat,age,gender,tumor_status,
                pathologic_T,pathologic_N,pathologic_M,pathologic_stage,
                clinical_T,clinical_N,clinical_M,clinical_stage,
                riskScore) %>% 
  column_to_rownames("TCGA_ID")

rt <- clinical %>% 
  dplyr::select(riskScore,everything())

####
age <- rt[,c(1:3,4)]
age <- age[!(apply(age, 1, function(x) any(x == "unknown"))), ]
age$age <- gsub("<=60", "age<=60", age$age)
age$age <- gsub(">60", "age>60", age$age)
##########age<=60
age <- age %>% 
  filter(age=="age<=60") %>% 
  dplyr::select(-age)
median_score <- median(age$riskScore)
age$riskScore <- ifelse(age$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(age$futime, age$fustat) 
####以这个KM函数为主！####
##2026-04-02
####以这个KM函数为主！####
####以这个KM函数为主！####
####以这个KM函数为主！####
your.km.plot <- function(genes, data, plot_title = NULL) {
  print(genes)
  
  # 取分组变量
  group <- data[[genes]]
  group <- factor(group, levels = c("low", "high"))
  
  # 复制一份数据并加入group
  data2 <- data
  data2$group <- group
  
  # 生存对象
  your.surv <- Surv(data2$futime, data2$fustat)
  
  # 拟合
  fit <- survfit(your.surv ~ group, data = data2)
  
  # log-rank
  sdf <- survdiff(your.surv ~ group, data = data2, rho = 0)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  print(p.val)
  
  # 作图
  photo2 <- ggsurvplot(
    fit,
    data = data2,
    legend.title = genes,
    legend.labs = c("low", "high"),
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    xlab = "Time in years",
    xlim = c(0, max(data2$futime, na.rm = TRUE) + 1),
    break.time.by = 1,
    size = 1.5,
    palette = c("#2878b5", "#c82423")
  )
  
  # 标题要加在主图上，不是直接 + ggtitle 到 ggsurvplot对象
  if (!is.null(plot_title)) {
    photo2$plot <- photo2$plot + ggtitle(plot_title)
  }
  
# 风险表标题
photo2$table <- photo2$table + labs(title = "Number at risk")
  
  print(photo2)
  return(photo2)
}
##2.测试函数功能
your.km.plot("riskScore",data = age, plot_title = "Age <= 60")
library(export)
graph2ppt(file = "output/Figure/23_age_under60_risk_KM.ppt",width=6.3, height=5.5)

###
age <- rt[,c(1:3,4)]
age <- age[!(apply(age, 1, function(x) any(x == "unknown"))), ]
age$age <- gsub("<=60", "age<=60", age$age)
age$age <- gsub(">60", "age>60", age$age)
##########age>60
age <- age %>% 
  filter(age=="age>60") %>% 
  dplyr::select(-age)
median_score <- median(age$riskScore)
age$riskScore <- ifelse(age$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(age$futime, age$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- age[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(age$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("Age>60")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = age, plot_title = "Age > 60")
library(export)
 graph2ppt(file = "output/Figure/23_age_up60_risk_KM.ppt",width=6.3, height=5.5)


gender <- rt[,c(1:3,5)]
gender <- gender[!(apply(gender, 1, function(x) any(x == "unknown"))), ]
##########gender
gender <- gender %>%
  filter(gender=="female") %>%
  dplyr::select(-gender)
median_score <- median(gender$riskScore)
gender$riskScore <- ifelse(gender$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(gender$futime, gender$fustat)
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- gender[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high"))
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(gender$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")
#   )+ ggtitle("female")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = gender, plot_title = "female")
library(export)
 graph2ppt(file = "output/Figure/35_female_risk_KM.ppt",width=6.3, height=5.5)

####
gender <- rt[,c(1:3,5)]
gender <- gender[!(apply(gender, 1, function(x) any(x == "unknown"))), ]
##########gender
gender <- gender %>%
  filter(gender=="male") %>%
  dplyr::select(-gender)
median_score <- median(gender$riskScore)
gender$riskScore <- ifelse(gender$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(gender$futime, gender$fustat)
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- gender[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high"))
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(gender$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")
#   )+ ggtitle("male")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = gender, plot_title = "male")
library(export)
graph2ppt(file = "output/Figure/35_male_risk_KM.ppt",width=6.3, height=5.5)


tumor_status <- rt[,c(1:3,6)]
tumor_status <- tumor_status[!(apply(tumor_status, 1, function(x) any(x == "unknown"))), ]
##########tumor_status
tumor_status <- tumor_status %>%
  filter(tumor_status=="Tumor free") %>%
  dplyr::select(-tumor_status)
median_score <- median(tumor_status$riskScore)
tumor_status$riskScore <- ifelse(tumor_status$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(tumor_status$futime, tumor_status$fustat)
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- tumor_status[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high"))
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(tumor_status$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")
#   )+ ggtitle("Tumor free")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = tumor_status, plot_title = "Tumor free")
library(export)
graph2ppt(file = "output/Figure/35_tumor_free_risk_KM.ppt",width=6.3, height=5.5)

####
tumor_status <- rt[,c(1:3,6)]
tumor_status <- tumor_status[!(apply(tumor_status, 1, function(x) any(x == "unknown"))), ]
##########tumor_status
tumor_status <- tumor_status %>%
  filter(tumor_status=="With tumor") %>%
  dplyr::select(-tumor_status)
median_score <- median(tumor_status$riskScore)
tumor_status$riskScore <- ifelse(tumor_status$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(tumor_status$futime, tumor_status$fustat)
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- tumor_status[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high"))
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(tumor_status$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")
#   )+ ggtitle("With tumor")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = tumor_status, plot_title = "With tumor")
library(export)
 graph2ppt(file = "output/Figure/35_with_tumor_risk_KM.ppt",width=6.3, height=5.5)














####
pathologic_T <- rt[,c(1:3,7)]
pathologic_T <- pathologic_T[!(apply(pathologic_T, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "T1/T2" 替换为 "T1+T2"，"T3/T4" 替换为 "T3+T4"
pathologic_T$pathologic_T <- gsub("T1", "I", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("T2", "I", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("T3", "II", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("T4", "II", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("II", "T3+T4", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("I", "T1+T2", pathologic_T$pathologic_T)
##########pathologic_T
pathologic_T <- pathologic_T %>% 
  filter(pathologic_T=="T1+T2") %>% 
  dplyr::select(-pathologic_T)
median_score <- median(pathologic_T$riskScore)
pathologic_T$riskScore <- ifelse(pathologic_T$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_T$futime, pathologic_T$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_T[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_T$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("pathologic_T1T2")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_T, plot_title = "pathologic_T1T2")
library(export)
graph2ppt(file = "output/Figure/27_pT1T2_risk_KM.ppt",width=6.3, height=5.5)

####
pathologic_T <- rt[,c(1:3,7)]
pathologic_T <- pathologic_T[!(apply(pathologic_T, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "T1/T2" 替换为 "T1+T2"，"T3/T4" 替换为 "T3+T4"
pathologic_T$pathologic_T <- gsub("T1", "I", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("T2", "I", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("T3", "II", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("T4", "II", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("II", "T3+T4", pathologic_T$pathologic_T)
pathologic_T$pathologic_T <- gsub("I", "T1+T2", pathologic_T$pathologic_T)
##########pathologic_T
pathologic_T <- pathologic_T %>% 
  filter(pathologic_T=="T3+T4") %>% 
  dplyr::select(-pathologic_T)
median_score <- median(pathologic_T$riskScore)
pathologic_T$riskScore <- ifelse(pathologic_T$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_T$futime, pathologic_T$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_T[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_T$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("pathologic_T3T4")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_T, plot_title = "pathologic_T3T4")
library(export)
 graph2ppt(file = "output/Figure/28_pT3T4_risk_KM.ppt",width=6.3, height=5.5)


####
pathologic_N <- rt[,c(1:3,8)]
pathologic_N <- pathologic_N[!(apply(pathologic_N, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "N0" 替换为 "N-"，"N1/2" 替换为 "N+"
pathologic_N$pathologic_N <- gsub("N0", "N-", pathologic_N$pathologic_N)
pathologic_N$pathologic_N <- gsub("N1", "N+", pathologic_N$pathologic_N)
pathologic_N$pathologic_N <- gsub("N2", "N+", pathologic_N$pathologic_N)
pathologic_N$pathologic_N <- gsub("N3", "N+", pathologic_N$pathologic_N)
##########pathologic_N
pathologic_N <- pathologic_N %>% 
  filter(pathologic_N=="N-") %>% 
  dplyr::select(-pathologic_N)
median_score <- median(pathologic_N$riskScore)
pathologic_N$riskScore <- ifelse(pathologic_N$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_N$futime, pathologic_N$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_N[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_N$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("pathologic_N-")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_N, plot_title = "pathologic_N-")
library(export)
 graph2ppt(file = "output/Figure/29_pN-_risk_KM.ppt",width=6.3, height=5.5)

####
pathologic_N <- rt[,c(1:3,8)]
pathologic_N <- pathologic_N[!(apply(pathologic_N, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "N0" 替换为 "N-"，"N1/2" 替换为 "N+"
pathologic_N$pathologic_N <- gsub("N0", "N-", pathologic_N$pathologic_N)
pathologic_N$pathologic_N <- gsub("N1", "N+", pathologic_N$pathologic_N)
pathologic_N$pathologic_N <- gsub("N2", "N+", pathologic_N$pathologic_N)
pathologic_N$pathologic_N <- gsub("N3", "N+", pathologic_N$pathologic_N)
##########pathologic_N
pathologic_N <- pathologic_N %>% 
  filter(pathologic_N=="N+") %>% 
  dplyr::select(-pathologic_N)
median_score <- median(pathologic_N$riskScore)
pathologic_N$riskScore <- ifelse(pathologic_N$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_N$futime, pathologic_N$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_N[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_N$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("pathologic_N+")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_N, plot_title = "pathologic_N+")
library(export)
 graph2ppt(file = "output/Figure/30_pN+_risk_KM.ppt",width=6.3, height=5.5)


####
pathologic_M <- rt[,c(1:3,9)]
pathologic_M <- pathologic_M[!(apply(pathologic_M, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "M0" 替换为 "M-"，"M1" 替换为 "M+"
pathologic_M$pathologic_M <- gsub("M0", "M-", pathologic_M$pathologic_M)
pathologic_M$pathologic_M <- gsub("M1", "M+", pathologic_M$pathologic_M)
##########pathologic_M
pathologic_M <- pathologic_M %>% 
  filter(pathologic_M=="M-") %>% 
  dplyr::select(-pathologic_M)
median_score <- median(pathologic_M$riskScore)
pathologic_M$riskScore <- ifelse(pathologic_M$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_M$futime, pathologic_M$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_M[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_M$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("pathologic_M-")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_M, plot_title = "pathologic_M-")
library(export)
 graph2ppt(file = "output/Figure/31_pM-_risk_KM.ppt",width=6.3, height=5.5)

######M+的数据太少了，跑不出来####
# pathologic_M <- rt[,c(1:3,9)]
# pathologic_M <- pathologic_M[!(apply(pathologic_M, 1, function(x) any(x == "unknown"))), ]
# # 将数据框中某一列的 "M0" 替换为 "M-"，"M1" 替换为 "M+"
# pathologic_M$pathologic_M <- gsub("M0", "M-", pathologic_M$pathologic_M)
# pathologic_M$pathologic_M <- gsub("M1", "M+", pathologic_M$pathologic_M)
# ##########pathologic_M
# pathologic_M <- pathologic_M %>%
#   filter(pathologic_M=="M+") %>%
#   dplyr::select(-pathologic_M)
# median_score <- median(pathologic_M$riskScore)
# pathologic_M$riskScore <- ifelse(pathologic_M$riskScore <= median_score, "low", "high")
# ##1.写函数
# genes <- "riskScore"
# your.surv <- Surv(pathologic_M$futime, pathologic_M$fustat)
# # your.km.plot <- function(genes,data){
# #   print(genes)
# #   group <- pathologic_M[,genes] #分组
# #   survival_dat <- data.frame(group = group)
# #   group <- factor(group, levels = c("low", "high"))
# #   fit <- survfit(your.surv ~ group)
# #   sdf <- survdiff(your.surv ~ group,rho=0)
# #   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
# #   p.val
# #   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
# #                         legend.title = genes,#定义图例的名称
# #                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
# #                         #legend = "top",#图例位置
# #                         pval = T, #在图上添加log rank检验的p值
# #                         #pval.method = TRUE,#添加p值的检验方法
# #                         conf.int = F,#添加置信区间
# #                         risk.table = TRUE, #在图下方添加风险表
# #                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
# #                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
# #                         #linetype = "strata", #改变不同组别的生存曲线的线型
# #                         #surv.median.line = "hv", #标注出中位生存时间
# #                         xlab = "Time in years", #x轴标题
# #                         xlim = c(0,max(pathologic_M$futime)+1), #展示x轴的范围
# #                         break.time.by = 2, #x轴间隔
# #                         size = 1.5, #线条大小
# #                         #ggtheme = theme_bw(), #为图形添加网格
# #                         palette = c("#2878b5", "#c82423")#图形颜色风格
# #   ) + ggtitle("pathologic_M+")
# # 
# #   photo2  #看一下图
# #   # 修改图例
# #   # 修改风险表的图例名称
# #   photo2$table <- photo2$table + labs(
# #     title = "Number at risk")
# #   # Changing the font size, style and color of photo2
# #   # survival curves, risk table
# #   photo2  #再看一下图
# # }
# ##2.测试函数功能
# your.km.plot("riskScore",data = pathologic_M, plot_title = "pathologic_M+")
# library(export)
# # graph2ppt(file = "output/Figure/32_pM+_risk_KM.ppt",width=6.3, height=5.5)
# 

####
pathologic_stage <- rt[,c(1:3,10)]
pathologic_stage <- pathologic_stage[!(apply(pathologic_stage, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "I/II" 替换为 "I+II"，"III/IV" 替换为 "III+IV"
pathologic_stage$pathologic_stage <- gsub("Stage IV", "34", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("Stage III", "34", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("Stage II", "12", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("Stage I", "12", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("12", "I+II", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("34", "III+IV", pathologic_stage$pathologic_stage)
##########pathologic_stage
pathologic_stage <- pathologic_stage %>% 
  filter(pathologic_stage=="I+II") %>% 
  dplyr::select(-pathologic_stage)
median_score <- median(pathologic_stage$riskScore)
pathologic_stage$riskScore <- ifelse(pathologic_stage$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_stage$futime, pathologic_stage$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_stage[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_stage$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("pathologic_StageI+II")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_stage, plot_title = "pathologic_StageI+II")
library(export)
 graph2ppt(file = "output/Figure/33_pstage12_risk_KM.ppt",width=6.3, height=5.5)

####
pathologic_stage <- rt[,c(1:3,10)]
pathologic_stage <- pathologic_stage[!(apply(pathologic_stage, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "I/II" 替换为 "I+II"，"III/IV" 替换为 "III+IV"
pathologic_stage$pathologic_stage <- gsub("Stage IV", "34", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("Stage III", "34", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("Stage II", "12", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("Stage I", "12", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("12", "I+II", pathologic_stage$pathologic_stage)
pathologic_stage$pathologic_stage <- gsub("34", "III+IV", pathologic_stage$pathologic_stage)
##########pathologic_stage
pathologic_stage <- pathologic_stage %>% 
  filter(pathologic_stage=="III+IV") %>% 
  dplyr::select(-pathologic_stage)
median_score <- median(pathologic_stage$riskScore)
pathologic_stage$riskScore <- ifelse(pathologic_stage$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(pathologic_stage$futime, pathologic_stage$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- pathologic_stage[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(pathologic_stage$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   )  + ggtitle("pathologic_StageIII+IV")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = pathologic_stage, plot_title = "pathologic_StageIII+IV")
library(export)
graph2ppt(file = "output/Figure/34_pstage34_risk_KM.ppt",width=6.3, height=5.5)



















####
clinical_T <- rt[,c(1:3,11)]
clinical_T <- clinical_T[!(apply(clinical_T, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "T1/T2" 替换为 "T1+T2"，"T3/T4" 替换为 "T3+T4"
clinical_T$clinical_T <- gsub("T1", "I", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("T2", "I", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("T3", "II", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("T4", "II", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("II", "T3+T4", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("I", "T1+T2", clinical_T$clinical_T)
##########clinical_T
clinical_T <- clinical_T %>% 
  filter(clinical_T=="T1+T2") %>% 
  dplyr::select(-clinical_T)
median_score <- median(clinical_T$riskScore)
clinical_T$riskScore <- ifelse(clinical_T$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_T$futime, clinical_T$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_T[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_T$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_T1T2")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_T, plot_title = "clinical_T1T2")
library(export)
 graph2ppt(file = "output/Figure/27_cT1T2_risk_KM.ppt",width=6.3, height=5.5)

####
clinical_T <- rt[,c(1:3,11)]
clinical_T <- clinical_T[!(apply(clinical_T, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "T1/T2" 替换为 "T1+T2"，"T3/T4" 替换为 "T3+T4"
clinical_T$clinical_T <- gsub("T1", "I", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("T2", "I", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("T3", "II", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("T4", "II", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("II", "T3+T4", clinical_T$clinical_T)
clinical_T$clinical_T <- gsub("I", "T1+T2", clinical_T$clinical_T)
##########clinical_T
clinical_T <- clinical_T %>% 
  filter(clinical_T=="T3+T4") %>% 
  dplyr::select(-clinical_T)
median_score <- median(clinical_T$riskScore)
clinical_T$riskScore <- ifelse(clinical_T$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_T$futime, clinical_T$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_T[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_T$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_T3T4")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_T, plot_title = "clinical_T3T4")
library(export)
graph2ppt(file = "output/Figure/28_cT3T4_risk_KM.ppt",width=6.3, height=5.5)


####
clinical_N <- rt[,c(1:3,12)]
clinical_N <- clinical_N[!(apply(clinical_N, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "N0" 替换为 "N-"，"N1/2" 替换为 "N+"
clinical_N$clinical_N <- gsub("N0", "N-", clinical_N$clinical_N)
clinical_N$clinical_N <- gsub("N1", "N+", clinical_N$clinical_N)
clinical_N$clinical_N <- gsub("N2", "N+", clinical_N$clinical_N)
clinical_N$clinical_N <- gsub("N3", "N+", clinical_N$clinical_N)
##########clinical_N
clinical_N <- clinical_N %>% 
  filter(clinical_N=="N-") %>% 
  dplyr::select(-clinical_N)
median_score <- median(clinical_N$riskScore)
clinical_N$riskScore <- ifelse(clinical_N$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_N$futime, clinical_N$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_N[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_N$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_N-")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_N, plot_title = "clinical_N-")
library(export)
graph2ppt(file = "output/Figure/29_cN-_risk_KM.ppt",width=6.3, height=5.5)

####
clinical_N <- rt[,c(1:3,12)]
clinical_N <- clinical_N[!(apply(clinical_N, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "N0" 替换为 "N-"，"N1/2" 替换为 "N+"
clinical_N$clinical_N <- gsub("N0", "N-", clinical_N$clinical_N)
clinical_N$clinical_N <- gsub("N1", "N+", clinical_N$clinical_N)
clinical_N$clinical_N <- gsub("N2", "N+", clinical_N$clinical_N)
clinical_N$clinical_N <- gsub("N3", "N+", clinical_N$clinical_N)
##########clinical_N
clinical_N <- clinical_N %>% 
  filter(clinical_N=="N+") %>% 
  dplyr::select(-clinical_N)
median_score <- median(clinical_N$riskScore)
clinical_N$riskScore <- ifelse(clinical_N$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_N$futime, clinical_N$fustat) 
# # your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_N[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_N$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_N+")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_N, plot_title = "clinical_N+")
library(export)
graph2ppt(file = "output/Figure/30_cN+_risk_KM.ppt",width=6.3, height=5.5)


####
clinical_M <- rt[,c(1:3,13)]
clinical_M <- clinical_M[!(apply(clinical_M, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "M0" 替换为 "M-"，"M1" 替换为 "M+"
clinical_M$clinical_M <- gsub("M0", "M-", clinical_M$clinical_M)
clinical_M$clinical_M <- gsub("M1", "M+", clinical_M$clinical_M)
##########clinical_M
clinical_M <- clinical_M %>% 
  filter(clinical_M=="M-") %>% 
  dplyr::select(-clinical_M)
median_score <- median(clinical_M$riskScore)
clinical_M$riskScore <- ifelse(clinical_M$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_M$futime, clinical_M$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_M[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_M$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_M-")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_M, plot_title = "clinical_M-")
library(export)
graph2ppt(file = "output/Figure/31_cM-_risk_KM.ppt",width=6.3, height=5.5)

####
clinical_M <- rt[,c(1:3,13)]
clinical_M <- clinical_M[!(apply(clinical_M, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "M0" 替换为 "M-"，"M1" 替换为 "M+"
clinical_M$clinical_M <- gsub("M0", "M-", clinical_M$clinical_M)
clinical_M$clinical_M <- gsub("M1", "M+", clinical_M$clinical_M)
##########clinical_M
clinical_M <- clinical_M %>% 
  filter(clinical_M=="M+") %>% 
  dplyr::select(-clinical_M)
median_score <- median(clinical_M$riskScore)
clinical_M$riskScore <- ifelse(clinical_M$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_M$futime, clinical_M$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_M[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_M$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_M+")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_M, plot_title = "clinical_M+")
library(export)
graph2ppt(file = "output/Figure/32_cM+_risk_KM.ppt",width=6.3, height=5.5)


####
clinical_stage <- rt[,c(1:3,14)]
clinical_stage <- clinical_stage[!(apply(clinical_stage, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "I/II" 替换为 "I+II"，"III/IV" 替换为 "III+IV"
clinical_stage$clinical_stage <- gsub("Stage IV", "34", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("Stage III", "34", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("Stage II", "12", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("Stage I", "12", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("12", "I+II", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("34", "III+IV", clinical_stage$clinical_stage)
##########clinical_stage
clinical_stage <- clinical_stage %>% 
  filter(clinical_stage=="I+II") %>% 
  dplyr::select(-clinical_stage)
median_score <- median(clinical_stage$riskScore)
clinical_stage$riskScore <- ifelse(clinical_stage$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_stage$futime, clinical_stage$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_stage[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_stage$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) + ggtitle("clinical_StageI+II")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_stage, plot_title = "clinical_StageI+II")
library(export)
graph2ppt(file = "output/Figure/33_cstage12_risk_KM.ppt",width=6.3, height=5.5)

####
clinical_stage <- rt[,c(1:3,14)]
clinical_stage <- clinical_stage[!(apply(clinical_stage, 1, function(x) any(x == "unknown"))), ]
# 将数据框中某一列的 "I/II" 替换为 "I+II"，"III/IV" 替换为 "III+IV"
clinical_stage$clinical_stage <- gsub("Stage IV", "34", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("Stage III", "34", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("Stage II", "12", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("Stage I", "12", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("12", "I+II", clinical_stage$clinical_stage)
clinical_stage$clinical_stage <- gsub("34", "III+IV", clinical_stage$clinical_stage)
##########clinical_stage
clinical_stage <- clinical_stage %>% 
  filter(clinical_stage=="III+IV") %>% 
  dplyr::select(-clinical_stage)
median_score <- median(clinical_stage$riskScore)
clinical_stage$riskScore <- ifelse(clinical_stage$riskScore <= median_score, "low", "high")
##1.写函数
genes <- "riskScore"
your.surv <- Surv(clinical_stage$futime, clinical_stage$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- clinical_stage[,genes] #分组
#   survival_dat <- data.frame(group = group)
#   group <- factor(group, levels = c("low", "high")) 
#   fit <- survfit(your.surv ~ group)
#   sdf <- survdiff(your.surv ~ group,rho=0)
#   p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
#   p.val
#   photo2 <-  ggsurvplot(fit,data = survival_dat, #这里很关键，不然会报错
#                         legend.title = genes,#定义图例的名称
#                         legend.labs = c("low","high"), #所以上面要因子化分组顺序
#                         #legend = "top",#图例位置
#                         pval = T, #在图上添加log rank检验的p值
#                         #pval.method = TRUE,#添加p值的检验方法
#                         conf.int = F,#添加置信区间
#                         risk.table = TRUE, #在图下方添加风险表
#                         #risk.table.col = "strata", #根据数据分组为风险表添加颜色
#                         risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
#                         #linetype = "strata", #改变不同组别的生存曲线的线型
#                         #surv.median.line = "hv", #标注出中位生存时间
#                         xlab = "Time in years", #x轴标题
#                         xlim = c(0,max(clinical_stage$futime)+1), #展示x轴的范围
#                         break.time.by = 2, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   )  + ggtitle("clinical_StageIII+IV")
#   
#   photo2  #看一下图
#   # 修改图例
#   # 修改风险表的图例名称 
#   photo2$table <- photo2$table + labs(
#     title = "Number at risk")
#   # Changing the font size, style and color of photo2
#   # survival curves, risk table
#   photo2  #再看一下图
# }
##2.测试函数功能
your.km.plot("riskScore",data = clinical_stage, plot_title = "clinical_StageIII+IV")
library(export)
graph2ppt(file = "output/Figure/34_cstage34_risk_KM.ppt",width=6.3, height=5.5)

