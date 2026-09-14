############################################绘制entire_KM生存图14_out2_out3_out4_out5
############################################绘制entire_KM生存图14_out2_out3_out4_out5
############################################绘制entire_KM生存图14_out2_out3_out4_out5
##清空数据
rm(list = ls())
library(tidyr)
library(dplyr)
library(tibble)
library(survminer)
library(survival)
library(ggplot2)
library(ggprism)
library(tidyverse)
library(timeROC)
library(biomaRt)
library(openxlsx)


#绘制test生存曲线
rt=read.table("output/riskoutside5.txt",header=T,sep="\t",check.names=F,row.names=1)

##1.写函数
genes <- "risk"
your.surv <- Surv(rt$futime, rt$fustat) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- rt[,genes] #分组
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
                        conf.int = F,#添加置信区间
                        risk.table = TRUE, #在图下方添加风险表
                        #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                        risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                        #linetype = "strata", #改变不同组别的生存曲线的线型
                        #surv.median.line = "hv", #标注出中位生存时间
                        xlab = "Time in years", #x轴标题
                        xlim = c(0,max(rt$futime)+1), #展示x轴的范围
                        break.time.by = 2, #x轴间隔
                        size = 1.5, #线条大小
                        #ggtheme = theme_bw(), #为图形添加网格
                        palette = c("#2878b5", "#c82423")#图形颜色风格
  ) 
  
  photo2  #看一下图
  # 修改图例
  # 修改风险表的图例名称 
  photo2$table <- photo2$table + labs(
    title = "Number at risk")
  # Changing the font size, style and color of photo2
  # survival curves, risk table
  photo2  #再看一下图
}
##2.测试函数功能
your.km.plot("risk",data = rt)

library(export)
graph2ppt(file = "output/Figure/999_im1_km.ppt",width=5.4, height=6.4)

load(file = "resource/IMvigor210CoreBiologies.Rdata")
metadata <- phenoData %>% 
  dplyr::select(binaryResponse) 

rt <- rt %>% 
  rownames_to_column("ID")
metadata <- metadata %>% 
  rownames_to_column("ID")
data <- merge(metadata,rt,by="ID")

data <- na.omit(data)
colnames(data)[2] <- "response"
#############################################################疗效差异表达
response <- data[,c("response","riskScore")]
response <- response %>% 
  dplyr::select(riskScore,response)
response <- na.omit(response)
colnames(response) <- c("riskScore","clinical")
##因子化，以便排序
response$clinical <- factor(response$clinical, levels = c("CR/PR", "SD/PD"))
##比较P值，选取比较对象
comparisons <- list(c("CR/PR", "SD/PD"))

e <- ggplot(response, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 4) +
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

graph2ppt(file = "output/Figure/999_im1_response1_compare.ppt",width=4.85, height=5.8)


df <- data %>% 
  dplyr::select(response,risk)

df <- data %>% 
  dplyr::select(response,risk)

# 计算每组的频数
df_summary <- df %>%
  group_by(risk, response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(risk) %>%
  mutate(percentage = count / sum(count) * 100)

# 构建列联表（contingency table）
contingency_table <- table(df$risk, df$response)

# 计算期望频数
expected_counts <- chisq.test(contingency_table)$expected

# 选择合适的检验方法
if (all(expected_counts >= 5)) {
  chi_test <- chisq.test(contingency_table)  # 期望频数都 ≥5，使用卡方检验
} else {
  chi_test <- fisher.test(contingency_table)  # 期望频数 <5，使用 Fisher 精确检验
}

# 提取 P 值
p_value <- chi_test$p.value  

# 确定显著性标记
sig_label <- ifelse(p_value < 0.0001, "****",  # p < 0.0001 -> ****
                    ifelse(p_value < 0.001, "***",  # p < 0.001 -> ***
                           ifelse(p_value < 0.01, "**",  # p < 0.01 -> **
                                  ifelse(p_value < 0.05, "*", "n.s."))))  # p < 0.05 -> *，否则 n.s.

# 绘制堆叠柱状图
p <- ggplot(df_summary, aes(x = risk, y = percentage, fill = response)) +
  geom_bar(stat = "identity", width = 0.6) +  # 堆叠柱状图
  geom_text(aes(label = paste0(round(percentage, 2), "%")), 
            position = position_stack(vjust = 0.5), size = 7) +  # 添加百分比标签
  scale_fill_manual(values = c("#008102", "#a00003")) +  # 指定颜色
  labs(x = "", y = "Percentage", fill = "") +  # 去掉坐标轴标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),  # 调整 X 轴字体
    legend.position = "top"  # 将图例放在顶部
  ) +
  geom_signif(
    comparisons = list(c("high", "low")),  # 比较 high risk 和 low risk
    annotations = sig_label,  # 显著性标记
    y_position = 110,  # 星号位置
    tip_length = 0.1,  # 线条长度
    size = 1.2  # 线条粗细
  ) + scale_y_continuous(breaks = seq(0, 100, by = 20)) +  # 设置 Y 轴间隔
  theme(axis.text.y = element_text(size = 12, face = "bold"))  # 调整 Y 轴字体

print(p)
graph2ppt(file = "output/Figure/999_im1_response2_compare.ppt",width=4.85, height=5.8)


#######R-LOOP里面的6没跑出来！！！####
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6
# ##清空数据
# rm(list = ls())
# 
# #绘制test生存曲线
# rt=read.table("output/riskoutside6.txt",header=T,sep="\t",check.names=F,row.names=1)
# 
# ##1.写函数
# genes <- "risk"
# your.surv <- Surv(rt$futime, rt$fustat) 
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- rt[,genes] #分组
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
#                         xlim = c(0,max(rt$futime)+1), #展示x轴的范围
#                         break.time.by = 1, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   ) 
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
# ##2.测试函数功能
# your.km.plot("risk",data = rt)
# 
# library(export)
# # graph2ppt(file = "output/Figure/999_im2_km.ppt",width=5.4, height=6.4)
# 
# library(openxlsx)
# cli <- read.xlsx("resource/GSE78220_cli.xlsx")
# metadata <- cli %>% 
#   dplyr::select(Patient.ID,irRECIST) 
# colnames(metadata)[1] <- "id"
# colnames(metadata)[2] <- "response"
# 
# 
# rt <- rt %>% 
#   rownames_to_column("ID")
# 
# ID <- read.xlsx("resource/GSE78220_ID.xlsx")
# 
# rt <- merge(ID,rt,by="ID")
# 
# data <- merge(metadata,rt,by="id")
# data <- na.omit(data)
# 
# data <- data %>%
#   mutate(response = case_when(
#     response %in% c("Complete Response", "Partial Response") ~ "CR/PR",
#     response == "Progressive Disease" ~ "SD/PD",
#     TRUE ~ response  # 其他值保持不变
#   ))
# 
# #############################################################疗效差异表达
# response <- data[,c("response","riskScore")]
# response <- response %>% 
#   dplyr::select(riskScore,response)
# response <- na.omit(response)
# colnames(response) <- c("riskScore","clinical")
# ##因子化，以便排序
# response$clinical <- factor(response$clinical, levels = c("CR/PR", "SD/PD"))
# ##比较P值，选取比较对象
# comparisons <- list(c("CR/PR", "SD/PD"))
# 
# e <- ggplot(response, aes(x = clinical, y = riskScore))
# e + geom_violin(aes(fill = clinical), trim = FALSE) + 
#   geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 4) +
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
# 
# # graph2ppt(file = "output/Figure/999_im2_response1_compare.ppt",width=4.85, height=5.8)
# 
# 
# df <- data %>% 
#   dplyr::select(response,risk)
# 
# # 计算每组的频数
# df_summary <- df %>%
#   group_by(risk, response) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   group_by(risk) %>%
#   mutate(percentage = count / sum(count) * 100)
# 
# # 构建列联表（contingency table）
# contingency_table <- table(df$risk, df$response)
# 
# # 计算期望频数
# expected_counts <- chisq.test(contingency_table)$expected
# 
# # 选择合适的检验方法
# if (all(expected_counts >= 5)) {
#   chi_test <- chisq.test(contingency_table)  # 期望频数都 ≥5，使用卡方检验
# } else {
#   chi_test <- fisher.test(contingency_table)  # 期望频数 <5，使用 Fisher 精确检验
# }
# 
# # 提取 P 值
# p_value <- chi_test$p.value  
# 
# # 确定显著性标记
# sig_label <- ifelse(p_value < 0.0001, "****",  # p < 0.0001 -> ****
#                     ifelse(p_value < 0.001, "***",  # p < 0.001 -> ***
#                            ifelse(p_value < 0.01, "**",  # p < 0.01 -> **
#                                   ifelse(p_value < 0.05, "*", "n.s."))))  # p < 0.05 -> *，否则 n.s.
# 
# # 绘制堆叠柱状图
# p <- ggplot(df_summary, aes(x = risk, y = percentage, fill = response)) +
#   geom_bar(stat = "identity", width = 0.6) +  # 堆叠柱状图
#   geom_text(aes(label = paste0(round(percentage, 2), "%")), 
#             position = position_stack(vjust = 0.5), size = 7) +  # 添加百分比标签
#   scale_fill_manual(values = c("#008102", "#a00003")) +  # 指定颜色
#   labs(x = "", y = "Percentage", fill = "") +  # 去掉坐标轴标题
#   theme_minimal() +  # 使用简洁主题
#   theme(
#     axis.text.x = element_text(size = 15, face = "bold"),  # 调整 X 轴字体
#     legend.position = "top"  # 将图例放在顶部
#   ) +
#   geom_signif(
#     comparisons = list(c("high", "low")),  # 比较 high risk 和 low risk
#     annotations = sig_label,  # 显著性标记
#     y_position = 110,  # 星号位置
#     tip_length = 0.1,  # 线条长度
#     size = 1.2  # 线条粗细
#   ) + scale_y_continuous(breaks = seq(0, 100, by = 20)) +  # 设置 Y 轴间隔
#   theme(axis.text.y = element_text(size = 12, face = "bold"))  # 调整 Y 轴字体
# 
# print(p)
# # graph2ppt(file = "output/Figure/999_im2_response2_compare.ppt",width=4.85, height=5.8)
# 


############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6_out7
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6_out7
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6_out7
##清空数据
rm(list = ls())

#绘制test生存曲线
rt=read.table("output/riskoutside7.txt",header=T,sep="\t",check.names=F,row.names=1)

##1.写函数
genes <- "risk"
your.surv <- Surv(rt$futime, rt$fustat) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- rt[,genes] #分组
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
                        conf.int = F,#添加置信区间
                        risk.table = TRUE, #在图下方添加风险表
                        #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                        risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                        #linetype = "strata", #改变不同组别的生存曲线的线型
                        #surv.median.line = "hv", #标注出中位生存时间
                        xlab = "Time in years", #x轴标题
                        xlim = c(0,max(rt$futime)+1), #展示x轴的范围
                        break.time.by = 1, #x轴间隔
                        size = 1.5, #线条大小
                        #ggtheme = theme_bw(), #为图形添加网格
                        palette = c("#2878b5", "#c82423")#图形颜色风格
  ) 
  
  photo2  #看一下图
  # 修改图例
  # 修改风险表的图例名称 
  photo2$table <- photo2$table + labs(
    title = "Number at risk")
  # Changing the font size, style and color of photo2
  # survival curves, risk table
  photo2  #再看一下图
}
##2.测试函数功能
your.km.plot("risk",data = rt)

rlibrary(export)
graph2ppt(file = "output/Figure/999_im3_km.ppt",width=5.4, height=6.4)



cli <- read.xlsx("resource/GSE91061_cli.xlsx")
metadata <- cli %>% 
  dplyr::select(Patient,Response) 
colnames(metadata)[1] <- "ID"
colnames(metadata)[2] <- "response"

rt <- rt %>% 
  rownames_to_column("ID")
data <- merge(metadata,rt,by="ID")
data <- na.omit(data)
data <- data[-10,]

data <- data %>%
  mutate(response = case_when(
    response %in% c("CR", "PR") ~ "CR/PR",
    response %in% c("SD", "PD") ~ "SD/PD",
    TRUE ~ response  # 其他值保持不变
  ))

#############################################################疗效差异表达
response <- data[,c("response","riskScore")]
response <- response %>% 
  dplyr::select(riskScore,response)
response <- na.omit(response)
colnames(response) <- c("riskScore","clinical")
##因子化，以便排序
response$clinical <- factor(response$clinical, levels = c("CR/PR", "SD/PD"))
##比较P值，选取比较对象
comparisons <- list(c("CR/PR", "SD/PD"))

e <- ggplot(response, aes(x = clinical, y = riskScore))
e + geom_violin(aes(fill = clinical), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 4) +
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

graph2ppt(file = "output/Figure/999_im3_response1_compare.ppt",width=4.85, height=5.8)


df <- data %>% 
  dplyr::select(response,risk)

# 计算每组的频数
df_summary <- df %>%
  group_by(risk, response) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(risk) %>%
  mutate(percentage = count / sum(count) * 100)

# 构建列联表（contingency table）
contingency_table <- table(df$risk, df$response)

# 计算期望频数
expected_counts <- chisq.test(contingency_table)$expected

# 选择合适的检验方法
if (all(expected_counts >= 5)) {
  chi_test <- chisq.test(contingency_table)  # 期望频数都 ≥5，使用卡方检验
} else {
  chi_test <- fisher.test(contingency_table)  # 期望频数 <5，使用 Fisher 精确检验
}

# 提取 P 值
p_value <- chi_test$p.value  

# 确定显著性标记
sig_label <- ifelse(p_value < 0.0001, "****",  # p < 0.0001 -> ****
                    ifelse(p_value < 0.001, "***",  # p < 0.001 -> ***
                           ifelse(p_value < 0.01, "**",  # p < 0.01 -> **
                                  ifelse(p_value < 0.05, "*", "n.s."))))  # p < 0.05 -> *，否则 n.s.

# 绘制堆叠柱状图
p <- ggplot(df_summary, aes(x = risk, y = percentage, fill = response)) +
  geom_bar(stat = "identity", width = 0.6) +  # 堆叠柱状图
  geom_text(aes(label = paste0(round(percentage, 2), "%")), 
            position = position_stack(vjust = 0.5), size = 7) +  # 添加百分比标签
  scale_fill_manual(values = c("#008102", "#a00003")) +  # 指定颜色
  labs(x = "", y = "Percentage", fill = "") +  # 去掉坐标轴标题
  theme_minimal() +  # 使用简洁主题
  theme(
    axis.text.x = element_text(size = 15, face = "bold"),  # 调整 X 轴字体
    legend.position = "top"  # 将图例放在顶部
  ) +
  geom_signif(
    comparisons = list(c("high", "low")),  # 比较 high risk 和 low risk
    annotations = sig_label,  # 显著性标记
    y_position = 110,  # 星号位置
    tip_length = 0.1,  # 线条长度
    size = 1.2  # 线条粗细
  ) + scale_y_continuous(breaks = seq(0, 100, by = 20)) +  # 设置 Y 轴间隔
  theme(axis.text.y = element_text(size = 12, face = "bold"))  # 调整 Y 轴字体

print(p)
graph2ppt(file = "output/Figure/999_im3_response2_compare.ppt",width=4.85, height=5.8)

############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6_out7_out8
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6_out7_out8
############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6_out7_out8
# ##清空数据
# rm(list = ls())
# 
# #绘制test生存曲线
# rt=read.table("output/riskoutside8.txt",header=T,sep="\t",check.names=F,row.names=1)
# 
# ##1.写函数
# genes <- "risk"
# your.surv <- Surv(rt$futime, rt$fustat)
# your.km.plot <- function(genes,data){
#   print(genes)
#   group <- rt[,genes] #分组
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
#                         xlim = c(0,max(rt$futime)+1), #展示x轴的范围
#                         break.time.by = 1, #x轴间隔
#                         size = 1.5, #线条大小
#                         #ggtheme = theme_bw(), #为图形添加网格
#                         palette = c("#2878b5", "#c82423")#图形颜色风格
#   )
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
# ##2.测试函数功能
# your.km.plot("risk",data = rt)
# 
# library(export)
# # graph2ppt(file = "output/Figure/14_entire_km_out8.ppt",width=6.3, height=5.5)
# 
# 
# clinical <- data.frame(
#   Sample_ID = c(3, 20, 258, 325, 378, 488, 541, 573, 618, 658, 678, 700, 720, 756, 825, 830, 947, 990, 1017, 1066, 1079, 1104, 1145, 1155, 1164, 1203, 1208, 1250, 1297, 1322, 1327, 1337, 1352,
#                 1358, 1401, 1412, 1425, 1443, 1456, 1490, 1508, 1510, 1528, 1554, 1589, 1619, 1637, 1708, 1711, 1751, 1778, 1809, 1873, 1883, 1960, 2107, 2126, 2133, 2133, 2317),
#   Mutation_burden = c(176, 93, 53, 217, 117, 353, 180, 366, 141, 117, 489, 68, 325, 603, 194, 223, 401, 43, 82, 42, 101, 1100, 178, 273, 175, 137, 242, 382, 427, 182, 441, 428, 1391,
#                       328, 42, 202, 260, 314, 157, 169, 18, 273, 146, 38, 862, 13, 263, 266, 308, 327, 343, 156, 388, 245, 3424, 569, 285, 211, 524, 576),
#   Global_methylation_level = c(0.462218282, 0.497886714, 0.414779054, 0.357762176, 0.414296576, 0.434316952, 0.362546225, 0.410597593, 0.349218209, 0.403780292, 0.396234214, 0.428694023, 0.294703592, 0.280030877, 0.251712972, 0.479584284, 0.387089023, 0.325908608, 0.417538488, 0.395929718, 0.392672276, 0.293149031, 0.421131362, 0.335438797, 0.293888515, 0.35971363, 0.34974738, 0.327435945, 0.416016984, 0.443427433, 0.42678156, 0.329661704, 0.22054474,
#                                0.377759517, 0.436672231, 0.264443355, 0.303500801, 0.270158864, 0.407933594, 0.315326936, 0.415487136, 0.377090295, 0.364790201, 0.400425652, 0.328377787, 0.455254226, 0.394335199, 0.450911862, 0.421240392, 0.390149399, 0.371084292, 0.324772186, 0.307758327, 0.413997537, 0.413389522, 0.312709532, 0.415388193, 0.377043232, 0.381138966, 0.323139645),
#   Aneuploidy_level = c(9963319.531, 0, 144952144.3, 373215792.9, 102187138.2, 535751134.6, 349091677.2, 114018200.9, 263496773.3, 141510150.7, 191366702.2, 17433981.68, 307380376.1, 658988768.2, 784914239.2, 0, 587994366.69, 285011126.5, 286678444.4, 48979857.76, 12815975.36, 895819215, 5860926.472, 443253779.2, 520968472.1, 437223965.6, 927702954.8, 861780460, 68793176.27, 227732642, 99666614.32, 436123975.5, 553574112.3,
#                        649708156.1, 17494258.78, 545902993.4, 346192818.8, 836459992.2, 76492185.01, 403292217.7, 7542228.411, 465391056.4, 374128096.1, 5947504.686, 297587719.9, 7041421.93, 265842882.3, 2076055.089, 401532642.4, 225455169.6, 128523671.2, 400743320.7, 486476347.6, 199412472, 349177478.5, 716225322, 145048834.8, 235929294.4, 93298957.65, 495812412.5),
#   PFS = c(28.1, 1.133333333, 0.966666667, 3.866666667, 13.733333333, 1.6, 1.2, 1.5, 1.366666667, 1.066666667, 22.433333333, 1.366666667, 4.066666667, 1.5, 1.433333333, 2.933333333, 20.6, 3.166666667, 1.266666667, 0.766666667, 5.666666667, 0.966666667, 2.433333333, 1.133333333, 0.366666667, 0.833333333, 0.933333333, 0.933333333, 10.5, 1.1, 6.833333333, 2.166666667, 8.333333333,
#           8.566666667, 0.933333333, 0.1, 1.966666667, 1.766666667, 7.5, 0.9, 2.166666667, 10.8, 9.3, 1.233333333, 3.966666667, 2.733333333, 2.6, 5.6, "na", 4.833333333, 0.766666667, 1.466666667, 1.233333333, 4.366666667, 1.166666667, 6.466666667, 2.1, 1.233333333, 3.766666667, 0.766666667),
#   PD_Event_1_Censoring_0 = c(1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
#                              1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, "na", 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1),
#   Clinical_benefit = c("DCB", "NDB", "NDB", "na", "DCB", "NDB", "NDB", "NDB", "NDB", "NDB", "DCB", "NDB", "NDB", "NDB", "NDB", "NDB", "DCB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "DCB", "NDB", "DCB",
#                        "DCB", "NDB", "NDB", "NDB", "NDB", "DCB", "NDB", "NDB", "DCB", "DCB", "NDB", "NDB", "NDB", "NDB", "DCB", "na", "NDB", "NDB", "NDB", "NDB", "NDB", "NDB", "DCB", "NDB", "NDB", "NDB", "NDB")
# )
# 
# clinical <- clinical[-49,]
# clinical <- clinical[-4,]
# clinical <- clinical %>%
#   dplyr::select(Sample_ID,Clinical_benefit)
# library(GEOquery)
# geo_data0 <- getGEO(filename = "resource\\GSE135222_series_matrix.txt.gz", getGPL = F)
# cli <- pData(geo_data0)
# cli <- cli %>%
#   dplyr::select(title) %>%
#   rownames_to_column("ID")
# 
# cli <- cli %>%
#   mutate(Sample_ID = as.numeric(gsub("^NSCLC (.*)", "\\1", title)))
# 
# data <- merge(clinical,cli,by="Sample_ID")
# rt <- rt %>%
#   rownames_to_column("ID")
# data2 <- merge(data,rt,by="ID")
# 
# data2 <- data2 %>%
#   mutate(Clinical_benefit = case_when(
#     Clinical_benefit %in% c("NDB") ~ "non-responder",
#     Clinical_benefit %in% c("DCB") ~ "responder",
#     TRUE ~ Clinical_benefit  # 其他值保持不变
#   ))
# 
# colnames(data2)[3] <- "response"
# #############################################################疗效差异表达
# response <- data2[,c("response","riskScore")]
# response <- response %>%
#   dplyr::select(riskScore,response)
# response <- na.omit(response)
# colnames(response) <- c("riskScore","clinical")
# ##因子化，以便排序
# response$clinical <- factor(response$clinical, levels = c("responder", "non-responder"))
# ##比较P值，选取比较对象
# comparisons <- list(c("responder", "non-responder"))
# 
# e <- ggplot(response, aes(x = clinical, y = riskScore))
# e + geom_violin(aes(fill = clinical), trim = FALSE) +
#   geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
#   geom_jitter(aes(color = clinical), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 2.5) +
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
# 
# # graph2ppt(file = "output/Figure/99_im4_response1_compare.ppt",width=4.85, height=5.8)
# 
# 
# df <- data2 %>%
#   dplyr::select(response,risk)
# 
# # 计算每组的频数
# df_summary <- df %>%
#   group_by(risk, response) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   group_by(risk) %>%
#   mutate(percentage = count / sum(count) * 100)
# 
# df_summary$response <- factor(df_summary$response, levels = c("responder", "non-responder"))
# 
# # 构建列联表（contingency table）
# contingency_table <- table(df$risk, df$response)
# 
# # 计算期望频数
# expected_counts <- chisq.test(contingency_table)$expected
# 
# # 选择合适的检验方法
# if (all(expected_counts >= 5)) {
#   chi_test <- chisq.test(contingency_table)  # 期望频数都 ≥5，使用卡方检验
# } else {
#   chi_test <- fisher.test(contingency_table)  # 期望频数 <5，使用 Fisher 精确检验
# }
# 
# # 提取 P 值
# p_value <- chi_test$p.value
# 
# # 确定显著性标记
# sig_label <- ifelse(p_value < 0.0001, "****",  # p < 0.0001 -> ****
#                     ifelse(p_value < 0.001, "***",  # p < 0.001 -> ***
#                            ifelse(p_value < 0.01, "**",  # p < 0.01 -> **
#                                   ifelse(p_value < 0.05, "*", "n.s."))))  # p < 0.05 -> *，否则 n.s.
# 
# # 绘制堆叠柱状图
# p <- ggplot(df_summary, aes(x = risk, y = percentage, fill = response)) +
#   geom_bar(stat = "identity", width = 0.6) +  # 堆叠柱状图
#   geom_text(aes(label = paste0(round(percentage, 2), "%")),
#             position = position_stack(vjust = 0.5), size = 7) +  # 添加百分比标签
#   scale_fill_manual(values = c("#008102", "#a00003")) +  # 指定颜色
#   labs(x = "", y = "Percentage", fill = "") +  # 去掉坐标轴标题
#   theme_minimal() +  # 使用简洁主题
#   theme(
#     axis.text.x = element_text(size = 15, face = "bold"),  # 调整 X 轴字体
#     legend.position = "top"  # 将图例放在顶部
#   ) +
#   geom_signif(
#     comparisons = list(c("high", "low")),  # 比较 high risk 和 low risk
#     annotations = sig_label,  # 显著性标记
#     y_position = 110,  # 星号位置
#     tip_length = 0.1,  # 线条长度
#     size = 1.2  # 线条粗细
#   ) + scale_y_continuous(breaks = seq(0, 100, by = 20)) +  # 设置 Y 轴间隔
#   theme(axis.text.y = element_text(size = 12, face = "bold"))  # 调整 Y 轴字体
# 
# print(p)
# # graph2ppt(file = "output/Figure/99_im4_response2_compare.ppt",width=4.7, height=4.7)