# ##加载R包
# rm(list = ls())
# 
# library(tidyr)
# library(dplyr)
# library(tibble)
# library(DESeq2)
# library(tibble)
# 
# HNSC <- readRDS("resource/TCGA.rds/TCGA-HNSC.rds")
# 
# load(file = "resource/gene_id_tcga_v36.Rdata")
# 
# HNSC <- as.data.frame(HNSC)
# HNSC <- HNSC %>%
#   rownames_to_column("ID")
# HNSC <- HNSC %>%
#   column_to_rownames("ID") %>%
#   t() %>%
#   as.data.frame() %>%
#   mutate(newcolumn = rowMeans(.)) %>%
#   arrange(desc(newcolumn))
# 
# HNSC$sample <- substr(rownames(HNSC),14,15)
# HNSC$ID <- substr(rownames(HNSC),1,16)
# HNSC <- HNSC %>%
#   rownames_to_column("TCGA_ID")
# 
# HNSC <- HNSC %>%
#   filter(sample %in% c("01", "11")) %>%
#   dplyr::select(-TCGA_ID) %>%
#   dplyr::select(-newcolumn) %>%
#   dplyr::select(ID,sample,everything()) %>%
#   distinct(ID,.keep_all = T)
# 
# HNSC <- HNSC %>%
#   mutate(sample = ifelse(sample == "01", "Tumor",
#                          ifelse(sample == "11", "Normal", sample)))
# 
# metadata <- HNSC %>%
#   dplyr::select(ID,sample) %>%
#   distinct(ID,.keep_all = T)
# 
# colnames(metadata) <- c("sample","group")
# 
# 
# HNSC <- HNSC %>%
#   dplyr::select(-sample) %>%
#   column_to_rownames("ID") %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column("gene_id")
# 
# expression <- merge(gene_id,HNSC,by="gene_id")
# expression <- expression %>%
#   filter(gene_type=="protein_coding") %>%
#   dplyr::select(-gene_type) %>%
#   dplyr::select(-gene_id) %>%
#   mutate(newcolumn = rowMeans(.[,-1])) %>%
#   arrange(desc(newcolumn)) %>%
#   distinct(gene_name,.keep_all = T) %>%
#   dplyr::select(-newcolumn) %>%
#   column_to_rownames("gene_name")
# 
# ### 3.核心环节，构建dds对象，前面的操作都是铺垫
# ### 要记住四个参数
# ### 要有数据countData，这里就是exprSet
# ### 要有分组信息，在colData中，这里是metadata
# ### design部分是分组信息，格式是~group,没有分组设置为~1
# ### 第一列如果是基因名称，需要自动处理，设置参数tidy=TRUE
# ### 对象是个复合体
# library(DESeq2)
# 
# dds <- DESeqDataSetFromMatrix(
#   countData = expression,
#   colData = metadata,
#   design=~group)
# 
# nrow(dds)
# rownames(dds)
# 
# ## 筛选样本，counts函数提取表达量数据，取出count为0的,我认为一半表达不为0的更好
# dds <- dds[rowMedians(counts(dds))>0,]
# nrow(dds)
# 
# ### 最终VST标准化
# system.time(vsd <- vst(dds, blind = FALSE))
# exprSet_vst <- as.data.frame(assay(vsd))
# 
# colnames(metadata) <- c("TCGA_id","group")
# 
# 
# index <- as.character(metadata$TCGA_id)
# 
# exprSet_vst <- exprSet_vst[,index]
# 
# identical(colnames(exprSet_vst),metadata$TCGA_id)#判断是否完全一致，包括顺序
# 
# 
# data <- as.data.frame(t(exprSet_vst))
# data <- data %>%
#   rownames_to_column("TCGA_id") %>%
#   inner_join(metadata,by = "TCGA_id") %>%
#   dplyr::select(TCGA_id,group,everything())
# 
# 
# ####使用lapply + function模式来是计算，lapply 以及do.call 的介绍，lapply 三步走：
# ##第1,写出单次处理的function
# my.wilcox = function(x){
#   dd <- wilcox.test(data[,x] ~ group, data = data)
#   Normalposition <- grep("Normal",data$group)
#   Tumorposition <- grep("Tumor",data$group)
#   Normalcol <- data[,x][Normalposition]
#   Tumorcol<- data[,x][Tumorposition]
#   NormalMean = mean(Normalcol)
#   TumorMean = mean(Tumorcol)
#   logFC=TumorMean-NormalMean
#   data.frame(gene=x,
#              p.value=dd$p.value,
#              logFC
#   )
# }
# 
# 
# ##测试函数功能
# my.wilcox("RPGR")
# # ##lapply批量作用于函数，返回list，为了方便，我们先用100基因试试
# # lapplylist = lapply(colnames(data)[-c(1:2)][1:100],my.wilcox)
# 
# ##lapply批量作用于函数，返回list
# lapplylist = lapply(colnames(data)[-c(1:2)],my.wilcox)
# 
# ##第3步do.call 转换list
# wilcox_data <- do.call(rbind,lapplylist)
# 
# # library(dplyr)
# # wilcox_210_p <- wilcox_data %>%
# #   filter(p.value < 0.05) %>%
# #   filter(abs(logFC) > 0.585)
# 
# save(wilcox_data,file = "output/wilcox_data_rloop.R")
# 
# save(exprSet_vst,file = "resource/HNSC_exprSet_vst.R")
# 
# ##加载R包
# rm(list = ls())
# 
# library(tidyr)
# library(dplyr)
# library(tibble)
# library(DESeq2)
# library(tibble)
# 
# HNSC <- readRDS("resource/TCGA.rds/TCGA-HNSC.rds")
# 
# load(file = "resource/gene_id_tcga_v36.Rdata")
# 
# HNSC <- as.data.frame(HNSC)
# HNSC <- HNSC %>%
#   rownames_to_column("ID")
# HNSC <- HNSC %>%
#   column_to_rownames("ID") %>%
#   t() %>%
#   as.data.frame() %>%
#   mutate(newcolumn = rowMeans(.)) %>%
#   arrange(desc(newcolumn))
# 
# HNSC$sample <- substr(rownames(HNSC),14,15)
# HNSC$ID <- substr(rownames(HNSC),1,12)
# HNSC <- HNSC %>%
#   rownames_to_column("TCGA_ID")
# 
# HNSC <- HNSC %>%
#   filter(sample=="01") %>%
#   dplyr::select(-sample) %>%
#   dplyr::select(-TCGA_ID) %>%
#   dplyr::select(-newcolumn) %>%
#   dplyr::select(ID,everything()) %>%
#   distinct(ID,.keep_all = T)
# 
# HNSC <- HNSC %>%
#   column_to_rownames("ID") %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column("gene_id")
# 
# expression <- merge(gene_id,HNSC,by="gene_id")
# expression <- expression %>%
#   filter(gene_type=="protein_coding") %>%
#   dplyr::select(-gene_type) %>%
#   dplyr::select(-gene_id) %>%
#   mutate(newcolumn = rowMeans(.[,-1])) %>%
#   arrange(desc(newcolumn)) %>%
#   distinct(gene_name,.keep_all = T) %>%
#   dplyr::select(-newcolumn) %>%
#   column_to_rownames("gene_name")
# 
# exptDesign_TCGA <- data.frame(
#   condition = rep("tumor", length(colnames(expression))),
#   row.names = colnames(expression)
# )
# 
# dds <- DESeqDataSetFromMatrix(
#   countData = expression,
#   colData = exptDesign_TCGA,
#   design = ~ 1)
# 
# nrow(dds)
# rownames(dds)
# 
# ## 筛选样本，counts函数提取表达量数据，取出count为0的,我认为一半表达不为0的更好
# dds <- dds[rowMedians(counts(dds))>0,]
# nrow(dds)
# 
# ### 最终VST标准化
# system.time(vsd <- vst(dds, blind = FALSE))
# exprSet_vst <- as.data.frame(assay(vsd))
# 
# save(exprSet_vst,file = "resource/HNSC_exprSet_vst.R")
rm(list = ls())
library(dplyr)
library(data.table)
library(tidyr)
library(tibble)
library(GSVA)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)
library(ggpubr)
library(rms)
library(foreign)
library(utils)
library(limma)
library(estimate)
library(WGCNA)
library(DESeq2)
library(GSVA)
library(matrixStats)
library(cowplot)
library(ggExtra)
load('resource/cellMarker_ssGSEA.Rdata')
load('resource/exp_TCGA.Rdata')

exprSet_vst <- exp_TCGA %>% 
  filter(OS.time>0.081) %>% 
  dplyr::select(-OS) %>% 
  dplyr::select(-OS.time) %>% 
  column_to_rownames("ID") %>% 
  t() %>% 
  as.data.frame()

platformMap2 <- data.table::fread("resource/gene_list_of_R-loop_regulators.txt",data.table = F)
index0 <- platformMap2$Gene_name
gene_name <- rownames(exprSet_vst)
index <- intersect(index0,gene_name)

load(file = "output/wilcox_data_rloop.R")
wilcox_data_p <- wilcox_data %>% 
  filter(p.value<0.001)  
index <- intersect(index,wilcox_data_p$gene)

breg <- exprSet_vst[index,]
breg <- na.omit(breg)
index <- rownames(breg)

colnames(exp_TCGA) <- gsub("-","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\(","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\)","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\;","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("/","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("@","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub(" ","_",colnames(exp_TCGA))
index <- gsub("-", "_", index)

genes <- index
library(survival)
##批量生存分析
res <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(OS.time, OS)~', genes[i]))
  x = coxph(surv, data = exp_TCGA)
  x = summary(x)
  p.value=signif(x$wald["pvalue"], digits=2)
  HR =signif(x$coef[2], digits=8);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  CI <- paste0("(", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res[i,1] = genes[i]
  res[i,2] = HR
  res[i,3] = CI
  res[i,4] = p.value
}
names(res) <- c("ID","HR","95% CI","p.value")
res_immune <- res
res_immune_p <- res_immune %>%
  filter(p.value<0.05)

index <- res_immune_p$ID
index <- gsub("_", "-", index)



load('resource/cellMarker_ssGSEA.Rdata')
load('resource/exp_TCGA.Rdata')

exp_TCGA <- exp_TCGA %>% 
  filter(OS.time>0.081) %>% 
  dplyr::select(-OS) %>% 
  dplyr::select(-OS.time) %>% 
  column_to_rownames("ID") %>% 
  t() %>% 
  as.data.frame()
exprSet_vst <- exp_TCGA

# -------------------------------------------------------------------------

# gsva_data <- gsva(as.matrix(exprSet_vst),cellMarker, method = "ssgsea")
gsva_data <- gsva(gsvaParam(as.matrix(exprSet_vst), cellMarker))

# ss <- gsva_data

# platformMap2 <- data.table::fread("resource/gene_list_of_R-loop_regulators.txt",data.table = F)
# index <- platformMap2$Gene_name
# 
# # index <- c("IL10", "IRF4", "PRKAR1A", "GRAP", "ACP5", 
# #            "IL12A", "TFEC", "NRDC", "TSC22D3", "SYT11", "IL6", 
# #            "SMAD3", "CD86", "SLAMF1", "IL10RA", "OAS3", 
# #            "IRF5", "FCGR2B", "LAX1", "PAG1", "HLA-DMB", "OAS1", 
# #            "TNFAIP3", "FLNA", "MYC", "CD96")
# ####CSH1表达太低，去除
breg <- exprSet_vst[index,]
breg <- na.omit(breg)




# -------------------------------------------------------------------------

# dir.create('ConsensusCluster/')新增文件夹
## 一致性聚类
results = ConsensusClusterPlus(as.matrix(breg),
                               maxK=9,
                               reps=100,
                               pItem=0.8,
                               pFeature=1,
                               tmyPal = c('navy','darkred'),
                               title='output/Figure/ConsensusCluster/',
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="pdf")

icl <- calcICL(results,title = 'output/Figure/ConsensusCluster/',plot = 'pdf')

sample_subtypes <- results[[2]][["consensusClass"]]
table(sample_subtypes)


immune_data <- rbind(sample_subtypes, breg)
immune_data <- immune_data %>% 
  as.data.frame() %>% 
  t() %>% 
  as.data.frame()
colnames(immune_data)[1] <- "sample_subtypes"
immune_data$sample_subtypes <- as.character(immune_data$sample_subtypes)  # 先转换为字符，避免因子自动转换
immune_data$sample_subtypes[immune_data$sample_subtypes == "1"] <- "C1"
immune_data$sample_subtypes[immune_data$sample_subtypes == "2"] <- "C2"
# 将 ID 列转换为因子，并指定因子水平
immune_data$sample_subtypes <- factor(immune_data$sample_subtypes, levels = c("C1", "C2"))
colnames(immune_data)[1] <- "group"
immune_data <- immune_data %>% 
  rownames_to_column("ID")
load('resource/exp_TCGA.Rdata')

data <- merge(immune_data,exp_TCGA)
data <- data %>% 
  dplyr::select(ID,OS,OS.time,group,everything())


genes <- "group"
your.surv <- Surv(data$OS.time, data$OS)
your.km.plot <- function(genes,data){
  print(genes)
  group <- data[,genes] #分组
  survival_dat <- data.frame(group = group)
  group <- factor(group, levels = c("C1", "C2"))
  fit <- survfit(your.surv ~ group)
  sdf <- survdiff(your.surv ~ group,rho=0)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n)-1)
  p.val
  photo2 <-  ggsurvplot(fit,data = data, #这里很关键，不然会报错
                        legend.title = genes,#定义图例的名称
                        legend.labs = c("C1","C2"), #所以上面要因子化分组顺序
                        #legend = "top",#图例位置
                        pval = T, #在图上添加log rank检验的p值
                        #pval.method = TRUE,#添加p值的检验方法
                        conf.int = TRUE,#添加置信区间
                        risk.table = TRUE, #在图下方添加风险表
                        #risk.table.col = "strata", #根据数据分组为风险表添加颜色
                        risk.table.y.text = F,#风险表Y轴是否显示分组的名称,F为以线条展示分组
                        #linetype = "strata", #改变不同组别的生存曲线的线型
                        #surv.median.line = "hv", #标注出中位生存时间
                        xlab = "Time in years", #x轴标题
                        xlim = c(0,max(data$OS.time)+1), #展示x轴的范围
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
your.km.plot("group",data = data)

library(export)
graph2ppt(file = "output/Figure/00_immune_C1C2_KM.ppt",width=6.3, height=5.5)

rt <- exprSet_vst %>% 
  rownames_to_column("ID")
write.table(rt, file="output/uniq.symbol.txt", sep="\t", quote=F, col.names=T, row.names = F)

filterCommonGenes(input.f="output/uniq.symbol.txt", 
                  output.f="output/commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds="output/commonGenes.gct",
              output.ds="output/estimateScore.gct")

scores <- read.table("output/estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])
rownames(scores) <- gsub("\\.", "\\-", rownames(scores))
scores <- as.data.frame(scores[,1:3])
scores <- scores %>% 
  rownames_to_column("ID")

group <- data %>% 
  dplyr::select(ID,group)

differ <- merge(group,scores,by="ID")
differ <- differ %>% 
  column_to_rownames("ID")


differ <- differ %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "Score")

ggplot(data = differ,aes(x=gene,y=Score,fill=group))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(label = "p.signif")+
  scale_fill_manual(values = c("C1" = "#00BFC4", "C2" = "#F8766D"))

library(export)
graph2ppt(file = "output/Cibersort",width=5,height=4.5)

ESTIMATEScore <- differ %>% 
  filter(gene=="ESTIMATEScore") %>% 
  dplyr::select(group,Score) 

ESTIMATEScore$group <- factor(ESTIMATEScore$group, levels = c("C1", "C2"))
comparisons <- list(c("C1", "C2"))
e <- ggplot(ESTIMATEScore, aes(x = group, y = Score))
e + geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = group), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 2.5) +
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
graph2ppt(file = "output/Figure/00_ESTIMATEScore_compare.ppt",width=4.85, height=5.8)


ImmuneScore <- differ %>% 
  filter(gene=="ImmuneScore") %>% 
  dplyr::select(group,Score) 

ImmuneScore$group <- factor(ImmuneScore$group, levels = c("C1", "C2"))
comparisons <- list(c("C1", "C2"))
e <- ggplot(ImmuneScore, aes(x = group, y = Score))
e + geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = group), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 2.5) +
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
graph2ppt(file = "output/Figure/00_ImmuneScore_compare.ppt",width=4.85, height=5.8)


StromalScore <- differ %>% 
  filter(gene=="StromalScore") %>% 
  dplyr::select(group,Score) 

StromalScore$group <- factor(StromalScore$group, levels = c("C1", "C2"))
comparisons <- list(c("C1", "C2"))
e <- ggplot(StromalScore, aes(x = group, y = Score))
e + geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.5, fill = c("#4DBBD5", "#E64B35"), color = "black", alpha = 1, outlier.shape = NA)+
  geom_jitter(aes(color = group), shape = 16, position = position_jitter(0.2), alpha = 0.65, size = 2.5) +
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
graph2ppt(file = "output/Figure/00_StromalScore_compare.ppt",width=4.85, height=5.8)



box <- gsva_data %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("ID")


differ2 <-  merge(group,box,by="ID")
differ2 <- differ2 %>% 
  pivot_longer(cols=3:30,
               names_to= "celltype",
               values_to = "Score")

differ2$celltype <- factor(differ2$celltype, levels = c(
  "Activated B cell",
  "Activated CD4 T cell",
  "Activated CD8 T cell",
  "Activated dendritic cell",
  "Central memory CD4 T cell",
  "Central memory CD8 T cell",
  "CD56bright natural killer cell",
  "CD56dim natural killer cell",
  "Effector memeory CD4 T cell",
  "Effector memeory CD8 T cell",
  "Eosinophil",
  "Gamma delta T cell",
  "Immature  B cell",
  "Immature dendritic cell",
  "Macrophage",
  "Mast cell",
  "Memory B cell",
  "MDSC",
  "Monocyte",
  "Natural killer cell",
  "Natural killer T cell",
  "Neutrophil",
  "Plasmacytoid dendritic cell",
  "Regulatory T cell",
  "T follicular helper cell",
  "Type 1 T helper cell",
  "Type 17 T helper cell",
  "Type 2 T helper cell"
), ordered = F)

library(ggplot2)
library(ggpubr)
ggplot(data = differ2, aes(x = celltype, y = Score)) +
  geom_boxplot(aes(fill = group), position = position_dodge(1), width = 0.3, outlier.shape = NA) +
  geom_violin(aes(colour = group), position = position_dodge(1), scale = "width", fill = NA) + # 小提琴透明填充，边框颜色根据 risk
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    panel.grid = element_blank()
  ) +
  stat_compare_means(aes(group = group), label = "p.signif") +
  scale_fill_manual(values = c("C1" = "#00BFC4", "C2" = "#F8766D")) +  # 设置箱图填充颜色
  scale_colour_manual(values = c("C1" = "#00BFC4", "C2" = "#F8766D"))  # 设置小提琴边框颜色
graph2ppt(file = "output/Figure/00_immune_box_c1c2.ppt",width=12.5, height=5.1)




# ##########################################################复杂热图
# 
# clinical=read.table("resource/TCGA.HNSC.CLI.txt",header=T,sep="\t",check.names=F,row.names=1) 
# index <- group$ID
# clinical <- clinical[index,]
# clinical$age <- ifelse(clinical$age > 60, ">60", "<=60")
# clinical <- clinical %>% 
#   dplyr::select(gender,age,clinical_T,clinical_N,clinical_M,clinical_stage,
#                 lymphovascular_invasion,perineural_invasion,vital_status,tumor_status)
# clinical <- clinical %>% 
#   mutate(group=group$group)
# 
# 
# ##因子化
# clinical$gender <- factor(clinical$gender,levels = c("female","male"),order=T)
# clinical$age <- factor(clinical$age,levels = c("<=60",">60"),order=T)
# clinical$clinical_T <- factor(clinical$clinical_T,levels = c("T1","T2","T3","T4","unknown"),order=T)
# clinical$clinical_N <- factor(clinical$clinical_N,levels = c("N0","N1","N2","N3","unknown"),order=T)
# clinical$clinical_M <- factor(clinical$clinical_M,levels = c("M0","M1","unknown"),order=T)
# clinical$clinical_stage <- factor(clinical$clinical_stage,levels = c("Stage I","Stage II","Stage III","Stage IV","unknown"),order=T)
# clinical$lymphovascular_invasion <- factor(clinical$lymphovascular_invasion,levels = c("Negative","Positive","unknown"),order=T)
# clinical$perineural_invasion <- factor(clinical$perineural_invasion,levels = c("Negative","Positive","unknown"),order=T)
# clinical$vital_status <- factor(clinical$vital_status,levels = c("Alive","Dead"),order=T)
# clinical$tumor_status <- factor(clinical$tumor_status,levels = c("Tumor free","With tumor","unknown"),order=T)
# clinical$group <- factor(clinical$group,levels = c("C1","C2"),order=T)
# 
# ###########################################把免疫细胞表达矩阵搞出来
# immune_cell <- as.data.frame(t(gsva_data))
# immune_cell <- immune_cell %>% 
#   rownames_to_column("ID")
# ##################################两个矩阵是否顺序相同
# # identical(rownames(clinical),immune_cell$ID)
# clinical <- clinical[match(immune_cell$ID, rownames(clinical)),]
# identical(rownames(clinical),immune_cell$ID)
# 
# ssgsea_df <- column_to_rownames(immune_cell,"ID")
# 
# 
# 
# library(ComplexHeatmap)
# library(circlize)
# 
# 
# 
# columnAnno <- HeatmapAnnotation(gender = clinical$gender,
#                                 age = clinical$age,
#                                 clinical_T = clinical$clinical_T,
#                                 clinical_N = clinical$clinical_N,
#                                 clinical_M = clinical$clinical_M,
#                                 clinical_stage = clinical$clinical_stage,
#                                 lymphovascular_invasion = clinical$lymphovascular_invasion,
#                                 perineural_invasion = clinical$perineural_invasion,
#                                 tumor_status = clinical$tumor_status,
#                                 vital_status = clinical$vital_status,
#                                 group = clinical$group)
# 
# ####这几行代码的整体作用是先对 ssgsea_df 进行标准化和转置处理，
# ####然后对得到的标准化后的数据进行阈值处理，将数据的取值范围限定在 -2 到 2 之间。
# ####这种操作常见于生物信息学数据分析流程中，例如在对基因集富集分析结果进行进一步处理和可视化之前
# ####，通过标准化和阈值处理可以使数据更适合后续的分析和展示。
# # scaled_ssgsea <- scale(t(ssgsea_df))
# # scaled_ssgsea[scaled_ssgsea>2] <- 2
# # scaled_ssgsea[scaled_ssgsea< -2] <- -2
# scaled_ssgsea <- t(ssgsea_df)
# 
# 
# 
# # 
# # ComplexHeatmap::Heatmap(scaled_ssgsea, na_col = "white",show_column_names = F,
# #                         row_names_side = "left",name = "fraction",
# #                         column_order = c(rownames(ssgsea_df)[c(grep("c1",clinical$group),grep("c2",clinical$group))]),
# #                         column_split = clinical$group, column_title = NULL,
# #                         cluster_columns = F,
# #                         top_annotation = columnAnno)
# 
# 
# 
# 
# # 检查维度是否一致
# if (ncol(scaled_ssgsea) != nrow(clinical)) {
#   stop("scaled_ssgsea 的列数和 clinical 的行数不匹配，请检查数据。")
# }
# 
# # 生成 column_order
# column_order <- rownames(ssgsea_df)[order(factor(clinical$group, levels = c("c1", "c2")))]
# 
# 
# 
# # 绘制热图
# ComplexHeatmap::Heatmap(scaled_ssgsea, 
#                         na_col = "white",
#                         show_column_names = FALSE,
#                         row_names_side = "left",
#                         name = "fraction",
#                         column_order = column_order,
#                         column_split = clinical$group, 
#                         column_title = NULL,
#                         cluster_columns = FALSE,
#                         top_annotation = columnAnno)
# # graph2pdf(file = "output/Figure/immune_C1C2_heatmap",width=20, height=10)
# 
# 
# # 定义注释颜色
# anno_colors = list(
#   gender = c("female" = "#F4A99B", "male" = "#015493"),
#   age = c("<=60" = "#9193B4", ">60" = "#2F2D54"),
#   clinical_T = c("T1" = "#79ACAE", "T2" = "#618C95", "T3" = "#496C7C", "T4" = "#324C63", "unknown" = "#E3E3E3"),
#   clinical_N = c("N0" = "#A5B978", "N1" = "#859346", "N2" = "#656D14", "N3" = "#E1E3E3", "unknown" = "#E3E3E3"),
#   clinical_M = c("M0" = "#FBE3C0", "M1" = "#8A4B43", "unknown" = "#E3E3E3"),
#   clinical_stage = c("Stage I" = "#EEA5B4", "Stage II" = "#E57A86", "Stage III" = "#DC4F58", "Stage IV" = "#D4242A", "unknown" = "#E3E3E3"),
#   lymphovascular_invasion = c("Negative" = "#FEE4E8", "Positive" = "#CE8AD8", "unknown" = "#E3E3E3"),
#   perineural_invasion = c("Negative" = "#FAAE5F", "Positive" = "#972D36", "unknown" = "#E3E3E3"),
#   tumor_status = c("Tumor free" = "#00BFC4", "With tumor" = "#F8766D", "unknown" = "#E3E3E3"),
#   vital_status = c("Alive" = "#BABABA", "Dead" = "black"),
#   group = c("C1" = "#2878b5", "C2" = "#c82423")
# )
# 
# # 创建顶部注释并指定颜色
# columnAnno <- HeatmapAnnotation(gender = clinical$gender,
#                                 age = clinical$age,
#                                 clinical_T = clinical$clinical_T,
#                                 clinical_N = clinical$clinical_N,
#                                 clinical_M = clinical$clinical_M,
#                                 clinical_stage = clinical$clinical_stage,
#                                 lymphovascular_invasion = clinical$lymphovascular_invasion,
#                                 perineural_invasion = clinical$perineural_invasion,
#                                 tumor_status = clinical$tumor_status,
#                                 vital_status = clinical$vital_status,
#                                 group = clinical$group,
#                                 col = anno_colors)
# 
# # 绘制热图
# ComplexHeatmap::Heatmap(scaled_ssgsea, 
#                         na_col = "white",
#                         show_column_names = FALSE,
#                         row_names_side = "left",
#                         name = "fraction",
#                         column_order = column_order,
#                         column_split = clinical$group, 
#                         column_title = NULL,
#                         cluster_columns = FALSE,
#                         top_annotation = columnAnno)
# 
# 
# 
# # graph2pdf(file = "output/Figure/immune_C1C2_heatmap",width=20, height=10)
# # save(clinical,file = "output/clinical.Rdata")
# 
# 
# 
# 
# 
# 


##########################################################复杂热图

clinical=read.table("resource/TCGA.HNSC.CLI.txt",header=T,sep="\t",check.names=F,row.names=1) 
index <- group$ID
clinical <- clinical[index,]
clinical$age <- ifelse(clinical$age > 60, ">60", "<=60")
clinical <- clinical %>% 
  dplyr::select(gender,age,pathologic_T,pathologic_N,pathologic_stage,
                lymphovascular_invasion,perineural_invasion,vital_status,tumor_status)
clinical <- clinical %>% 
  mutate(group=group$group)


##因子化
clinical$gender <- factor(clinical$gender,levels = c("female","male"),order=T)
clinical$age <- factor(clinical$age,levels = c("<=60",">60"),order=T)
clinical$pathologic_T <- factor(clinical$pathologic_T,levels = c("T1","T2","T3","T4","unknown"),order=T)
clinical$pathologic_N <- factor(clinical$pathologic_N,levels = c("N0","N1","N2","N3","unknown"),order=T)
# clinical$clinical_M <- factor(clinical$clinical_M,levels = c("M0","M1","unknown"),order=T)
clinical$pathologic_stage <- factor(clinical$pathologic_stage,levels = c("Stage I","Stage II","Stage III","Stage IV","unknown"),order=T)
clinical$lymphovascular_invasion <- factor(clinical$lymphovascular_invasion,levels = c("Negative","Positive","unknown"),order=T)
clinical$perineural_invasion <- factor(clinical$perineural_invasion,levels = c("Negative","Positive","unknown"),order=T)
clinical$vital_status <- factor(clinical$vital_status,levels = c("Alive","Dead"),order=T)
clinical$tumor_status <- factor(clinical$tumor_status,levels = c("Tumor free","With tumor","unknown"),order=T)
clinical$group <- factor(clinical$group,levels = c("C1","C2"),order=T)

###########################################把免疫细胞表达矩阵搞出来
immune_cell <- as.data.frame(t(gsva_data))
immune_cell <- immune_cell %>% 
  rownames_to_column("ID")
##################################两个矩阵是否顺序相同
# identical(rownames(clinical),immune_cell$ID)
clinical <- clinical[match(immune_cell$ID, rownames(clinical)),]
identical(rownames(clinical),immune_cell$ID)

ssgsea_df <- column_to_rownames(immune_cell,"ID")



library(ComplexHeatmap)
library(circlize)



columnAnno <- HeatmapAnnotation(gender = clinical$gender,
                                age = clinical$age,
                                pathologic_T = clinical$pathologic_T,
                                pathologic_N = clinical$pathologic_N,
                                
                                pathologic_stage = clinical$pathologic_stage,
                                lymphovascular_invasion = clinical$lymphovascular_invasion,
                                perineural_invasion = clinical$perineural_invasion,
                                tumor_status = clinical$tumor_status,
                                vital_status = clinical$vital_status,
                                group = clinical$group)

####这几行代码的整体作用是先对 ssgsea_df 进行标准化和转置处理，
####然后对得到的标准化后的数据进行阈值处理，将数据的取值范围限定在 -2 到 2 之间。
####这种操作常见于生物信息学数据分析流程中，例如在对基因集富集分析结果进行进一步处理和可视化之前
####，通过标准化和阈值处理可以使数据更适合后续的分析和展示。
# scaled_ssgsea <- scale(t(ssgsea_df))
# scaled_ssgsea[scaled_ssgsea>2] <- 2
# scaled_ssgsea[scaled_ssgsea< -2] <- -2
scaled_ssgsea <- t(ssgsea_df)



# 
# ComplexHeatmap::Heatmap(scaled_ssgsea, na_col = "white",show_column_names = F,
#                         row_names_side = "left",name = "fraction",
#                         column_order = c(rownames(ssgsea_df)[c(grep("c1",clinical$group),grep("c2",clinical$group))]),
#                         column_split = clinical$group, column_title = NULL,
#                         cluster_columns = F,
#                         top_annotation = columnAnno)




# 检查维度是否一致
if (ncol(scaled_ssgsea) != nrow(clinical)) {
  stop("scaled_ssgsea 的列数和 clinical 的行数不匹配，请检查数据。")
}

# 生成 column_order
column_order <- rownames(ssgsea_df)[order(factor(clinical$group, levels = c("c1", "c2")))]



# 绘制热图
ComplexHeatmap::Heatmap(scaled_ssgsea, 
                        na_col = "white",
                        show_column_names = FALSE,
                        row_names_side = "left",
                        name = "fraction",
                        column_order = column_order,
                        column_split = clinical$group, 
                        column_title = NULL,
                        cluster_columns = FALSE,
                        top_annotation = columnAnno)
graph2pdf(file = "output/Figure/immune_C1C2_heatmap",width=20, height=10)


# 定义注释颜色
anno_colors = list(
  gender = c("female" = "#F4A99B", "male" = "#015493"),
  age = c("<=60" = "#9193B4", ">60" = "#2F2D54"),
  pathologic_T = c("T1" = "#79ACAE", "T2" = "#618C95", "T3" = "#496C7C", "T4" = "#324C63", "unknown" = "#E3E3E3"),
  pathologic_N = c("N0" = "#8ED873", "N1" = "#A5B978", "N2" = "#859346", "N3" = "#656D14", "unknown" = "#E3E3E3"),
  
  pathologic_stage = c("Stage I" = "#EEA5B4", "Stage II" = "#E57A86", "Stage III" = "#DC4F58", "Stage IV" = "#D4242A", "unknown" = "#E3E3E3"),
  lymphovascular_invasion = c("Negative" = "#FEE4E8", "Positive" = "#CE8AD8", "unknown" = "#E3E3E3"),
  perineural_invasion = c("Negative" = "#FAAE5F", "Positive" = "#972D36", "unknown" = "#E3E3E3"),
  tumor_status = c("Tumor free" = "#00BFC4", "With tumor" = "#F8766D", "unknown" = "#E3E3E3"),
  vital_status = c("Alive" = "#BABABA", "Dead" = "black"),
  group = c("C1" = "#2878b5", "C2" = "#c82423")
)

# 创建顶部注释并指定颜色
columnAnno <- HeatmapAnnotation(gender = clinical$gender,
                                age = clinical$age,
                                pathologic_T = clinical$pathologic_T,
                                pathologic_N = clinical$pathologic_N,
                                
                                pathologic_stage = clinical$pathologic_stage,
                                lymphovascular_invasion = clinical$lymphovascular_invasion,
                                perineural_invasion = clinical$perineural_invasion,
                                tumor_status = clinical$tumor_status,
                                vital_status = clinical$vital_status,
                                group = clinical$group,
                                col = anno_colors)

# 绘制热图
ComplexHeatmap::Heatmap(scaled_ssgsea, 
                        na_col = "white",
                        show_column_names = FALSE,
                        row_names_side = "left",
                        name = "fraction",
                        column_order = column_order,
                        column_split = clinical$group, 
                        column_title = NULL,
                        cluster_columns = FALSE,
                        top_annotation = columnAnno)



graph2pdf(file = "output/Figure/00_immune_C1C2_heatmap",width=20, height=12)
save(clinical,file = "output/clinical.Rdata")


