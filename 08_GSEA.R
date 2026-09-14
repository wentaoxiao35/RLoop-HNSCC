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
library(WGCNA)
library(DESeq2)

##############################################绘制GSEA 47
##############################################绘制GSEA 47
##############################################绘制GSEA 47
# rm(list = ls())
# load(file = "resource/gene_id_tcga_v36.Rdata")
# HNSC <- readRDS("resource/TCGA.rds/TCGA-HNSC.rds")
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
# exprSet_vst <- exprSet_vst %>% 
#   t() %>% 
#   as.data.frame()
# 
# exprSet_vst <- exprSet_vst %>% 
#   rownames_to_column("ID")

load('resource/exp_TCGA.Rdata')
exp_TCGA <- exp_TCGA %>% 
  filter(OS.time>0.081) %>% 
  dplyr::select(-OS) %>% 
  dplyr::select(-OS.time)
exprSet_vst <- exp_TCGA

risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)
risk <- risk %>%
  rownames_to_column("ID") %>%
  dplyr::select(ID,riskScore)

exprSet_vst <- merge(exprSet_vst,risk,by="ID")
exprSet_vst <- exprSet_vst %>% 
  column_to_rownames("ID") %>% 
  dplyr::select(-riskScore)

sample <- exprSet_vst %>% 
  rownames_to_column("ID") %>% 
  dplyr::select(ID)

exprSet_vst <- exprSet_vst %>% 
  t() %>% 
  as.data.frame()

risk <- merge(sample,risk,by="ID")
median_value <- median(risk$riskScore, na.rm = TRUE)
risk$riskScore <- ifelse(risk$riskScore <= median_value, "low", "high")

metadata <- risk
colnames(metadata)[2] <- "group"
colnames(metadata)[1] <- "TCGA_ID"

exprSet_1 <- exprSet_vst %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("TCGA_ID")

###数据处理
data <- inner_join(metadata,exprSet_1,by="TCGA_ID")
data$group <- as.factor(data$group)
table(data$group)
# save(data,file = "resource/data_vst.Rdata")
#画图
library(ggpubr)
ggboxplot(
  data, x = "group", y = "TP53",
  color = "group", palette = c("#00AFBB", "#E7B800"),
  add = "jitter"
)+
  stat_compare_means(method = "wilcox.test")

####使用lapply + function模式来是计算，lapply 以及do.call 的介绍，lapply 三步走：
##第1,写出单次处理的function
my.wilcox = function(x){
  dd <- wilcox.test(data[,x] ~ group, data = data)
  lowposition <- grep("low",data$group)
  highposition <- grep("high",data$group)
  lowcol <- data[,x][lowposition]
  highcol<- data[,x][highposition]
  lowMean = mean(lowcol)
  highMean = mean(highcol)
  logFC=highMean-lowMean
  data.frame(gene=x,
             p.value=dd$p.value,
             logFC
  )
}

##测试函数功能
my.wilcox("RPGR")
# ##lapply批量作用于函数，返回list，为了方便，我们先用100基因试试
# lapplylist = lapply(colnames(data)[-c(1:2)][1:100],my.wilcox)

##lapply批量作用于函数，返回list
lapplylist = lapply(colnames(data)[-c(1:2)],my.wilcox)

##第3步do.call 转换list
wilcox_data <- do.call(rbind,lapplylist)
save(wilcox_data,file = "output/wilcox_data_logfc.Rdata")



rm(list = ls())
library(clusterProfiler)
load(file = "output/wilcox_data_logfc.Rdata")
##############################################################
gene_df <- wilcox_data
## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$gene
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
#################################################################
### GSEA 变化在于gene set
#################################################################
### 1.hallmarks gene set
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("resource/h.all.v2023.1.Hs.symbols.gmt")
### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =hallmarks)

yd <- as.data.frame(y)
### 看整体分布
library(ggplot2)
enrichplot::dotplot(y,showCategory=8,
        split=".sign",
        font.size = 9,
        label_format = 60)+facet_grid(~.sign)

library(export)
graph2ppt(file = "output/Figure/GSEA_Hallmark.pptx",width=9.5, height=7.5)

hallmarks <- read.gmt("resource/c5.go.v2024.1.Hs.symbols.gmt")
### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =hallmarks)

yd <- as.data.frame(y)
### 看整体分布
library(ggplot2)
dotplot(y,showCategory=12,
        split=".sign",
        font.size = 9,
        label_format = 60)+facet_grid(~.sign)

library(export)
graph2ppt(file = "output/Figure/GSEA_GO.pptx",width=9.5, height=7.5)
y1 <- data.frame(y)
y1 <- filter(y,NES>0)
dotplot(y1,showCategory=15,
        split=".sign",
        font.size = 9,
        label_format = 40)+facet_grid(~.sign)
yd1 <- as.data.frame(y1)
yd1 <- yd1 %>% 
  arrange(desc(NES)) 

index <- yd1$ID[1:10]
library(enrichplot)
gseaplot2(y, geneSetID = index)
graph2pdf(file = "output/Figure/GSEA_GO_up10",width=16, height=12)


y2 <- data.frame(y)
y2 <- filter(y,NES<0)
dotplot(y2,showCategory=15,
        split=".sign",
        font.size = 9,
        label_format = 40)+facet_grid(~.sign)
yd2 <- as.data.frame(y2)
yd2 <- yd2 %>% 
  arrange(NES) 

index <- yd2$ID[1:10]
library(enrichplot)
gseaplot2(y, geneSetID = index)
graph2pdf(file = "output/Figure/GSEA_GO_down10",width=16, height=12)




# df <- ggplot2::fortify(y, showCategory = 30, split=".sign")
### https://github.com/YuLab-SMU/enrichplot/blob/master/R/method-fortify.R

# ### 选择需要呈现的来作图
# library(enrichplot)
# 
# pathway.id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
# gseaplot2(y,
#           color = "blue",
#           geneSetID = pathway.id,
#           pvalue_table = T)
# # graph2ppt(file = "output/HALLMARK_TNFA_SIGNALING_VIA_NFKB.pptx",width=11.5, height=6)
# 
# # gseaplot2(y,10,color = "red",pvalue_table = T)
# # gseaplot2(y, geneSetID = 1:3)


####
### 1.hallmarks gene set
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("resource/h.all.v2023.1.Hs.symbols.gmt")
### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =hallmarks)

library(enrichplot)
gseaplot2(y, geneSetID = c("HALLMARK_ALLOGRAFT_REJECTION",
                           "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                           "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                           "HALLMARK_IL2_STAT5_SIGNALING",
                           "HALLMARK_KRAS_SIGNALING_DN",
                           "HALLMARK_BILE_ACID_METABOLISM"))
library(export)
# graph2ppt(file = "output/Figure/GSEA_Hallmark_down.pptx",width=9.9, height=9.2)
graph2pdf(file = "output/Figure/GSEA_Hallmark_down.pdf",width=9.9, height=9.2)

gseaplot2(y, geneSetID = c("HALLMARK_MYC_TARGETS_V1",
                           "HALLMARK_MTORC1_SIGNALING",
                           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                           "HALLMARK_GLYCOLYSIS",
                           "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                           "HALLMARK_MYOGENESIS",
                           "HALLMARK_APICAL_JUNCTION",
                           "HALLMARK_HYPOXIA"))
library(export)
# graph2ppt(file = "output/Figure/GSEA_Hallmark_up.pptx",width=9.9, height=9.2)
graph2pdf(file = "output/Figure/GSEA_Hallmark_up.pdf",width=9.9, height=9.2)







# ###如果结果都不好，试试找差异基因再做分析
# ###################差异分析后GO,hallmark
# rm(list = ls())
# library(clusterProfiler)
# load(file = "output/wilcox_data_logfc.Rdata")
# p_wilcox <- wilcox_data %>% 
#   filter(p.value < 0.05) %>% 
#   filter(abs(logFC) >0.264)
# library(clusterProfiler)
# gene <- as.character(p_wilcox[,1])
# gene = bitr(gene,
#             fromType="SYMBOL",
#             toType="ENTREZID",
#             OrgDb="org.Hs.eg.db")
# head(gene)
# ### 尝试新的GO聚类的方法
# go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
# # save(go,file = "output/go.Rdata")
# library(ggplot2)
# p <- barplot(go, label_format = 60,split="ONTOLOGY") +
#   facet_grid(ONTOLOGY~., scale="free")
# p
# # graph2ppt(file="output/diff_go.pptx",width=8.2, height=6.8)
# 
# ######################################################################
# ### 我们知道了富集分析的原理，那么只要有个基因集，我们就可以用来富集
# ### 第一，输入差异基因
# gene <- gene$SYMBOL
# ### 第二，基因集(之后会介绍这些基因集合的来源)
# hallmark <- read.gmt("resource/h.all.v2023.1.Hs.symbols.gmt")
# ### 通用的富集分析
# y <- enricher(gene,TERM2GENE =hallmark)
# test <- as.data.frame(y)
# ### 画图
# dotplot(y)
# barplot(y)
# # graph2ppt(file="output/diff_hallmark.pptx",width=8.2, height=6.8)