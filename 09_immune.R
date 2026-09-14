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
# install.packages("estimate", repos="http://R-Forge.R-project.org")
# install.packages("estimate_1.0.13.tar.gz",repos = NULL,type = "source")
# options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# if(!require("data.table")) install.packages("data.table",update = F,ask = F)
# if(!require("GSVA")) BiocManager::install("GSVA",update = F,ask = F)
# if(!require("matrixStats")) BiocManager::install("matrixStats",update = F,ask = F)
############################################绘制immune 48
############################################绘制immune 48
############################################绘制immune 48
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

data <- exprSet_vst %>% 
  rownames_to_column("ID")



rm(metadata)
rm(risk)
rm(median_value)
rm(sample)


write.table(data, file="output/uniq.symbol.txt", sep="\t", quote=F, col.names=T, row.names = F)

####
filterCommonGenes(input.f="output/uniq.symbol.txt", 
                  output.f="output/commonGenes.gct", 
                  id="GeneSymbol")
####
estimateScore(input.ds="output/commonGenes.gct",
              output.ds="output/estimateScore.gct")

####
scores <- read.table("output/estimateScore.gct", skip=2, header=T, check.names=F)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])
rownames(scores) <- gsub("\\.", "\\-", rownames(scores))
scores <- as.data.frame(scores[,1:3])

risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)
risk <- risk %>% 
  rownames_to_column("id") %>%
  dplyr::select(id,riskScore)
colnames(risk)[1] <- "ID"

scores <- scores %>% 
  rownames_to_column("ID")

data <- merge(risk,scores,by="ID")

data <- data %>% 
  column_to_rownames("ID")

# 计算riskScore的中位数
median_riskScore <- median(data$riskScore, na.rm = TRUE)
# 创建新列risk，如果riskScore大于或等于中位数，则为'high'，否则为'low'
data$risk <- ifelse(data$riskScore >= median_riskScore, 'high', 'low')
data$risk <- factor(data$risk, levels = c('low', 'high'))
data <- data %>% 
  dplyr::select(riskScore,risk,everything())

####
data <- data[,-1]
data <- data %>% 
  pivot_longer(cols=-1,
               names_to= "gene",
               values_to = "Score")
ggplot(data = data,aes(x=gene,y=Score,fill=risk))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(label = "p.format")+
  scale_fill_manual(values = c("low" = "#2878b5", "high" = "#c82423"))
library(export)
graph2ppt(file = "output/Figure/999_Cibersort.pptx",width=5,height=4.5)


##加载免疫细胞cell marker
load(file = "resource/cellMarker_ssGSEA.Rdata")
##加载表达矩阵
data <- exprSet_vst
data <- as.matrix(data)


# ####老版本
# gsva_data <- gsva(data, cellMarker, method = "ssgsea")

####新版本
gsva_data <- gsva(gsvaParam(as.matrix(data), cellMarker))
#Error: [matrixStats (>= 1.2.0)] useNames = NA is defunct. 
#Instead, specify either useNames = TRUE or useNames = FALSE. See also ?
#matrixStats::matrixStats.options
#remotes::install_version("matrixStats", version="1.1.0") 
#restart your session and run previous scripts
### 简单作图看一下
pheatmap(gsva_data, #热图的数据
         cluster_rows = F,#行聚类
         cluster_cols = F,#列聚类，可以看出样本之间的区分度
         ##annotation_col =annotation_col, #标注样本分类
         annotation_legend=TRUE, # 显示注释
         show_rownames = T,# 显示行名
         show_colnames = F,# 显示行名
         ##scale = "row", #以行来标准化
         color =colorRampPalette(c("blue", "white","red"))(100),#调色
         #filename = "heatmap_F.pdf",#是否保存
         cellwidth = 1, cellheight = 10,# 格子比例
         fontsize = 10)


heatdata <- gsva_data %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("ID")
heatdata <- merge(risk,heatdata,by="ID")


rt <- heatdata %>% 
  column_to_rownames("ID") %>% 
  dplyr::select(-riskScore) %>% 
  t() %>% 
  as.data.frame()

Type <- heatdata %>% 
  column_to_rownames("ID") %>%
  dplyr::select(riskScore)
# 计算riskScore的中位数
median_riskScore <- median(Type$riskScore, na.rm = TRUE)
# 创建新列risk，如果riskScore大于或等于中位数，则为'high'，否则为'low'
Type$risk <- ifelse(Type$riskScore >= median_riskScore, 'high', 'low')
Type$risk <- factor(Type$risk, levels = c('low', 'high'))
Type <- Type %>% 
  dplyr::select(risk,riskScore)

##顺序
var="riskScore" 
Type=Type[order(Type[,var]),]
rt=rt[,row.names(Type)]

###调色
ann_colors <- list(risk = c(low = "#2878b5", high = "#c82423"),
                   riskScore = colorRampPalette(c("#2878b5", "white", "#c82423"))(100)
)


# breaksList <- seq(-2.5, 2.5, length.out = 100)
pheatmap(rt, annotation=Type, 
         color =colorRampPalette(c("blue", "white","red"))(100),
         # breaks = breaksList,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale="row",
         show_colnames=FALSE,
         fontsize=8,
         fontsize_row=7,
         fontsize_col=3,
         cellwidth = 0.8,
         cellheight = 15,
         annotation_colors = ann_colors)

library(export)
graph2ppt(file = "output/Figure/999_immune_heatmap.ppt",width=11, height=10)
save(heatdata,file = "output/ssGSVA.Rdata")


############################画相关图~~~~箱图
# 计算riskScore的中位数
median_riskScore <- median(heatdata$riskScore, na.rm = TRUE)
# 创建新列risk，如果riskScore大于或等于中位数，则为'high'，否则为'low'
heatdata$risk <- ifelse(heatdata$riskScore >= median_riskScore, 'high', 'low')
heatdata$risk <- factor(heatdata$risk, levels = c('low', 'high'))
heatdata <- heatdata %>% 
  dplyr::select(riskScore,risk,everything()) %>% 
  column_to_rownames("ID")

####
dd1 <- heatdata %>% 
  pivot_longer(cols=3:30,
               names_to= "celltype",
               values_to = "Score")

dd1$celltype <- factor(dd1$celltype,levels = c(
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
  "Type 2 T helper cell"),ordered = F)
library(ggplot2)
library(ggpubr)
ggplot(data = dd1, aes(x = celltype, y = Score)) +
  geom_boxplot(aes(fill = risk), position = position_dodge(1), width = 0.3, outlier.shape = NA) +
  geom_violin(aes(colour = risk), position = position_dodge(1), scale = "width", fill = NA) + # 小提琴透明填充，边框颜色根据 risk
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    panel.grid = element_blank()
  ) +
  stat_compare_means(aes(group = risk), label = "p.signif") +
  scale_fill_manual(values = c("low" = "#2878b5", "high" = "#c82423")) +  # 设置箱图填充颜色
  scale_colour_manual(values = c("low" = "#2878b5", "high" = "#c82423"))  # 设置小提琴边框颜色


library(export)
graph2ppt(file = "output/Figure/999_ssGSEA",width=12.5,height=5.1)

###############免疫检查点
ICB <- data.table::fread("resource/ICB.txt",data.table = F)
ICB_data <- exprSet_vst %>%
  rownames_to_column("gene")
ICB_data <- merge(ICB,ICB_data,by="gene")
ICB_data <- ICB_data %>%
  column_to_rownames("gene") %>%
  t() %>%
  as.data.frame()
ICB_risk <- heatdata %>%
  dplyr::select(riskScore,risk)
ICB_final <- merge(ICB_risk,ICB_data,by = "row.names", all = F)
ICB_final <- ICB_final %>%
  column_to_rownames("Row.names")

####
dd2 <- ICB_final %>%
  pivot_longer(cols=3:ncol(ICB_final),
               names_to= "ICB",
               values_to = "expression")

ICB_name <- colnames(ICB_data)

dd2$ICB <- factor(dd2$ICB,levels = ICB_name,ordered = F)

library(ggplot2)
library(ggpubr)
ggplot(data =dd2, aes(x = ICB, y = expression))+
  geom_boxplot(aes(fill = risk),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = risk),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=risk), label = "p.signif")+
  scale_fill_manual(values = c("low" = "#00BFC4", "high" = "#F8766D"))
library(export)
graph2ppt(file = "output/Figure/ssGSEA",width=19,height=7)


####只画阳性结果ICB
####
ICB_index <- colnames(ICB_final)
ICB_index <- ICB_index[!ICB_index %in% c("risk", "riskScore")]
###再去掉阴性结果
ICB_index <- ICB_index[!ICB_index %in% c("BTN2A1",
                                         "CD274",
                                         "CD70",
                                         "CD80",
                                         "CD86",
                                         "HLA-A",
                                         "HLA-B",
                                         "HLA-C",
                                         "HLA-E",
                                         "HLA-F",
                                         "HLA-G",
                                         "TDO2",
                                         "VTCN1",
                                         "PDCD1LG2"
                                         )]


ICB_data <- ICB_data[,ICB_index]
ICB_final <- ICB_final %>%
  dplyr::select(1:2)
ICB_final <- merge(ICB_risk,ICB_data,by = "row.names", all = F)
ICB_final <- ICB_final %>%
  column_to_rownames("Row.names")

dd2 <- ICB_final %>%
  pivot_longer(cols=3:ncol(ICB_final),
               names_to= "ICB",
               values_to = "expression")

ICB_name <- colnames(ICB_data)

dd2$ICB <- factor(dd2$ICB,levels = ICB_name,ordered = F)

library(ggplot2)
library(ggpubr)
ggplot(data =dd2, aes(x = ICB, y = expression))+
  geom_boxplot(aes(fill = risk),position = position_dodge(1),width=.7,outlier.shape = NA)+
  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=risk), label = "p.signif")+
  scale_fill_manual(values = c("low" = "#00BFC4", "high" = "#F8766D"))
library(export)
graph2ppt(file = "output/Figure/999_ICB",width=16,height=5.0)


# ###########################气泡图
# #####免疫细胞相关性分析
# ####
# exprSet1 <- as.data.frame(t(heatdata))
# gene <- "riskScore"
# ####
# batch_cor <- function(gene){
#   y <- as.numeric(exprSet1[gene,])
#   rownames <- rownames(exprSet1)[3:30]
#   do.call(rbind,future_lapply(rownames, function(x){
#     dd  <- cor.test(as.numeric(exprSet1[x,]),y,method="spearman")
#     data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
#   }))
# }
# library(future.apply)
# # plan(multiprocess)
# system.time(dd <- batch_cor(gene))
# 
# ####
# library(ggpubr)
# library(reshape)
# write.table(dd,file="output/dd.xls",sep="\t",row.names=F,quote=F)
# bc <- read.table("output/dd.xls",header = T,sep="\t",check.names=F)
# bc <- bc %>% 
#   arrange(cor)
# cor_index <- bc$mRNAs
# bc$Cor <- abs(bc$cor)
# library(tibble)
# library(dplyr)
# 
# # bc$mRNAs <- factor(bc$mRNAs,levels = c("Effector memeory CD4 T cell",
# #                                        "Type 2 T helper cell",
# #                                        "Activated CD4 T cell",
# #                                        "Immature dendritic cell",
# #                                        "Effector memeory CD8 T cell",
# #                                        "Immature  B cell",
# #                                        "Effector memeory CD8 T cell",
# #                                        "Activated CD8 T cell",
# #                                        "Type 1 T helper cell",
# #                                        "Macrophage",
# #                                        "Activated B cell",
# #                                        "Gamma delta T cell",
# #                                        "Natural killer T cell",
# #                                        "Activated CD4 T cell",
# #                                        "CD56bright natural killer cell",
# #                                        "Natural killer cell",
# #                                        "Central memory CD8 T cell",
# #                                        "Effector memeory CD4 T cell",
# #                                        "Immature dendritic cell",
# #                                        "Mast cell",
# #                                        "Plasmacytoid dendritic cell",
# #                                        "Type 17 T helper cell",
# #                                        "Monocyte",
# #                                        "Neutrophil",
# #                                        "Eosinophil",
# #                                        "CD56dim natural killer cell",
# #                                        "Type 2 T helper cell",
# #                                        "Memory B cell"),ordered = F)
# bc$mRNAs <- factor(bc$mRNAs,levels = cor_index,ordered = F)
# ggplot(bc, aes(mRNAs, cor),rotate = T) + 
#   geom_segment(aes(xend=mRNAs, yend = 0),linetype = "dashed",colour = "orange") +
#   geom_point(aes(color = p.value, size = Cor)) +
#   scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
#   scale_color_continuous(low="brown1", high="cyan", guide=guide_colorbar(reverse=TRUE))+
#   scale_size_continuous(range=c(2, 8)) +
#   theme_minimal() + 
#   xlab("Cor") +
#   ylab(NULL)+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))
# library(export)
# # graph2ppt(file = "output/immune_correlation",width=12,height=5)
# 
# ############################################绘制risk immune 相关性散点图
# ############################################绘制risk immune 相关性散点图
# ############################################绘制risk immune 相关性散点图
# rm(list = ls())
# load(file = "output/ssGSVA.Rdata")
# data <- heatdata %>% 
#   column_to_rownames("ID")
# colnames(data) <- gsub(" ","_",colnames(data))
# 
# ### 测试数据和相关性分析的方法
# ggplot(data,aes(riskScore,Activated_CD8_T_cell))+
#   geom_point(col="#45b2cb",size=3,alpha=0.7,stroke=1)+
#   geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=1.5,col="#45b2cb")+
#   stat_cor(method = "pearson", digits = 3, size=5)+
#   theme_bw()+
#   theme(
#     plot.title = element_text(hjust = 0.5),
#     plot.margin = margin(1, 1, 1, 1, "cm"),
#     axis.title = element_text(size = 15),  # 调整坐标轴标题的字体大小
#     axis.text = element_text(size = 15)    # 调整坐标轴刻度标签的字体大小
#   )
# library(export)
# graph2ppt(file = "output/Activated_CD8_T_cell",width=6,height=5.4)



rm(list = ls())
#install.packages("scales")
#install.packages("ggtext")

#加载
library(limma)
library(scales)
library(ggplot2)
library(ggtext)



#读取免疫细胞浸润文件
immune=read.csv("resource/infiltration_estimation_for_tcga.csv", header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune=avereps(immune)

#读取矩阵文件
risk=read.table("output/riskoutside1.txt", header=T, sep="\t", check.names=F,row.names = 1)


#对风险文件和免疫细胞浸润文件取交集，得到交集样品
sameSample=intersect(row.names(risk), row.names(immune))
#以"A1BG"这个基因的表达水平为例
data=risk[sameSample, "riskScore"]
immune=immune[sameSample,]

#与免疫细胞的相关性分析
x=as.numeric(data)
outTab=data.frame()
for(i in colnames(immune)){
  y=as.numeric(immune[,i])
  corT=cor.test(x, y, method="pearson")
  cor=corT$estimate
  pvalue=corT$p.value
  if(pvalue<0.05){
    outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
  }
}

#导出结果
write.table(file="output/corResult.txt", outTab, sep="\t", quote=F, row.names=F)

#气泡图
corResult=read.table("output/corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(hue_pal()(length(levels(b$Software))),table(b$Software))     #定义颜色
# pdf(file="cor.pdf", width=10, height=10)       #保存图片
ggplot(data=b, aes(x=cor, y=immune, color=Software))+
  labs(x="Correlation coefficient",y="Immune cell")+
  geom_point(size=4.1)+
  theme(panel.background=element_rect(fill="white",size=1,color="black"),
        panel.grid=element_line(color="grey75",size=0.5),
        axis.ticks = element_line(size=0.5),
        axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))
# dev.off()
library(export)
graph2pdf(file = "output/999_mul_immune",width=7.5, height=14)



rm(list = ls())
#加载
library(ComplexHeatmap) 
library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
library(limma)
library(scales)
library(ggplot2)
library(ggtext)



#读取免疫细胞浸润文件
immune=read.csv("resource/infiltration_estimation_for_tcga.csv", header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
colnames(immune) = gsub("_","--",colnames(immune))
rownames(immune)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
#dim(immune)
immune=avereps(immune)
#dim(immune)

#读入
risk=read.table("output/riskoutside1.txt",header=T,sep="\t",row.names=1,check.names=F)

#取交集样本
sameSample=intersect(rownames(immune),rownames(risk))
immune=immune[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]

#创建type
immMethod <- sub("^.*--", "", colnames(immune))
type <- data.frame(immMethod, row.names = colnames(immune))
names(type) <- "Methods"

#标准化
standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

# 数据整理
annCol <- data.frame(RiskType = risk$risk,
                     row.names = rownames(risk),
                     stringsAsFactors = F)
annRow <- data.frame(row.names = colnames(immune),
                     Methods = factor(immMethod,levels = unique(immMethod)),
                     stringsAsFactors = F)

annColors <- list("RiskType" = c("high" = "#AA1701","low" = "#0048A1"))

indata <- t(immune)
indata <- indata[,colSums(indata) > 0]
plotdata <- standarize.fun(indata,halfwidth = 2)
samorder <- rownames(risk[order(risk$riskScore),])
plotdata1 <- plotdata[rownames(annRow[which(annRow$Methods == "TIMER"),,drop = F]),]
plotdata2 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT"),,drop = F]),]
plotdata3 <- plotdata[rownames(annRow[which(annRow$Methods == "CIBERSORT-ABS"),,drop = F]),]
plotdata4 <- plotdata[rownames(annRow[which(annRow$Methods == "QUANTISEQ"),,drop = F]),]
plotdata5 <- plotdata[rownames(annRow[which(annRow$Methods == "MCPCOUNTER"),,drop = F]),]
plotdata6 <- plotdata[rownames(annRow[which(annRow$Methods == "XCELL"),,drop = F]),]
plotdata7 <- plotdata[rownames(annRow[which(annRow$Methods == "EPIC"),,drop = F]),]

#画图
hm1 <- pheatmap(mat = as.matrix(plotdata1[,samorder]),
                border_color = NA,
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                annotation_col = annCol[samorder,,drop = F],
                annotation_colors = annColors,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "TIMER")

hm2 <- pheatmap(mat = as.matrix(plotdata2[,samorder]),
                border_color = NA,
                color = greenred(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "CIBERSORT")

hm3 <- pheatmap(mat = as.matrix(plotdata3[,samorder]),
                border_color = NA,
                color = blueyellow(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "CIBERSORT-ABS")

hm4 <- pheatmap(mat = as.matrix(plotdata4[,samorder]),
                border_color = NA,
                color = bluered(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "QUANTISEQ")

hm5 <- pheatmap(mat = as.matrix(plotdata5[,samorder]),
                border_color = NA,
                color = inferno(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "MCPCOUNTER")

hm6 <- pheatmap(mat = as.matrix(plotdata6[,samorder]),
                border_color = NA,
                color = viridis(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "XCELL")

hm7 <- pheatmap(mat = as.matrix(plotdata7[,samorder]),
                border_color = NA,
                color = magma(64), 
                cluster_rows = F,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = F,
                cellwidth = 0.8,
                cellheight = 10,
                gaps_col = table(annCol$RiskType)[2],
                name = "EPIC")

pdf("immune.pdf", width = 12,height = 24)
draw(hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5 %v% hm6 %v% hm7, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
invisible(dev.off())

draw(hm1 %v% hm2 %v% hm3 %v% hm4 %v% hm5 %v% hm6 %v% hm7, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")

