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
library(clusterProfiler)
library(org.Hs.eg.db)

##############################################绘制WGCNA 46
##############################################绘制WGCNA 46
##############################################绘制WGCNA 46
##http://www.bio-info-trainee.com/7674.html
#WGCNA
# rm(list = ls())
# load(file = "resource/gene_id_tcga_v36.Rdata")
# HNSC <- readRDS("resource/TCGA.rds/TCGA-HNSC.rds")
# HNSC <- as.data.frame(HNSC)
# 
# HNSC <- HNSC %>% 
#   rownames_to_column("ID")
# 
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

exprSet_vst <- exprSet_vst %>% 
  t() %>% 
  as.data.frame()

datExpr0 <- exprSet_vst

###################################################################################
### 1.goodSamplesGenes 的输入条件
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

### 如果没有达标就需要筛选

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    #printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###过滤后的最新表达数据
datExpr <- datExpr0 %>% 
  t() %>% 
  as.data.frame()

###再重新与risk取交集，得到最终risk
risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)
risk <- risk %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,riskScore)
sample <- datExpr %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id)
risk <- merge(sample,risk,by="id")
median_value <- median(risk$riskScore, na.rm = TRUE)
risk$riskScore <- ifelse(risk$riskScore <= median_value, "low", "high")
colnames(risk)[1] <- "ID"
colnames(risk)[2] <- "condition"


### 设置多线程
# enableWGCNAThreads()


### 软阈值的确定，
### 使用函数pickSoftThreshold
### 人为给定一些值
powers = c(c(1:10), seq(from = 12, to=30, by=2))
length(powers)
### 这个power的长度，决定了最终返回结果的行数目

# Call the network topology analysis function
### 模块这里都是作用于表达量数据，没有性状数据的事情
### 稍微耗时间。
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType="signed", verbose = 5)

test <- sft$fitIndices
### 这时候结果已经出来，但是可以通过画图再来看

### 画图呈现结果
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

library(export)
graph2ppt(file = "output/Figure/46_Scale_independence.ppt",width=8.5, height=6.2)

### 一步法网络构建以及模块发现
### power = 12，说明理由
##  maxBlockSize = 5000, 根据自己的电脑可以调整，16GB内存20000，32GB内存30000
### TOMType 
### 这一步是所有分析的基石
### 防止报错
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 12, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 100,
                       networkType="signed",
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "yourname", 
                       verbose = 3)


cor<-stats::cor



### 返回的net是个列表，Module Epigene 表达数据也在里面
table(net$colors)
### 0表示的是不在模块内的基因

### 画图
# open a graphics window
sizeGrWindow(12, 9)
### 标签转为颜色
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plot(net$dendrograms[[1]])
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

library(export)
# graph2ppt(file = "output/Figure/47_Cluster_dendrogram.ppt",width=9, height=7.5)
graph2pdf(file = "output/Figure/47_Cluster_dendrogram.pdf",width=9, height=7.5)
### 保存数据
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
geneTree = net$dendrograms[[1]]
# Recalculate MEs with color labels
test <- moduleEigengenes(datExpr, moduleColors)
### 提取结果
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
### 重新排序，改变的是列的位置 ，没有实际意义，主要用于绘图
MEs = orderMEs(MEs0)
#############################################################################
### 很重要的步骤，很简单的操作
### 求相关性，datTraits 的要求是什么？？

datTraits <- risk
datTraits$condition <- factor(datTraits$condition,levels = c("low","high"))
design=model.matrix(~0+ datTraits$condition)
design <- as.data.frame(design)
colnames(design)=levels(datTraits$condition)

moduleTraitCor = cor(MEs, design, use = "p")
### 算p值
nSamples <- nrow(datExpr)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
### 神奇操作！！
textMatrix =  paste(signif(moduleTraitCor, 2), "(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


### 如何把基因数目显示出来
dd1 <- data.frame(ME=names(MEs),color=substring(names(MEs),3))
dd2 <- data.frame(table(moduleColors))
colnames(dd2) <- c("color","num")
## 合并
dd3 <- merge(dd1,dd2,by="color")
rownames(dd3) <- dd3$ME
dd3 <- dd3[names(MEs),]

### c(bottom, left, top, right)
par(mar = c(6, 5, 3, 3))
ynumbers <- paste0(names(MEs),paste0("(",dd3$num,")"))
ynumbers
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = ynumbers,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### 最重要的图已经获取了
#######################################################################################
library(export)
graph2ppt(file = "output/Figure/48_Module_trait_relationships.ppt",width=5.5, height=6.6)
### 提取感兴趣基因作图
#######################################################################################

### 提取感兴趣基因作图
### 肯定是从表达矩阵中提取
module = "yellow"
moduleGenes = moduleColors==module
datExpr <- as.data.frame(datExpr)
dd <- datExpr[,moduleGenes]

library(dplyr)
library(tibble)
library(tidyr)
test <- dd[1:10,1:10]
dd <- dd %>% 
  rownames_to_column("sample") %>% 
  mutate(group = datTraits$condition) %>% 
  dplyr::select(group,everything()) %>% 
  pivot_longer(cols = 3:ncol(.),
               names_to = "gene",
               values_to = "expression") 

# library(ggplot2)
# ggplot(dd,aes(x = group,y = expression)) +
#   geom_boxplot(aes(fill=group))

###############################################################################
### 探索这个模块中的基因

design
mylove = as.data.frame(design$high)
names(mylove) = "mylove"

## 从第三位开始选取
modNames = substring(names(MEs), 3)

### 计算基因的相关性
### module membership MM
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

### Gene Significance GS
geneTraitSignificance = as.data.frame(cor(datExpr, mylove$mylove, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(mylove), sep="")
names(GSPvalue) = paste("p.GS.", names(mylove), sep="");

### 举例子
module = "yellow"
column = match(module, modNames)
table(moduleColors)
### 提取模块内的基因
moduleGenes = moduleColors==module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
### 本质上就是一个散点图
MM <- geneModuleMembership[moduleGenes, column]
GS <- geneTraitSignificance[moduleGenes, 1]

verboseScatterplot(MM,GS,
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


library(export)
 graph2ppt(file = "output/Figure/Module_membership_yellow",width=6, height=6.1)


###############################################################################
### 自己动手来找hub gene
MM  <- geneModuleMembership[moduleGenes, column]
MMP <- MMPvalue[moduleGenes, column]
GS <-  geneTraitSignificance[moduleGenes, 1]
GSP <- GSPvalue[moduleGenes, 1]

mydata_yellow <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
# mydata_black <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
# mydata_greenyellow <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
# mydata_darkgreen <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
# mydata_tan <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
# mydata_green <- data.frame(moduleGenes=colnames(datExpr)[moduleGenes],MM,MMP,GS,GSP)
mydata <- rbind(mydata_yellow)



###############################################################################
### KEGG分析
gene <- as.character(mydata$moduleGenes)
# #基因名称转换，返回的是数据框
# library(clusterProfiler)
# library(org.Hs.eg.db)
# gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# head(gene)
# # write.csv(gene,file = "output/total.csv",row.names = F)
# # 
# # 
#####################################################################
#####################################################################
### KEGG分析
### 第一，输入差异基因
gene <- gene
### 第二，基因集(之后会介绍这些基因集合的来源)
kegg <- read.gmt("resource/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
### 通用的富集分析
y <- enricher(gene,TERM2GENE =kegg)
test <- as.data.frame(y)
### 画图
dotplot(y)
barplot(y)
library(export)
graph2ppt(file = "output/Figure/00_WGCNA_KEGG",width=7.8,height=6.4)
### cc分析
gene <- gene
### 第二，基因集(之后会介绍这些基因集合的来源)
cc <- read.gmt("resource/c5.go.cc.v2023.1.Hs.symbols.gmt")
### 通用的富集分析
y <- enricher(gene,TERM2GENE =cc)
test <- as.data.frame(y)
### 画图
dotplot(y)
barplot(y)

### mf分析
gene <- gene
### 第二，基因集(之后会介绍这些基因集合的来源)
mf <- read.gmt("resource/c5.go.mf.v2023.1.Hs.symbols.gmt")
### 通用的富集分析
y <- enricher(gene,TERM2GENE =mf)
test <- as.data.frame(y)
### 画图
dotplot(y)
barplot(y)

### bp分析
gene <- gene
### 第二，基因集(之后会介绍这些基因集合的来源)
bp <- read.gmt("resource/c5.go.bp.v2023.1.Hs.symbols.gmt")
### 通用的富集分析
y <- enricher(gene,TERM2GENE =bp)
test <- as.data.frame(y)
### 画图
dotplot(y)
barplot(y)

### GO分析
gene <- gene
### 第二，基因集(之后会介绍这些基因集合的来源)
go <- read.gmt("resource/c5.go.v2024.1.Hs.symbols.gmt")
### 通用的富集分析
y <- enricher(gene,TERM2GENE =go)
test <- as.data.frame(y)
### 画图
dotplot(y)
barplot(y)
library(export)
graph2ppt(file = "output/Figure/00_WGCNA_GO",width=7.8,height=6.4)




# ### 尝试新的GO聚类的方法
# gene = bitr(gene,
#             fromType="SYMBOL",
#             toType="ENTREZID",
#             OrgDb="org.Hs.eg.db")
# go <- enrichGO(gene = gene$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all")
# library(ggplot2)
# p <- barplot(go, label_format = 60,split="ONTOLOGY") +
#   facet_grid(ONTOLOGY~., scale="free")
# p
# # library(export)
# # graph2ppt(file = "output/WGCNA_go",width=8.4, height=7.8)