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

##############################################wilcox C1 C2 差异表达基因
load(file = "output/clinical.Rdata")
load(file = "resource/exp_TCGA.Rdata")
clinical <- clinical %>% 
  rownames_to_column("ID")%>% 
  dplyr::select(ID,group)

data <- merge(clinical,exp_TCGA,by="ID")
data <- data %>% 
  dplyr::select(-OS) %>% 
  dplyr::select(-OS.time)


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
  lowposition <- grep("C1",data$group)
  highposition <- grep("C2",data$group)
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

##lapply批量作用于函数，返回list
lapplylist = lapply(colnames(data)[-c(1:2)],my.wilcox)

##第3步do.call 转换list
immune_wilcox_data <- do.call(rbind,lapplylist)
save(immune_wilcox_data,file = "output/immune_wilcox_data_logfc.Rdata")



rm(list = ls())
library(clusterProfiler)
load(file = "output/immune_wilcox_data_logfc.Rdata")
##############################################################
gene_df <- immune_wilcox_data
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
GO <- read.gmt("resource/c5.go.v2024.1.Hs.symbols.gmt")
### 主程序GSEA
y <- GSEA(geneList,TERM2GENE =GO)

yd <- as.data.frame(y)
### 看整体分布
library(ggplot2)
dotplot(y,showCategory=12,
        split=".sign",
        font.size = 9,
        label_format = 60)+facet_grid(~.sign)

library(export)
graph2ppt(file = "output/GSEA_Hallmark.pptx",width=9.8, height=6.4)
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
# graph2pdf(file = "output/Figure/immune_C1C2_up10",width=16, height=12)


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
# graph2pdf(file = "output/Figure/immune_C1C2_down10",width=16, height=12)



differ <- immune_wilcox_data %>% 
  filter(p.value < 0.001) %>% 
  filter(abs(logFC) > 1)



dataset <- gene_df

# 设置p_value和logFC的阈值
cut_off_pvalue = 0.001  #统计显著性
cut_off_logFC = 1      #差异倍数值

# 根据阈值参数，上调基因设置为‘up’，下调基因设置为‘Down’，无差异设置为‘Stable’，并保存到change列中
dataset$change = ifelse(dataset$p.value < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                        ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
head(dataset)



p <- ggplot(
  # 数据、映射、颜色
  dataset, aes(x = logFC, y = -log10(p.value), colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # 辅助线
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # 坐标轴
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p <- p + scale_x_continuous(limits=c(-5,5),breaks=seq(-4,4,2))

p


library(export)
# graph2ppt(file = "output/Figure/00_immune_c1c2differ_heatmap",width=7.8,height=6.4)



index <- differ$gene



load(file = "resource/exp_TCGA.Rdata")
exp_TCGA <- exp_TCGA %>% 
  select(ID, OS, OS.time,all_of(index))

colnames(exp_TCGA)[2] <- "fustat"
colnames(exp_TCGA)[3] <- "futime"

##处理一下，不然会报错
colnames(exp_TCGA) <- gsub("-","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\(","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\)","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\;","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("/","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("@","_",colnames(exp_TCGA))

genes <- colnames(exp_TCGA)[4:ncol(exp_TCGA)]
library(survival)
##批量生存分析
res <- data.frame()
for (i in 1:length(genes)) {
  print(i)
  surv =as.formula(paste('Surv(futime, fustat)~', genes[i]))
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
res <- res %>%
  filter(p.value<0.05)

##########因为前面替换了符号，所以要找到这些基因把名字改回去
gene_ <- res %>%
  filter(grepl("_", ID))

gene <- res$ID
gene <- gsub("_", "-", gene)

save(gene,file = "output/immune_C1C2differ_cox_gene.Rdata")

###############################################################################
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
# graph2ppt(file = "output/Figure/00_immune_c1c2differ_cox_kegg",width=7.8,height=6.4)
### cc分析
gene <- gene
### 第二，基因集(之后会介绍这些基因集合的来源)
GO <- read.gmt("resource/c5.go.v2024.1.Hs.symbols.gmt")
### 通用的富集分析
y <- enricher(gene,TERM2GENE =GO)
test <- as.data.frame(y)
### 画图
dotplot(y)
barplot(y)
library(export)
# graph2ppt(file = "output/Figure/00_immune_c1c2differ_cox_GO",width=7.8,height=6.4)
