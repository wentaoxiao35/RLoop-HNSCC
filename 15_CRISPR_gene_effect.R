##清空数据
rm(list = ls())

##加载R包
library(tidyr)
library(dplyr)
library(tibble)

CRISPRGeneEffect = data.table::fread(file = "resource/CRISPRGeneEffect.csv",data.table = F)
colnames(CRISPRGeneEffect)[1] <- c("ModelID")

cell = data.table::fread(file = "resource/Model.csv",data.table = F)

HNSC <- cell %>%
  filter(OncotreePrimaryDisease=="Head and Neck Squamous Cell Carcinoma") %>%
  dplyr::select(ModelID,CellLineName)

HNSC_PLEG <- merge(HNSC,CRISPRGeneEffect)

HNSC_PLEG <- HNSC_PLEG %>%
  column_to_rownames("ModelID") %>%
  dplyr::select(-CellLineName) %>%
  t() %>%
  as.data.frame()

HNSC_PLEG <- HNSC_PLEG %>%
  rownames_to_column("Gene")

HNSC_PLEG$Gene <- gsub("\\s*\\(.*", "", HNSC_PLEG$Gene)

mygene <- HNSC_PLEG %>% 
  filter(Gene == "OLR1") %>% 
  column_to_rownames("Gene") %>% 
  t() %>% 
  as.data.frame()

mygene$cell_line <- rownames(mygene)
index <- rownames(mygene)
cellname <- cell[cell$ModelID %in% index,] 
cellname <- dplyr::select(cellname,ModelID,CellLineName)

colnames(mygene)[2] <- "ModelID"

mygene <- inner_join(mygene,cellname,by="ModelID")
mygene <- dplyr::select(mygene,-ModelID)

library(ggplot2)
# 为了确保fdps是按照从大到小排序的，我们可以先对数据框进行排序
mygene <- mygene[order(mygene$OLR1), ]  # 注意这里的负号，它确保了降序排序

# 接下来，使用ggplot2绘制柱状图
p <- ggplot(data = mygene, aes(x = reorder(CellLineName, OLR1), y = OLR1, fill = CellLineName)) +
  geom_bar(stat = "identity", color = "black") +  # 使用identity来直接绘制y的值，并添加黑色边框
  xlab("CellLineName") + ylab("Gene Expression (OLR1)") +
  ggtitle("Gene Expression by CellLineName (OLR1)") +  # 可选：添加图表标题
  theme_minimal() +  # 使用一个更简洁的主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 可选：将x轴标签旋转45度，以避免重叠
p

p <- ggplot(data = mygene, aes(x = reorder(CellLineName, OLR1), y = OLR1, fill = "all_same")) +
  geom_bar(stat = "identity", color = "black", fill = "#4ABBAD") +  # 使用identity来直接绘制y的值，并添加黑色边框，同时设置所有柱子的填充颜色为淡红色
  xlab("Cell Line") + ylab("CERES score of OLR1") +
  ggtitle("CRISPRGeneEffect of OLR1") +  # 可选：添加图表标题
  theme_minimal() +  # 使用一个更简洁的主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 可选：将x轴标签旋转45度，以避免重叠。
        legend.position = "none")  # 如果不想显示图例（因为你没有基于不同类别填充颜色），可以添加这一行来隐藏图例
# 显示图表
print(p)
library(export)
# graph2ppt(file = "output/Figure/CRISPRGeneEffect of OLR1.ppt",width=14, height=4)
