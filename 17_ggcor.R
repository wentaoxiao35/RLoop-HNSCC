rm(list = ls())
# devtools::install_github("hannet91/ggcor")
library(devtools)
library(ggcor)
library(ggplot2)
library(dplyr)
library(GSVA)
library(tibble)

data1=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)
data1 <- data1 %>%
  dplyr::select(-futime) %>%
  dplyr::select(-fustat) %>%
  dplyr::select(-risk)

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
gsva_data <- gsva(gsvaParam(as.matrix(exprSet_vst), cellMarker))

data2 <- as.data.frame(gsva_data)
data2 <- data2 %>%
  t() %>%
  as.data.frame()

data1 <- data1 %>%
  mutate(riskScore = riskScore + 10)


library(vegan)
library(dplyr)
# data("varechem")
# data("varespec")
# set.seed(20191224)

mantel <- fortify_mantel(data1, data2, 
                         spec.select = list(
                           SERPINE1 = 1:1,
                           HMGA2 = 2:2,
                           DKK1 = 3:3,
                           KDM5D = 4:4,
                           SYCP2 = 5:5,
                           ZNF831 = 6:6,
                           riskScore = 7:7)) %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf), 
                 labels = c("<0.2", "0.2-0.4",">0.4"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                       labels = c("<0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE))

quickcor(data2, type = "lower") + 
  geom_square() + 
  add_link(mantel, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", 
                       midpoint = 0.4) + # 负相关是红色，正相关是蓝色
  coord_cartesian(clip = "off")+  # 允许超出绘图区域
  theme(legend.position = "right")  # 确保图例仍然位于右侧


w i# quickcor(data2, type = "upper") + geom_square() + 
#   add_link(mantel, mapping = aes(colour = p.value, size = r),
#            diag.label = TRUE) +
#   scale_size_manual(values = c(0.5, 1.5, 3)) +
#   scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +  #指定Mantel相关性值的显著性程度与颜色之间的顺序
#   add_diag_label() + remove_axis("x")


# quickcor(data2, type = "lower") + 
#   geom_square() + 
#   add_link(mantel, mapping = aes(colour = p.value, size = r),
#            diag.label = TRUE) +
#   scale_size_manual(values = c(0.5, 1.5, 3)) +
#   scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
#   coord_cartesian(clip = "off")+  # 允许超出绘图区域
#   theme(legend.position = "right")  # 确保图例仍然位于右侧



mantel <- fortify_mantel(data1, data2, 
                         spec.select = list(riskScore = 7:7)) %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.1, 0.3, Inf), 
                 labels = c("<0.2", "0.1-0.3",">0.4"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                       labels = c("<0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE))

quickcor(data2, type = "lower") + 
  geom_square() + 
  add_link(mantel, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", 
                       midpoint = 0.4) + # 负相关是红色，正相关是蓝色
  coord_cartesian(clip = "off")+  # 允许超出绘图区域
  theme(legend.position = "right")  # 确保图例仍然位于右侧
library(export)
graph2ppt(file = "output/Figure/999_Rho_cell-cell.pptx",width=9.5, height=7.5)
