#https://mp.weixin.qq.com/s/t7MR_5mGXEoIqml7lPvaOw
#https://mp.weixin.qq.com/s/wEoAqlRR4uN-2LGq_0X9-Q
rm(list = ls())
# install.packages("TCGAplot_4.0.0.zip",repos=NULL)#安装
# install.packages("fmsb")   #安装fmsb包
# install.packages("forestplot")
# install.packages("psych")
# BiocManager::install("ComplexHeatmap")
# install.packages("tinyarray")


install.packages(c("devtools", "usethis", "roxygen2"))
pkgdir <- "/Users/xiaowentao/Desktop/TCGAplot/TCGAplot/"   # 改成你的实际路径
list.files(pkgdir)

devtools::build(pkgdir)

library(TCGAplot)

#####泛癌分析
pan_boxplot('SERPINE1',palette="lancet",legend="right")+theme(axis.text.x = element_text(angle = 45, hjust = 1))
library(export)
graph2ppt(file = "output/Figure/pan_cancer.ppt",width=13, height=4)
pan_paired_boxplot('SERPINE1',palette="lancet",legend="right")
# graph2ppt(file = "output/Figure/pan_cancer_pair.ppt",width=11, height=4)


graph2ppt(file = "output/Figure/pan_cancer_pair.ppt",width=13, height=4)
gene_TMB_radar('SERPINE1',method = "pearson")
gene_MSI_radar('SERPINE1',method = "pearson")
# 
# 
# ####基因与免疫的关系
# gene_checkpoint_heatmap('OLR1',method="pearson")
# gene_chemokine_heatmap('OLR1',method="pearson")
# gene_receptor_heatmap('OLR1',method="pearson")
# gene_immustimulator_heatmap('OLR1',method="pearson")
# gene_immuinhibitor_heatmap('OLR1',method="pearson")
# gene_immucell_heatmap('OLR1',method="pearson")
# gene_immunescore_heatmap('OLR1',method="pearson")
# gene_immunescore_triangle('OLR1',method="pearson")

####Cox回归分析
pan_forest("SERPINE1")


graph2ppt(file = "output/Figure/pan_cox.ppt",width=8, height=9)

