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
library(biomaRt)

############################################绘制train森林图01
############################################绘制train森林图01
############################################绘制train森林图01
#读取train输入文件             
rt=read.table("output/lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)  

#使用train组构建cox模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#绘制森林图
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 1, 
         refLabel = "reference", 
         noDigits = 2)

library(export)
graph2ppt(file = "output/Figure/01_train_forest.ppt",width=7.6, height=5.3)

############################################绘制coef图02
############################################绘制coef图02
############################################绘制coef图02
##清空数据
rm(list = ls())
#读取train输入文件             
rt=read.table("output//lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)  

#使用train组构建cox模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型相关信息
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

rt <- as.data.frame(outTab)
rt <- rt %>%
  dplyr::select(id,coef)
colnames(rt)[1] <- "gene"
# 创建一个新变量group，根据coef的正负值分为Protective和Risk
rt$group <- ifelse(rt$coef >= 0, "Risk", "Protective")

bar <- rt %>% 
  arrange(coef)

bar$gene <- factor(bar$gene,levels = bar$gene)
bar$coef <- as.numeric(bar$coef)

ggplot(data = bar,aes(x = gene,y = coef,fill = group)) +
  geom_col() +
  xlab('') + ylab('') +
  # 主题
  theme_prism(border = T)

# 绘图
p <- ggplot(data = bar,aes(x = gene,y = coef,fill = group)) +
  geom_col() +
  xlab('Gene') + ylab('LASSO Cox coefficient') +
  # 主题
  theme_prism(border = T) +
  # 填充颜色
  scale_fill_manual(values = c('Risk'= '#c82423','Protective'='#037F77')) +
  # 翻转坐标轴
  coord_flip() + ylim(-0.4,0.4)
p
library(export)
graph2ppt(file = "output/Figure/02_lasso_cox_coef.ppt",width=7.3, height=5.3)

############################################绘制train_KM生存图03
############################################绘制train_KM生存图03
############################################绘制train_KM生存图03
##清空数据
rm(list = ls())

#绘制train生存曲线
rt=read.table("output/riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)
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
 graph2ppt(file = "output/Figure/03_train_km.ppt",width=6.3, height=5.5)

############################################绘制train_热图04
############################################绘制train_热图04
############################################绘制train_热图04
##清空数据
rm(list = ls())
library(pheatmap)
#绘制train风险热图
rt=read.table("output/riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
annotation$type <- factor(annotation$type, levels = c("low", "high"))
type_colors <- c("low" = "#2878b5", "high" = "#c82423")
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border=FALSE,#去掉网格线
         fontsize_row=11,
         fontsize_col=3,
         scale = "row",
         show_colnames = F,
         cellwidth = 1,#调整格子宽度
         cellheight = 30,#调整格子长度
         color = colorRampPalette(c("#2878b5","white","#c82423"))(100),
         annotation_colors = list(type = type_colors))

library(export)
 graph2ppt(file = "output/Figure/04_train_heatmap.ppt",width=8.8, height=5.8)

############################################绘制train_风险图1 05
############################################绘制train_风险图1 05
############################################绘制train_风险图1 05
##清空数据
rm(list = ls())
rt=read.table("output/riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>100]=100
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk Score",
     col=c(rep("#2878b5",lowLength),
           rep("#c82423",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/05_train_risk.ppt",width=7.3, height=5)

############################################绘制train_风险图2生存状态 06
############################################绘制train_风险图2生存状态 06
############################################绘制train_风险图2生存状态 06
##清空数据
rm(list = ls())
#绘制train生存状态图
rt=read.table("output/riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)        
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="#c82423"
color[color==0]="#2878b5"
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/06_train_status.ppt",width=7.3, height=5)

############################################绘制train_roc 07
############################################绘制train_roc 07
############################################绘制train_roc 07
##清空数据
rm(list = ls())
library(timeROC)
# 读取数据
rt = read.table("output/riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)

var="riskScore"

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt[,var], cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green",'blue','red'),lwd=2,bty = 'n')
library(export)
 graph2ppt(file = "output/Figure/07_train_roc.ppt",width=5, height=5)

############################################绘制train_gene_location 08
############################################绘制train_gene_location 08
############################################绘制train_gene_location 08
##清空数据
rm(list = ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(biomaRt)

# # 获取ensembl hsp基因组全信息
# hsp <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                         dataset = "hsapiens_gene_ensembl",
#                         host = "https://www.ensembl.org")
# save(hsp,file = "resource/hsp.Rdata")
# load(file = "resource/hsp.Rdata")
#
# # 提取设置
# attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position")
#
# # 提取
# hsp_info <- getBM(attributes = attributes, mart = hsp)
# head(hsp_info)
# data.frame(table(hsp_info$chromosome_name))
# save(hsp_info,file = "resource/human_gene.Rdata")

rm(list = ls())
load(file = "resource/human_gene.Rdata")
library(survival)
library(survminer)
library(RCircos)
#读取train输入文件
rt=read.table("output/lassoSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)

#使用train组构建cox模型
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型相关信息
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
index <- outTab[,"id"]

# #将没有hgnc_symbol名称的snoRNA重新命名
# hsp_info[22477, 2] <- "ACA24"
# hsp_info[59382, 2] <- "ACA61"
# hsp_info[28646, 2] <- "U3"
# hsp_info[63527, 2] <- "SNORD19B"

# index <- c("MRPL10","POLR2C","PSMC1","TRIM37")

data <- subset(hsp_info,hgnc_symbol %in% index)
#减去重名
# data <- data[-2,]
# data <- data[-1,] 减去重名snoRNA

data <- data[,-1]
colnames(data) <- c("Gene","Chromosome","chromStart","chromEnd")
data <- dplyr::select(data,Chromosome,chromStart,chromEnd,Gene)
chr <- "chr"
data$Chromosome <- paste(chr,data$Chromosome)
data$Chromosome <- gsub(" ","",data$Chromosome)
data$Chromosome <- as.factor(data$Chromosome)
data$Gene <- as.factor(data$Gene)
geneInfor <- data

data(UCSC.HG19.Human.CytoBandIdeogram);# 导入内置人类染色体数据
RCircos.Set.Core.Components(cyto.info=UCSC.HG19.Human.CytoBandIdeogram,
                            chr.exclude=NULL, tracks.inside=10, tracks.outside=0);

RCircos.List.Plot.Parameters()# 列出所有绘图参数

params <- RCircos.Get.Plot.Parameters()
params$text.size <- 1.1  # 增大字体大小，默认值为0.8，可以根据需要调整
RCircos.Reset.Plot.Parameters(params)

RCircos.Set.Plot.Area();# height和width指定生成图片的长和宽，compress指定生成的图片是否需要压缩

RCircos.Chromosome.Ideogram.Plot();# 绘制染色体图形

RCircos.Gene.Connector.Plot(geneInfor,
                            track.num=1,
                            side='in');# 绘图

RCircos.Gene.Name.Plot(geneInfor,
                       name.col=4,
                       track.num=2,
                       side='in')# 绘图
library(export)
#graph2ppt(file = "output/Figure/08_train_gene_location.ppt",width=5, height=5)
graph2jpg(file = "output/Figure/08_train_gene_location.jpg",width=5, height=5)

############################################绘制test_KM生存图09
############################################绘制test_KM生存图09
############################################绘制test_KM生存图09
##清空数据
rm(list = ls())

#绘制test生存曲线
rt=read.table("output/riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)

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
 graph2ppt(file = "output/Figure/09_test_km.ppt",width=6.3, height=5.5)

############################################绘制test_热图10
############################################绘制test_热图10
############################################绘制test_热图10
##清空数据
rm(list = ls())
library(pheatmap)
#绘制test风险热图
rt=read.table("output/riskTest.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
annotation$type <- factor(annotation$type, levels = c("low", "high"))
type_colors <- c("low" = "#2878b5", "high" = "#c82423")
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation,
         cluster_rows = FALSE,#行聚类
         cluster_cols = FALSE,
         border=FALSE,#去掉网格线
         fontsize_row=8,
         fontsize_col=3,
         scale = "row",
         show_colnames = F,
         cellwidth = 1,#调整格子宽度
         cellheight = 40,#调整格子长度
         color = colorRampPalette(c("#2878b5","white","#c82423"))(100),
         annotation_colors = list(type = type_colors))

library(export)
 graph2ppt(file = "output/Figure/10_test_heatmap.ppt",width=8.8, height=5.8)

############################################绘制test_风险图1 11
############################################绘制test_风险图1 11
############################################绘制test_风险图1 11
##清空数据
rm(list = ls())
rt=read.table("output/riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>100]=100
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk Score",
     col=c(rep("#2878b5",lowLength),
           rep("#c82423",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/11_test_risk.ppt",width=7.3, height=5)

############################################绘制test_风险图2生存状态 12
############################################绘制test_风险图2生存状态 12
############################################绘制test_风险图2生存状态 12
##清空数据
rm(list = ls())
#绘制train生存状态图
rt=read.table("output/riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)        
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="#c82423"
color[color==0]="#2878b5"
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/12_test_status.ppt",width=7.3, height=5)

############################################绘制test_roc 13
############################################绘制test_roc 13
############################################绘制test_roc 13
##清空数据
rm(list = ls())
library(timeROC)
# 读取数据
rt = read.table("output/riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)

var="riskScore"

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt[,var], cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green",'blue','red'),lwd=2,bty = 'n')
library(export)
 graph2ppt(file = "output/Figure/13_test_roc.ppt",width=5, height=5)


############################################绘制entire_KM生存图14
############################################绘制entire_KM生存图14
############################################绘制entire_KM生存图14
##清空数据
rm(list = ls())

#绘制test生存曲线
rt=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)

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
 graph2ppt(file = "output/Figure/14_entire_km.ppt",width=6.3, height=5.5)

############################################绘制entire_热图15
############################################绘制entire_热图15
############################################绘制entire_热图15
##清空数据
rm(list = ls())
library(pheatmap)
#绘制entire风险热图
rt=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
annotation$type <- factor(annotation$type, levels = c("low", "high"))
type_colors <- c("low" = "#2878b5", "high" = "#c82423")
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation,
         cluster_rows = FALSE,#行聚类
         cluster_cols = FALSE,
         border=FALSE,#去掉网格线
         fontsize_row=8,
         fontsize_col=3,
         scale = "row",
         show_colnames = F,
         cellwidth = 0.5,#调整格子宽度###特别注意entire组是test/train的1/2，由1改为0.5，不然太长了报错
         cellheight = 40,#调整格子长度
         color = colorRampPalette(c("#2878b5","white","#c82423"))(100),
         annotation_colors = list(type = type_colors))

library(export)
 graph2ppt(file = "output/Figure/15_entire_heatmap.ppt",width=8.8, height=5.8)
# 
# ##图片元素太多PPT打不开
# pdf(file="output/Figure/15_entire_heatmap.pdf",
#     width = 8.8, #图片的宽度
#     height = 5.8, #图片的高度
# )
# dev.off()
############################################绘制entire_风险图1 16
############################################绘制entire_风险图1 16
############################################绘制entire_风险图1 16
##清空数据
rm(list = ls())
rt=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>100]=100
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk Score",
     col=c(rep("#2878b5",lowLength),
           rep("#c82423",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/16_entire_risk.ppt",width=7.3, height=5)

############################################绘制entire_风险图2生存状态 17
############################################绘制entire_风险图2生存状态 17
############################################绘制entire_风险图2生存状态 17
##清空数据
rm(list = ls())
#绘制entire生存状态图
rt=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)        
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="#c82423"
color[color==0]="#2878b5"
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/17_entire_status.ppt",width=7.3, height=5)

############################################绘制entire_roc 18
############################################绘制entire_roc 18
############################################绘制entire_roc 18
##清空数据
rm(list = ls())
library(timeROC)
# 读取数据
rt = read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)

var="riskScore"

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt[,var], cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green",'blue','red'),lwd=2,bty = 'n')
library(export)
 graph2ppt(file = "output/Figure/18_entire_roc.ppt",width=5, height=5)





############################################绘制entire_KM生存图14_out2
############################################绘制entire_KM生存图14_out2
############################################绘制entire_KM生存图14_out2
##清空数据
rm(list = ls())

#绘制test生存曲线
rt=read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1)

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
 graph2ppt(file = "output/Figure/14_entire_km_out2.ppt",width=6.3, height=5.5)

############################################绘制entire_热图15
############################################绘制entire_热图15
############################################绘制entire_热图15
##清空数据
rm(list = ls())
library(pheatmap)
#绘制entire风险热图
rt=read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
annotation$type <- factor(annotation$type, levels = c("low", "high"))
type_colors <- c("low" = "#2878b5", "high" = "#c82423")
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation,
         cluster_rows = FALSE,#行聚类
         cluster_cols = FALSE,
         border=FALSE,#去掉网格线
         fontsize_row=8,
         fontsize_col=3,
         scale = "row",
         show_colnames = F,
         cellwidth = 2.56,#调整格子宽度###特别注意entire组是test/train的1/2，由1改为0.5，不然太长了报错
         cellheight = 40,#调整格子长度
         color = colorRampPalette(c("#2878b5","white","#c82423"))(100),
         annotation_colors = list(type = type_colors))

library(export)
 graph2ppt(file = "output/Figure/15_entire_heatmap_out2.ppt",width=8.8, height=5.8)
# 
# ##图片元素太多PPT打不开
# pdf(file="output/Figure/15_entire_heatmap.pdf",
#     width = 8.8, #图片的宽度
#     height = 5.8, #图片的高度
# )
# dev.off()
############################################绘制entire_风险图1 16_out2
############################################绘制entire_风险图1 16_out2
############################################绘制entire_风险图1 16_out2
##清空数据
rm(list = ls())
rt=read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>100]=100
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk Score",
     col=c(rep("#2878b5",lowLength),
           rep("#c82423",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/16_entire_risk_out2.ppt",width=7.3, height=5)

############################################绘制entire_风险图2生存状态 17
############################################绘制entire_风险图2生存状态 17
############################################绘制entire_风险图2生存状态 17
##清空数据
rm(list = ls())
#绘制entire生存状态图
rt=read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1)        
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="#c82423"
color[color==0]="#2878b5"
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/17_entire_status_out2.ppt",width=7.3, height=5)

############################################绘制entire_roc 18
############################################绘制entire_roc 18
############################################绘制entire_roc 18
##清空数据
rm(list = ls())
library(timeROC)
# 读取数据
rt = read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1)

var="riskScore"

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt[,var], cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green",'blue','red'),lwd=2,bty = 'n')
library(export)
 graph2ppt(file = "output/Figure/18_entire_roc_out2.ppt",width=5, height=5)





############################################绘制entire_KM生存图14_out2_out3
############################################绘制entire_KM生存图14_out2_out3
############################################绘制entire_KM生存图14_out2_out3
##清空数据
rm(list = ls())

#绘制test生存曲线
rt=read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1)

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
graph2ppt(file = "output/Figure/14_entire_km_out3.ppt",width=6.3, height=5.5)

############################################绘制entire_热图15
############################################绘制entire_热图15
############################################绘制entire_热图15
##清空数据
rm(list = ls())
library(pheatmap)
#绘制entire风险热图
rt=read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1) 
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
annotation$type <- factor(annotation$type, levels = c("low", "high"))
type_colors <- c("low" = "#2878b5", "high" = "#c82423")
rownames(annotation)=rownames(rt)
pheatmap(rt1, 
         annotation=annotation,
         cluster_rows = FALSE,#行聚类
         cluster_cols = FALSE,
         border=FALSE,#去掉网格线
         fontsize_row=8,
         fontsize_col=3,
         scale = "row",
         show_colnames = F,
         cellwidth = 0.918,#调整格子宽度###特别注意entire组是test/train的1/2，由1改为0.5，不然太长了报错
         cellheight = 40,#调整格子长度
         color = colorRampPalette(c("#2878b5","white","#c82423"))(100),
         annotation_colors = list(type = type_colors))

library(export)
 graph2ppt(file = "output/Figure/15_entire_heatmap_out3.ppt",width=8.8, height=5.8)
# 
# ##图片元素太多PPT打不开
# pdf(file="output/Figure/15_entire_heatmap.pdf",
#     width = 8.8, #图片的宽度
#     height = 5.8, #图片的高度
# )
# dev.off()
############################################绘制entire_风险图1 16_out2
############################################绘制entire_风险图1 16_out2
############################################绘制entire_风险图1 16_out2
##清空数据
rm(list = ls())
rt=read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1)     
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>100]=100
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk Score",
     col=c(rep("#2878b5",lowLength),
           rep("#c82423",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/16_entire_risk_out3.ppt",width=7.3, height=5)

############################################绘制entire_风险图2生存状态 17
############################################绘制entire_风险图2生存状态 17
############################################绘制entire_风险图2生存状态 17
##清空数据
rm(list = ls())
#绘制entire生存状态图
rt=read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1)        
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$fustat)
color[color==1]="#c82423"
color[color==0]="#2878b5"
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
       pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))

library(export)
 graph2ppt(file = "output/Figure/17_entire_status_out3.ppt",width=7.3, height=5)

############################################绘制entire_roc 18
############################################绘制entire_roc 18
############################################绘制entire_roc 18
##清空数据
rm(list = ls())
library(timeROC)
# 读取数据
rt = read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1)

var="riskScore"

ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
               marker=rt[,var], cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)

plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green",'blue','red'),lwd=2,bty = 'n')
library(export)
 graph2ppt(file = "output/Figure/18_entire_roc_out3.ppt",width=5, height=5)


# ############################################绘制entire_KM生存图14_out2_out3_out4
# ############################################绘制entire_KM生存图14_out2_out3_out4
# ############################################绘制entire_KM生存图14_out2_out3_out4
# ##清空数据
# rm(list = ls())
# 
# #绘制test生存曲线
# rt=read.table("output/riskoutside4.txt",header=T,sep="\t",check.names=F,row.names=1)
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
#                         conf.int = TRUE,#添加置信区间
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
# # graph2ppt(file = "output/Figure/14_entire_km_out4.ppt",width=6.3, height=5.5)
# 
# ############################################绘制entire_热图15
# ############################################绘制entire_热图15
# ############################################绘制entire_热图15
# ##清空数据
# rm(list = ls())
# library(pheatmap)
# #绘制entire风险热图
# rt=read.table("output/riskoutside4.txt",header=T,sep="\t",check.names=F,row.names=1) 
# rt=rt[order(rt$riskScore),]
# rt1=rt[c(3:(ncol(rt)-2))]
# rt1=t(rt1)
# annotation=data.frame(type=rt[,ncol(rt)])
# annotation$type <- factor(annotation$type, levels = c("low", "high"))
# type_colors <- c("low" = "#2878b5", "high" = "#c82423")
# rownames(annotation)=rownames(rt)
# pheatmap(rt1, 
#          annotation=annotation,
#          cluster_rows = FALSE,#行聚类
#          cluster_cols = FALSE,
#          border=FALSE,#去掉网格线
#          fontsize_row=8,
#          fontsize_col=3,
#          scale = "row",
#          show_colnames = F,
#          cellwidth = 2,#调整格子宽度###特别注意entire组是test/train的1/2，由1改为0.5，不然太长了报错
#          cellheight = 20,#调整格子长度
#          color = colorRampPalette(c("#2878b5","white","#c82423"))(100),
#          annotation_colors = list(type = type_colors))
# 
# library(export)
# # graph2ppt(file = "output/Figure/15_entire_heatmap_out4.ppt",width=8.8, height=5.8)
# # 
# # ##图片元素太多PPT打不开
# # pdf(file="output/Figure/15_entire_heatmap.pdf",
# #     width = 8.8, #图片的宽度
# #     height = 5.8, #图片的高度
# # )
# # dev.off()
# ############################################绘制entire_风险图1 16_out2
# ############################################绘制entire_风险图1 16_out2
# ############################################绘制entire_风险图1 16_out2
# ##清空数据
# rm(list = ls())
# rt=read.table("output/riskoutside4.txt",header=T,sep="\t",check.names=F,row.names=1)     
# rt=rt[order(rt$riskScore),]
# riskClass=rt[,"risk"]
# lowLength=length(riskClass[riskClass=="low"])
# highLength=length(riskClass[riskClass=="high"])
# line=rt[,"riskScore"]
# line[line>100]=100
# plot(line,
#      type="p",
#      pch=20,
#      xlab="Patients (increasing risk socre)",
#      ylab="Risk Score",
#      col=c(rep("#2878b5",lowLength),
#            rep("#c82423",highLength)))
# trainMedianScore=median(rt$riskScore)
# abline(h=trainMedianScore,v=lowLength,lty=2)
# legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
#        pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))
# 
# library(export)
# # graph2ppt(file = "output/Figure/16_entire_risk_out4.ppt",width=7.3, height=5)
# 
# ############################################绘制entire_风险图2生存状态 17
# ############################################绘制entire_风险图2生存状态 17
# ############################################绘制entire_风险图2生存状态 17
# ##清空数据
# rm(list = ls())
# #绘制entire生存状态图
# rt=read.table("output/riskoutside4.txt",header=T,sep="\t",check.names=F,row.names=1)        
# rt=rt[order(rt$riskScore),]
# riskClass=rt[,"risk"]
# lowLength=length(riskClass[riskClass=="low"])
# highLength=length(riskClass[riskClass=="high"])
# color=as.vector(rt$fustat)
# color[color==1]="#c82423"
# color[color==0]="#2878b5"
# plot(rt$futime,
#      pch=19,
#      xlab="Patients (increasing risk socre)",
#      ylab="Survival time (years)",
#      col=color)
# abline(v=lowLength,lty=2)
# legend("topleft", legend=c("Low Risk", "High Risk"), bty="n", 
#        pch=20, pt.cex=1.5, col=c("#2878b5", "#c82423"))
# 
# library(export)
# # graph2ppt(file = "output/Figure/17_entire_status_out4.ppt",width=7.3, height=5)
# 
# ############################################绘制entire_roc 18
# ############################################绘制entire_roc 18
# ############################################绘制entire_roc 18
# ##清空数据
# rm(list = ls())
# library(timeROC)
# # 读取数据
# rt = read.table("output/riskoutside4.txt",header=T,sep="\t",check.names=F,row.names=1)
# 
# var="riskScore"
# 
# ROC_rt=timeROC(T=rt$futime, delta=rt$fustat,
#                marker=rt[,var], cause=1,
#                weighting='aalen',
#                times=c(1,2,3), ROC=TRUE)
# 
# plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
# plot(ROC_rt,time=2,col='blue',add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=3,col='red',add=TRUE,title=FALSE,lwd=2)
# legend('bottomright',
#        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
#          paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
#          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
#        col=c("green",'blue','red'),lwd=2,bty = 'n')
# library(export)
# # graph2ppt(file = "output/Figure/18_entire_roc_out4.ppt",width=5, height=5)



############################################绘制entire_KM生存图14_out2_out3_out4_out5
############################################绘制entire_KM生存图14_out2_out3_out4_out5
############################################绘制entire_KM生存图14_out2_out3_out4_out5
##清空数据
rm(list = ls())

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
                        conf.int = TRUE,#添加置信区间
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

library(export)
 graph2ppt(file = "output/Figure/14_entire_km_out5.ppt",width=6.3, height=5.5)

# ############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6
# ############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6
# ############################################绘制entire_KM生存图14_out2_out3_out4_out5_out6
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
#                         conf.int = TRUE,#添加置信区间
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
# # graph2ppt(file = "output/Figure/14_entire_km_out6.ppt",width=6.3, height=5.5)

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
                        conf.int = TRUE,#添加置信区间
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

library(export)
 graph2ppt(file = "output/Figure/14_entire_km_out7.ppt",width=6.3, height=5.5)

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
#                         conf.int = TRUE,#添加置信区间
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
#  graph2ppt(file = "output/Figure/14_entire_km_out8.ppt",width=6.3, height=5.5)
# 
# 
# # library(TCGAplot)
# # pan_boxplot("PYGL",palette="lancet",legend="right",method="wilcox.test")
# # pan_forest("PYGL",adjust=F)
