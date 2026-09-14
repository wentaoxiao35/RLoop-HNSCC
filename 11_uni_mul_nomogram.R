###########################################绘制riskscore_uni_multivariate 21
###########################################绘制riskscore_uni_multivariate 21
###########################################绘制riskscore_uni_multivariate 21
##清空数据
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
rm(list = ls())
risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)     
clinical=read.table("resource/TCGA.HNSC.CLI.txt",header=T,sep="\t",check.names=F,row.names=1) 
index <- rownames(risk) 
clinical_new <- clinical[index,]
clinical_new <- clinical_new %>% 
  rownames_to_column("TCGA_ID")
risk <- risk %>% 
  rownames_to_column("TCGA_ID")
clinical_final <- merge(risk,clinical_new,by="TCGA_ID")
clinical <- clinical_final %>% 
  dplyr::select(TCGA_ID,futime,fustat,gender,age,
                pathologic_T,pathologic_N,
                pathologic_stage,lymphovascular_invasion,perineural_invasion,riskScore) %>% 
  column_to_rownames("TCGA_ID")


clinical <- clinical[!(apply(clinical, 1, function(x) any(x == "unknown"))), ]

clinical$gender[clinical$gender == "female"] <- 0
clinical$gender[clinical$gender == "male"] <- 1
clinical$pathologic_T[clinical$pathologic_T == "T1"] <- 1
clinical$pathologic_T[clinical$pathologic_T == "T2"] <- 1
clinical$pathologic_T[clinical$pathologic_T == "T3"] <- 2
clinical$pathologic_T[clinical$pathologic_T == "T4"] <- 2
clinical$pathologic_N[clinical$pathologic_N == "N0"] <- 1
clinical$pathologic_N[clinical$pathologic_N == "N1"] <- 1
clinical$pathologic_N[clinical$pathologic_N == "N2"] <- 2
clinical$pathologic_N[clinical$pathologic_N == "N3"] <- 2
# clinical$clinical_M[clinical$clinical_M == "M0"] <- 0
# clinical$clinical_M[clinical$clinical_M == "M1"] <- 1
clinical$pathologic_stage[clinical$pathologic_stage == "Stage I"] <- 1
clinical$pathologic_stage[clinical$pathologic_stage == "Stage II"] <- 2
clinical$pathologic_stage[clinical$pathologic_stage == "Stage III"] <- 3
clinical$pathologic_stage[clinical$pathologic_stage == "Stage IV"] <- 4
clinical$lymphovascular_invasion[clinical$lymphovascular_invasion == "Negative"] <- 0
clinical$lymphovascular_invasion[clinical$lymphovascular_invasion == "Positive"] <- 1
clinical$perineural_invasion[clinical$perineural_invasion == "Negative"] <- 0
clinical$perineural_invasion[clinical$perineural_invasion == "Positive"] <- 1
clinical[] <- lapply(clinical, as.numeric)

rt <- clinical 

univariatedata <- rt
uniTab=data.frame()

for(i in colnames(univariatedata[,3:ncol(univariatedata)])){
  cox <- coxph(Surv(futime, fustat) ~ univariatedata[,i], data = univariatedata)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file = "output/riskoutside_uniCox.txt",sep="\t",row.names=F,quote=F)
uniTab=uniTab[as.numeric(as.character(uniTab[,"pvalue"])) <0.05,]

multivariatedata=univariatedata[,c("futime","fustat",as.vector(uniTab[,"id"]))]
multiCox=coxph(Surv(futime, fustat) ~ ., data = multivariatedata)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file = "output/riskoutside_multiCox.txt",sep="\t",row.names=F,quote=F)

bioForest=function(coxFile=null,forestFile=null,height=null,forestCol=null){
  univariatedata <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(univariatedata)
  hr <- sprintf("%.3f",univariatedata$"HR")
  hrLow  <- sprintf("%.3f",univariatedata$"HR.95L")
  hrHigh <- sprintf("%.3f",univariatedata$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(univariatedata$pvalue<0.001, "<0.001", sprintf("%.3f", univariatedata$pvalue))
  n <- nrow(univariatedata)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  ####
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  ####
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
}
####
bioForest(coxFile="output/riskoutside_uniCox.txt",forestFile="output/uniForest.pdf", forestCol=c("green","green"))
library(export)
graph2ppt(file = "output/Figure/21_univariate.ppt",width=8.5, height=4.5)
####
bioForest(coxFile="output/riskoutside_multiCox.txt",forestFile="output/multiForest.pdf",forestCol=c("red","red"))
library(export)
graph2ppt(file = "output/Figure/22_multivariate.ppt",width=8.5, height=4.5)


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
library(pec)

############################################绘制nomogram
############################################绘制nomogram
############################################绘制nomogram
##清空数据
rm(list = ls())
risk=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)     
clinical=read.table("resource/TCGA.HNSC.CLI.txt",header=T,sep="\t",check.names=F,row.names=1) 
index <- rownames(risk) 
clinical_new <- clinical[index,]
clinical_new <- clinical_new %>% 
  rownames_to_column("TCGA_ID")
risk <- risk %>% 
  rownames_to_column("TCGA_ID")
clinical_final <- merge(risk,clinical_new,by="TCGA_ID")
clinical <- clinical_final %>% 
  dplyr::select(TCGA_ID,futime,fustat,gender,age,
                pathologic_T,pathologic_N,
                pathologic_stage,lymphovascular_invasion,perineural_invasion,riskScore) %>% 
  column_to_rownames("TCGA_ID")


rows_with_unknown <- apply(clinical, 1, function(x) any(x == "unknown"))

# 过滤掉包含"unknown"的行
clinical <- clinical[!rows_with_unknown, ]

# ##因子化
# clinical$gender <- factor(clinical$gender,levels = c("female","male"),order=T)
clinical$pathologic_T <- factor(clinical$pathologic_T,levels = c("T1","T2","T3","T4"),order=T)
clinical$pathologic_N <- factor(clinical$pathologic_N,levels = c("N0","N1","N2","N3"),order=T)
# clinical$pathologic_M <- factor(clinical$pathologic_M,levels = c("M0","M1"),order=T)
clinical$pathologic_stage <- factor(clinical$pathologic_stage,levels = c("Stage I","Stage II","Stage III","Stage IV"),order=T)
# clinical$lymphatic_invasion <- factor(clinical$lymphatic_invasion,levels = c("NO","YES"),order=T)
# clinical$venous_invasion <- factor(clinical$venous_invasion,levels = c("NO","YES"),order=T)
##因子化
# clinical$gender <- factor(clinical$gender,levels = c("female","male"),order=T)
# clinical$clinical_T <- factor(clinical$clinical_T,levels = c("T1","T2","T3","T4"),order=T)
# clinical$clinical_N <- factor(clinical$clinical_N,levels = c("N0","N1","N2","N3"),order=T)
# clinical$clinical_M <- factor(clinical$clinical_M,levels = c("M0","M1"),order=T)
# clinical$clinical_stage <- factor(clinical$clinical_stage,levels = c("I","II","III","IV"),order=T)
clinical$lymphovascular_invasion <- factor(clinical$lymphovascular_invasion,levels = c("Negative","Positive"),order=T)
clinical$perineural_invasion <- factor(clinical$perineural_invasion,levels = c("Negative","Positive"),order=T)
# clinical$risk <- factor(clinical$risk,levels = c("low","high"),order=T)

##因子化也不能算
# ####如果要输出的nomogram变量不是数字的话，把下面的映射杠掉！！！########
# # 映射表，先不映射？
# map_pathologic_T <- setNames(1:4, c("T1", "T2", "T3", "T4"))
# map_pathologic_N <- setNames(1:4, c("N0", "N1", "N2","N3"))
# # map_clinical_M <- setNames(1:2, c("M0", "M1"))
# map_pathologic_stage <- setNames(1:4, c("Stage I","Stage II","Stage III","Stage IV"))
# map_invasion <- setNames(0:1, c("Negative", "Positive"))
# 
# # 应用映射
# clinical$pathologic_T <- map_pathologic_T[clinical$pathologic_T]
# clinical$pathologic_N <- map_pathologic_N[clinical$pathologic_N]
# # clinical$clinical_M <- map_clinical_M[clinical$clinical_M]
# clinical$pathologic_stage <- map_pathologic_stage[clinical$pathologic_stage]
# clinical$lymphovascular_invasion <- map_invasion[clinical$lymphovascular_invasion]
# clinical$perineural_invasion <- map_invasion[clinical$perineural_invasion]
# 
# str(clinical)
# clinical$fustat <- as.numeric(clinical$fustat)
# clinical$age <- as.numeric(clinical$age)
# # clinical$clinical_stage <- as.numeric(clinical$clinical_stage)
# clinical$lymphovascular_invasion <- as.numeric(clinical$lymphovascular_invasion)
# clinical$perineural_invasion <- as.numeric(clinical$perineural_invasion)
# clinical$pathologic_T <- as.numeric(clinical$pathologic_T)
# clinical$pathologic_N <- as.numeric(clinical$pathologic_N)
# # clinical$clinical_M <- as.numeric(clinical$clinical_M)
# clinical$pathologic_stage <- as.numeric(clinical$pathologic_stage)
# ####
# clinical$fustat <- as.numeric(clinical$fustat)
# clinical$futime <- as.numeric(clinical$futime)
# clinical$age <- as.numeric(clinical$age)
# clinical$riskScore <- as.numeric(clinical$riskScore)

####如果要输出的nomogram变量不是数字的话，把上面的杠掉#########










str(clinical)

## 用cph()函数，基于cox比例风险模型构建nomogram
ddist <- datadist(clinical) #将数据打包
options(datadist="ddist")
options(contrasts = c("contr.treatment", "contr.treatment"))
units(clinical$futime) <- "year" #定义时间的单位

cox_cph <- cph(Surv(futime, fustat) ~ pathologic_T+pathologic_N+pathologic_stage+
                 lymphovascular_invasion+perineural_invasion+riskScore,
               x=T, y=T, surv=T,
               data=clinical)
#cox_cph <- step(cox_cph,direction = "both")
#?step()

surv <- Survival(cox_cph)

cox_cph_nomogram<- nomogram(cox_cph,fun=list(function(x) surv(1, x), function(x) surv(3, x),function(x) surv(5, x)),##调整相应时间即可
                            lp= F,
                            funlabel=c('1-Year Survival','3-Year survival','5-Year Survival'),
                            maxscale=100,
                            fun.at=c('0.9','0.7','0.5','0.3','0.1'))


par(mar = c(5, 12, 3, 2))   # 下、左、上、右边距

###label.every = 2
plot(
  cox_cph_nomogram,
  lplabel = "Linear Predictor",
  xfrac = 0.4, # 左侧标签距离坐标轴的距离
  tcl = -0.2, # 刻度长短和方向 
  lmgp = 0.2,# 坐标轴标签距离坐标轴远近
  points.label = "Points",
  total.points.label = "Total Points",
  cap.labels = FALSE,
  cex.var = 0.8,# 左侧标签字体大小
  cex.axis = 0.9,# 坐标轴字体大小
  col.grid = gray(c(0.85, 0.95)),
  label.every = 2
)


# 

library(export)
graph2ppt(file = "output/Figure/Nomogram_tcga.ppt",width=8.3, height=6.7)


##计算nomogram分数，用于评价
Nomogram<-predict(cox_cph,clinical,type="lp")#!!!type="lp",是他没错
##加入到数据中的最后一列
clinical$Nomogram <- Nomogram

# ##计算nomogram分数，用于评价
# univariatedata$Nomogram <- univariatedata$age * ifelse(univariatedata$age == 1, 0, 33.49963) +
#   univariatedata$histological_type * ifelse(univariatedata$histological_type == 1, 0, 72.70534) +
#   univariatedata$stage * ifelse(univariatedata$stage == 1, 0, ifelse(univariatedata$stage == 2, 33.33333, ifelse(univariatedata$stage == 3, 66.66667, 100))) +
#   univariatedata$POSTN * ifelse(univariatedata$POSTN == 1, 0, 69.07288)


## 1.C-index评价模型
# ①  cph()函数构建的模型，其C-index计算
validate(cox_cph, method="boot", B=1000, dxy=T)   #重复模拟1000次
rcorrcens(Surv(futime,fustat) ~ predict(cox_cph), data = clinical)
## 该模型中，C-index为 1-0.191=0.809 ;  或者Dxy =-0.618, C-Index=|-0.418/2|+ 0.5 =0.809.

##################################################### 2.Calibration curve验证模型
library(rms)
dd <- datadist(clinical)
options(datadist="dd")
####
time=1
f_1 <- cph(Surv(futime, fustat) ~ pathologic_T+pathologic_N+pathologic_stage+
             lymphovascular_invasion+perineural_invasion+riskScore,
           x=T, y=T, surv=T, data=clinical, time.inc=time)
P_1 <- calibrate(f_1, cmethod="KM", method="boot", u=time, m=92, B=1000)
plot(P_1,
     add = F,
     subtitles = F,
     cex.subtitles = 0.8,
     lwd = 2,
     lty = 1,
     errbar.col = "red",
     xlim = c(0,1),
     ylim = c(0,1),
     xlab="Nomogram-Predicted Probability of 1, 3, 5-Year OS",
     ylab="Actual 1, 3, 5-Year OS(proportion)",
     col="red",
     sub=F)
abline(0, 1, col = "black")
####
time=3
f_3 <- cph(Surv(futime, fustat) ~ pathologic_T+pathologic_N+pathologic_stage+
             lymphovascular_invasion+perineural_invasion+riskScore,
           x=T, y=T, surv=T, data=clinical, time.inc=time)
P_3 <- calibrate(f_3, cmethod="KM", method="boot", u=time, m=92, B=1000)
plot(P_3,
     add = T,
     subtitles = F,
     cex.subtitles = 0.8,
     lwd = 2,
     lty = 1,
     errbar.col = "orange",
     xlim = c(0,1),
     ylim = c(0,1),
     xlab="Nomogram-Predicted Probability of 1, 3, 5-Year OS",
     ylab="Actual 1, 3, 5-Year OS(proportion)",
     col="orange",
     sub=F)
abline(0, 1, col = "black")
####
time=5
f_5 <- cph(Surv(futime, fustat) ~ pathologic_T+pathologic_N+pathologic_stage+
             lymphovascular_invasion+perineural_invasion+riskScore,
           x=T, y=T, surv=T, data=clinical, time.inc=time)
P_5 <- calibrate(f_5, cmethod="KM", method="boot", u=time, m=92, B=1000)
plot(P_5,
     add = T,
     subtitles = F,
     cex.subtitles = 0.8,
     lwd = 2,
     lty = 1,
     errbar.col = "blue",
     xlim = c(0,1),
     ylim = c(0,1),
     xlab="Nomogram-Predicted Probability of 1, 3, 5-Year OS",
     ylab="Actual 1, 3, 5-Year OS(proportion)",
     col="blue",
     sub=F)
abline(0, 1, col = "black")
####
####
legend("bottomright",legend = c("1 year","3 year","5 year"), 
       col = c("red","orange","blue"),lwd = 2)
abline(0, 1, col = "black")
library(export)
graph2ppt(file = "output/Figure/calibration_tcga.ppt",width=6.7, height=6.6)

##################################################################C-index time
index_data <- clinical %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,futime,fustat,pathologic_T,pathologic_N,pathologic_stage,
                  lymphovascular_invasion,perineural_invasion,riskScore,Nomogram)

colnames(index_data)[2:3] <- c("time","status")

data <- index_data
library(riskRegression)
library(pec)
ddist <- datadist(data)
options(datadist = "ddist")


models = list(
  T_stage=cph(Surv(time,status)~pathologic_T,data=data,x=TRUE,y=TRUE,surv = T),
  N_stage=cph(Surv(time,status)~pathologic_N,data=data,x=TRUE,y=TRUE,surv = T),
  stage=cph(Surv(time,status)~pathologic_stage,data=data,x=TRUE,y=TRUE,surv = T),
  lymphovascular_invasion=cph(Surv(time,status)~lymphovascular_invasion,data=data,x=TRUE,y=TRUE,surv = T),
  perineural_invasion=cph(Surv(time,status)~perineural_invasion,data=data,x=TRUE,y=TRUE,surv = T),
  Riskscore=cph(Surv(time,status)~riskScore,data=data,x=TRUE,y=TRUE,surv = T), 
  Nomogram=cph(Surv(time,status)~Nomogram,data=data,x=TRUE,y=TRUE,surv = T)
)

times <- c(1,2,3,4,5,6,7,8,9,10)
cindex<- cindex(models,
                formula=Surv(time,status)~1,
                eval.times = times,
                data=data)
plot(cindex)


#用ggplot2画，但这次是多个模型,cindex$AppCindex里面是3个向量
cindex$AppCindex

cindex_df <- data.frame(
  Time = times,
  do.call(cbind,cindex$AppCindex)
)
cindex_df


dat = pivot_longer(cindex_df,cols = 2:ncol(cindex_df),
                   names_to = "model",
                   values_to = "cindex")
head(dat)

dat$model <- factor(dat$model, levels = c("Nomogram",
                                          "Riskscore", 
                                          "T_stage",
                                          "N_stage",
                                          "stage",
                                          "lymphovascular_invasion",
                                          "perineural_invasion"))

library(ggplot2)
ggplot(dat, aes(x = Time, y = cindex)) +
  geom_line(aes(color = model),linewidth = 1.5) + 
  scale_color_brewer(palette = "Set1")+
  ylim(0,1)+
  labs(title = "Time-dependent C-index", x = "Time (years)", y = "C-index") + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 20, by = 2))

library(export)
graph2ppt(file = "output/Figure/Nomogram+C_index.ppt",width=7, height=5)


####################################################################ROC
var="Nomogram"

ROC_rt=timeROC(T=clinical$futime, delta=clinical$fustat,
               marker=clinical[,var], cause=1,
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
graph2ppt(file = "output/Figure/nomogram_roc.ppt",width=5, height=5)

#########################################################Nomogarm KM
median_score <- median(clinical$Nomogram)
clinical$Nomogram <- ifelse(clinical$Nomogram <= median_score, "low", "high")
##1.写函数
genes <- "Nomogram"
your.surv <- Surv(clinical$futime, clinical$fustat) 
your.km.plot <- function(genes,data){
  print(genes)
  group <- clinical[,genes] #分组
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
                        xlim = c(0,max(clinical$futime)+1), #展示x轴的范围
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
your.km.plot("Nomogram",data = clinical)
library(export)
graph2ppt(file = "output/Figure/nomogarm_KM.ppt",width=6.3, height=5.5)
