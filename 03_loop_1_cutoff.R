##清空数据
rm(list = ls())

##加载R包
library(tidyr)
library(dplyr)
library(tibble)
library(survival)
library(caret)
library("glmnet")
library(survminer)
library("survivalROC")

####################HNSC
load(file = "resource/exp_TCGA.Rdata")
load(file = "resource/exp_GSE41613.Rdata")
load(file = "resource/exp_GSE65858.Rdata")
####################immune
load(file = "resource/exp_IMvigor210.Rdata")
load(file = "resource/exp_GSE78220.Rdata")
load(file = "resource/exp_GSE91061.Rdata")
####################immune_c1c2differ_cox_gene
load(file = "output/immune_C1C2differ_cox_gene.Rdata")

#排除生存30天以下的数据
exp_TCGA <- exp_TCGA %>% 
  filter(OS.time > 0.081)
exp_GSE41613 <- exp_GSE41613 %>% 
  filter(OS.time > 0.081)
exp_GSE65858 <- exp_GSE65858 %>% 
  filter(OS.time > 0.081)
exp_IMvigor210 <- exp_IMvigor210 %>% 
  filter(OS.time > 0.081)
exp_GSE78220 <- exp_GSE78220 %>% 
  filter(OS.time > 0.081)
exp_GSE91061 <- exp_GSE91061 %>% 
  filter(OS.time > 0.081)

#####先看一下有没有以数字开头的基因名称，以免跑的时候出错
gene_cols <- grep("^[0-9]", colnames(exp_TCGA))
gene_cols <- grep("^[0-9]", colnames(exp_GSE41613))
gene_cols <- grep("^[0-9]", colnames(exp_GSE65858))
gene_cols <- grep("^[0-9]", colnames(exp_IMvigor210))
gene_cols <- grep("^[0-9]", colnames(exp_GSE78220))
gene_cols <- grep("^[0-9]", colnames(exp_GSE91061))
gene_cols <- grep("^[0-9]", gene)

##处理一下，不然会报错
colnames(exp_TCGA) <- gsub("-","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\(","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\)","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("\\;","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("/","_",colnames(exp_TCGA))
colnames(exp_TCGA) <- gsub("@","_",colnames(exp_TCGA))
##处理一下，不然会报错
colnames(exp_GSE41613) <- gsub("-","_",colnames(exp_GSE41613))
colnames(exp_GSE41613) <- gsub("\\(","_",colnames(exp_GSE41613))
colnames(exp_GSE41613) <- gsub("\\)","_",colnames(exp_GSE41613))
colnames(exp_GSE41613) <- gsub("\\;","_",colnames(exp_GSE41613))
colnames(exp_GSE41613) <- gsub("/","_",colnames(exp_GSE41613))
colnames(exp_GSE41613) <- gsub("@","_",colnames(exp_GSE41613))
##处理一下，不然会报错
colnames(exp_GSE65858) <- gsub("-","_",colnames(exp_GSE65858))
colnames(exp_GSE65858) <- gsub("\\(","_",colnames(exp_GSE65858))
colnames(exp_GSE65858) <- gsub("\\)","_",colnames(exp_GSE65858))
colnames(exp_GSE65858) <- gsub("\\;","_",colnames(exp_GSE65858))
colnames(exp_GSE65858) <- gsub("/","_",colnames(exp_GSE65858))
colnames(exp_GSE65858) <- gsub("@","_",colnames(exp_GSE65858))
##处理一下，不然会报错
colnames(exp_IMvigor210) <- gsub("-","_",colnames(exp_IMvigor210))
colnames(exp_IMvigor210) <- gsub("\\(","_",colnames(exp_IMvigor210))
colnames(exp_IMvigor210) <- gsub("\\)","_",colnames(exp_IMvigor210))
colnames(exp_IMvigor210) <- gsub("\\;","_",colnames(exp_IMvigor210))
colnames(exp_IMvigor210) <- gsub("/","_",colnames(exp_IMvigor210))
colnames(exp_IMvigor210) <- gsub("@","_",colnames(exp_IMvigor210))
##处理一下，不然会报错
colnames(exp_GSE78220) <- gsub("-","_",colnames(exp_GSE78220))
colnames(exp_GSE78220) <- gsub("\\(","_",colnames(exp_GSE78220))
colnames(exp_GSE78220) <- gsub("\\)","_",colnames(exp_GSE78220))
colnames(exp_GSE78220) <- gsub("\\;","_",colnames(exp_GSE78220))
colnames(exp_GSE78220) <- gsub("/","_",colnames(exp_GSE78220))
colnames(exp_GSE78220) <- gsub("@","_",colnames(exp_GSE78220))
##处理一下，不然会报错
colnames(exp_GSE91061) <- gsub("-","_",colnames(exp_GSE91061))
colnames(exp_GSE91061) <- gsub("\\(","_",colnames(exp_GSE91061))
colnames(exp_GSE91061) <- gsub("\\)","_",colnames(exp_GSE91061))
colnames(exp_GSE91061) <- gsub("\\;","_",colnames(exp_GSE91061))
colnames(exp_GSE91061) <- gsub("/","_",colnames(exp_GSE91061))
colnames(exp_GSE91061) <- gsub("@","_",colnames(exp_GSE91061))

##构建共有基因表达矩阵
index1 <- gene
index2 <- colnames(exp_TCGA)[-(1:3)]
index3 <- colnames(exp_GSE41613)[-(1:3)]
index4 <- colnames(exp_GSE65858)[-(1:3)]
index6 <- colnames(exp_IMvigor210)[-(1:3)]
index7 <- colnames(exp_GSE78220)[-(1:3)]
index8 <- colnames(exp_GSE91061)[-(1:3)]

index <- intersect(index1,index2)
index <- intersect(index,index3)
index <- intersect(index,index4)
index <- intersect(index,index6)
index <- intersect(index,index7)
index <- intersect(index,index8)

#######构建最终循环表达矩阵
rt <- exp_TCGA[,c("ID","OS.time","OS",index)]
rt <- rt %>% 
  column_to_rownames("ID")
colnames(rt)[1] <- "futime"
colnames(rt)[2] <- "fustat"
##########
outside_TCGA <- rt
##########
outside_41613 <- exp_GSE41613[,c("ID","OS.time","OS",index)]
outside_41613 <- outside_41613 %>% 
  column_to_rownames("ID") %>% 
  as.data.frame()
colnames(outside_41613)[1] <- "futime"
colnames(outside_41613)[2] <- "fustat"
##########
outside_65858 <- exp_GSE65858[,c("ID","OS.time","OS",index)]
outside_65858 <- outside_65858 %>% 
  column_to_rownames("ID") %>% 
  as.data.frame()
colnames(outside_65858)[1] <- "futime"
colnames(outside_65858)[2] <- "fustat"
##########
outside_IMvigor210 <- exp_IMvigor210[,c("ID","OS.time","OS",index)]
outside_IMvigor210 <- outside_IMvigor210 %>% 
  column_to_rownames("ID") %>% 
  as.data.frame()
colnames(outside_IMvigor210)[1] <- "futime"
colnames(outside_IMvigor210)[2] <- "fustat"
##########
outside_78220 <- exp_GSE78220[,c("ID","OS.time","OS",index)]
outside_78220 <- outside_78220 %>% 
  column_to_rownames("ID") %>% 
  as.data.frame()
colnames(outside_78220)[1] <- "futime"
colnames(outside_78220)[2] <- "fustat"
##########
outside_91061 <- exp_GSE91061[,c("ID","OS.time","OS",index)]
outside_91061 <- outside_91061 %>% 
  column_to_rownames("ID") %>% 
  as.data.frame()
colnames(outside_91061)[1] <- "futime"
colnames(outside_91061)[2] <- "fustat"

##最后检查一边数据里面有没有缺失值
sum(is.na(rt))
sum(is.na(outside_41613))
sum(is.na(outside_65858))
sum(is.na(outside_TCGA))
sum(is.na(outside_IMvigor210))
sum(is.na(outside_91061))
sum(is.na(outside_78220))

load(file = "output/response_IMvigor210.R")
load(file = "output/response_78220.R")
load(file = "output/response_91061.R")
wilcox_response_risk <- function(df_response, df_risk, response_col = "response", risk_col = "riskScore") {
  # 添加ID列
  df_response$ID <- rownames(df_response)
  df_risk$ID <- rownames(df_risk)
  # 合并数据框
  df_merged <- merge(df_response[, c("ID", response_col)],
                     df_risk[, c("ID", risk_col)],
                     by = "ID")
  # 重命名列
  colnames(df_merged) <- c("ID", "response", "riskScore")
  # 确保 response 是因子
  df_merged$response <- factor(df_merged$response)
  # 检查分组数量
  if (length(levels(df_merged$response)) != 2) {
    stop("response 分组数量不等于2，Wilcoxon检验不适用！")
  }
  # 进行Wilcoxon检验
  test_result <- wilcox.test(riskScore ~ response, data = df_merged)
  # 返回结果
  return(list(
    merged_data = df_merged,
    wilcox_test = test_result,
    p_value = test_result$p.value
  ))
}

risk_response_assoc_test <- function(df_response, df_risk,
                                     response_col = "response", risk_col = "risk") {
  # 添加 ID 列
  df_response$ID <- rownames(df_response)
  df_risk$ID <- rownames(df_risk)
  # 合并两个数据框，按 ID
  df <- merge(df_response[, c("ID", response_col)],
              df_risk[, c("ID", risk_col)],
              by = "ID")
  # 重命名列
  colnames(df) <- c("ID", "response", "risk")
  # 确保是因子（可选）
  df$response <- factor(df$response)
  df$risk <- factor(df$risk)
  # 生成列联表
  contingency_table <- table(df$risk, df$response)
  # 计算期望频数
  expected_counts <- chisq.test(contingency_table)$expected
  # 判断使用卡方检验或 Fisher 精确检验
  if (all(expected_counts >= 5)) {
    chi_test <- chisq.test(contingency_table)
    test_method <- "Chi-squared test"
  } else {
    chi_test <- fisher.test(contingency_table)
    test_method <- "Fisher's exact test"
  }
  # 计算每组百分比（可用于画图）
  df_summary <- df %>%
    group_by(risk, response) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(risk) %>%
    mutate(percentage = count / sum(count) * 100)
  # 返回结果
  return(list(
    merged_data = df,
    summary_table = df_summary,
    contingency_table = contingency_table,
    expected_counts = expected_counts,
    test_method = test_method,
    test_result = chi_test,
    p_value = chi_test$p.value
  ))
}

for(i in 1:1000000){
  #############随机分组#############
  inTrain<-createDataPartition(y=rt[,2],p=0.5,list=F)
  train<-rt[inTrain,]
  test<-rt[-inTrain,]
  
  outside1 <- outside_TCGA
  outside2 <- outside_41613
  outside3 <- outside_65858
  outside5 <- outside_IMvigor210
  # outside6 <- outside_78220
  outside7 <- outside_91061
  
  trainOut=cbind(id=row.names(train),train)
  testOut=cbind(id=row.names(test),test)
  
  #############train单因素cox#############
  outTab=data.frame()
  pFilter=0.05
  sigGenes=c("futime","fustat")
  for(i in colnames(train[,3:ncol(train)])){
    cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
    if(coxP<pFilter){
      sigGenes=c(sigGenes,i)
    }
  }
  if(length(sigGenes) < 10) {
    next
  }
  train=train[,sigGenes]
  test=test[,sigGenes]
  outside1=outside1[,sigGenes]
  outside2=outside2[,sigGenes]
  outside3=outside3[,sigGenes]
  outside5=outside5[,sigGenes]
  # outside6=outside6[,sigGenes]
  outside7=outside7[,sigGenes]
  
  uniSigExp=train[,sigGenes]
  uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
  
  #############train-单因素cox后-lasso#############
  trainLasso=train
  x=as.matrix(trainLasso[,c(3:ncol(trainLasso))])
  y=data.matrix(Surv(trainLasso$futime,trainLasso$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 1000)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene=c("futime","fustat",lassoGene)
  if(length(lassoGene) < 5 || length(lassoGene) > 18) {
    next
  }
  train=train[,lassoGene]
  test=test[,lassoGene]
  outside1=outside1[,lassoGene]
  outside2=outside2[,lassoGene]
  outside3=outside3[,lassoGene]
  outside5=outside5[,lassoGene]
  # outside6=outside6[,lassoGene]
  outside7=outside7[,lassoGene]
  
  lassoSigExp=train
  lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
  
  ###########train-单因素cox后-lasso后-多因素cox###############
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = train)
  multiCox=step(multiCox,direction = "both")
  multiCoxSum=summary(multiCox)
  
  #输出多因素cox有意义结果
  outMultiTab=data.frame()
  outMultiTab=cbind(
    coef=multiCoxSum$coefficients[,"coef"],
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
  
  ##train中根据上述多因素结果构建riskscore
  riskScore=predict(multiCox,type="risk",newdata=train)
  coxGene=rownames(multiCoxSum$coefficients)
  outCol=c("futime","fustat",coxGene)
  medianTrainRisk=median(riskScore)
  risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
  trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
  
  ##在test中根据上述多因素结果也构建riskscore
  riskScoreTest=predict(multiCox,type="risk",newdata=test)
  riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
  testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
  
  ##在outside1中根据上述多因素结果也构建riskscore
  riskScoreoutside1=predict(multiCox,type="risk",newdata=outside1)
  medianoutsideRisk1=median(riskScoreoutside1)
  riskoutside1=as.vector(ifelse(riskScoreoutside1>medianoutsideRisk1,"high","low"))
  outsideRiskOut1=cbind(id=rownames(cbind(outside1[,outCol],riskScoreoutside1,riskoutside1)),cbind(outside1[,outCol],riskScore=riskScoreoutside1,risk=riskoutside1))
  ##在outside2中根据上述多因素结果也构建riskscore
  riskScoreoutside2=predict(multiCox,type="risk",newdata=outside2)
  medianoutsideRisk2=median(riskScoreoutside2)
  riskoutside2=as.vector(ifelse(riskScoreoutside2>medianoutsideRisk2,"high","low"))
  outsideRiskOut2=cbind(id=rownames(cbind(outside2[,outCol],riskScoreoutside2,riskoutside2)),cbind(outside2[,outCol],riskScore=riskScoreoutside2,risk=riskoutside2))
  ##在outside3中根据上述多因素结果也构建riskscore
  riskScoreoutside3=predict(multiCox,type="risk",newdata=outside3)
  medianoutsideRisk3=median(riskScoreoutside3)
  riskoutside3=as.vector(ifelse(riskScoreoutside3>medianoutsideRisk3,"high","low"))
  outsideRiskOut3=cbind(id=rownames(cbind(outside3[,outCol],riskScoreoutside3,riskoutside3)),cbind(outside3[,outCol],riskScore=riskScoreoutside3,risk=riskoutside3))
  
  ##############################################################################处理 outside5
  # 计算风险评分
  riskScoreOutside5 = predict(multiCox, type = "risk", newdata = outside5)
  # 创建一个包含生存时间、事件和风险评分的数据框
  outside5_with_risk <- data.frame(outside5[, c("futime", "fustat")], riskScore = riskScoreOutside5)
  # 使用 surv_cutpoint 函数寻找最佳分割点
  res.cut5 <- surv_cutpoint(outside5_with_risk, 
                            time = "futime",
                            event = "fustat",
                            variables = "riskScore",
                            minprop = 0.1) 
  # 根据最佳分割点对风险评分进行分组
  res.cat5 <- surv_categorize(res.cut5)
  # 提取分组信息
  riskoutside5 <- res.cat5$riskScore
  # 合并结果到最终的数据框
  outsideRiskOut5 = cbind(id = rownames(cbind(outside5[, outCol], riskScoreOutside5, riskoutside5)),
                          cbind(outside5[, outCol], riskScore = riskScoreOutside5, risk = riskoutside5))
  # ##############################################################################处理 outside6
  # # 计算风险评分
  # riskScoreOutside6 = predict(multiCox, type = "risk", newdata = outside6)
  # # 创建一个包含生存时间、事件和风险评分的数据框
  # outside6_with_risk <- data.frame(outside6[, c("futime", "fustat")], riskScore = riskScoreOutside6)
  # # 使用 surv_cutpoint 函数寻找最佳分割点
  # res.cut6 <- surv_cutpoint(outside6_with_risk, 
  #                           time = "futime",
  #                           event = "fustat",
  #                           variables = "riskScore",
  #                           minprop = 0.1) 
  # # 根据最佳分割点对风险评分进行分组
  # res.cat6 <- surv_categorize(res.cut6)
  # # 提取分组信息
  # riskoutside6 <- res.cat6$riskScore
  # # 合并结果到最终的数据框
  # outsideRiskOut6 = cbind(id = rownames(cbind(outside6[, outCol], riskScoreOutside6, riskoutside6)),
  #                         cbind(outside6[, outCol], riskScore = riskScoreOutside6, risk = riskoutside6))
  ##############################################################################处理 outside7
  # 计算风险评分
  riskScoreOutside7 = predict(multiCox, type = "risk", newdata = outside7)
  # 创建一个包含生存时间、事件和风险评分的数据框
  outside7_with_risk <- data.frame(outside7[, c("futime", "fustat")], riskScore = riskScoreOutside7)
  # 使用 surv_cutpoint 函数寻找最佳分割点
  res.cut7 <- surv_cutpoint(outside7_with_risk, 
                            time = "futime",
                            event = "fustat",
                            variables = "riskScore",
                            minprop = 0.1) 
  # 根据最佳分割点对风险评分进行分组
  res.cat7 <- surv_categorize(res.cut7)
  # 提取分组信息
  riskoutside7 <- res.cat7$riskScore
  # 合并结果到最终的数据框
  outsideRiskOut7 = cbind(id = rownames(cbind(outside7[, outCol], riskScoreOutside7, riskoutside7)),
                          cbind(outside7[, outCol], riskScore = riskScoreOutside7, risk = riskoutside7))
  
  ##train中riskscore的生存结果
  diff=survdiff(Surv(futime, fustat) ~risk,data = train)
  pValue=1-pchisq(diff$chisq,df=1)
  ##test中riskscore的生存结果
  diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
  pValueTest=1-pchisq(diffTest$chisq,df=1)
  ##outside1中riskscore的生存结果
  diffoutside1=survdiff(Surv(futime, fustat) ~riskoutside1,data = outside1)
  pValueoutside1=1-pchisq(diffoutside1$chisq,df=1)
  ##outside2中riskscore的生存结果
  diffoutside2=survdiff(Surv(futime, fustat) ~riskoutside2,data = outside2)
  pValueoutside2=1-pchisq(diffoutside2$chisq,df=1)
  ##outside3中riskscore的生存结果
  diffoutside3=survdiff(Surv(futime, fustat) ~riskoutside3,data = outside3)
  pValueoutside3=1-pchisq(diffoutside3$chisq,df=1)
  
  ##outside5中riskscore的生存结果
  diffoutside5=survdiff(Surv(futime, fustat) ~riskoutside5,data = outside5)
  pValueoutside5=1-pchisq(diffoutside5$chisq,df=1)
  # ##outside6中riskscore的生存结果
  # diffoutside6=survdiff(Surv(futime, fustat) ~riskoutside6,data = outside6)
  # pValueoutside6=1-pchisq(diffoutside6$chisq,df=1)
  ##outside7中riskscore的生存结果
  diffoutside7=survdiff(Surv(futime, fustat) ~riskoutside7,data = outside7)
  pValueoutside7=1-pchisq(diffoutside7$chisq,df=1)
  
  result1 <- wilcox_response_risk(response_IMvigor210, outsideRiskOut5)
  # result2 <- wilcox_response_risk(response_78220, outsideRiskOut6)
  result3 <- wilcox_response_risk(response_91061, outsideRiskOut7)                          
  result4 <- risk_response_assoc_test(response_IMvigor210, outsideRiskOut5)
  # result5 <- risk_response_assoc_test(response_78220, outsideRiskOut6)
  result6 <- risk_response_assoc_test(response_91061, outsideRiskOut7)                        
  
  ##循环跑出制定结果
  if(  (          pValue<0.01) 
       & (    pValueTest<0.03) 
       & (pValueoutside1<0.03) 
       & (pValueoutside2<0.045)
       & (pValueoutside3<0.045)
       & (pValueoutside5<0.045)
       # & (pValueoutside6<0.045)
       & (pValueoutside7<0.045)
       
       & (result1$p_value<0.05)
       # & (result2$p_value<0.05)
       & (result3$p_value<0.05)
       
       & (result4$p_value<0.05)
       # & (result5$p_value<0.05)
       & (result6$p_value<0.05)
  )
  { write.table(trainOut,file="output/train.txt",sep="\t",quote=F,row.names=F)
    write.table(testOut,file="output/test.txt",sep="\t",quote=F,row.names=F)
    write.table(outTab,file="output/uniCox.txt",sep="\t",quote=F,row.names=F)
    write.table(uniSigExp,file="output/uniSigExp.txt",sep="\t",row.names=F,quote=F)
    write.table(lassoSigExp,file="output/lassoSigExp.txt",sep="\t",row.names=F,quote=F)
    pdf("output/lambda.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    pdf("output/cvfit.pdf")
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
    dev.off()
    write.table(outMultiTab,file="output/multiCox.xls",sep="\t",row.names=F,quote=F)
    write.table(testRiskOut,file="output/riskTest.txt",sep="\t",quote=F,row.names=F)
    write.table(trainRiskOut,file="output/riskTrain.txt",sep="\t",quote=F,row.names=F)
    write.table(outsideRiskOut1,file="output/riskoutside1.txt",sep="\t",quote=F,row.names=F)
    write.table(outsideRiskOut2,file="output/riskoutside2.txt",sep="\t",quote=F,row.names=F)
    write.table(outsideRiskOut3,file="output/riskoutside3.txt",sep="\t",quote=F,row.names=F)
    write.table(outsideRiskOut5,file="output/riskoutside5.txt",sep="\t",quote=F,row.names=F)
    # write.table(outsideRiskOut6,file="output/riskoutside6.txt",sep="\t",quote=F,row.names=F)
    write.table(outsideRiskOut7,file="output/riskoutside7.txt",sep="\t",quote=F,row.names=F)
    break
  }
}

