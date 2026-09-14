rm(list = ls())
##加载R包
library(tidyr)
library(dplyr)
library(tibble)
##清空数据riskTest riskTest riskTest riskTest riskTest
##清空数据riskTest riskTest riskTest riskTest riskTest
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1) 
riskTest=read.table("output/riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)
##找出coef
cor <- multiCox %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskTest
rx=riskTest
exp1 <- rx %>% 
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore) 
rx <- rx %>% 
  rownames_to_column("id")
exp1 <- exp1 %>% 
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskTest <- rx %>% 
  column_to_rownames("id") 

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskTest <- swap_columns(riskTest, "riskscore", "riskScore")

riskTest <- riskTest %>% 
  dplyr::select(-riskScore) %>% 
  dplyr::rename(riskScore = riskscore) %>% 
  dplyr::select(-risk)

median_value <- median(riskTest$riskScore, na.rm = TRUE)
riskTest$risk <- ifelse(riskTest$riskScore <= median_value, "low", "high")
write.table(riskTest,file="output/riskTest.txt",sep="\t",quote=F,row.names=T)

##清空数据riskTrain riskTrain riskTrain riskTrain riskTrain
##清空数据riskTrain riskTrain riskTrain riskTrain riskTrain
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1) 
riskTrain=read.table("output/riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1) 
##找出coef
cor <- multiCox %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskTrain
rx=riskTrain
exp1 <- rx %>% 
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore) 
rx <- rx %>% 
  rownames_to_column("id")
exp1 <- exp1 %>% 
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskTrain <- rx %>% 
  column_to_rownames("id") 

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskTrain <- swap_columns(riskTrain, "riskscore", "riskScore")

riskTrain <- riskTrain %>% 
  dplyr::select(-riskScore) %>% 
  dplyr::rename(riskScore = riskscore) %>% 
  dplyr::select(-risk)

median_value <- median(riskTrain$riskScore, na.rm = TRUE)
riskTrain$risk <- ifelse(riskTrain$riskScore <= median_value, "low", "high")

write.table(riskTrain,file="output/riskTrain.txt",sep="\t",quote=F,row.names=T)


##清空数据riskoutside1 riskoutside1 riskoutside1 riskoutside1 riskoutside1
##清空数据riskoutside1 riskoutside1 riskoutside1 riskoutside1 riskoutside1
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1) 
riskoutside1=read.table("output/riskoutside1.txt",header=T,sep="\t",check.names=F,row.names=1)
##找出coef
cor <- multiCox %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskoutside
rx=riskoutside1
exp1 <- rx %>% 
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore) 
rx <- rx %>% 
  rownames_to_column("id")
exp1 <- exp1 %>% 
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskoutside1 <- rx %>% 
  column_to_rownames("id") 

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskoutside1 <- swap_columns(riskoutside1, "riskscore", "riskScore")

riskoutside1 <- riskoutside1 %>% 
  dplyr::select(-riskScore) %>% 
  dplyr::rename(riskScore = riskscore) %>% 
  dplyr::select(-risk)

median_value <- median(riskoutside1$riskScore, na.rm = TRUE)
riskoutside1$risk <- ifelse(riskoutside1$riskScore <= median_value, "low", "high")

write.table(riskoutside1,file="output/riskoutside1.txt",sep="\t",quote=F,row.names=T)



##清空数据riskoutside2 riskoutside2 riskoutside2 riskoutside2 riskoutside2
##清空数据riskoutside2 riskoutside2 riskoutside2 riskoutside2 riskoutside2
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1) 
riskoutside2=read.table("output/riskoutside2.txt",header=T,sep="\t",check.names=F,row.names=1)
##找出coef
cor <- multiCox %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskoutside2
rx=riskoutside2
exp1 <- rx %>% 
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore) 
rx <- rx %>% 
  rownames_to_column("id")
exp1 <- exp1 %>% 
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskoutside2 <- rx %>% 
  column_to_rownames("id") 

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskoutside2 <- swap_columns(riskoutside2, "riskscore", "riskScore")

riskoutside2 <- riskoutside2 %>% 
  dplyr::select(-riskScore) %>% 
  dplyr::rename(riskScore = riskscore) %>% 
  dplyr::select(-risk)

median_value <- median(riskoutside2$riskScore, na.rm = TRUE)
riskoutside2$risk <- ifelse(riskoutside2$riskScore <= median_value, "low", "high")

write.table(riskoutside2,file="output/riskoutside2.txt",sep="\t",quote=F,row.names=T)


##清空数据riskoutside3 riskoutside3 riskoutside3 riskoutside3 riskoutside3
##清空数据riskoutside3 riskoutside3 riskoutside3 riskoutside3 riskoutside3
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1) 
riskoutside3=read.table("output/riskoutside3.txt",header=T,sep="\t",check.names=F,row.names=1)
##找出coef
cor <- multiCox %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskoutside3
rx=riskoutside3
exp1 <- rx %>% 
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore) 
rx <- rx %>% 
  rownames_to_column("id")
exp1 <- exp1 %>% 
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskoutside3 <- rx %>% 
  column_to_rownames("id") 

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskoutside3 <- swap_columns(riskoutside3, "riskscore", "riskScore")

riskoutside3 <- riskoutside3 %>% 
  dplyr::select(-riskScore) %>% 
  dplyr::rename(riskScore = riskscore) %>% 
  dplyr::select(-risk)

median_value <- median(riskoutside3$riskScore, na.rm = TRUE)
riskoutside3$risk <- ifelse(riskoutside3$riskScore <= median_value, "low", "high")

write.table(riskoutside3,file="output/riskoutside3.txt",sep="\t",quote=F,row.names=T)


# ##清空数据riskoutside4 riskoutside4 riskoutside4 riskoutside4 riskoutside4
# ##清空数据riskoutside4 riskoutside4 riskoutside4 riskoutside4 riskoutside4
# rm(list = ls())
# ###重新计算riskscore
# multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1) 
# riskoutside4=read.table("output/riskoutside4.txt",header=T,sep="\t",check.names=F,row.names=1)
# ##找出coef
# cor <- multiCox %>% 
#   rownames_to_column("id") %>% 
#   dplyr::select(id,coef)
# cor$coef <- as.numeric(cor$coef)
# cor_list <- setNames(cor$coef, cor$id)
# 
# ##riskoutside4
# rx=riskoutside4
# exp1 <- rx %>% 
#   dplyr::select(-c(futime,fustat,risk,riskScore))
# 
# # 初始化 riskscore 向量
# riskscores <- numeric(nrow(exp1))
# # 计算每个样本的 riskscore
# for (i in 1:nrow(exp1)) {
#   # 对当前样本的每个基因值乘以对应的 cor 值，并求和
#   riskscores[i] <- sum(exp1[i, ] * cor_list)
# }
# # 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
# exp1$riskscore <- riskscores
# exp1 <- exp1 %>%
#   dplyr::select(riskscore) 
# rx <- rx %>% 
#   rownames_to_column("id")
# exp1 <- exp1 %>% 
#   rownames_to_column("id")
# rx <- merge(rx, exp1, by = "id")
# riskoutside4 <- rx %>% 
#   column_to_rownames("id") 
# 
# ##调换列名，其他顺序不变
# swap_columns <- function(df, col1, col2) {
#   # 检查提供的列名是否都在数据框中
#   if (!col1 %in% names(df) || !col2 %in% names(df)) {
#     stop("One or both column names do not exist in the dataframe.")
#   }
#   # 获取列的位置
#   col1_index <- which(names(df) == col1)
#   col2_index <- which(names(df) == col2)
#   # 创建一个列索引的向量
#   column_indices <- seq_along(df)
#   # 交换两个列索引
#   column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
#   # 使用更新的索引向量重新排列列
#   df <- df[, column_indices]
#   return(df)
# }
# riskoutside4 <- swap_columns(riskoutside4, "riskscore", "riskScore")
# 
# riskoutside4 <- riskoutside4 %>% 
#   dplyr::select(-riskScore) %>% 
#   dplyr::rename(riskScore = riskscore) %>% 
#   dplyr::select(-risk)
# 
# median_value <- median(riskoutside4$riskScore, na.rm = TRUE)
# riskoutside4$risk <- ifelse(riskoutside4$riskScore <= median_value, "low", "high")
# 
# write.table(riskoutside4,file="output/riskoutside4.txt",sep="\t",quote=F,row.names=T)


##清空数据riskoutside5 riskoutside5 riskoutside5 riskoutside5 riskoutside5
##清空数据riskoutside5 riskoutside5 riskoutside5 riskoutside5 riskoutside5
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1)
riskoutside5=read.table("output/riskoutside5.txt",header=T,sep="\t",check.names=F,row.names=1)
##找出coef
cor <- multiCox %>%
  rownames_to_column("id") %>%
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskoutside5
rx=riskoutside5
exp1 <- rx %>%
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore)
rx <- rx %>%
  rownames_to_column("id")
exp1 <- exp1 %>%
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskoutside5 <- rx %>%
  column_to_rownames("id")

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskoutside5 <- swap_columns(riskoutside5, "riskscore", "riskScore")

# riskoutside5 <- riskoutside5 %>%
#   dplyr::select(-riskScore) %>%
#   dplyr::rename(riskScore = riskscore) %>%
#   dplyr::select(-risk)

riskoutside5 <- riskoutside5 %>%
  dplyr::select(-riskScore) %>%
  dplyr::rename(riskScore = riskscore)

# median_value <- median(riskoutside5$riskScore, na.rm = TRUE)
# riskoutside5$risk <- ifelse(riskoutside5$riskScore <= median_value, "low", "high")

write.table(riskoutside5,file="output/riskoutside5.txt",sep="\t",quote=F,row.names=T)


##清空数据riskoutside6 riskoutside6 riskoutside6 riskoutside6 riskoutside6
##清空数据riskoutside6 riskoutside6 riskoutside6 riskoutside6 riskoutside6
# rm(list = ls())
# ###重新计算riskscore
# multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1)
# riskoutside6=read.table("output/riskoutside6.txt",header=T,sep="\t",check.names=F,row.names=1)
# ##找出coef
# cor <- multiCox %>%
#   rownames_to_column("id") %>%
#   dplyr::select(id,coef)
# cor$coef <- as.numeric(cor$coef)
# cor_list <- setNames(cor$coef, cor$id)
# 
# 
# ##riskoutside6
# rx=riskoutside6
# exp1 <- rx %>%
#   dplyr::select(-c(futime,fustat,risk,riskScore))
# 
# # 初始化 riskscore 向量
# riskscores <- numeric(nrow(exp1))
# # 计算每个样本的 riskscore
# for (i in 1:nrow(exp1)) {
#   # 对当前样本的每个基因值乘以对应的 cor 值，并求和
#   riskscores[i] <- sum(exp1[i, ] * cor_list)
# }
# # 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
# exp1$riskscore <- riskscores
# exp1 <- exp1 %>%
#   dplyr::select(riskscore)
# rx <- rx %>%
#   rownames_to_column("id")
# exp1 <- exp1 %>%
#   rownames_to_column("id")
# rx <- merge(rx, exp1, by = "id")
# riskoutside6 <- rx %>%
#   column_to_rownames("id")
# 
# ##调换列名，其他顺序不变
# swap_columns <- function(df, col1, col2) {
#   # 检查提供的列名是否都在数据框中
#   if (!col1 %in% names(df) || !col2 %in% names(df)) {
#     stop("One or both column names do not exist in the dataframe.")
#   }
#   # 获取列的位置
#   col1_index <- which(names(df) == col1)
#   col2_index <- which(names(df) == col2)
#   # 创建一个列索引的向量
#   column_indices <- seq_along(df)
#   # 交换两个列索引
#   column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
#   # 使用更新的索引向量重新排列列
#   df <- df[, column_indices]
#   return(df)
# }
# riskoutside6 <- swap_columns(riskoutside6, "riskscore", "riskScore")
# 
# # riskoutside6 <- riskoutside6 %>%
# #   dplyr::select(-riskScore) %>%
# #   dplyr::rename(riskScore = riskscore) %>%
# #   dplyr::select(-risk)
# 
# riskoutside6 <- riskoutside6 %>%
#   dplyr::select(-riskScore) %>%
#   dplyr::rename(riskScore = riskscore)
# 
# # median_value <- median(riskoutside6$riskScore, na.rm = TRUE)
# # riskoutside6$risk <- ifelse(riskoutside6$riskScore <= median_value, "low", "high")
# 
# write.table(riskoutside6,file="output/riskoutside6.txt",sep="\t",quote=F,row.names=T)
# 

##清空数据riskoutside7 riskoutside7 riskoutside7 riskoutside7 riskoutside7
##清空数据riskoutside7 riskoutside7 riskoutside7 riskoutside7 riskoutside7
rm(list = ls())
###重新计算riskscore
multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1)
riskoutside7=read.table("output/riskoutside7.txt",header=T,sep="\t",check.names=F,row.names=1)
##找出coef
cor <- multiCox %>%
  rownames_to_column("id") %>%
  dplyr::select(id,coef)
cor$coef <- as.numeric(cor$coef)
cor_list <- setNames(cor$coef, cor$id)

##riskoutside7
rx=riskoutside7
exp1 <- rx %>%
  dplyr::select(-c(futime,fustat,risk,riskScore))

# 初始化 riskscore 向量
riskscores <- numeric(nrow(exp1))
# 计算每个样本的 riskscore
for (i in 1:nrow(exp1)) {
  # 对当前样本的每个基因值乘以对应的 cor 值，并求和
  riskscores[i] <- sum(exp1[i, ] * cor_list)
}
# 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
exp1$riskscore <- riskscores
exp1 <- exp1 %>%
  dplyr::select(riskscore)
rx <- rx %>%
  rownames_to_column("id")
exp1 <- exp1 %>%
  rownames_to_column("id")
rx <- merge(rx, exp1, by = "id")
riskoutside7 <- rx %>%
  column_to_rownames("id")

##调换列名，其他顺序不变
swap_columns <- function(df, col1, col2) {
  # 检查提供的列名是否都在数据框中
  if (!col1 %in% names(df) || !col2 %in% names(df)) {
    stop("One or both column names do not exist in the dataframe.")
  }
  # 获取列的位置
  col1_index <- which(names(df) == col1)
  col2_index <- which(names(df) == col2)
  # 创建一个列索引的向量
  column_indices <- seq_along(df)
  # 交换两个列索引
  column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
  # 使用更新的索引向量重新排列列
  df <- df[, column_indices]
  return(df)
}
riskoutside7 <- swap_columns(riskoutside7, "riskscore", "riskScore")

# riskoutside7 <- riskoutside7 %>%
#   dplyr::select(-riskScore) %>%
#   dplyr::rename(riskScore = riskscore) %>%
#   dplyr::select(-risk)

riskoutside7 <- riskoutside7 %>%
  dplyr::select(-riskScore) %>%
  dplyr::rename(riskScore = riskscore) 

# median_value <- median(riskoutside7$riskScore, na.rm = TRUE)
# riskoutside7$risk <- ifelse(riskoutside7$riskScore <= median_value, "low", "high")

write.table(riskoutside7,file="output/riskoutside7.txt",sep="\t",quote=F,row.names=T)


# ##清空数据riskoutside8 riskoutside8 riskoutside8 riskoutside8 riskoutside8
# ##清空数据riskoutside8 riskoutside8 riskoutside8 riskoutside8 riskoutside8
# rm(list = ls())
# ###重新计算riskscore
# multiCox=read.table("output/multiCox.xls",header=T,sep="\t",check.names=F,row.names=1)
# riskoutside8=read.table("output/riskoutside8.txt",header=T,sep="\t",check.names=F,row.names=1)
# ##找出coef
# cor <- multiCox %>%
#   rownames_to_column("id") %>%
#   dplyr::select(id,coef)
# cor$coef <- as.numeric(cor$coef)
# cor_list <- setNames(cor$coef, cor$id)
# 
# ##riskoutside8
# rx=riskoutside8
# exp1 <- rx %>%
#   dplyr::select(-c(futime,fustat,risk,riskScore))
# 
# # 初始化 riskscore 向量
# riskscores <- numeric(nrow(exp1))
# # 计算每个样本的 riskscore
# for (i in 1:nrow(exp1)) {
#   # 对当前样本的每个基因值乘以对应的 cor 值，并求和
#   riskscores[i] <- sum(exp1[i, ] * cor_list)
# }
# # 将 riskscore 向量转换为数据框的一列，并添加到 exp 数据框中
# exp1$riskscore <- riskscores
# exp1 <- exp1 %>%
#   dplyr::select(riskscore)
# rx <- rx %>%
#   rownames_to_column("id")
# exp1 <- exp1 %>%
#   rownames_to_column("id")
# rx <- merge(rx, exp1, by = "id")
# riskoutside8 <- rx %>%
#   column_to_rownames("id")
# 
# ##调换列名，其他顺序不变
# swap_columns <- function(df, col1, col2) {
#   # 检查提供的列名是否都在数据框中
#   if (!col1 %in% names(df) || !col2 %in% names(df)) {
#     stop("One or both column names do not exist in the dataframe.")
#   }
#   # 获取列的位置
#   col1_index <- which(names(df) == col1)
#   col2_index <- which(names(df) == col2)
#   # 创建一个列索引的向量
#   column_indices <- seq_along(df)
#   # 交换两个列索引
#   column_indices[c(col1_index, col2_index)] <- column_indices[c(col2_index, col1_index)]
#   # 使用更新的索引向量重新排列列
#   df <- df[, column_indices]
#   return(df)
# }
# riskoutside8 <- swap_columns(riskoutside8, "riskscore", "riskScore")
# 
# # riskoutside8 <- riskoutside8 %>%
# #   dplyr::select(-riskScore) %>%
# #   dplyr::rename(riskScore = riskscore) %>%
# #   dplyr::select(-risk)
# 
# riskoutside8 <- riskoutside8 %>%
#   dplyr::select(-riskScore) %>%
#   dplyr::rename(riskScore = riskscore) 
# 
# # median_value <- median(riskoutside8$riskScore, na.rm = TRUE)
# # riskoutside8$risk <- ifelse(riskoutside8$riskScore <= median_value, "low", "high")
# 
# write.table(riskoutside8,file="output/riskoutside8.txt",sep="\t",quote=F,row.names=T)
# 
