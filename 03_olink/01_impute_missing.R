library(readxl)
library(missForest)
library(doParallel)

# ===== All data =====

dt <- read_xlsx("03_olink/data/olink_all_data.xlsx")
dt <- as.matrix(dt)
sample_names <- dt[,1]
dt <- dt[,-1]

dt <- apply(dt, 2, as.numeric)
rownames(dt) <- sample_names

# Examine missingness
par(mfrow = c(1, 2))   # 1 row, 2 columns
hist(
  rowMeans(is.na(dt)),
  breaks = 10,
  col = "skyblue",
  border = "white",
  main = "Missingness per Row",
  xlab = "Proportion Missing",
  las = 1
)
abline(v = .45,
       col = "red",
       lwd = 2,
       lty = 2)
hist(
  colMeans(is.na(dt)),
  breaks = 30,
  col = "salmon",
  border = "white",
  main = "Missingness per Column",
  xlab = "Proportion Missing",
  las = 1
)
abline(v = .4,
       col = "red",
       lwd = 2,
       lty = 2)
par(mfrow = c(1, 1))   # reset

# Drop samples with high missingness
dt <- dt[rowMeans(is.na(dt))<.45,]

# Drop features that are completely missing
dt <- dt[,colMeans(is.na(dt))<.4]

registerDoParallel(cores=10)
mf_result <- missForest(dt, parallelize = "variables")
# $OOBerror
# NRMSE 
# 0.1258091

write.csv(data.frame(pct_missing = rowMeans(is.na(dt))), "03_olink/data/olink_all_data_imputed_pct_missing.csv")
write.csv(mf_result$ximp, "03_olink/data/olink_all_data_imputed.csv")


# ===== QC data =====

dt <- read_xlsx("03_olink/data/olink_all_data_qc.xlsx")
dt <- as.matrix(dt)
sample_names <- dt[,1]
dt <- dt[,-1]

dt <- apply(dt, 2, as.numeric)
rownames(dt) <- sample_names

# Examine missingness
par(mfrow = c(1, 2))   # 1 row, 2 columns
hist(
  rowMeans(is.na(dt)),
  breaks = 30,
  col = "skyblue",
  border = "white",
  main = "Missingness per Row",
  xlab = "Proportion Missing",
  las = 1
)
abline(v = .74,
       col = "red",
       lwd = 2,
       lty = 2)
hist(
  colMeans(is.na(dt)),
  breaks = 30,
  col = "salmon",
  border = "white",
  main = "Missingness per Column",
  xlab = "Proportion Missing",
  las = 1
)
abline(v = .5,
       col = "red",
       lwd = 2,
       lty = 2)
par(mfrow = c(1, 1))   # reset

# Drop samples with high missingness
dt <- dt[rowMeans(is.na(dt))<.74,]

# Drop features that are completely missing
dt <- dt[,colMeans(is.na(dt))<.5]

registerDoParallel(cores=10)
mf_result <- missForest(dt, parallelize = "variables")
# $OOBerror
# NRMSE 
# 0.07655388

write.csv(data.frame(pct_missing = rowMeans(is.na(dt))), "03_olink/data/olink_all_data_imputed_pct_qc_missing.csv")
write.csv(mf_result$ximp, "03_olink/data/olink_all_data_qc_imputed.csv")