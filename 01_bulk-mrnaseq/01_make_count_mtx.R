path <- "/Users/plezar/Documents/TIME-Lab/SpaceU87_Nov2024/featureCounts/"
files <- list.files(path, pattern = "COUNTS_")

merge_out <- function (files) {
  df <- read.delim(paste0(path, files), header= T)
  as.matrix(df[,2])
}

results <- lapply(files, merge_out)
names(results) <- gsub("_counts.txt", "", gsub("COUNTS_", "", files))

results <- bind_cols(results)
results <- as.data.frame(results)
rownames(results) <- read.delim(paste0(path, files[1]), header= T)[,1]

write.csv(results, "/Users/plezar/Documents/TIME-Lab/SpaceU87_Nov2024/count_mtx.csv")
