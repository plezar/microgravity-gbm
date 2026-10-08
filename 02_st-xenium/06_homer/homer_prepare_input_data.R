top_genes <- read.csv("/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/results/xenium_cNMF_k4_topgenes.csv", row.names = 1)
write.table(data.frame(Acc = top_genes$X1),
            "02_st-xenium/results/homer_in/GEP_1.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = top_genes$X2),
            "02_st-xenium/results/homer_in/GEP_2.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = top_genes$X3),
            "02_st-xenium/results/homer_in/GEP_3.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(data.frame(Acc = top_genes$X4),
            "02_st-xenium/results/homer_in/GEP_4.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)