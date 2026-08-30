# =========================
# 0. Libraries
# =========================
library(data.table)
library(readr)
getwd()

# =========================
# 1. Load + process data
# =========================

microarray_DEG <- readRDS(snakemake@input[["array"]])
rna_DEG        <- readRDS(snakemake@input[["rna"]])

# =========================
# 1. DEG processing function

# =========================

process_expr_dt <- function(list_dt,
                            logfc_cut = 0,
                            padj_cut = 0.05,
                            padj_col = NULL,
                            min_deg = 50) {
  
  list_dt <- lapply(list_dt, as.data.table)
  
  if (is.null(padj_col)) {
    padj_col <- intersect(c("adj.P.Val", "BH"), names(list_dt[[1]]))[1]
  }
  
  up <- sapply(list_dt, function(x)
    x[get(padj_col) < padj_cut & logFC > 1, .N]
  )
  
  down <- sapply(list_dt, function(x)
    x[get(padj_col) < padj_cut & logFC < -1, .N]
  )
  
  keep <- up >= min_deg & down >= min_deg
  list_dt <- list_dt[keep]
  
  names(list_dt) <- names(keep)[keep]
  
  lapply(list_dt, function(x) {
    x[, logFC := fifelse(
      get(padj_col) < padj_cut & abs(logFC) > logfc_cut,
      as.character(logFC),
      "NS"
    )]
    x[, .(TAIR_ID, logFC)]
  })
}

microarray_DEG <- process_expr_dt(microarray_DEG, logfc_cut = 0, padj_cut = 0.05)
rna_DEG <- process_expr_dt(rna_DEG, logfc_cut = 0, padj_cut = 0.05)

# =========================
# 3. Load Arabidopsis genes
# =========================

ara_genes <- readRDS("resource/Arabidopsis_genes.rds")
ara_genes <- ara_genes[ ,.(TAIR_ID)]
  
# =========================
# 4. Merge DEG tables
# =========================
merge_deg_list <- function(deg_list, ara_genes) {
  
    dt_list <- lapply(names(deg_list), function(name) {
    dt <- as.data.table(deg_list[[name]])
    setnames(dt, "logFC", name)
    dt
  })
  
  Reduce(
    function(x, y) merge(x, y, by = "TAIR_ID", all.x = TRUE, sort = FALSE),
    c(list(ara_genes), dt_list)
  )
}

ma_all  <- merge_deg_list(microarray_DEG, ara_genes)
rna_all <- merge_deg_list(rna_DEG, ara_genes)

# =========================
# 5. Convert NS/NA → numeric
# =========================

convert_na_to_LE <- function(dt) {
  cols <- setdiff(names(dt), "TAIR_ID")
  
  dt[, (cols) := lapply(.SD, function(x) {
    ifelse(is.na(x), "LE", x)
  }), .SDcols = cols]
  
  dt
} 

rna_all <- convert_na_to_LE(rna_all) #run only for rna data!

setkey(ma_all, TAIR_ID)
setkey(rna_all, TAIR_ID)

allDEG <- ma_all[rna_all] #NA in microarray data means that expression of these genes are missing in this type of data
all_data <- copy(allDEG)

convert_to_numeric <- function(dt) {
  cols <- setdiff(names(dt), "TAIR_ID")
  
  dt[, (cols) := lapply(.SD, function(x) {
    
    out <- rep(0, length(x))
    
    valid <- !(is.na(x) | x %in% c("LE", "NS"))
    
    out[valid] <- as.numeric(x[valid])
    
    out
  }), .SDcols = cols]
  
  dt
}

all_data <- convert_to_numeric(all_data)
all_data_da <- copy(all_data)

all_data_da <- convert_to_numeric(all_data_da)
cols <- setdiff("TAIR_ID", names(all_data_da))
index <- which(rowSums(all_data_da[ ,-1] == 0) == length(names(all_data_da)) - 1)
all_data_da <- all_data_da[-index, ]


# =========================
# 6. Order columns by time
# =========================
order_by_time <- function(dt) {
  cols <- setdiff(names(dt), "TAIR_ID")
  time <- as.numeric(sub("h_.*", "", cols))
  ordered_cols <- cols[order(time)]
  setcolorder(dt, c("TAIR_ID", ordered_cols))
  dt
}

allDEG <- order_by_time(allDEG)
all_data_da <- order_by_time(all_data_da)

# =========================
# 7. Clustering
# =========================
scale_data <- scale(all_data_da[, -1], center = TRUE, scale = TRUE)
scale_dist <- dist(t(scale_data))
scale_hclust <- hclust(scale_dist, method = "ward.D2")

# IMPORTANT:
# plot(scale_hclust)
# cutree(scale_hclust, 5)
# Choose k manually based on dendrogram. Ater you choose run cutree(scale_hclust, k) where k is the number of desired clusters
# Then look which cluster that datasets you want to keep assigned to and write the special rules to expract them

# =========================
# 8. Time Helper
# =========================
extract_time <- function(x) {
  as.numeric(sub("h_.*", "", x))
}

# =========================
# 9. Cluster rules
# =========================

# STUDY-SPECIFIC RULES
rules <- list(
  ER = quote(time < 12 | sample == "24h_4C_T87_cells_GSE31837"),
  LR = quote(time %in% c(12, 24)),
  VR = quote(time >= 48)
)

cluster_rule_map <- c(
  "1" = "ER",
  "4" = "LR",
  "5" = "VR"
)

filter_by_cluster_rules <- function(all_data, hclust_obj, k, rules, cluster_rule_map) {
  
  clusters <- cutree(hclust_obj, k = k)
  samples <- names(clusters)
  time <- extract_time(samples)
  
  df <- data.table(sample = samples, cluster = clusters, time = time)
  df[, keep := FALSE]
  
  for (cl in unique(df$cluster)) {
    
    rule_name <- cluster_rule_map[as.character(cl)]
    if (is.na(rule_name)) next
    
    rule <- rules[[rule_name]]
    
    df[cluster == cl, keep := eval(rule)]
  }
  
  keep_samples <- df[keep == TRUE, sample]
  
  result <- all_data[, c("TAIR_ID", keep_samples), with = FALSE]
  
  list(data = result, decisions = df)
}

res <- filter_by_cluster_rules(
  all_data_da,
  scale_hclust,
  k = 5,
  rules,
  cluster_rule_map
)

all_data_filtered <- res$data
res$decisions[grep("24h.*GSE31837", sample), cluster := 4]

filter_logFC <- function(dt, logfc_cut){
  cols <- setdiff(names(dt), "TAIR_ID")
  
  dt[, (cols) := lapply(.SD, function(x) {
    fifelse(x > logfc_cut | x < -logfc_cut, x, 0)
  }), .SDcols = cols]
  dt
}

all_data_filtered <- filter_logFC(all_data_filtered, logfc_cut = 1)

# =========================
# 10. Cluster labeling
# =========================
cluster_map <- c("1" = "ER", "4" = "LR", "5" = "VR")

res$decisions[, cluster_label := cluster_map[as.character(cluster)]]
decision_table <- res$decisions[keep == TRUE]

cluster_samples <- split(
  decision_table$sample,
  decision_table$cluster_label
)

allDEG <- allDEG[ ,c("TAIR_ID", decision_table$sample), with = F]

# =========================
# 11. Count DEGs
# =========================
get_cluster_counts <- function(all_data, cluster_samples) {
  
  lapply(cluster_samples, function(samples) {
    
    sub <- all_data[, ..samples]
    
    list(
      up   = rowSums(sub > 1, na.rm = TRUE),
      down = rowSums(sub < -1, na.rm = TRUE)
    )
  })
}

results_list <- get_cluster_counts(all_data_filtered, cluster_samples)

# =========================
# 12. Final count table
# =========================
final_dt <- data.table(
  TAIR_ID = all_data_filtered$TAIR_ID,
  upER   = results_list$ER$up,
  upLR   = results_list$LR$up,
  upVR   = results_list$VR$up,
  downER = results_list$ER$down,
  downLR = results_list$LR$down,
  downVR = results_list$VR$down
)

# =========================
# 13. Binomial thresholds
# =========================
get_k_threshold <- function(n, p = 0.05, alpha = 0.01) {
  
  ks <- 0:n
  
  pvals <- sapply(ks, function(k) {
    binom.test(k, n, p = p, alternative = "greater")$p.value
  })
  
  min(ks[pvals < alpha])
}

cluster_sizes <- sapply(cluster_samples, length)
k_thresh <- sapply(cluster_sizes, get_k_threshold)

# =========================
# 14. Filter robust DEGs
# =========================
robust_dt <- final_dt[
  (upER   >= k_thresh["ER"] | downER >= k_thresh["ER"]) |
    (upLR   >= k_thresh["LR"] | downLR >= k_thresh["LR"]) |
    (upVR   >= k_thresh["VR"] | downVR >= k_thresh["VR"])
]

# =========================
# 15. Save results
# =========================

fwrite(allDEG, snakemake@output[["all"]], sep = ",")
fwrite(robust_dt, snakemake@output[["final"]], sep = ",")

write_robustDEG <- function(dt, k_cut, out) {
  
  cols <- setdiff(names(dt), "TAIR_ID")
  
  invisible(lapply(seq_along(cols), function(i) {
    
    col <- cols[i]
    key <- gsub("up|down", "", col)
    k <- k_cut[[key]]
    
    fwrite(
      dt[get(col) >= k, .(TAIR_ID)],
      out[i], 
      sep = ","
    )
  }))
}

write_robustDEG(robust_dt, k_thresh, snakemake@output[["deg"]])
write.table(t(k_thresh), snakemake@output[["threshold"]], quote = F, row.names = F)
