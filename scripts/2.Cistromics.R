library(data.table)

# =========================
# STEP 2. TF enrichment (CisCross)
# =========================

# -------------------------
# 1. Load binomial thresholds
# -------------------------
# k_thresh contains ER, LR, VR thresholds (binomial-based)

k_thresh_dt <- fread(snakemake@input[["threshold"]])

expand_k_thresh <- function(k_thresh) {
  # Expand ER/LR/VR thresholds to up/down directions
  c(
    upER   = k_thresh$ER,
    downER = k_thresh$ER,
    
    upLR   = k_thresh$LR,
    downLR = k_thresh$LR,
    
    upVR   = k_thresh$VR,
    downVR = k_thresh$VR
  )
}

k_thresh_expanded <- expand_k_thresh(k_thresh_dt)

# -------------------------
# 2. Load DEG count table
# -------------------------
deg_counts_dt <- fread(snakemake@input[["count_table"]])

# -------------------------
# 3. Extract robust DEG sets
# -------------------------
get_deg_sets <- function(count_table, k_thresh_vec) {
  
  deg_sets <- lapply(names(k_thresh_vec), function(cond) {
    count_table[get(cond) >= k_thresh_vec[[cond]], TAIR_ID]
  })
  
  names(deg_sets) <- names(k_thresh_vec)
  
  return(deg_sets)
}

deg_sets <- get_deg_sets(deg_counts_dt, k_thresh_expanded)

# -------------------------
# 4. Load CisCross results
# -------------------------
load_ciscross <- function(path) {
  
  files <- list.files(path = path, pattern = "\\.csv$", full.names = TRUE)
  
  tf_list <- lapply(files, fread)
  
  names(tf_list) <- gsub("\\.csv$", "", basename(files))
  
  return(tf_list)
}

tf_raw_list <- load_ciscross(snakemake@input[["ciscross_dir"]])

# -------------------------
# 5. Filter significant TF targets
# -------------------------
# NOTE: Column names are dataset-specific (spaces!)
# We standardize them here

tf_filtered_list <- lapply(tf_raw_list, function(dt) {
  
  # Filter by adjusted p-value
  dt <- dt[dt[["adjusted p_value"]] < 0.05, ]
  
  # Standardize column name
  setnames(dt, "TAIR ID", "TAIR_ID")
  
  # Keep only relevant columns
  dt[, .(TAIR_ID, adjusted_p_value = `adjusted p_value`)]
})

# -------------------------
# 6. Ensure correct naming
# -------------------------
# IMPORTANT:
# Names must match DEG sets exactly

names(tf_filtered_list) <- c("downER","upER","downLR","upLR","downVR","upVR")
names(tf_raw_list)      <- names(tf_filtered_list)

# -------------------------
# 7. Define TF-DEG mapping
# -------------------------
# Study-specific rule:
# TFs must overlap both up- and down-regulated DEGs within a phase

tf_deg_map <- list(
  upER   = c("upER", "downER"),
  downER = c("downER", "upER"),
  
  upLR   = c("upLR", "downLR"),
  downLR = c("downLR", "upLR"),
  
  upVR   = c("upVR", "downVR"),
  downVR = c("downVR", "upVR")
)

# -------------------------
# 8. Compute TF-DEG overlap
# -------------------------
get_tf_deg_overlap <- function(tf_sets, deg_sets, tf_deg_map) {
  
  tf_overlap <- lapply(names(tf_deg_map), function(name) {
    
    tf_dt <- tf_sets[[name]]
    deg_groups <- tf_deg_map[[name]]
    
    # Combine DEG sets
    deg_genes <- unlist(deg_sets[deg_groups], use.names = FALSE)
    
    # Keep TFs targeting DEG genes
    tf_hits <- tf_dt[TAIR_ID %in% deg_genes]
    
    unique(tf_hits$TAIR_ID)
  })
  
  names(tf_overlap) <- names(tf_deg_map)
  
  return(tf_overlap)
}

tf_deg_overlap <- get_tf_deg_overlap(tf_filtered_list, deg_sets, tf_deg_map)

# -------------------------
# 9. Collect candidate TFs
# -------------------------
tf_candidates <- unique(unlist(tf_deg_overlap, use.names = FALSE))

# -------------------------
# 10. Extract TF statistics
# -------------------------
tf_candidate_stats <- lapply(tf_raw_list, function(dt) {
  
  dt <- dt[dt$`TAIR ID` %in% tf_candidates, .(`TAIR ID`, `adjusted p_value`)]
  
  dt <- unique(dt, by = "TAIR ID")
  
  dt[order(`TAIR ID`)]
})

# -------------------------
# 11. Build final TF table
# -------------------------
tf_tables <- Map(function(dt, nm) {
  
  dt <- as.data.table(dt)
  
  setnames(dt, c("TAIR ID", "adjusted p_value"), c("TAIR_ID", nm))
  
  dt
  
}, tf_candidate_stats, names(tf_candidate_stats))

tf_final_table <- Reduce(function(x, y) {
  merge(x, y, by = "TAIR_ID", all = TRUE)
}, tf_tables)

# -------------------------
# 12. Save result
# -------------------------
fwrite(tf_final_table, snakemake@output[["robustTF_pval"]])