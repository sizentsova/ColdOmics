library(rvest)
library(stringr)
library(data.table)
library(tools)

# =========================
# 1. Load CisCross HTML/TXT reports
# =========================

ciscross_dir <- snakemake@input[["ciscross_dir"]]

report_files <- list.files(
  ciscross_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

# Read HTML reports
report_html <- lapply(report_files, read_html)

# Convert reports into line-by-line text
report_lines <- lapply(report_html, function(x) {
  
  html_text(x) |>
    strsplit("\r\n") |>
    unlist()
})

names(report_lines) <- tools::file_path_sans_ext(basename(report_files))

# =========================
# 2. Helper function
# Extract table between markers
# =========================

extract_ciscross_table <- function(lines,
                                   start_pattern,
                                   end_pattern,
                                   first_only = FALSE) {
  
  start_idx <- grep(start_pattern, lines)
  end_idx   <- grep(end_pattern, lines)
  
  if (first_only) {
    start_idx <- head(start_idx, 1)
    end_idx   <- head(end_idx, 1)
  }
  
  if (length(start_idx) == 0 || length(end_idx) == 0) {
    return(NULL)
  }
  
  section_lines <- lines[start_idx:(end_idx - 1)]
  
  dt <- fread(
    text = paste(section_lines, collapse = "\n"),
    sep = "\t",
    header = TRUE
  )
  
  # Clean column names
  setnames(dt, gsub("\\.", " ", names(dt)))
  
  return(dt)
}

# =========================
# 3. Extract UP enrichment tables
# =========================

up_tables <- lapply(report_lines, function(lines) {
  
  extract_ciscross_table(
    lines,
    start_pattern = "^TF Name.*UP",
    end_pattern   = "^DONE UP"
  )
})

# =========================
# 4. Extract DOWN enrichment tables
# =========================

down_tables <- lapply(report_lines, function(lines) {
  
  extract_ciscross_table(
    lines,
    start_pattern = "^TF Name.*DOWN",
    end_pattern   = "^DONE DOWN"
  )
})

# =========================
# 4. Make CisCross summary table
# =========================


tf_list <- list(
  
  downER = down_tables$ciscross_ER_1500p_0.05,
  upER   = up_tables$ciscross_ER_1500p_0.05,
  
  downLR = down_tables$ciscross_LR_1500p_0.05,
  upLR   = up_tables$ciscross_LR_1500p_0.05,
  
  downVR = down_tables$ciscross_VR_1500p_0.05,
  upVR   = up_tables$ciscross_VR_1500p_0.05
)

tf_list <- lapply(tf_list, function(dt) {
  
  # Standardize column name
  setnames(dt, "TAIR ID", "TAIR_ID")
  
  # Keep only relevant columns
  dt[, .(TAIR_ID, `adjusted p_value`)]
})

names(tf_list) <- c("downER","upER","downLR","upLR","downVR","upVR") # IMPORTANT:Names must match DEG sets exactly

# -------------------------
# 5. Load binomial thresholds
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
# 6. Load DEG count table
# -------------------------
deg_counts_dt <- fread(snakemake@input[["count_table"]])

# -------------------------
# 7. Extract robust DEG sets
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
# 8. Define TF-DEG mapping
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
# 9. Compute TF-DEG overlap
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

tf_deg_overlap <- get_tf_deg_overlap(tf_list, deg_sets, tf_deg_map)

# -------------------------
# 10. Collect candidate TFs
# -------------------------
tf_candidates <- unique(unlist(tf_deg_overlap, use.names = FALSE))

# -------------------------
# 11. Extract TF statistics
# -------------------------
tf_candidate_stats <- lapply(tf_list, function(dt) {
  
  dt <- dt[dt$`TAIR_ID` %in% tf_candidates, .(`TAIR_ID`, `adjusted p_value`)]
  
  dt <- unique(dt, by = "TAIR_ID")
  
  dt[order(`TAIR_ID`)]
})

# -------------------------
# 11. Build final TF table
# -------------------------
tf_tables <- Map(function(dt, nm) {
  
  dt <- as.data.table(dt)
  
  setnames(dt, c("adjusted p_value"), c(nm))
  
  dt
  
}, tf_candidate_stats, names(tf_candidate_stats))

tf_final_table <- Reduce(function(x, y) {
  merge(x, y, by = "TAIR_ID", all = TRUE)
}, tf_tables)

# -------------------------
# 12. Save result
# -------------------------
fwrite(tf_final_table, snakemake@output[["tf_enrichment"]])

