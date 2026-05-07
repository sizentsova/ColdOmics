
library(data.table)

# =========================
# 1. Load data
# =========================

tf_cis   <- fread(snakemake@input[["tf_cis"]])
tf_reg   <- fread(snakemake@input[["tf_reg"]])
tf_gen   <- fread(snakemake@input[["tf_gen"]])
tf_clima <- fread(snakemake@input[["tf_clima"]])

# =========================
# 2. Filter climatomics
# =========================

tf_clima <- tf_clima[
  (p_upper < 0.05 | p_lower < 0.05) &
    label != "ALL_TF"
]

# =========================
# 3. Standardize column names
# =========================

setnames(tf_reg, "TAIR ID", "TAIR_ID")

# =========================
# 4. Extract TF sets
# =========================

cis_ids   <- unique(tf_cis$TAIR_ID)
reg_ids   <- unique(tf_reg$TAIR_ID)
gen_ids   <- unique(tf_gen$TAIR_ID)
clima_ids <- unique(tf_clima$label)

# =========================
# 5. Build unified TF table
# =========================

all_tf <- unique(c(cis_ids, reg_ids, gen_ids, clima_ids))

tf_summary <- data.table(TAIR_ID = all_tf)

tf_summary[, transcriptomics := TAIR_ID %in% cis_ids]
tf_summary[, cistromics      := TAIR_ID %in% cis_ids]
tf_summary[, regulomics      := TAIR_ID %in% reg_ids]
tf_summary[, genomics        := TAIR_ID %in% gen_ids]
tf_summary[, climatomics     := TAIR_ID %in% clima_ids]

# =========================
# 6. Add gene names (BioMart)
# =========================

ara_genes <- readRDS("resource/Arabidopsis_genes.rds")
ara_genes <- ara_genes[ ,.(TAIR_ID, GeneName)]

# =========================
# 7. Merge annotation
# =========================

final_tf_table <- merge(
  ara_genes,
  tf_summary,
  by = "TAIR_ID",
  all = TRUE
)

# =========================
# 9. Save
# =========================

fwrite(final_tf_table, snakemake@output[["tf_int"]])