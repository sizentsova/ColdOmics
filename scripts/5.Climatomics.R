library(data.table)
library(readxl)

# -----------------------------
# 1. Load data
# -----------------------------
mutation_dt <- fread(snakemake@input[["alt"]])
                     
temperature_dt <- as.data.table(
  read_xlsx(snakemake@input[["CHELSA_temperature_data"]], sheet = 1)
)

setnames(
  temperature_dt,
  old = c("EcotypeID", "Winter_tas_min", "Spring_tas_min", "Summer_tas_min", "Autumn_tas_min"),
  new = c("EcotypeID", "winter", "spring", "summer", "autumn")
)

#REMOVE temperature anomalous ecotypes
bad_ids <- c("7183", "8424")

mutation_dt <- mutation_dt[!EcotypeID %in% bad_ids]
temperature_dt <- temperature_dt[!EcotypeID %in% bad_ids]

mutation_dt <- mutation_dt[EcotypeID %in% temperature_dt$EcotypeID] #remove ecotypes of bad quality

# -----------------------------
# 2. Long format (seasons only)
# -----------------------------
long_dt <- melt(
  temperature_dt,
  id.vars = "EcotypeID",
  measure.vars = c("winter", "spring", "summer", "autumn"),
  variable.name = "season",
  value.name = "temperature"
)

# -----------------------------
# 3. Permutation test
# -----------------------------
perm_test <- function(temp, grp, B = 9999) {
  obs <- mean(temp[grp == "cold_TF"]) - mean(temp[grp == "control"])
  
  null <- replicate(B, {
    shuffled <- sample(grp)
    mean(temp[shuffled == "cold_TF"]) - mean(temp[shuffled == "control"])
  })
  
  list(
    obs = obs,
    p_upper = (sum(null >= obs) + 1) / (B + 1),
    p_lower = (sum(null <= obs) + 1) / (B + 1),
    se_null = sd(null) / sqrt(B),
    mean_control = mean(temp[grp == "control"]),
    mean_cold = mean(temp[grp == "cold_TF"])
  )
}

# -----------------------------
# 4. Helper function to run MC
# -----------------------------
run_mc <- function(dt, cold_ids, label, level_type) {
  
  tmp <- copy(dt)
  
  tmp[, group := fifelse(
    EcotypeID %in% cold_ids, "cold_TF", "control")
  ]
  
  # keep only relevant ecotypes
  tmp <- tmp[!is.na(group)]
  
  tmp[, {
    res <- perm_test(temperature, group)
    
    .(
      mean_control = res$mean_control,
      mean_cold_TF = res$mean_cold,
      observed_difference = res$obs,
      p_upper = res$p_upper,
      p_lower = res$p_lower,
      se_null = res$se_null,
      label = label,
      level = level_type
    )
    
  }, by = season]
}

# -----------------------------
# 5. GLOBAL analysis (all TFs)
# -----------------------------
all_cold_ids <- unique(mutation_dt$EcotypeID)

global_results <- run_mc(long_dt, all_cold_ids, label = "ALL_TF", level_type = "global")

# -----------------------------
# 6. GENE-wise analysis
# -----------------------------
gene_list <- split(mutation_dt, by = "TAIR_ID", keep.by = FALSE)

gene_results <- rbindlist(lapply(names(gene_list), function(gene)  {
  
  cold_ids <- unique(gene_list[[gene]][["EcotypeID"]])
  
  run_mc(long_dt, cold_ids,label = gene, level_type = "gene")
  
}))


# -----------------------------
# 7. Combine everything
# -----------------------------
final_results <- rbindlist(list(global_results, gene_results))

fwrite(final_results, snakemake@output[["mc"]])
