

suppressPackageStartupMessages({
  library('vcfR')
  library('data.table')
})


# =========================================================
# Load VCF and BED annotation
# =========================================================

vcf <- read.vcfR(snakemake@input[["vcf"]])  # (rename high_variants → vcf)
bed <- read.table(snakemake@input[["bed"]])

setnames(bed, c("chr", "start", "end", "TAIR_ID"))

# FIX section of VCF (variant-level metadata)
vcf_fix <- as.data.table(vcf@fix)   # (rename fix → vcf_fix)

# =========================================================
# 1. Parse INFO field (split by "|")
# =========================================================

info_split <- strsplit(vcf_fix$INFO, "\\|")  # (rename info_list → info_split)

# remove empty tokens
info_split <- lapply(info_split, function(x) x[x != ""])

# =========================================================
# 2. Find HIGH annotation positions in INFO
# =========================================================

high_idx <- lapply(info_split, function(x) which(x == "HIGH"))  # (rename index → high_idx)

# =========================================================
# 3. Expand region around HIGH hits (x-1, x, x+2)
# =========================================================

expanded_idx <- lapply(high_idx, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(integer(0))
  c(rbind(x - 1, x, x + 2))
})  # (rename index_ls → expanded_idx)

# =========================================================
# 4. Filter INFO per variant using expanded index
# =========================================================

info_filtered <- Map(function(info, idx) {
  info <- info[idx]
  unique(info)
}, info_split, expanded_idx)

# overwrite INFO column with filtered version
vcf_fix[, INFO := vapply(info_filtered, paste, collapse = "|", FUN.VALUE = character(1))]

# =========================================================
# 5. Build annotation table (variant metadata)
# =========================================================

variant_dt <- data.table(   # (rename dt → variant_dt)
  TAIR_ID = vapply(info_filtered, function(x) x[3], character(1)),
  effect  = vapply(info_filtered, function(x) x[1], character(1)),
  variant_index = seq_along(info_filtered)   # (rename ID_variant → variant_index)
)

# keep only variants present in BED region
variant_dt <- variant_dt[TAIR_ID %in% bed$TAIR_ID]

# =========================================================
# 6. Subset VCF and genotype matrix
# =========================================================

vcf_fix <- vcf_fix[variant_dt$variant_index]
gt.matrix <- extract.gt(vcf)[variant_dt$variant_index, ]

# =========================================================
# 7. Filter variants by ALT frequency (≤ 5%)
# =========================================================

alt_count <- rowSums(gt.matrix == "1|1", na.rm = TRUE)
n_samples <- ncol(gt.matrix)

keep <- alt_count <= floor(0.05 * n_samples)

vcf_fix <- vcf_fix[keep]
gt.matrix <- gt.matrix[keep, ]
variant_dt <- variant_dt[keep]

# =========================================================
# 8. Extract genotype groups per variant
# =========================================================

samples_alt <- apply(gt.matrix == "1|1", 1, function(x) colnames(gt.matrix)[which(x)])  # (samples_11)
samples_ref <- apply(gt.matrix == "0|0", 1, function(x) colnames(gt.matrix)[which(x)])  # (samples_00)
samples_na  <- apply(is.na(gt.matrix), 1, function(x) colnames(gt.matrix)[which(x)])

# genotype counts
n_alt <- lengths(samples_alt)
n_ref <- lengths(samples_ref)
n_na  <- lengths(samples_na)

# =========================================================
# 9. Build variant coordinate table
# =========================================================

variant_parts <- tstrsplit(rownames(gt.matrix), "_", fixed = TRUE)

af_dt <- data.table(   # (rename AF_df → af_dt)
  variant_name = rownames(gt.matrix),
  chr = variant_parts[[1]],
  pos = variant_parts[[2]],
  variant_index = variant_parts[[3]],   # (rename ID_variant → variant_id)
  alt_count = n_alt,
  ref_count = n_ref,
  na_count  = n_na
)

# =========================================================
# 10. Convert genotype lists into long format
# =========================================================

make_genotype_dt <- function(sample_list, gt_label, af_dt, vcf_fix) {
  rbindlist(lapply(seq_along(sample_list), function(i) {
    
    if (length(sample_list[[i]]) == 0) return(NULL)
    
    data.table(
      EcotypeID = sample_list[[i]],   # (rename ID → sample_id)
      chr = af_dt$chr[i],
      pos = af_dt$pos[i],
      variant_index = af_dt$variant_index[i],
      REF = vcf_fix$REF[i],
      ALT = vcf_fix$ALT[i],
      genotype = gt_label
    )
  }), fill = TRUE)
}

# build genotype tables
gen_alt <- make_genotype_dt(samples_alt, "1|1", af_dt, vcf_fix)
gen_ref <- make_genotype_dt(samples_ref, "0|0", af_dt, vcf_fix)
gen_na <- make_genotype_dt(samples_na, "NA", af_dt, vcf_fix)

# =========================================================
# 11. Merge with variant metadata
# =========================================================

variant_dt[, variant_index := as.character(variant_index)]

alt_final <- merge(variant_dt, gen_alt, by = "variant_index", sort = FALSE)
ref_final <- merge(variant_dt, gen_ref, by = "variant_index", sort = FALSE)
na_final  <- merge(variant_dt, gen_na,  by = "variant_index", sort = FALSE)

write.csv(alt_final, snakemake@output[["alt"]], row.names = F)
write.csv(ref_final, snakemake@output[["ref"]], row.names = F)
write.csv(na_final, snakemake@output[["na"]], row.names = F)
