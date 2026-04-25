library(biomaRt)
library(data.table)

robustTF <- fread(snakemake@input[["tf_enrichment"]])

# Connect to Ensembl Plants
mart <- useMart(
  biomart = "plants_mart",
  dataset = "athaliana_eg_gene",
  host = "https://plants.ensembl.org"
)

# Retrieve gene coordinates
robustTF_coord <- getBM(
  attributes = c("chromosome_name","start_position","end_position",
                 "ensembl_gene_id"),
  filters = c("ensembl_gene_id"),
  values = robustTF$TAIR_ID,
  mart = mart
)

# Write BED file
write.table(
  robustTF_coord,
  file = snakemake@output[["tf_bed"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
