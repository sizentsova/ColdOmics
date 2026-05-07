
library(data.table)

robustTF <- fread(snakemake@input[["tf_enrichment"]])

#load Arabidopsis genes with coordinates

ara_genes <- readRDS("resource/Arabidopsis_genes.rds")

#Overlap Arabidopsis genes with robust TF

robustTF_coord <- ara_genes[ara_genes$TAIR_ID %in% robustTF$TAIR_ID, ]
robustTF_coord <- robustTF_coord[ ,.(chromosome_name, start_position, end_position, TAIR_ID)]

# Write BED file
write.table(
  robustTF_coord,
  file = snakemake@output[["tf_bed"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
