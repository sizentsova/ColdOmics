# 2a.CisCross_local.R -- offline replacement for the CisCross web service.
#
# Runs the vendored CisCross enrichment engine (resource/ciscross/) on the
# robust DEG lists from step 1 and writes one report per response phase in the
# exact tab-separated layout the web service produces, so that
# scripts/2.Cistromics.R and scripts/3.Regulomics.R consume it unchanged.
#
# Emitted sections (the only three ColdOmics parses):
#   TF Name ... UP                          -> DONE UP          (2.Cistromics.R)
#   TF Name ... DOWN                        -> DONE DOWN        (2.Cistromics.R)
#   TF Name ... Target type of TF regulation-> DONE TF_targets  (3.Regulomics.R)

library(data.table)

source(snakemake@input[["engine"]])

idx        <- load_index(snakemake@input[["index"]])
fdr_cutoff <- as.numeric(snakemake@params[["fdr"]])
phases     <- c("ER", "LR", "VR")

# Optional custom background (e.g. expressed genes only); "" = whole genome.
bg_file    <- snakemake@params[["background"]]
background <- NULL
if (is.character(bg_file) && length(bg_file) == 1L && nzchar(bg_file)) {
  if (!file.exists(bg_file))
    stop("cistromics.background file not found: ", bg_file)
  background <- unique(toupper(trimws(readLines(bg_file, warn = FALSE))))
  background <- background[nzchar(background)]
}

# -------------------------------------------------------------------------
# Input
# -------------------------------------------------------------------------
read_deg <- function(path) {
  dt <- fread(path, header = TRUE)
  if (!nrow(dt)) return(character(0))
  g <- unique(toupper(trimws(dt[[1]])))
  g[nzchar(g) & grepl("^AT[1-5CM]G[0-9]{5}$", g)]
}

# -------------------------------------------------------------------------
# Enrichment -> one row per TF, ranked
# -------------------------------------------------------------------------
# collapse_best_per_tf() keeps each TF's most significant peak set, matching the
# web service's TF-level table. The BH correction inside ciscross_enrich() runs
# across all peak sets, so it is re-applied across the retained TFs to give a
# TF-level adjusted p-value comparable to the web report.
enrich_tfs <- function(genes) {
  if (!length(genes)) return(data.table())
  res <- ciscross_enrich(genes, idx, background = background, correction = "BH")
  tf  <- collapse_best_per_tf(res)
  tf[, adj_p := p.adjust(p_value, "BH")]
  tf[adj_p <= fdr_cutoff][order(adj_p, p_value)]
}

# -------------------------------------------------------------------------
# TF -> TF regulatory network (the regulome)
# -------------------------------------------------------------------------
# This section is NOT a TF -> all-target-genes list. In the web report the
# nodes are exactly the enriched TFs that are themselves differentially
# expressed in this phase -- the ones the TF_Regulator section annotates as
# transcriptionally regulated (UA/DR/... rather than NTR) -- and every target
# is another node. Edge A -> B when B's promoter is bound by A's best peak set;
# self-loops (a TF binding its own promoter) are kept, as in the web report.
#
# 3.Regulomics.R counts rows per TF, so this is an out-degree: getting the node
# set wrong changes Connectivity_index by orders of magnitude. Reference values
# from resource/ciscross_example/ciscross_ER_1500p_0.05.txt: 16 nodes,
# 77 edges, mean out-degree 4.812.
tf_network <- function(tf_dt, deg_up, deg_down) {
  if (!nrow(tf_dt)) return(data.table())
  nodes <- tf_dt[tf_agi %chin% union(deg_up, deg_down)]
  if (!nrow(nodes)) return(data.table())
  nodes[, direction := fifelse(tf_agi %chin% deg_up, "U", "D")]

  edges <- rbindlist(lapply(seq_len(nrow(nodes)), function(i) {
    bound <- toupper(idx$incidence[[nodes$peakset_id[i]]])
    hit   <- nodes[tf_agi %chin% bound]
    if (!nrow(hit)) return(NULL)
    data.table(tf_name     = nodes$tf_name[i],
               tf_agi      = nodes$tf_agi[i],
               tf_dir      = nodes$direction[i],
               target_name = hit$tf_name,
               target_agi  = hit$tf_agi,
               target_dir  = hit$direction)
  }))
  if (!nrow(edges)) return(edges)
  unique(edges, by = c("tf_agi", "target_agi"))
}

# -------------------------------------------------------------------------
# Report writer (tab-separated, CRLF, as produced by the web service)
# -------------------------------------------------------------------------
fmt_p <- function(x) formatC(x, format = "E", digits = 5)

enrich_block <- function(tf_dt, direction) {
  head <- paste("TF Name", "TAIR ID", "TF Family", "p_value",
                "adjusted p_value", direction, sep = "\t")
  rows <- if (nrow(tf_dt))
    sprintf("%s\t%s\t%s\t%s\t%s\t", tf_dt$tf_name, tf_dt$tf_agi, tf_dt$family,
            fmt_p(tf_dt$p_value), fmt_p(tf_dt$adj_p))
  else character(0)
  c(head, rows, paste0("DONE ", direction, "\t\t\t\t"))
}

targets_block <- function(edges) {
  head <- paste("TF Name", "TAIR ID", "Type of TF regulation", "TF Target Name",
                "Target TAIR ID", "Target type of TF regulation", sep = "\t")
  # The regulation-type columns carry the node's DEG direction (U/D). The web
  # service appends an activator/repressor call (UA, UR, ...) that cannot be
  # derived from binding data alone, so it is not reproduced. No ColdOmics
  # script reads these columns.
  rows <- if (nrow(edges))
    sprintf("%s\t%s\t%s\t%s\t%s\t%s", edges$tf_name, edges$tf_agi,
            edges$tf_dir, edges$target_name, edges$target_agi,
            edges$target_dir)
  else character(0)
  c(head, rows, "DONE TF_targets\t\t\t\t")
}

write_report <- function(path, up_dt, down_dt, edges) {
  lines <- c(enrich_block(up_dt,   "UP"),
             enrich_block(down_dt, "DOWN"),
             targets_block(edges))
  con <- file(path, open = "wb")           # binary: exact CRLF on every OS
  on.exit(close(con))
  writeLines(lines, con, sep = "\r\n")
}

# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------
outputs <- unlist(snakemake@output)

for (ph in phases) {
  up_genes   <- read_deg(snakemake@input[[paste0("up", ph)]])
  down_genes <- read_deg(snakemake@input[[paste0("down", ph)]])

  up_tf   <- enrich_tfs(up_genes)
  down_tf <- enrich_tfs(down_genes)

  # The web report carries one network per phase covering the whole submission,
  # so nodes are drawn from TFs enriched in either direction.
  both_tf <- unique(rbindlist(list(up_tf, down_tf)), by = "tf_agi")
  edges   <- tf_network(both_tf, up_genes, down_genes)

  out <- grep(sprintf("ciscross_%s_", ph), outputs, value = TRUE)
  if (length(out) != 1L)
    stop("Expected exactly one output file for phase ", ph, "; got ",
         length(out))

  write_report(out, up_tf, down_tf, edges)

  n_nodes <- if (nrow(edges)) uniqueN(edges$tf_agi) else 0L
  message(sprintf(
    "%s: %d up-TFs, %d down-TFs | regulome %d nodes, %d edges (%d/%d DEGs)",
    ph, nrow(up_tf), nrow(down_tf), n_nodes, nrow(edges),
    length(up_genes), length(down_genes)))
}
