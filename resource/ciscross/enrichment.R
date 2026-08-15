# enrichment.R -- CisCross enrichment engine (R).
#
# Given a gene list and a precomputed incidence index (see build_index.R), test
# each TF peak set for over-representation of binding in the promoters of the
# input genes vs. the rest of the genome, using a one-sided Fisher's exact test,
# then correct across all peak sets (Benjamini-Hochberg FDR by default).
#
#   Foreground F = input genes present in the universe
#   Background   = universe \ F
#   Per peak set (bound = genes whose promoter overlaps a peak):
#            bound   not-bound
#     F        a          b        a = |F & bound|
#     Bg       c          d        c = |Bg & bound|
#     Fisher's exact test, alternative = "greater"

suppressPackageStartupMessages(library(data.table))

load_index <- function(path) readRDS(path)

# genes:       character vector of AGI ids (input list = foreground)
# index:       list from build_index.R / load_index()
# background:  optional character vector of AGI ids to use as the background
#              instead of the whole genome. The test universe becomes
#              foreground union background; ids not in the index universe are
#              dropped. NULL (default) = whole-genome background (the paper's
#              default and what the web service uses by default).
# Returns data.table sorted by ascending FDR.
ciscross_enrich <- function(genes, index, background = NULL,
                            correction = c("BH", "bonferroni", "none"),
                            alternative = "greater") {
  correction <- match.arg(correction)
  genes <- unique(toupper(trimws(genes)))
  universe <- toupper(index$universe)

  fg      <- intersect(genes, universe)
  dropped <- setdiff(genes, universe)
  n_fg    <- length(fg)
  if (n_fg == 0L) stop("None of the input genes are in the gene universe.")

  custom_bg <- !is.null(background)
  if (custom_bg) {
    bg_in <- intersect(unique(toupper(trimws(background))), universe)
    uni   <- union(fg, bg_in)          # restricted universe = foreground + bg
  } else {
    uni   <- universe                   # whole genome
  }
  n_uni <- length(uni)
  n_bg  <- n_uni - n_fg
  if (n_bg <= 0L)
    stop("Background is empty after intersecting the gene universe.")
  fg_set <- fg  # character vector; membership via %chin%

  res <- rbindlist(lapply(names(index$incidence), function(psid) {
    bound <- toupper(index$incidence[[psid]])
    a <- sum(fg_set %chin% bound)
    # genes bound within the test universe. With whole-genome background the
    # incidence is already a subset of the universe, so length(bound) is exact.
    n_bound <- if (custom_bg) sum(bound %chin% uni) else length(bound)
    b <- n_fg - a
    c <- n_bound - a
    d <- n_bg - c
    p <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                     alternative = alternative)$p.value
    # fold enrichment: (a/n_fg) / (n_bound/n_uni)
    fold <- if (n_bound > 0 && n_fg > 0) (a / n_fg) / (n_bound / n_uni) else NA_real_
    data.table(peakset_id = psid, a = a, b = b, c = c, d = d,
               fold_enrichment = fold, p_value = p)
  }))

  res[, FDR := switch(correction,
                      BH         = p.adjust(p_value, "BH"),
                      bonferroni = p.adjust(p_value, "bonferroni"),
                      none       = p_value)]

  res <- merge(res, index$meta, by = "peakset_id", all.x = TRUE, sort = FALSE)
  setorder(res, FDR, p_value)
  setcolorder(res, c("peakset_id", "tf_name", "tf_agi", "family", "dna_source",
                     "a", "b", "c", "d", "n_peaks", "fold_enrichment",
                     "p_value", "FDR"))
  attr(res, "n_foreground")     <- n_fg
  attr(res, "n_dropped")        <- length(dropped)
  attr(res, "dropped")          <- dropped
  attr(res, "n_universe")       <- n_uni
  attr(res, "n_background")     <- n_bg
  attr(res, "background_custom") <- custom_bg
  res[]
}

# Collapse the per-peak-set result to one row per TF (its most significant peak
# set), adding a `rank` column. This is the view comparable to the web service's
# TF-level table, since it merges col/colamp/replicate tracks of the same TF.
collapse_best_per_tf <- function(res) {
  by_tf <- res[order(FDR, p_value)][, .SD[1], by = tf_agi]
  setorder(by_tf, FDR, p_value)
  by_tf[, rank := .I]
  setcolorder(by_tf, c("rank", "tf_name", "tf_agi", "family", "peakset_id"))
  by_tf[]
}

# gene <-> enriched-TF membership (CisCross-Light style): for the top peak sets,
# which of the input genes have a bound promoter.
ciscross_light <- function(genes, index, res, fdr_cutoff = 0.05) {
  genes <- unique(toupper(trimws(genes)))
  top <- res[FDR <= fdr_cutoff]
  if (nrow(top) == 0L) return(data.table())
  rbindlist(lapply(top$peakset_id, function(psid) {
    bound <- toupper(index$incidence[[psid]])
    hit <- intersect(genes, bound)
    if (!length(hit)) return(NULL)
    data.table(gene_id = hit, peakset_id = psid,
               tf_name = top[peakset_id == psid, tf_name][1],
               tf_agi  = top[peakset_id == psid, tf_agi][1])
  }))
}
