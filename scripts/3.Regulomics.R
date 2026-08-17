
library(rvest)
library(stringr)
library(data.table)
library(tools)

# =========================
# 1. Load CisCross HTML/TXT reports
# =========================

report_files <- unname(unlist(snakemake@input[["reports"]]))

# Read HTML reports
report_html <- lapply(report_files, read_html)

# Convert reports into line-by-line text
report_lines <- lapply(report_html, function(x) {
  
  html_text(x) |>
    strsplit("\r\n") |>
    unlist()
})

report_names <- tools::file_path_sans_ext(basename(report_files))
report_phases <- sub("^ciscross_(ER|LR|VR)_.*$", "\\1", report_names)
if (any(!report_phases %in% c("ER", "LR", "VR")) || anyDuplicated(report_phases)) {
  stop("Expected exactly one CisCross report for each of ER, LR and VR")
}
names(report_lines) <- report_phases

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
# 3. Extract regulation tables
# =========================

regulation_tables <- lapply(report_lines, function(lines) {
  
  extract_ciscross_table(
    lines,
    start_pattern = "^TF Name.*Target type of TF regulation",
    end_pattern   = "^DONE TF_targets",
    first_only    = TRUE
  )
})

# =========================
# 4. Name tables by file
# =========================

names(regulation_tables) <- report_phases

# =========================
# 5. Combine + annotate
# =========================

regulome_dt <- rbindlist(
  Map(function(dt, label) {
    dt[, Response_type := label]
    dt
  }, regulation_tables, names(regulation_tables)),
  use.names = TRUE,
  fill = TRUE
)

# =========================
# 6. Compute connectivity
# =========================

connectivity_dt <- regulome_dt[
                   , .(Number_of_nodes = .N),
                   by = .(`TAIR ID`, Response_type)
][order(Response_type, -Number_of_nodes)]

# =========================
# 7. Summary statistics
# =========================

connectivity_summary <- connectivity_dt[
  , .(Connectivity_index = mean(Number_of_nodes)),
  by = Response_type
]

# =========================
# 8. Save results
# =========================

fwrite(connectivity_dt, snakemake@output[["connectivity_per_TF"]])
fwrite(connectivity_summary, snakemake@output[["connectivity_network"]])
