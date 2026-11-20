# Run from the release directory with Rscript run_all.R [output-directory].
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_arg)) {
  script <- sub("^--file=", "", file_arg[[1]])
} else {
  # The lab's Windows launcher supplies the script path as its first trailing argument.
  script <- args[[1]]
  args <- args[-1]
}
root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
output <- if (length(args)) args[[1]] else file.path(root, "results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
protected <- normalizePath(file.path(root, c("inputs", "reference", "analysis")), winslash = "/")
if (any(output == protected | startsWith(paste0(output, "/"), paste0(protected, "/")))) {
  stop("Choose an output folder outside analysis, inputs and reference.")
}
required <- c("data.table", "limma", "sandwich", "igraph")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages(library(data.table))
source(file.path(root, "analysis", "_common.R"))
for (step in list.files(file.path(root, "analysis"), pattern = "^[0-9].*[.]R$")) {
  message("Running ", step)
  source(file.path(root, "analysis", step), local = new.env())
}
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
write_table(do.call(rbind, comparisons), "validation.tsv")
message("Reproduction checks passed. Results: ", output)
