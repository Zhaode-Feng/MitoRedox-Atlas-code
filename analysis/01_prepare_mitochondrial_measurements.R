# Chapter 1 / Fig. 1: median centre log intensities and adjust matched PTM features.
frozen <- readRDS(file.path(inputs, "quantified_measurements.rds"))
samples <- frozen$sample_info
layers <- data.table(
  Layer = c("Protein abundance", "Unadjusted SNO", "Unadjusted SSG",
            "Protein-adjusted SNO", "Protein-adjusted SSG"),
  Type = c("Protein", "SNO", "SSG", "SNO", "SSG"),
  Matrix = c("protein_center", "sno_center", "ssg_center", "sno_corrected", "ssg_corrected")
)
meta <- list(Protein = frozen$protein_meta, SNO = frozen$sno_meta, SSG = frozen$ssg_meta)
matrices <- setNames(lapply(seq_len(nrow(layers)), function(i) {
  x <- frozen[[layers$Matrix[i]]]
  stopifnot(nrow(x) == nrow(meta[[layers$Type[i]]]))
  rownames(x) <- meta[[layers$Type[i]]]$FeatureID
  x[, samples$Sample, drop = FALSE]
}), layers$Layer)
reconstructed <- lapply(c("Protein", "SNO", "SSG"), function(assay) {
  raw <- frozen[[paste0(tolower(assay), "_raw")]][, samples$Sample, drop = FALSE]
  normalised <- log2(raw)
  centres <- apply(normalised, 2L, median, na.rm = TRUE)
  centres[!is.finite(centres)] <- 0
  normalised <- sweep(normalised, 2L, centres, "-")
  normalised[, frozen$exclusions[[assay]]] <- NA_real_
  rownames(normalised) <- meta[[assay]]$FeatureID
  normalised
})
for (i in 2:3) {
  parent_rows <- match(meta[[layers$Type[i]]]$Accession, meta$Protein$Accession)
  reconstructed[[i + 2L]] <- reconstructed[[i]] - reconstructed[[1L]][parent_rows, , drop = FALSE]
}
for (i in seq_along(matrices)) {
  check_numbers(paste("Preprocessing", names(matrices)[i]), as.vector(reconstructed[[i]]),
                as.vector(matrices[[i]]), 1e-10)
}
matrices[] <- reconstructed
tissues <- c("Brain", "Heart", "Liver", "Mus")
stopifnot(nrow(samples) == 32L, !anyDuplicated(samples$Sample))
for (i in seq_len(nrow(layers))) {
  metadata <- meta[[layers$Type[i]]]
  stopifnot(!anyDuplicated(metadata$FeatureID), "IsMito" %in% names(metadata))
}
count <- rbindlist(lapply(seq_len(3L), function(i) {
  m <- matrices[[i]]
  mito <- meta[[layers$Type[i]]]$IsMito
  data.table(Assay = layers$Type[i], AllFeatures = nrow(m), Mitochondrial = sum(mito))
}))
check_numbers("Fig. 1B mitochondrial feature counts", count$Mitochondrial, c(1256, 3118, 2913), 0)
write_table(count, "01_resource/mitochondrial_feature_counts.tsv")
coverage <- rbindlist(lapply(seq_len(3L), function(i) {
  rbindlist(lapply(tissues, function(tissue) {
    observed <- rowSums(is.finite(matrices[[i]][, samples$Tissue == tissue, drop = FALSE])) > 0L
    data.table(Assay = layers$Type[i], Tissue = if (tissue == "Mus") "Muscle" else tissue,
               `All detected features` = sum(observed),
               `High-confidence mitochondrial features` = sum(observed & meta[[layers$Type[i]]]$IsMito))
  }))
}))
compare_by("Fig. 1 tissue coverage", coverage, read_ref("Tissue_coverage.csv"),
           c("Assay", "Tissue"), c("All detected features", "High-confidence mitochondrial features"), 0)
write_table(coverage, "01_resource/tissue_coverage.tsv")
assign("samples", samples, envir = parent.env(environment()))
assign("layers", layers, envir = parent.env(environment()))
assign("meta", meta, envir = parent.env(environment()))
assign("matrices", matrices, envir = parent.env(environment()))
assign("tissues", tissues, envir = parent.env(environment()))
