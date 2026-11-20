# Chapters 2 and 3 / Figs. 2–3: fit the original limma families, then match display keys.
fit_contrast <- function(matrix, group_a, n_a, n_b, eligible_rest) {
  a <- rowSums(is.finite(matrix[, group_a, drop = FALSE]))
  b <- rowSums(is.finite(matrix[, !group_a, drop = FALSE]))
  eligible <- a >= n_a & b >= n_b & eligible_rest
  design <- cbind(Focal = as.numeric(group_a), Reference = as.numeric(!group_a))
  adjusted <- effect <- p <- rep(NA_real_, nrow(matrix))
  if (sum(eligible) >= 2L) {
    fit <- limma::lmFit(matrix[eligible, , drop = FALSE], design)
    fit <- limma::eBayes(limma::contrasts.fit(fit, c(1, -1)))
    effect[eligible] <- fit$coefficients[, 1]
    p[eligible] <- fit$p.value[, 1]
    adjusted[eligible] <- p.adjust(p[eligible], method = "BH")
  }
  data.table(Feature_ID = rownames(matrix), Testable = eligible,
             N_A = a, N_B = b, log2_fold_change = effect, P_value = p, BH_adjusted_P = adjusted)
}

for (analysis in c("Tissue", "Ageing")) {
  expected <- read_ref(paste0(analysis, "_differential_results.csv"))
  tissue_fit <- analysis == "Tissue"
  pieces <- vector("list", length(tissues) * nrow(layers))
  for (i in seq_len(nrow(layers))) {
    layer <- layers$Layer[i]
    annot <- meta[[layers$Type[i]]]
    universe <- if (tissue_fit || startsWith(layer, "Protein-adjusted")) annot$IsMito else rep(TRUE, nrow(annot))
    x <- matrices[[layer]][universe, , drop = FALSE]
    for (j in seq_along(tissues)) {
      tissue <- tissues[j]
      in_tissue <- samples$Tissue == tissue
      if (tissue_fit) {
        other_tissues <- setdiff(tissues, tissue)
        covered_rest <- rowSums(vapply(other_tissues, function(t) {
          rowSums(is.finite(x[, samples$Tissue == t, drop = FALSE])) > 0L
        }, logical(nrow(x))))
        result <- fit_contrast(x, in_tissue, 4L, 6L, covered_rest >= 2L)
      } else {
        result <- fit_contrast(x[, in_tissue, drop = FALSE],
                               samples$Age[in_tissue] == "Old", 2L, 2L, TRUE)
      }
      result[, `:=`(Tissue = if (tissue == "Mus") "Skeletal muscle" else tissue, Measurement_layer = layer)]
      pieces[[(i - 1L) * length(tissues) + j]] <- result
    }
  }
  full <- rbindlist(pieces)
  keys <- c("Tissue", "Measurement_layer", "Feature_ID")
  columns <- c("log2_fold_change", "P_value", "BH_adjusted_P")
  compare_by(analysis, full, expected, keys, columns)
  fitted <- full[expected[, ..keys], on = keys]
  fitted[, Meets_nominal_P_and_effect_threshold := Testable &
           P_value < .05 & abs(log2_fold_change) >= .5]
  check_numbers(paste(analysis, "joint threshold"),
                as.integer(fitted$Meets_nominal_P_and_effect_threshold),
                as.integer(expected$Meets_nominal_P_and_effect_threshold), 0)
  write_table(fitted, paste0("02_differential/", tolower(analysis), "_submission_contrasts.tsv"))
  if (!tissue_fit) {
    n <- fitted[Meets_nominal_P_and_effect_threshold &
                  Measurement_layer %in% c("Protein-adjusted SNO", "Protein-adjusted SSG"),
                .(N = uniqueN(Feature_ID)), by = Measurement_layer][order(Measurement_layer)]
    check_numbers("Fig. 3C unique age-associated SNO and SSG features", n$N, c(410, 319), 0)
  }
}

# Fig. 3E keeps a separate cohort of testable, matched observations.
membership <- read_input("ageing/mitochondrial_membership.tsv")
matched <- read_input("ageing/matched_effects.tsv")
joined <- merge(matched, membership, by = c("Tissue", "PTM_type", "FeatureID"))
stopifnot(nrow(joined) == nrow(matched))
responsive <- function(p, effect) is.finite(p) & p < .05 & is.finite(effect) & abs(effect) >= .5
retained <- joined[Mito_HighConf & matched_testable &
                     (responsive(protein_P, protein_logFC) |
                        responsive(raw_ptm_P, raw_ptm_logFC) |
                        responsive(corrected_ptm_P, corrected_ptm_logFC))]
retained[, class := fcase(delta_corrected_ptm >= .5, "RELATIVE_PTM_GAIN",
                          delta_corrected_ptm <= -.5, "RELATIVE_PTM_LOSS",
                          default = "PROTEOME_ALIGNED")]
composition <- retained[, .N, by = .(Tissue, PTM_type, class)]
composition[, denominator := sum(N), by = .(Tissue, PTM_type)]
composition[, fraction := N / denominator]
reference <- read_ref("Ageing_matched_layer_composition.csv")
reference[, Tissue := sub("^Muscle$|^Skeletal muscle$", "Mus", Tissue)]
compare_by("Fig. 3E matched three-layer composition", composition, reference,
           c("Tissue", "PTM_type", "class"), c("N", "denominator", "fraction"))
write_table(composition, "02_differential/matched_layer_composition.tsv")
