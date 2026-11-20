# Chapter 5 / Fig. 5 and S11: sequence and structural associations.
universe <- read_input("sequence/site_inclusion_audit.tsv")
scope <- paste(samples$Tissue, samples$Age, sep = "_")
site_scores <- rbindlist(lapply(c("SNO", "SSG"), function(assay) {
  reference <- universe[Omics == assay]
  x <- matrices[[paste("Protein-adjusted", assay)]]
  keys <- sub("_M1$", "", rownames(x))
  keep <- grepl("_M1$", rownames(x)) & keys %in% reference$ExactSiteID
  grouped <- as.data.table(x[keep, , drop = FALSE])
  grouped[, ExactSiteID := keys[keep]]
  aggregate <- grouped[, lapply(.SD, function(v) median(v[is.finite(v)])),
                        by = ExactSiteID, .SDcols = samples$Sample]
  m <- as.matrix(aggregate[, samples$Sample, with = FALSE])
  stopifnot(setequal(reference$ExactSiteID, aggregate$ExactSiteID))
  totals <- rowSums(is.finite(m))
  per_group <- sapply(unique(scope), function(group) rowSums(is.finite(m[, scope == group, drop = FALSE])))
  eligible <- totals >= 16L & rowSums(per_group >= 2L) == 8L
  # Scale in the full high-confidence M1 universe before applying coverage and sequence filters.
  for (j in seq_len(ncol(m))) {
    v <- m[, j]
    good <- is.finite(v)
    if (sum(good) < 5L) { m[, j] <- NA_real_; next }
    center <- median(v[good])
    width <- mad(v[good], center = center, constant = 1.4826)
    if (!is.finite(width) || width <= 1e-8) width <- sd(v[good])
    m[, j] <- if (is.finite(width) && width > 1e-8) (v - center) / width else NA_real_
  }
  scores <- rowMeans(sapply(unique(scope), function(group) rowMeans(m[, scope == group, drop = FALSE],
                                                                    na.rm = TRUE)))
  chosen <- eligible & reference[match(aggregate$ExactSiteID, ExactSiteID), SequenceEligible]
  data.table(Omics = assay, ExactSiteID = aggregate$ExactSiteID[chosen],
             FiniteSamples = totals[chosen], SiteScore = scores[chosen])
}))
canonical <- read_input("sequence/canonical_site_scores.tsv")
stopifnot(nrow(site_scores) == nrow(canonical))
compare_by("Fig. 5 canonical SiteScore", site_scores, canonical, c("Omics", "ExactSiteID"),
           c("FiniteSamples", "SiteScore"))
write_table(site_scores, "05_sequence/site_scores.tsv")

sequence_env <- new.env()
source(file.path(root, "analysis", "_sequence_models.R"), local = sequence_env)
out_sequence <- rbindlist(list(
  sequence_env$boot[PTM == "SSG_minus_SNO" & B == 5000L],
  sequence_env$hc[PTM %in% c("SNO", "SSG")]
), fill = TRUE)
context <- new.env()
source(file.path(root, "analysis", "_sequence_context_models.R"), local = context)
eisenberg <- rbindlist(lapply(c("SNO", "SSG"), function(assay) {
  d <- copy(context$d[Omics == assay])
  d[, Feature_z := zscore(Eisenberg_PM6)]
  ans <- hc1_fit(d, SiteScore ~ Feature_z, "Feature_z", "Accession")
  setnames(ans, c("Estimate", "N_sites"), c("Beta", "N_records"))
  cbind(data.table(PTM = assay, Feature = "Eisenberg_PM6"), ans)
}))
out_sequence <- rbindlist(list(out_sequence, eisenberg), fill = TRUE)
expected <- read_ref("Figure5CDE_sequence_models.csv")
compare_by("Fig. 5C–E sequence effects", out_sequence, expected, c("PTM", "Feature"),
           c("Beta", "CI_low", "CI_high", "P", "N_records", "N_proteins"))
write_table(out_sequence, "05_sequence/submission_sequence_effects.tsv")

expected_context <- read_ref("FigureS11ABCD_context_models.csv")
model_columns <- c("Estimate", "SE", "CI_low", "CI_high", "P", "N_sites", "N_proteins")
compare_by("S11A within- and between-protein sequences", context$A,
           expected_context[Panels == "Supplementary Figure 11A"],
           c("Omics", "Feature", "Component"), model_columns)
compare_by("S11C compartment GRAVY", context$D[Feature == "GRAVY_PM6"],
           expected_context[Panels == "Supplementary Figure 11C" & Feature == "GRAVY_PM6"],
           c("Omics", "Feature", "Scope"), model_columns)

expected_positions <- read_ref("Figure5A_position_amino_acid_models.csv")
seq_data <- sequence_env$d
position_grid <- CJ(Omics = c("SNO", "SSG"), Position = c(-6:-1, 1:6),
                    AminoAcid = strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1L]])
positions <- rbindlist(lapply(seq_len(nrow(position_grid)), function(i) {
  r <- position_grid[i]
  d <- copy(seq_data[Omics == r$Omics])
  index <- 7L + r$Position
  d[, Present := as.integer(substr(Sequence13, index, index) == r$AminoAcid)]
  support <- d[, .(N = .N, G = uniqueN(Accession)), by = Present]
  if (nrow(support) < 2L || any(support$N < 10L) || any(support$G < 8L)) return(NULL)
  result <- hc1_fit(d, SiteScore ~ Present, "Present", "Accession")
  cbind(r[, .(Omics, Position, AminoAcid)], result)
}))
positions[, FDR_primary := p.adjust(P, "BH"), by = Omics]
positions[, FDR_global := p.adjust(P, "BH")]
stopifnot(nrow(positions) == nrow(expected_positions[is.finite(Estimate)]))
compare_by("Fig. 5A position effects", positions, expected_positions[is.finite(Estimate)],
           c("Omics", "Position", "AminoAcid"), c(model_columns, "FDR_primary", "FDR_global"))
write_table(positions, "05_sequence/position_models.tsv")

structure <- read_input("local_context/structure_features.tsv")
paired <- read_input("local_context/shared_structure_sites.tsv")
fit_structure <- function(d, feature, outcome) {
  d <- copy(d)
  d[, `:=`(Feature_z = zscore(get(feature)), pLDDT_z = zscore(SitePLDDT))]
  hc1_fit(d, reformulate(c("Feature_z", "pLDDT_z"), outcome), "Feature_z", "Accession")
}
struct_ref <- rbindlist(list(
  read_ref("Figure5BFG_structural_models.csv"),
  expected_context[Panels == "Supplementary Figure 11B"]
), fill = TRUE)
keys <- c("Panels", "Omics", "Feature")
struct_out <- rbindlist(lapply(seq_len(nrow(struct_ref)), function(i) {
  r <- struct_ref[i]
  if (r$AnalysisID == "CommonSitePreference") {
    result <- fit_structure(paired[PrimaryPLDDT70 == TRUE], r$Feature, "PTMPreference")
  } else {
    d <- structure[Omics == r$Omics & PrimaryPLDDT70 == TRUE & is.finite(SiteScore)]
    result <- fit_structure(d, r$Feature, "SiteScore")
  }
  cbind(r[, ..keys], result)
}))
compare_by("Fig. 5B,F,G and S11B structure models", struct_out, struct_ref, keys, model_columns)
write_table(struct_out, "05_sequence/structural_models.tsv")

# Compartment-specific structural slopes share the pooled within-PTM scaling.
structure[, Compartment := mapply(context$lookup, Accession, Gene, USE.NAMES = FALSE)]
structure_scope <- structure[PrimaryPLDDT70 == TRUE &
                             Compartment %in% c("Matrix", "MIM", "MOM", "IMS")]
compartments <- structure_scope[, .(N = .N, G = uniqueN(Accession)), by = .(Omics, Compartment)]
linear_features <- c("CysRelativeSASA", "CANeighbors_8A", "SpatialHydropathy_8A")
compartment_fits <- rbindlist(lapply(c("SNO", "SSG"), function(assay) {
  eligible <- compartments[Omics == assay & N >= 30L & G >= 20L, Compartment]
  d <- copy(structure_scope[Omics == assay & Compartment %in% eligible])
  d[, `:=`(pLDDT_z = zscore(SitePLDDT), Compartment = factor(Compartment))]
  rbindlist(lapply(linear_features, function(feature) {
    d[, Feature_z := zscore(get(feature))]
    full <- lm(SiteScore ~ Feature_z * Compartment + pLDDT_z, data = d)
    v <- sandwich::vcovCL(full, cluster = d$Accession, type = "HC1", cadjust = TRUE)
    index <- grep("Feature_z:Compartment|Compartment.*:Feature_z", names(coef(full)))
    b <- coef(full)[index]
    heterogeneity <- pchisq(drop(crossprod(b, solve(v[index, index, drop = FALSE], b))),
                            df = length(index), lower.tail = FALSE)
    rbindlist(lapply(eligible, function(compartment) {
      result <- hc1_fit(d[Compartment == compartment],
                        SiteScore ~ Feature_z + pLDDT_z, "Feature_z", "Accession")
      cbind(data.table(Omics = assay, Scope = compartment, Feature = feature,
                       HeterogeneityP = heterogeneity), result)
    }))
  }))
}))
cite <- c("Omics", "Scope", "Feature")
compartment_out <- rbindlist(list(compartment_fits, context$D[Feature == "GRAVY_PM6"]), fill = TRUE)
compartment_out[, FDR_primary := p.adjust(P, "BH"), by = Omics]
compartment_out[, FDR_global := p.adjust(P, "BH")]
compare_by("S11C compartment slopes and Wald tests", compartment_out,
           expected_context[Panels == "Supplementary Figure 11C"], cite,
           c(model_columns, "HeterogeneityP", "FDR_primary", "FDR_global"))
write_table(compartment_out, "05_sequence/compartment_models.tsv")

external <- read_input("local_context/macrophage_ssg.tsv")
external[, Outcome_z := zscore(log1p(SSG_Occ))]
features <- c("GRAVY_PM6", "RSA_num", "pCoil_num")
model <- rbindlist(lapply(features, function(feature) {
  d <- if (feature == "GRAVY_PM6") copy(external[StrictSequence13 == TRUE]) else copy(external)
  if (feature == "GRAVY_PM6") {
    kd <- c(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-.4,H=-3.2,I=4.5,K=-3.9,
            L=3.8,M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-.8,T=-.7,V=4.2,W=-.9,Y=-1.3)
    d[, GRAVY_PM6 := vapply(strsplit(Sequence13, "", fixed = TRUE),
                           function(aa) mean(kd[aa[-7L]]), numeric(1))]
  }
  d[, Feature_z := zscore(get(feature))]
  result <- hc1_fit(d, Outcome_z ~ Feature_z + factor(CoarseCompartment), "Feature_z", "Protein")
  cbind(data.table(Feature = feature), result)
}))
model[, FDR_primary := p.adjust(P, "BH")]
compare_by("S11D macrophage SSG", model, expected_context[Panels == "Supplementary Figure 11D"],
           "Feature", c(model_columns, "FDR_primary"))
write_table(model, "05_sequence/external_models.tsv")
