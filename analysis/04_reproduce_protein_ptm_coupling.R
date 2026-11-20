# Chapter 4 / Fig. 4: accession-level resampling of matched measurements.
data <- read_input("coupling/matched_observations.tsv")
age_effect <- read_input("coupling/age_effect_input.tsv")
expected <- read_ref("Figure4_protein_coupling.csv")
data[, c("residual_Protein", "residual_SNO", "residual_SSG") :=
       lapply(.SD, function(x) x - mean(x)),
     by = .(Tissue, Accession, Age), .SDcols = c("Protein", "SNO", "SSG")]
rho <- function(x) c(cor(x[, 1L], x[, 2L], method = "spearman"),
                      cor(x[, 1L], x[, 3L], method = "spearman"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
results <- lapply(which(expected$Analysis != "age_effect"), function(i) {
  row <- expected[i]
  tissue <- sub("^Muscle$|^Skeletal muscle$", "Mus", row$Tissue)
  d <- data[Tissue == tissue & Age == row$Age]
  columns <- c("Protein", "SNO", "SSG")
  if (row$Analysis == "residual") columns <- paste0("residual_", columns)
  ids <- sort(unique(d$Accession), method = "radix")
  blocks <- split(seq_len(nrow(d)), factor(d$Accession, levels = ids))
  x <- as.matrix(d[, ..columns])
  observed <- rho(x)
  stopifnot(length(ids) == row$n_proteins)
  set.seed(row$seed)
  boot <- replicate(5000L, {
    picked <- sample.int(length(ids), length(ids), replace = TRUE)
    correlations <- rho(x[unlist(blocks[picked], use.names = FALSE), , drop = FALSE])
    correlations[1L] - correlations[2L]
  })
  ci <- quantile(boot, c(.025, .975), type = 7, names = FALSE)
  result <- row[, .(Analysis, Tissue, Age)]
  result[, `:=`(rho_SNO = observed[1], rho_SSG = observed[2],
                 delta_rho = observed[1] - observed[2], CI_low = ci[1], CI_high = ci[2])]
  result
})
# Panels E–F retain the original fixed-rank multinomial bootstrap.
# One RNG stream is consumed in the order of the archived tissue input.
columns <- c("Protein", "SNO", "SSG")
reconstructed <- data[, lapply(.SD, function(x) mean(x[Age == "Old"]) - mean(x[Age == "Young"])),
                      by = .(Tissue, Accession), .SDcols = columns]
setnames(reconstructed, columns, c("P_age", "SNO_age", "SSG_age"))
compare_by("Fig. 4E age input effects", reconstructed, age_effect, c("Tissue", "Accession"),
           c("P_age", "SNO_age", "SSG_age"))
set.seed(20260914L)
age_results <- lapply(unique(age_effect$Tissue), function(tissue) {
  d <- age_effect[Tissue == tissue]
  x <- as.matrix(d[, .(P_age, SNO_age, SSG_age)])
  observed <- rho(x)
  ranks <- apply(x, 2L, rank)
  n <- nrow(ranks)
  a <- ranks[, 1L]; b <- ranks[, 2L]; c <- ranks[, 3L]
  moments <- cbind(a, b, c, a * a, b * b, c * c, a * b, a * c)
  draws <- rmultinom(5000L, n, rep(1 / n, n))
  s <- crossprod(moments, draws)
  va <- s[4L, ] - s[1L, ]^2 / n
  vb <- s[5L, ] - s[2L, ]^2 / n
  vc <- s[6L, ] - s[3L, ]^2 / n
  db <- (s[7L, ] - s[1L, ] * s[2L, ] / n) / sqrt(va * vb) -
        (s[8L, ] - s[1L, ] * s[3L, ] / n) / sqrt(va * vc)
  ci <- quantile(db, c(.025, .975), type = 7, names = FALSE)
  data.table(Analysis = "age_effect", Tissue = tissue, Age = NA_character_,
             rho_SNO = observed[1L], rho_SSG = observed[2L],
             delta_rho = observed[1L] - observed[2L], CI_low = ci[1L], CI_high = ci[2L])
})
out <- rbindlist(c(results, age_results))
expected[is.na(Age) | Age == "", Age := NA_character_]
compare_by("Fig. 4 paired Spearman and bootstrap", out, expected,
           c("Analysis", "Tissue", "Age"), c("rho_SNO", "rho_SSG", "delta_rho", "CI_low", "CI_high"))
write_table(out, "04_coupling/protein_ptm_spearman.tsv")

site <- read_input("coupling/shared_site_samples.tsv")
site[, c("SNO_residual", "SSG_residual") := lapply(.SD, function(x) x - mean(x)),
     by = .(Tissue, SiteKey, Age), .SDcols = c("SNO", "SSG")]
shared <- site[, .(baseline = cor(SNO, SSG, method = "spearman"),
                   residual = cor(SNO_residual, SSG_residual, method = "spearman")), by = Tissue]
shared <- melt(shared, id.vars = "Tissue", variable.name = "Analysis", value.name = "rho")
ref <- read_ref("FigureS10_shared_site_coupling.csv")
ref[, Tissue := sub("^Muscle$|^Skeletal muscle$", "Mus", Tissue)]
compare_by("S10 sample-level shared-site Spearman", shared, ref[Analysis %in% c("baseline", "residual")],
           c("Tissue", "Analysis"), "rho")
write_table(shared, "04_coupling/shared_site_spearman.tsv")

site_age <- read_input("coupling/shared_site_age_effects.tsv")
age_shared <- site_age[, {
  fit <- cor.test(SNO_logFC, SSG_logFC, method = "spearman", exact = FALSE)
  .(N = .N, SpearmanRho = unname(fit$estimate), PValue = fit$p.value,
    SNO_MedianLogFC = median(SNO_logFC), SSG_MedianLogFC = median(SSG_logFC))
}, by = Tissue]
age_shared[, FDR_BH := p.adjust(PValue, "BH")]
compare_by("S10C shared-site age effects", age_shared, ref[Analysis == "age_effect"], "Tissue",
           c("N", "SpearmanRho", "PValue", "SNO_MedianLogFC", "SSG_MedianLogFC", "FDR_BH"))
write_table(age_shared, "04_coupling/shared_site_age_effects.tsv")
