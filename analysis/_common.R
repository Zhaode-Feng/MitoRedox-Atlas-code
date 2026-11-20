inputs <- file.path(root, "inputs")
refs <- file.path(root, "reference")
comparisons <- list()

read_input <- function(path, ...) fread(file.path(inputs, path), ...)
read_ref <- function(path, ...) fread(file.path(refs, path), ...)
write_table <- function(data, path) {
  file <- file.path(output, path)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  fwrite(data, file, sep = "\t", na = "NA")
}
check_numbers <- function(label, actual, expected, tolerance = 1e-8) {
  stopifnot(length(actual) == length(expected))
  both <- is.finite(actual) & is.finite(expected)
  delta <- if (any(both)) max(abs(actual[both] - expected[both])) else 0
  finite_mismatch <- sum(is.finite(actual) != is.finite(expected))
  pass <- finite_mismatch == 0 && delta <= tolerance
  comparisons[[length(comparisons) + 1L]] <<- data.table(
    Check = label, Values = length(actual), Compared = sum(both),
    MaxAbsoluteDifference = delta, MissingnessMismatches = finite_mismatch,
    Tolerance = tolerance, Pass = pass
  )
  if (!pass) stop(label, " differs: max error ", signif(delta, 8),
                  ", finite-state mismatches ", finite_mismatch)
}
compare_by <- function(label, actual, expected, keys, columns, tolerance = 1e-8) {
  stopifnot(!anyDuplicated(actual[, ..keys]), !anyDuplicated(expected[, ..keys]))
  aligned <- actual[expected, on = keys, which = TRUE]
  if (anyNA(aligned)) stop(label, ": expected rows missing: ", sum(is.na(aligned)))
  for (column in columns) {
    check_numbers(paste(label, column, sep = " / "), actual[[column]][aligned],
                  expected[[column]], tolerance)
  }
}
zscore <- function(x) as.numeric(scale(x))
hc1_fit <- function(data, formula, term, cluster) {
  used <- data[complete.cases(data[, unique(c(all.vars(formula), cluster)), with = FALSE])]
  fit <- lm(formula, data = used)
  groups <- uniqueN(used[[cluster]])
  variance <- sandwich::vcovCL(fit, cluster = used[[cluster]], type = "HC1", cadjust = TRUE)
  estimate <- unname(coef(fit)[term])
  se <- sqrt(variance[term, term])
  width <- qt(.975, groups - 1L) * se
  data.table(Estimate = estimate, SE = se, CI_low = estimate - width, CI_high = estimate + width,
             P = 2 * pt(-abs(estimate / se), groups - 1L),
             N_sites = nrow(used), N_proteins = groups)
}
