# Original sequence regression and resampling logic. Paths adapted for this release.
suppressPackageStartupMessages(library(sandwich))
out <- file.path(output, "05_sequence")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read <- function(n) read_input(file.path('sequence', n), na.strings = c('NA', ''))
write <- function(x, n) write_table(x, file.path('05_sequence', n))
d <- read('sequence_master.tsv')
canonical <- read('canonical_site_scores.tsv')
stopifnot(!anyDuplicated(d[, .(Omics, ExactSiteID)]),
          !anyDuplicated(canonical[, .(Omics, ExactSiteID)]),
          nrow(d) == 1736L, all(!is.na(d$Accession)),
          all(nchar(d$Sequence13) == 13L), all(substr(d$Sequence13, 7, 7) == 'C'),
          all(grepl('^[ACDEFGHIKLMNPQRSTVWY]{13}$', d$Sequence13)))
check <- merge(d[, .(Omics, ExactSiteID, SiteScore)],
               canonical[, .(Omics, ExactSiteID, Canonical = SiteScore)],
               by = c('Omics', 'ExactSiteID'), all = TRUE)
stopifnot(nrow(check) == nrow(d), all(is.finite(check$Canonical)),
          max(abs(check$SiteScore - check$Canonical)) < 1e-10)
kd <- c(A=1.8,C=2.5,D=-3.5,E=-3.5,F=2.8,G=-.4,H=-3.2,I=4.5,K=-3.9,
        L=3.8,M=1.9,N=-3.5,P=-1.6,Q=-3.5,R=-4.5,S=-.8,T=-.7,V=4.2,W=-.9,Y=-1.3)
bases <- c('GRAVY', 'HydrophobicFraction', 'PolarFraction', 'ChargedFraction', 'ProFraction')
calc <- function(s, w) {
  aa <- strsplit(s, '', fixed=TRUE)[[1]][c((7-w):6, 8:(7+w))]
  c(mean(kd[aa]), mean(kd[aa] > 1), mean(kd[aa] < 0),
    mean(aa %in% c('D','E','K','R')), mean(aa == 'P'))
}
feature_qc <- list()
for(w in c(4L, 6L)) {
  mat <- t(vapply(d$Sequence13, calc, numeric(5), w=w))
  for(j in seq_along(bases)) {
    f <- paste0(bases[j], '_PM', w)
    if(f %in% names(d)) {
      err <- max(abs(d[[f]] - mat[,j]))
      stopifnot(err < 1e-10)
      feature_qc[[length(feature_qc)+1]] <- data.table(Check=f, MaxAbsError=err, Pass=TRUE)
    }
    set(d, j=f, value=mat[,j])
  }
}
d[, SequencePM6 := paste0(substr(Sequence13,1,6), substr(Sequence13,8,13))]
stopifnot(all(abs(d$HydrophobicFraction_PM6 + d$PolarFraction_PM6 - 1) < 1e-12),
          all(nchar(d$SequencePM6) == 12L), all(is.finite(d$SiteScore)))
write(d, 'sequence_pm6_site_features.tsv')
write(rbindlist(feature_qc), 'feature_reconstruction_qc.tsv')

hc_fit <- function(dat, f, interaction=FALSE) {
  mu <- mean(dat[[f]]); sdv <- sd(dat[[f]])
  stopifnot(is.finite(sdv), sdv > 0)
  dat <- copy(dat); dat[, x := (get(f)-mu)/sdv]
  if(interaction) dat[, PTM := factor(Omics, levels=c('SNO','SSG'))]
  fit <- lm(if(interaction) SiteScore ~ x*PTM else SiteScore ~ x, data=dat)
  term <- if(interaction) 'x:PTMSSG' else 'x'
  G <- uniqueN(dat$Accession)
  V <- vcovCL(fit, cluster=dat$Accession, type='HC1', cadjust=TRUE)
  beta <- unname(coef(fit)[term]); se <- sqrt(V[term,term]); crit <- qt(.975,G-1)
  result <- data.table(Feature=f, Model=if(interaction) 'SiteScore ~ Feature_z * PTM' else 'SiteScore ~ Feature_z',
    PTM=if(interaction) 'SSG_minus_SNO' else unique(dat$Omics), Outcome='corrected_Overall_SiteScore',
    Beta=beta, SE=se, CI_low=beta-crit*se, CI_high=beta+crit*se,
    P=2*pt(-abs(beta/se), G-1), DF=G-1L, N_records=nrow(dat), N_proteins=G,
    FeatureMean=mu, FeatureSD=sdv, Scale=if(interaction) 'pooled_PTM_records' else 'within_PTM',
    CI_method='accession_cluster_HC1_t', B=NA_integer_, ValidBootstrap=NA_integer_)
  list(data=dat, fit=fit, term=term, result=result)
}

# First reproduce the actual legacy coefficient/HC1 tables and Run019 interactions.
legacy <- rbindlist(list(read('historical_ssg.tsv'),read('historical_sno.tsv')))
legacy <- legacy[Scope == 'Overall']
audit <- list()
for(p in c('SNO','SSG')) for(f in c(paste0(bases,'_PM4'),'GRAVY_PM6')) {
  r <- hc_fit(d[Omics==p],f)$result
  old <- legacy[Omics==p & Feature==f]
  stopifnot(nrow(old)==1)
  errs <- c(abs(r$Beta-old$Estimate),abs(r$CI_low-old$CI_low),abs(r$CI_high-old$CI_high))
  stopifnot(max(errs)<1e-8)
  audit[[length(audit)+1]] <- data.table(PTM=p, Feature=f, Check='legacy_HC1_beta_and_CI', MaxAbsError=max(errs), Pass=TRUE)
}
ri <- read('run019_interaction_full.tsv')
for(f in paste0(bases,'_PM4')) {
  r <- hc_fit(d,f,TRUE)$result; old <- ri[Feature==f]
  err <- abs(r$Beta-old$BetaInteraction)
  stopifnot(nrow(old)==1, err<1e-8)
  audit[[length(audit)+1]] <- data.table(PTM='SSG_minus_SNO',Feature=f,Check='Run019_interaction_point',MaxAbsError=err,Pass=TRUE)
}
write(rbindlist(audit), 'legacy_replication_qc.tsv')

# Shared protein resampling for all five fixed features within a branch.
# Duplicate draws enter through integer cluster weights; feature z-scales stay fixed.
boot_branch <- function(dat, interaction, seed, branch) {
  clusters <- unique(dat$Accession); G <- length(clusters); B <- 5000L
  stopifnot(G>=10L, nrow(dat)>=if(interaction)40L else 20L)
  set.seed(seed); counts <- rmultinom(B, G, rep(1/G,G))
  saveRDS(list(Accession=clusters, counts=counts, seed=seed),
          file.path(out,paste0('bootstrap_counts_',branch,'.rds')))
  rows <- list(); draws <- list(); hc <- list()
  for(f in paste0(bases,'_PM6')) {
    model <- hc_fit(dat,f,interaction)
    hc[[length(hc)+1]] <- model$result
    X <- model.matrix(model$fit); y <- model$data$SiteScore; k <- ncol(X)
    ix <- match(model$term,colnames(X))
    sx <- lapply(clusters,function(a) {
      i <- which(model$data$Accession==a)
      list(A=as.numeric(crossprod(X[i,,drop=FALSE])), b=as.numeric(crossprod(X[i,,drop=FALSE],y[i])))
    })
    AA <- do.call(cbind,lapply(sx,`[[`,'A')); BB <- do.call(cbind,lapply(sx,`[[`,'b'))
    sol <- function(cnt) unname(solve(matrix(AA %*% cnt,k,k),as.numeric(BB %*% cnt))[ix])
    stopifnot(abs(sol(rep(1,G))-model$result$Beta)<1e-10)
    values <- vapply(seq_len(B),function(i) tryCatch(sol(counts[,i]),error=function(e)NA_real_),numeric(1))
    stopifnot(sum(is.finite(values))>=.95*B)
    for(n in if(interaction)c(1000L,5000L) else 5000L) {
      z <- values[seq_len(n)]; q <- quantile(z,c(.025,.975),na.rm=TRUE,type=7)
      r <- copy(model$result)
      r[, `:=`(SE=NA_real_,P=NA_real_,DF=NA_integer_,CI_low=q[1],CI_high=q[2],
               CI_method='accession_cluster_percentile_bootstrap',B=n,ValidBootstrap=sum(is.finite(z)),
               Seed=seed,BootstrapMedian=median(z,na.rm=TRUE))]
      rows[[length(rows)+1]] <- r
    }
    draws[[length(draws)+1]] <- data.table(PTM=branch,Feature=f,Iteration=seq_len(B),Beta=values)
  }
  list(boot=rbindlist(rows),hc=rbindlist(hc),draws=rbindlist(draws))
}
branches <- list(
  boot_branch(d[Omics=='SNO'],FALSE,20260928L+100L,'SNO'),
  boot_branch(d[Omics=='SSG'],FALSE,20260928L+200L,'SSG'),
  boot_branch(d,TRUE,20260928L+300L,'SSG_minus_SNO'))
hc <- rbindlist(lapply(branches,`[[`,'hc'))
boot <- rbindlist(lapply(branches,`[[`,'boot'))
draws <- rbindlist(lapply(branches,`[[`,'draws'))
write(hc,'all_pm6_HC1.tsv'); write(boot,'all_pm6_bootstrap.tsv')
write(draws,'bootstrap_coefficients.tsv')
primary <- rbindlist(list(hc[PTM!='SSG_minus_SNO'],boot[PTM=='SSG_minus_SNO' & B==1000L]),fill=TRUE)
primary[, Window := 'Cys_PM6_excluding_central_Cys']
write(primary,'primary_forest_results.tsv')
write(boot[B==5000L],'bootstrap5000_review_results.tsv')
write(primary[,.(PTM,Feature,N_records,N_proteins,FeatureMean,FeatureSD,Scale)],'model_scale_and_coverage.tsv')
coverage <- d[,.(N_records=.N,N_proteins=uniqueN(Accession),MinFiniteSamples=min(FiniteSamples)),by=Omics]
stopifnot(coverage[Omics=='SNO',N_records]==920L,coverage[Omics=='SSG',N_records]==816L)
write(coverage,'coverage.tsv')
recheck <- boot[PTM=='SSG_minus_SNO',.(Feature,B,CI_low,CI_high)]
write(recheck,'interaction_B1000_B5000_comparison.tsv')
write(data.table(R=as.character(getRversion()),data_table=as.character(packageVersion('data.table')),
                 sandwich=as.character(packageVersion('sandwich')),MasterSeed=20260928L), 'environment.tsv')
cat('ANALYSIS_COMPLETE: all five PM6 features; legacy replication and bootstrap QC passed\n')
