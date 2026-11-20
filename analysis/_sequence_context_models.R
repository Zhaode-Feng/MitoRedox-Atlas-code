# Original sequence regression and resampling logic. Paths adapted for this release.
suppressPackageStartupMessages(library(sandwich))
read <- function(n) read_input(file.path('sequence', n), na.strings = c('NA', ''))
write <- function(x, n) write_table(x, file.path('05_sequence', n))
d <- data.table::fread(file.path(output, '05_sequence/sequence_pm6_site_features.tsv'))
oldA <- read('historical_mundlak.tsv');oldD <- read('historical_compartment.tsv')
stopifnot(nrow(d)==1736L,!anyDuplicated(d[,.(Omics,ExactSiteID)]),
          all(nchar(d$Sequence13)==13L),all(substr(d$Sequence13,7,7)=='C'))
ei <- c(A=.62,C=.29,D=-.90,E=-.74,F=1.19,G=.48,H=-.40,I=1.38,K=-1.50,
        L=1.06,M=.64,N=-.78,P=.12,Q=-.85,R=-2.53,S=-.18,T=-.05,V=1.08,W=.81,Y=.26)
eisenberg <- function(s,w) mean(ei[strsplit(s,'',fixed=TRUE)[[1]][c((7-w):6,8:(7+w))]])
stopifnot(max(abs(vapply(d$Sequence13,eisenberg,numeric(1),w=4L)-d$Eisenberg_PM4))<1e-10)
d[,Eisenberg_PM6:=vapply(Sequence13,eisenberg,numeric(1),w=6L)]
bases <- c('GRAVY','HydrophobicFraction','PolarFraction','ChargedFraction','ProFraction','Eisenberg')
features <- paste0(bases,'_PM6')
stopifnot(all(vapply(d[,..features],function(x)all(is.finite(x)),logical(1))))
write(d,'sequence_pm6_six_features.tsv')
audits <- list();scales <- list();members <- list()

cluster_fit <- function(dat,formula,terms) {
  G <- uniqueN(dat$Accession)
  stopifnot(nrow(dat)>=20L,G>=8L)
  fit <- lm(formula,data=dat)
  stopifnot(fit$rank==length(coef(fit)))
  V <- vcovCL(fit,cluster=dat$Accession,type='HC1',cadjust=TRUE)
  beta <- coef(fit)[terms];se <- sqrt(diag(V)[terms]);crit <- qt(.975,G-1L)
  list(fit=fit,V=V,result=data.table(TermCoefficient=terms,Estimate=as.numeric(beta),SE=as.numeric(se),
     CI_low=as.numeric(beta-crit*se),CI_high=as.numeric(beta+crit*se),
     P=as.numeric(2*pt(-abs(beta/se),G-1L)),ClusterDF=G-1L,
     N_sites=nrow(dat),N_proteins=G,CI_method='accession_cluster_HC1_t'))
}
mundlak <- function(p,f,record=FALSE) {
  z <- copy(d[Omics==p & is.finite(SiteScore) & is.finite(get(f))])
  z[,ProteinSiteCount:=.N,by=Accession];z <- z[ProteinSiteCount>=2L]
  mu <- mean(z[[f]]);sdv <- sd(z[[f]])
  z[,Feature_z:=(get(f)-mu)/sdv]
  z[,Feature_between:=mean(Feature_z),by=Accession]
  z[,Feature_within:=Feature_z-Feature_between]
  ans <- cluster_fit(z,SiteScore~Feature_within+Feature_between,c('Feature_within','Feature_between'))$result
  ans[,`:=`(Omics=p,Feature=f,Component=c('Within_protein','Between_protein'),
             Model='SiteScore ~ Feature_within + Feature_between',FeatureMean=mu,FeatureSD=sdv,
             EstimateOrigin='new_PM6')]
  if(record){
    scales[[length(scales)+1L]] <<- data.table(Panel='A',Omics=p,Feature=f,Center=mu,SD=sdv,ScalePopulation='PTM_multi_site_proteins',N=nrow(z))
    members[[length(members)+1L]] <<- z[,.(Panel='A',Omics,ExactSiteID,Accession,Feature=f,Feature_z,Feature_within,Feature_between)]
  }
  ans
}
for(p in c('SNO','SSG'))for(f in unique(oldA$Feature)){
  z <- mundlak(p,f);ref <- oldA[Omics==p & Feature==f][match(z$Component,Component)]
  err <- max(abs(as.matrix(z[,.(Estimate,SE,CI_low,CI_high,P)])-as.matrix(ref[,.(Estimate,SE,CI_low,CI_high,P)])))
  stopifnot(nrow(ref)==2L,err<1e-8,all(z$N_sites==ref$N_sites),all(z$N_proteins==ref$N_proteins))
  audits[[length(audits)+1L]] <- data.table(Check='Mundlak_legacy_reproduction',Omics=p,Feature=f,MaxAbsError=err,Pass=TRUE)
}
A <- rbindlist(lapply(c('SNO','SSG'),function(p)rbindlist(lapply(features,function(f)mundlak(p,f,TRUE)))))
A[,FDR_primary:=p.adjust(P,'BH'),by=.(Omics,Component)]
A[,FDR_global:=p.adjust(P,'BH')]
write(A,'FigureS12_A_mundlak_PM6.tsv')

# Rebuild the historical unambiguous compartment map (accession first, gene fallback).
annotation <- read('mitocarta_submitochondrial_map.tsv')
split_tokens <- function(x) unique(toupper(trimws(unlist(strsplit(ifelse(is.na(x),'',x),'[|;, ]+')))))
normalise_loc <- function(x) {
  z <- toupper(trimws(unlist(strsplit(ifelse(is.na(x),'',x),'[|;,/]+'))))
  z <- ifelse(grepl('MIM|INNER',z),'MIM',ifelse(grepl('MOM|OUTER',z),'MOM',
       ifelse(grepl('IMS|INTERMEMBRANE',z),'IMS',ifelse(grepl('MATRIX',z),'Matrix',NA))))
  unique(z[!is.na(z)])
}
acc_map <- list();gene_map <- list()
for(i in seq_len(nrow(annotation))){
  loc <- normalise_loc(annotation$MitoCarta3.0_SubMitoLocalization[i]);if(!length(loc))next
  for(a in split_tokens(annotation$UniProt[i]))acc_map[[a]] <- unique(c(acc_map[[a]],loc))
  for(g in split_tokens(annotation$Symbol[i]))gene_map[[g]] <- unique(c(gene_map[[g]],loc))
}
lookup <- function(acc,gene){
  z <- acc_map[[toupper(acc)]];if(is.null(z))z <- gene_map[[toupper(gene)]]
  if(is.null(z)||!length(z))return(NA_character_)
  if(length(z)>1L)return('Ambiguous')
  z
}
d[,Compartment:=mapply(lookup,Accession,Gene,USE.NAMES=FALSE)]
write(d[,.(Omics,ExactSiteID,Accession,Gene,Compartment)],'compartment_assignment.tsv')
eligible <- d[Compartment %in% c('Matrix','MIM','MOM','IMS'),.(N_sites=.N,N_proteins=uniqueN(Accession)),by=.(Omics,Compartment)]
eligible[,Eligible:=N_sites>=30L & N_proteins>=20L]
write(eligible,'compartment_eligibility.tsv')
comp_model <- function(p,f,record=FALSE){
  good <- eligible[Omics==p & Eligible==TRUE,Compartment]
  z <- copy(d[Omics==p & Compartment %in% good])
  mu <- mean(z[[f]]);sdv <- sd(z[[f]])
  z[,Feature_z:=(get(f)-mu)/sdv]
  z[,Compartment:=factor(Compartment)]
  full <- cluster_fit(z,SiteScore~Feature_z*Compartment,'Feature_z')
  ix <- grep('Feature_z:Compartment|Compartment.*:Feature_z',names(coef(full$fit)))
  b <- coef(full$fit)[ix];vv <- full$V[ix,ix,drop=FALSE]
  hp <- pchisq(as.numeric(t(b)%*%solve(vv,b)),df=length(ix),lower.tail=FALSE)
  ans <- rbindlist(lapply(c('Matrix','MIM','MOM','IMS')[c('Matrix','MIM','MOM','IMS')%in%good],function(cc){
    a <- cluster_fit(z[Compartment==cc],SiteScore~Feature_z,'Feature_z')$result
    a[,`:=`(Omics=p,Scope=cc,Feature=f,Model='SiteScore ~ Feature_z',Covariates='None',
             HeterogeneityP=hp,EstimateOrigin='new_PM6',FeatureMean=mu,FeatureSD=sdv)]
    a
  }))
  if(record){
    scales[[length(scales)+1L]] <<- data.table(Panel='D',Omics=p,Feature=f,Center=mu,SD=sdv,ScalePopulation='PTM_pooled_eligible_compartments',N=nrow(z))
    members[[length(members)+1L]] <<- z[,.(Panel='D',Omics,ExactSiteID,Accession,Feature=f,Scope=as.character(Compartment),Feature_z)]
  }
  ans
}
for(p in c('SNO','SSG')){
  z <- comp_model(p,'GRAVY_PM4');ref <- oldD[Omics==p & Feature=='GRAVY_PM4'][match(z$Scope,Scope)]
  err <- max(abs(as.matrix(z[,.(Estimate,SE,CI_low,CI_high,P,HeterogeneityP)])-as.matrix(ref[,.(Estimate,SE,CI_low,CI_high,P,HeterogeneityP)])))
  stopifnot(nrow(ref)==3L,err<1e-8,all(z$N_sites==ref$N_sites),all(z$N_proteins==ref$N_proteins))
  audits[[length(audits)+1L]] <- data.table(Check='Compartment_GRAVY_PM4_reproduction',Omics=p,Feature='GRAVY_PM4',MaxAbsError=err,Pass=TRUE)
}
newG <- rbindlist(lapply(c('SNO','SSG'),function(p)comp_model(p,'GRAVY_PM6',TRUE)))
retained <- copy(oldD[Feature!='GRAVY_PM4'])
retained[,`:=`(EstimateOrigin='unchanged_frozen_3D',CI_method='accession_cluster_HC1_t',ClusterDF=N_proteins-1L)]
keep <- c('Omics','Scope','Feature','Model','Covariates','Estimate','SE','CI_low','CI_high','P','N_sites','N_proteins','HeterogeneityP','EstimateOrigin','CI_method','ClusterDF')
D <- rbindlist(list(newG,retained[,..keep]),fill=TRUE)
D[,FDR_primary:=p.adjust(P,'BH'),by=Omics]
D[,FDR_global:=p.adjust(P,'BH')]
hetero <- unique(D[,.(Omics,Feature,HeterogeneityP)])
hetero[,HeterogeneityFDR:=p.adjust(HeterogeneityP,'BH'),by=Omics]
D <- merge(D,hetero[,.(Omics,Feature,HeterogeneityFDR)],by=c('Omics','Feature'),sort=FALSE)
D[,DirectionConsistent:=uniqueN(sign(Estimate))==1L,by=.(Omics,Feature)]
D[,CompartmentConsistent:=DirectionConsistent & HeterogeneityFDR>=.05]
stopifnot(nrow(A)==24L,nrow(D)==24L,all(is.finite(A$CI_low)),all(is.finite(D$CI_low)))
preserved <- D[EstimateOrigin=='unchanged_frozen_3D']
joined <- merge(preserved,oldD,by=c('Omics','Scope','Feature'),suffixes=c('_new','_old'))
stopifnot(nrow(joined)==18L)
for(f in c('Estimate','SE','CI_low','CI_high','P','N_sites','N_proteins','HeterogeneityP'))stopifnot(max(abs(joined[[paste0(f,'_new')]]-joined[[paste0(f,'_old')]]))==0)
write(D,'FigureS12_D_compartment_PM6.tsv')
write(newG,'compartment_GRAVY_PM6_new_results.tsv')
write(rbindlist(audits),'legacy_reproduction_qc.tsv')
write(rbindlist(scales),'standardization_audit.tsv')
write(rbindlist(members,fill=TRUE),'model_membership_and_predictors.tsv')
write(A[,.(N_sites=unique(N_sites),N_proteins=unique(N_proteins)),by=Omics],'mundlak_coverage.tsv')
write(data.table(R=as.character(getRversion()),data_table=as.character(packageVersion('data.table')),
                 sandwich=as.character(packageVersion('sandwich'))),'environment.tsv')
cat('ANALYSIS_COMPLETE: A=24 PM6 component estimates; D=6 new PM6 + 18 preserved 3D estimates\n')
