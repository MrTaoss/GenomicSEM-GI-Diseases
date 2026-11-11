# Genomic SEM analysis of common chronic gastrointestinal diseases
# Associated manuscript: Resolving the latent genetic architecture and spatial
# enrichment patterns of common gastrointestinal diseases
#
# R 4.5.2 | GenomicSEM 0.0.5 | lavaan 0.6-21
# Inputs: GRCh37 harmonized GWAS files, munged HapMap3 summary statistics,
#         1000 Genomes Phase 3 EUR LD scores, and the GenomicSEM SNP reference.
# Edit project_dir to point to the prepared analysis directory.
# This script covers factor modeling and factor GWAS; upstream GWAS QC and
# downstream FUMA, Q_SNP local-LD filtering, and functional analyses are separate.

suppressPackageStartupMessages({
  library(GenomicSEM)
  library(data.table)
  library(Matrix)
})

project_dir <- "path/to/project"
munge_dir   <- file.path(project_dir, "munged")
gwas_dir    <- file.path(project_dir, "harmonized_gwas")
ref_dir     <- file.path(project_dir, "reference")
ld_dir      <- paste0(file.path(ref_dir, "eur_w_ld_chr"), "/")
snp_ref     <- file.path(ref_dir, "reference.1000G.maf.0.005.txt")
out_dir     <- file.path(project_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

traits9  <- c("CD", "UC", "GERD", "PUD", "GASTD", "IBS", "CONST", "DIVD", "HEM")
traits10 <- append(traits9, "CED", after = 2)

# Sample prevalence is 0.5 for traits analyzed with effective sample size.
prev <- data.table(
  trait = traits10,
  sample = c(0.5, 0.5, 0.5, 175991/412728, 16413/438625,
             28941/407065, 0.5, 57028/500186, 0.5, 0.5),
  population = c(0.001165, 0.001986, 0.007, 0.1412, 0.10,
                 0.205, 0.0908, 0.1034, 0.20, 0.2603)
)

run_ldsc <- function(traits, label, select = NULL) {
  p <- prev[match(traits, prev$trait)]
  args <- list(
    traits = file.path(munge_dir, paste0(traits, ".sumstats.gz")),
    sample.prev = p$sample,
    population.prev = p$population,
    ld = ld_dir, wld = ld_dir,
    trait.names = traits,
    n.blocks = 200,
    ldsc.log = file.path(out_dir, paste0("ldsc_", label)),
    stand = TRUE
  )
  if (!is.null(select)) args$select <- select
  fit <- do.call(GenomicSEM::ldsc, args)
  dimnames(fit$S) <- list(traits, traits)
  dimnames(fit$I) <- list(traits, traits)
  dimnames(fit$S_Stand) <- list(traits, traits)
  lower <- which(lower.tri(matrix(FALSE, length(traits), length(traits)), diag = TRUE), arr.ind = TRUE)
  param_names <- paste0(traits[lower[, "row"]], "__", traits[lower[, "col"]])
  dimnames(fit$V) <- list(param_names, param_names)
  dimnames(fit$V_Stand) <- list(param_names, param_names)
  saveRDS(fit, file.path(out_dir, paste0("ldsc_", label, ".rds")))
  fit
}

# 1. Initial ten-phenotype structure and primary nine-phenotype structure.
ldsc10 <- run_ldsc(traits10, "10traits")
ldsc9  <- run_ldsc(traits9, "9traits")

run_efa <- function(fit, label) {
  pa_dir <- file.path(out_dir, paste0("parallel_analysis_", label))
  dir.create(pa_dir, showWarnings = FALSE)
  old_wd <- getwd()
  setwd(pa_dir)
  set.seed(20260729)
  capture.output(
    GenomicSEM::paLDSC(S = fit$S_Stand, V = fit$V_Stand,
                       r = 500, p = 0.95, save.pdf = TRUE,
                       diag = FALSE, fa = FALSE),
    file = "results.txt"
  )
  setwd(old_wd)
  efa <- lapply(1:4, function(k)
    stats::factanal(covmat = fit$S, factors = k,
                    rotation = if (k == 1L) "none" else "promax")
  )
  names(efa) <- paste0("EFA", 1:4)
  saveRDS(efa, file.path(out_dir, paste0("efa_", label, ".rds")))
  efa
}

efa10 <- run_efa(ldsc10, "10traits")
efa9  <- run_efa(ldsc9, "9traits")

# 2. Confirmatory models; FIBD comprises CD and UC.
models <- list(
  one_factor = '
    FGI =~ CD + UC + GERD + PUD + GASTD + IBS + CONST + DIVD + HEM
  ',
  two_factor = '
    FGI =~ GERD + PUD + GASTD + IBS + CONST + DIVD + HEM
    FIBD =~ CD + UC
    FGI ~~ FIBD
  ',
  upper_only = '
    FUPPER =~ GERD + PUD + GASTD + DIVD
    FBOWEL =~ IBS + CONST + HEM
    FIBD =~ CD + UC
    FUPPER ~~ FBOWEL
    FUPPER ~~ FIBD
    FBOWEL ~~ FIBD
  ',
  bowel_only = '
    FUPPER =~ GERD + PUD + GASTD
    FBOWEL =~ IBS + CONST + HEM + DIVD
    FIBD =~ CD + UC
    FUPPER ~~ FBOWEL
    FUPPER ~~ FIBD
    FBOWEL ~~ FIBD
  ',
  cross_loading = '
    FUPPER =~ GERD + PUD + GASTD + DIVD
    FBOWEL =~ IBS + CONST + HEM + DIVD
    FIBD =~ CD + UC
    FUPPER ~~ FBOWEL
    FUPPER ~~ FIBD
    FBOWEL ~~ FIBD
  ',
  equal_ibd = '
    FUPPER =~ GERD + PUD + GASTD + DIVD
    FBOWEL =~ IBS + CONST + HEM + DIVD
    FIBD =~ ibd*CD + ibd*UC
    FUPPER ~~ FBOWEL
    FUPPER ~~ FIBD
    FBOWEL ~~ FIBD
  '
)

fit_model <- function(covstruc, syntax, imp_cov = TRUE) {
  GenomicSEM::usermodel(
    covstruc = covstruc, estimation = "DWLS", model = syntax,
    CFIcalc = TRUE, std.lv = TRUE, imp_cov = imp_cov,
    fix_resid = FALSE, Q_Factor = FALSE
  )
}

full_names <- c("one_factor", "two_factor", "upper_only",
                "bowel_only", "cross_loading", "equal_ibd")
full_cfa <- lapply(models[full_names], function(x) fit_model(ldsc9, x))
saveRDS(full_cfa, file.path(out_dir, "cfa_full.rds"))

# 3. Chromosome-split validation: EFA on odd chromosomes, CFA on even.
odd  <- run_ldsc(traits9, "odd", select = "ODD")
even <- run_ldsc(traits9, "even", select = "EVEN")
odd_S <- odd$S
if (min(eigen(odd_S, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
  odd_S <- as.matrix(Matrix::nearPD(odd_S, corr = FALSE)$mat)
  dimnames(odd_S) <- dimnames(odd$S)
}
odd_efa <- lapply(1:4, function(k)
  stats::factanal(covmat = odd_S, factors = k,
                  rotation = if (k == 1L) "none" else "promax")
)
names(odd_efa) <- paste0("EFA", 1:4)
saveRDS(odd_efa, file.path(out_dir, "efa_odd.rds"))

even_names <- c("two_factor", "upper_only", "bowel_only", "cross_loading")
even_cfa <- lapply(models[even_names], function(x)
  fit_model(even, x, imp_cov = FALSE)
)
saveRDS(even_cfa, file.path(out_dir, "cfa_even.rds"))

fit_summary <- function(x, dataset) {
  rbindlist(lapply(names(x), function(nm) {
    z <- as.data.table(as.data.frame(x[[nm]]$modelfit))
    z[, `:=`(dataset = dataset, model = nm)]
    z
  }), fill = TRUE)
}
fwrite(rbindlist(list(fit_summary(full_cfa, "full"),
                     fit_summary(even_cfa, "even")), fill = TRUE),
       file.path(out_dir, "cfa_model_fit.tsv"), sep = "\t")

nested_chisq <- function(fits, restricted, unrestricted) {
  a <- as.data.frame(fits[[restricted]]$modelfit)
  b <- as.data.frame(fits[[unrestricted]]$modelfit)
  dchisq <- as.numeric(a$chisq[1] - b$chisq[1])
  ddf <- as.numeric(a$df[1] - b$df[1])
  data.table(restricted = restricted, unrestricted = unrestricted,
             delta_chisq = dchisq, delta_df = ddf,
             p = pchisq(dchisq, df = ddf, lower.tail = FALSE))
}

comparisons <- rbindlist(list(
  cbind(dataset = "full", nested_chisq(full_cfa, "upper_only", "cross_loading")),
  cbind(dataset = "full", nested_chisq(full_cfa, "bowel_only", "cross_loading")),
  cbind(dataset = "full", nested_chisq(full_cfa, "equal_ibd", "cross_loading")),
  cbind(dataset = "even", nested_chisq(even_cfa, "upper_only", "cross_loading")),
  cbind(dataset = "even", nested_chisq(even_cfa, "bowel_only", "cross_loading"))
), fill = TRUE)
fwrite(comparisons, file.path(out_dir, "cfa_nested_comparisons.tsv"), sep = "\t")

# 4. Factor GWAS, with unit-loading identification (GERD, IBS, CD).
snp_data <- GenomicSEM::sumstats(
  files = file.path(gwas_dir, paste0(traits9, ".txt")),
  ref = snp_ref,
  trait.names = traits9,
  se.logit = rep(TRUE, length(traits9)),
  OLS = NULL, linprob = NULL, N = NULL, betas = NULL,
  info.filter = 0.60, maf.filter = 0.01,
  keep.indel = FALSE, parallel = TRUE, cores = 2,
  ambig = FALSE, direct.filter = FALSE
)

factor_model <- '
  FUPPER =~ 1*GERD + PUD + GASTD + DIVD
  FBOWEL =~ 1*IBS + CONST + HEM + DIVD
  FIBD =~ 1*CD + UC
  FUPPER ~~ FBOWEL
  FUPPER ~~ FIBD
  FBOWEL ~~ FIBD
  FUPPER ~ SNP
  FBOWEL ~ SNP
  FIBD ~ SNP
'

# Run by chromosome, as in the source analysis.
for (chr in 1:22) {
  gwas_chr <- GenomicSEM::userGWAS(
    covstruc = ldsc9,
    SNPs = snp_data[snp_data$CHR == chr, , drop = FALSE],
    estimation = "DWLS",
    model = factor_model,
    printwarn = TRUE,
    sub = c("FUPPER~SNP", "FBOWEL~SNP", "FIBD~SNP"),
    cores = 4, toler = FALSE, SNPSE = FALSE,
    parallel = TRUE, GC = "standard", MPI = FALSE,
    smooth_check = TRUE, TWAS = FALSE,
    std.lv = FALSE, fix_measurement = TRUE,
    Q_SNP = TRUE, analytic = FALSE
  )
  saveRDS(gwas_chr, file.path(out_dir, sprintf("factor_gwas_chr%02d.rds", chr)))
}

# The saved userGWAS results are unfiltered. The original downstream scripts
# excluded SNP-level technical failures and Z_smooth > 1.96, followed by
# factor-specific Q_SNP P < 5e-8 filtering and EUR-LD proxy removal
# (r2 >= 0.10 within +/- 1 Mb) before FUMA and downstream analyses.
