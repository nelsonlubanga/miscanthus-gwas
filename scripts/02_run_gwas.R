# Run FarmCPU, BLINK and MLM (GAPIT) for every trait.
# Each trait runs in its own folder under results/gwas/<trait>/ so traits can
# run in parallel. Usage:
#   Rscript scripts/02_run_gwas.R            # all traits
#   Rscript scripts/02_run_gwas.R AUC_15     # selected traits

suppressPackageStartupMessages({
  library(data.table)
  library(GAPIT)
  library(parallel)
})
source("R/config.R")
source("R/gapit_cleanup.R")

GD <- readRDS(file.path(paths$processed, "GD.rds"))
GM <- as.data.frame(fread(file.path(paths$processed, "GM.csv")))
pheno <- as.data.frame(fread(file.path(paths$processed, "pheno_matched.csv")))

traits <- commandArgs(trailingOnly = TRUE)
if (length(traits) == 0) traits <- names(pheno)[-1]
stopifnot(all(traits %in% names(pheno)))

repo_root <- normalizePath(".")
gwas_root <- file.path(repo_root, paths$gwas)
dir.create(gwas_root, recursive = TRUE, showWarnings = FALSE)

run_trait <- function(trait) {
  out_dir <- file.path(gwas_root, trait)
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  log_con <- file("gapit.log", open = "wt")
  sink(log_con); sink(log_con, type = "message")
  on.exit({
    sink(type = "message"); sink(); close(log_con)
    setwd(repo_root)
  })

  Y <- pheno[!is.na(pheno[[trait]]), c("geno", trait)]
  res <- tryCatch({
    GAPIT(Y = Y, GD = GD, GM = GM,
          PCA.total = params$pca_total,
          model = params$models,
          SNP.MAF = 0,          # MAF already filtered in step 01
          file.output = TRUE)
    tidy_gapit_dir(out_dir, file.path(repo_root, "results/population_structure"))
    "ok"
  }, error = function(e) paste("error:", conditionMessage(e)))
  data.frame(trait = trait, n = nrow(Y), status = res)
}

t0 <- Sys.time()
status <- mclapply(traits, run_trait, mc.cores = min(params$n_cores, length(traits)),
                   mc.preschedule = FALSE)
status <- rbindlist(lapply(status, function(x)
  if (inherits(x, "try-error")) data.frame(trait = NA, n = NA, status = as.character(x)) else x))
print(status)
cat("Elapsed:", format(Sys.time() - t0), "\n")
fwrite(status, file.path(gwas_root, sprintf("run_status_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))))
