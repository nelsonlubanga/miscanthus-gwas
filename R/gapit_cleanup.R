# GAPIT writes ~50 files per trait. Keep the per-model p-value tables, PVE,
# QQ plots and the run log; the genotype PCA (identical for every trait) is
# copied once to results/population_structure/. Everything else is removed.

keep_patterns <- c(
  "^GAPIT\\.Association\\.GWAS_Results\\..*\\(NYC\\)\\.csv$",
  "^GAPIT\\.Association\\.PVE\\..*\\(NYC\\)\\.csv$",
  "^GAPIT\\.Association\\.QQ\\..*\\(NYC\\)\\.pdf$",
  "^gapit\\.log$"
)
pca_files <- c("GAPIT.Genotype.PCA.csv", "GAPIT.Genotype.PCA_2D.pdf",
               "GAPIT.Genotype.PCA_eigenvalues.csv", "GAPIT.Genotype.PCA_eigenValue.pdf")

tidy_gapit_dir <- function(trait_dir, pca_dir = "results/population_structure") {
  dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
  for (f in pca_files) {
    src <- file.path(trait_dir, f)
    if (file.exists(src) && !file.exists(file.path(pca_dir, f))) file.copy(src, pca_dir)
  }
  files <- list.files(trait_dir)
  keep <- Reduce(`|`, lapply(keep_patterns, grepl, x = files))
  unlink(file.path(trait_dir, files[!keep]))
  invisible(files[!keep])
}
