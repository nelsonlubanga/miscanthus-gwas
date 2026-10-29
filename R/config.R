# Shared settings for the Miscanthus (ABR33) GWAS pipeline.
# All paths are relative to the repository root; run scripts from there.

paths <- list(
  pheno      = "data/raw/pheno.final.txt",
  geno_hmp   = "data/raw/geno.final.hmp.txt",
  gff        = "data/reference/Msinensis_497_v7.1.gene.gff3.gz",
  annotation = "data/reference/Msinensis_497_v7.1.P14.annotation_info.txt.gz",
  processed  = "data/processed",
  gwas       = "results/gwas",
  consensus  = "results/consensus",
  candidates = "results/candidate_genes"
)

params <- list(
  models          = c("FarmCPU", "BLINK", "MLM"),
  pca_total       = 3,
  maf_min         = 0.05,
  snp_missing_max = 0.50,   # drop SNPs with >50% missing calls (GBS data)
  alpha           = 0.05,   # Bonferroni: alpha / n_SNPs
  min_models      = 2,      # consensus = significant in >= this many models
  window_bp       = 10000,  # +/- candidate-gene search window (LD decay, 03_LD_decay)
  n_cores         = 4       # traits run in parallel
)
