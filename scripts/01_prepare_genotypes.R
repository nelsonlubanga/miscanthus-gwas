# Convert the HapMap genotypes to 0/1/2 (minor-allele dosage), filter SNPs,
# and align genotypes with phenotypes. Outputs GAPIT-ready GD/GM objects.

suppressPackageStartupMessages(library(data.table))
source("R/config.R")
dir.create(paths$processed, recursive = TRUE, showWarnings = FALSE)

hmp <- fread(paths$geno_hmp, sep = "\t", colClasses = "character")
pheno <- fread(paths$pheno)
setnames(pheno, 1, "geno")

taxa_cols <- names(hmp)[12:ncol(hmp)]
cat("HapMap:", nrow(hmp), "SNPs x", length(taxa_cols), "taxa\n")

# Chromosome-anchored biallelic SNPs only
hmp <- hmp[!grepl("SCAFFOLD", chrom, ignore.case = TRUE)]
alleles <- tstrsplit(hmp$alleles, "/", fixed = TRUE)
keep <- lengths(strsplit(hmp$alleles, "/", fixed = TRUE)) == 2 &
  nchar(alleles[[1]]) == 1 & nchar(alleles[[2]]) == 1
hmp <- hmp[keep]
cat("After removing scaffolds / non-biallelic:", nrow(hmp), "SNPs\n")

# Taxa present in both datasets, in phenotype order
taxa <- intersect(pheno$geno, taxa_cols)
cat("Taxa with phenotype and genotype:", length(taxa), "\n")
pheno <- pheno[match(taxa, geno)]

iupac_het <- c("R", "Y", "S", "W", "K", "M")
calls <- as.matrix(hmp[, ..taxa])
a1 <- substr(hmp$alleles, 1, 1)
a2 <- substr(hmp$alleles, 3, 3)

# Dosage of allele 2 (0 = a1/a1, 1 = het, 2 = a2/a2)
dos <- matrix(NA_real_, nrow(calls), ncol(calls))
dos[calls == a1] <- 0
dos[calls == a2] <- 2
dos[calls %in% iupac_het] <- 1

miss <- rowMeans(is.na(dos))
p2 <- rowMeans(dos, na.rm = TRUE) / 2
maf <- pmin(p2, 1 - p2)

snp_ok <- miss <= params$snp_missing_max & !is.na(maf) & maf >= params$maf_min
cat(sprintf("SNPs failing missingness > %.2f: %d\n", params$snp_missing_max,
            sum(miss > params$snp_missing_max)))
cat(sprintf("SNPs failing MAF < %.2f: %d\n", params$maf_min,
            sum(miss <= params$snp_missing_max & (is.na(maf) | maf < params$maf_min))))

dos <- dos[snp_ok, , drop = FALSE]
p2 <- p2[snp_ok]
hmp <- hmp[snp_ok]

# Recode so 2 = minor-allele homozygote; impute missing with the SNP mean dosage
flip <- p2 > 0.5
dos[flip, ] <- 2 - dos[flip, ]
na_idx <- which(is.na(dos), arr.ind = TRUE)
dos[na_idx] <- (2 * pmin(p2, 1 - p2))[na_idx[, 1]]

GM <- data.frame(SNP = hmp[[1]],
                 Chromosome = as.integer(hmp$chrom),
                 Position = as.integer(hmp$pos))
GD <- data.frame(taxa = taxa, t(dos), check.names = FALSE)
names(GD)[-1] <- GM$SNP

cat("Final marker set:", nrow(GM), "SNPs x", length(taxa), "taxa\n")
cat(sprintf("Bonferroni threshold (alpha = %.2f): %.3e (-log10 = %.2f)\n",
            params$alpha, params$alpha / nrow(GM),
            -log10(params$alpha / nrow(GM))))

saveRDS(GD, file.path(paths$processed, "GD.rds"))
fwrite(GM, file.path(paths$processed, "GM.csv"))
fwrite(pheno, file.path(paths$processed, "pheno_matched.csv"))
fwrite(data.table(step = c("taxa", "snps_final", "bonferroni_p"),
                  value = c(length(taxa), nrow(GM), params$alpha / nrow(GM))),
       file.path(paths$processed, "marker_qc_summary.csv"))
