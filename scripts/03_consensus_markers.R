# Collect GAPIT results, apply the Bonferroni threshold, and keep markers that
# are significant for the same trait in at least `min_models` of the models.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("R/config.R")
dir.create(paths$consensus, recursive = TRUE, showWarnings = FALSE)

GM <- fread(file.path(paths$processed, "GM.csv"))
n_snps <- nrow(GM)
bonf <- params$alpha / n_snps
cat(sprintf("Bonferroni threshold: %.3e (%d SNPs)\n", bonf, n_snps))

# GAPIT 4 writes "(NYC)" (final) and "(Kansas)" copies; use the NYC files
files <- list.files(paths$gwas, pattern = "^GAPIT\\.Association\\.GWAS_Results\\..*\\(NYC\\)\\.csv$",
                    recursive = TRUE, full.names = TRUE)
pattern <- "^GAPIT\\.Association\\.GWAS_Results\\.([A-Za-z]+)\\.(.+)\\(NYC\\)\\.csv$"

all_res <- rbindlist(lapply(files, function(f) {
  b <- basename(f)
  d <- fread(f)
  setnames(d, names(d)[tolower(names(d)) == "effect"], "Effect")
  d[, .(SNP, Chr, Pos, P.value, MAF, Effect,
        model = sub(pattern, "\\1", b), trait = sub(pattern, "\\2", b))]
}))
all_res <- all_res[model %in% params$models]
cat("Result sets loaded:", uniqueN(all_res[, .(model, trait)]), "trait x model\n")

# Genomic inflation factor per trait x model
lambda <- all_res[, .(lambda_gc = median(qchisq(P.value, 1, lower.tail = FALSE), na.rm = TRUE) /
                        qchisq(0.5, 1)), by = .(trait, model)]
fwrite(lambda, file.path(paths$consensus, "genomic_inflation.csv"))

sig <- all_res[P.value < bonf][order(trait, Chr, Pos, model)]
sig[, neg_log10_p := -log10(P.value)]
fwrite(sig, file.path(paths$consensus, "significant_markers_all_models.csv"))

# Consensus: same SNP and trait significant in >= min_models models
consensus <- sig[, .(n_models = uniqueN(model),
                     models = paste(sort(unique(model)), collapse = "+"),
                     min_p = min(P.value),
                     p_FarmCPU = P.value[model == "FarmCPU"][1],
                     p_BLINK = P.value[model == "BLINK"][1],
                     p_MLM = P.value[model == "MLM"][1],
                     effect_FarmCPU = Effect[model == "FarmCPU"][1],
                     effect_BLINK = Effect[model == "BLINK"][1],
                     effect_MLM = Effect[model == "MLM"][1],
                     MAF = MAF[1]),
                 by = .(trait, SNP, Chr, Pos)][n_models >= params$min_models]
setorder(consensus, Chr, Pos, trait)
fwrite(consensus, file.path(paths$consensus, "consensus_markers.csv"))

summary_tab <- merge(
  dcast(sig[, .N, by = .(trait, model)], trait ~ model, value.var = "N", fill = 0),
  consensus[, .(n_consensus = .N), by = trait], by = "trait", all.x = TRUE)
summary_tab[is.na(n_consensus), n_consensus := 0]
fwrite(summary_tab, file.path(paths$consensus, "significant_counts_by_trait.csv"))

cat("\nSignificant markers per model:\n"); print(sig[, .N, by = model])
cat("\nConsensus marker-trait associations:", nrow(consensus),
    "| unique SNPs:", uniqueN(consensus$SNP), "\n")
print(consensus[, .(trait, SNP, models, min_p = signif(min_p, 3), MAF = round(MAF, 3))])

# Manhattan plot of consensus hits across traits (best p over models)
chr_len <- GM[, .(len = max(Position)), by = Chromosome][order(Chromosome)]
chr_len[, offset := cumsum(as.numeric(shift(len, fill = 0)))]
best <- all_res[, .(P = min(P.value)), by = .(SNP, Chr, Pos)]
best <- merge(best, chr_len, by.x = "Chr", by.y = "Chromosome")
best[, bp_cum := Pos + offset]
axis_df <- best[, .(center = (min(bp_cum) + max(bp_cum)) / 2), by = Chr]
cons_pts <- merge(unique(consensus[, .(SNP)]), best, by = "SNP")

p <- ggplot(best, aes(bp_cum, -log10(P))) +
  geom_point(aes(colour = factor(Chr %% 2)), size = 0.6, show.legend = FALSE) +
  geom_point(data = cons_pts, colour = "#D6604D", size = 2) +
  geom_hline(yintercept = -log10(bonf), linetype = "dashed") +
  scale_colour_manual(values = c("grey35", "grey65")) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$Chr, expand = c(0.01, 0)) +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))),
       title = "Best p-value per SNP across all traits and models",
       subtitle = sprintf("Dashed line: Bonferroni (%.2e); red: significant in >= %d models",
                          bonf, params$min_models)) +
  theme_classic(base_size = 12)
ggsave(file.path(paths$consensus, "manhattan_consensus.png"), p, width = 11, height = 4.5, dpi = 300)
