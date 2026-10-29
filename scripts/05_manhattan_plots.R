# One Manhattan figure per trait: FarmCPU, BLINK and MLM stacked, Bonferroni
# threshold, consensus SNPs highlighted and labelled with their nearest gene.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("R/config.R")
out_dir <- "results/manhattan"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GM <- fread(file.path(paths$processed, "GM.csv"))
bonf <- params$alpha / nrow(GM)
consensus <- fread(file.path(paths$consensus, "consensus_markers.csv"))
nearest <- fread(file.path(paths$candidates, "candidate_genes_by_marker.csv"))[
  , .(SNP, nearest_gene, distance_bp)]

# Cumulative genome coordinates
chr_len <- GM[, .(len = max(Position)), by = Chromosome][order(Chromosome)]
chr_len[, offset := cumsum(as.numeric(shift(len, fill = 0)))]
axis_df <- merge(GM, chr_len, by = "Chromosome")[
  , .(center = (min(Position + offset) + max(Position + offset)) / 2), by = Chromosome]

files <- list.files(paths$gwas, pattern = "^GAPIT\\.Association\\.GWAS_Results\\..*\\(NYC\\)\\.csv$",
                    recursive = TRUE, full.names = TRUE)
pattern <- "^GAPIT\\.Association\\.GWAS_Results\\.([A-Za-z]+)\\.(.+)\\(NYC\\)\\.csv$"
res <- rbindlist(lapply(files, function(f) {
  d <- fread(f, select = c("SNP", "Chr", "Pos", "P.value"))
  d[, `:=`(model = sub(pattern, "\\1", basename(f)), trait = sub(pattern, "\\2", basename(f)))]
}))
res <- res[model %in% params$models]
res <- merge(res, chr_len[, .(Chromosome, offset)], by.x = "Chr", by.y = "Chromosome")
res[, `:=`(bp_cum = Pos + offset, logp = -log10(P.value),
           model = factor(model, levels = params$models))]


plot_trait <- function(tr) {
  d <- res[trait == tr]
  cons_snps <- consensus[trait == tr, SNP]
  d[, class := fifelse(SNP %in% cons_snps, "consensus",
                fifelse(P.value < bonf, "single", "ns"))]
  # Consensus SNPs: filled + labelled where significant in that model, hollow otherwise
  d[class == "consensus" & P.value >= bonf, class := "consensus_ns"]
  lab <- merge(d[class == "consensus"], nearest, by = "SNP", all.x = TRUE)
  lab[, label := fifelse(is.na(nearest_gene) | nearest_gene == "", paste0(SNP, "\n(no gene within 10 kb)"),
                         paste0(SNP, "\n", nearest_gene))]
  # Stagger labels of nearby SNPs (within 6% of the genome) in the same panel
  genome_len <- max(res$bp_cum)
  setorder(lab, model, bp_cum)
  lab[, cluster := cumsum(c(TRUE, diff(bp_cum) > 0.06 * genome_len)), by = model]
  lab[, `:=`(k = seq_len(.N) - 1L, top = max(logp)), by = .(model, cluster)]
  ymax <- max(d$logp, -log10(bonf))
  step <- ymax * 0.16
  lab[, label_y := top + ymax * 0.06 + k * step]
  ymax <- max(ymax, lab$label_y + step * 0.9) * 1.03

  ggplot(d, aes(bp_cum, logp)) +
    geom_point(data = d[class == "ns"], aes(colour = factor(Chr %% 2)), size = 0.9,
               show.legend = FALSE) +
    scale_colour_manual(values = c("#2166AC", "#67A9CF")) +
    geom_hline(yintercept = -log10(bonf), linetype = "dashed", linewidth = 0.6,
               colour = "#D7191C") +
    geom_point(data = d[class == "single"], shape = 21, size = 2.4, fill = "#FF9900",
               colour = "black", stroke = 0.4) +
    geom_point(data = d[class == "consensus_ns"], shape = 21, size = 2.8, fill = "white",
               colour = "#E31A1C", stroke = 0.9) +
    geom_point(data = d[class == "consensus"], shape = 21, size = 3.2, fill = "#E31A1C",
               colour = "black", stroke = 0.5) +
    geom_segment(data = lab, aes(xend = bp_cum, yend = label_y - step * 0.05),
                 colour = "grey50", linewidth = 0.25) +
    geom_label(data = lab, aes(y = label_y, label = label), vjust = 0, size = 2.5,
               lineheight = 0.9, label.size = 0, fill = alpha("white", 0.85),
               label.padding = unit(0.1, "lines")) +
    facet_wrap(~model, ncol = 1, strip.position = "right") +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$Chromosome,
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    labs(title = tr,
         subtitle = sprintf(
           "Red dashed: Bonferroni p = %.2e  |  red: significant in >= %d models (%d SNP%s; hollow = not significant in this model)  |  orange: single-model hits",
           bonf, params$min_models, length(cons_snps), ifelse(length(cons_snps) == 1, "", "s")),
         x = "Chromosome", y = expression(-log[10](italic(p)))) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"),
          panel.spacing = unit(0.8, "lines"))
}

traits <- sort(unique(res$trait))
pdf_all <- file.path(out_dir, "manhattan_all_traits.pdf")
pdf(pdf_all, width = 11, height = 8)
for (tr in traits) {
  p <- plot_trait(tr)
  ggsave(file.path(out_dir, paste0("manhattan_", tr, ".png")), p,
         width = 11, height = 8, dpi = 300, device = ragg::agg_png)
  print(p)
}
invisible(dev.off())
cat("Wrote", length(traits), "trait figures and", pdf_all, "\n")
