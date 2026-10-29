# Candidate genes within +/- window_bp of each consensus marker, annotated with
# Phytozome M. sinensis v7.1 functional descriptions (Arabidopsis / rice hits).

suppressPackageStartupMessages(library(data.table))
source("R/config.R")
dir.create(paths$candidates, recursive = TRUE, showWarnings = FALSE)

consensus <- fread(file.path(paths$consensus, "consensus_markers.csv"))
if (nrow(consensus) == 0) stop("No consensus markers - nothing to annotate.")

# Gene models (GFF3 'gene' features)
gff <- fread(cmd = sprintf("gzip -dc '%s' | grep -v '^#' | awk -F'\\t' '$3==\"gene\"'", paths$gff),
             header = FALSE, sep = "\t",
             col.names = c("seqid", "source", "type", "start", "end", "score",
                           "strand", "phase", "attributes"))
genes <- gff[grepl("^Chr", seqid), .(Chr = as.integer(sub("^Chr0?", "", seqid)),
                 gene = sub(".*Name=([^;]+).*", "\\1", attributes),
                 gene_start = start, gene_end = end, strand)][!is.na(Chr)]
cat("Gene models on chromosomes:", nrow(genes), "\n")

# One SNP-level window per unique consensus marker
snps <- consensus[, .(traits = paste(sort(unique(trait)), collapse = "; "),
                      models = paste(sort(unique(unlist(strsplit(models, "+", fixed = TRUE)))),
                                     collapse = "+"),
                      min_p = min(min_p)),
                  by = .(SNP, Chr, Pos)]
snps[, `:=`(win_start = pmax(1L, Pos - params$window_bp), win_end = Pos + params$window_bp)]

setkey(genes, Chr, gene_start, gene_end)
hits <- foverlaps(snps, genes, by.x = c("Chr", "win_start", "win_end"),
                  by.y = c("Chr", "gene_start", "gene_end"), nomatch = NULL)
hits[, distance_bp := fifelse(Pos >= gene_start & Pos <= gene_end, 0L,
                       fifelse(Pos < gene_start, gene_start - Pos, Pos - gene_end))]
hits[, location := fcase(distance_bp == 0, "within gene",
                         (strand == "+" & Pos < gene_start) | (strand == "-" & Pos > gene_end), "upstream",
                         default = "downstream")]

# Functional annotation (per transcript -> keep primary transcript per locus)
anno <- fread(cmd = sprintf("gzip -dc '%s'", paths$annotation), sep = "\t", quote = "")
setnames(anno, sub("^#", "", names(anno)))
anno <- anno[order(locusName, transcriptName)][!duplicated(locusName)]
clean_defline <- function(x) {
  x <- sub("^\\(\\d+ of \\d+\\)\\s*", "", x)
  x <- sub("^[A-Za-z0-9:._/]+\\s+-\\s+", "", x)
  trimws(x)
}
anno <- anno[, .(gene = locusName, Pfam, Panther, GO, KO,
                 arabidopsis_hit = `Best-hit-arabi-name`,
                 arabidopsis_function = clean_defline(`Best-hit-arabi-defline`),
                 rice_hit = `Best-hit-rice-name`,
                 rice_function = clean_defline(`Best-hit-rice-defline`))]

cand <- merge(hits, anno, by = "gene", all.x = TRUE)
cand[, putative_function := fcoalesce(
  fifelse(arabidopsis_function %in% c("", NA), NA_character_, arabidopsis_function),
  fifelse(rice_function %in% c("", NA), NA_character_, rice_function),
  "Uncharacterised")]
setorder(cand, Chr, Pos, distance_bp)
cand <- cand[, .(SNP, Chr, Pos, traits, models, min_p, gene, gene_start, gene_end,
                 strand, distance_bp, location, putative_function,
                 arabidopsis_hit, arabidopsis_function, rice_hit, rice_function,
                 Pfam, Panther, GO, KO)]
fwrite(cand, file.path(paths$candidates, "candidate_genes.csv"))

# Compact table: nearest gene per marker plus all genes in window
nearest <- cand[, .SD[1], by = SNP]
per_snp <- merge(snps[, .(SNP, Chr, Pos, traits, models, min_p)],
                 cand[, .(n_genes_in_window = .N,
                          genes_in_window = paste(gene, collapse = "; ")), by = SNP],
                 by = "SNP", all.x = TRUE)
per_snp <- merge(per_snp, nearest[, .(SNP, nearest_gene = gene, distance_bp, location,
                                      putative_function)], by = "SNP", all.x = TRUE)
per_snp[is.na(n_genes_in_window), n_genes_in_window := 0L]
setorder(per_snp, Chr, Pos)
fwrite(per_snp, file.path(paths$candidates, "candidate_genes_by_marker.csv"))

cat(sprintf("Window: +/- %d kb\n", params$window_bp / 1000))
cat("Consensus SNPs:", nrow(snps), "| with >= 1 gene in window:",
    sum(per_snp$n_genes_in_window > 0), "| candidate genes:", uniqueN(cand$gene), "\n\n")
options(width = 200)
print(per_snp[, .(SNP, traits, models, min_p = signif(min_p, 3), nearest_gene,
                  distance_bp, location, putative_function = substr(putative_function, 1, 60))])
