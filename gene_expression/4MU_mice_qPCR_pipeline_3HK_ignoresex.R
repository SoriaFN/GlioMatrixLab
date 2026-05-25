# ============================================================
# Generic qPCR analysis pipeline — THREE housekeeping genes
# Sample ID format: Group_Age_Sex_Replicate (sex collapsed)
# ΔCt (vs geomean of 3 HKs) → ΔΔCt → log2FC → limma stats → plots
#
# Input requirements:
#   - Excel file, one sheet
#   - Row 1 = header
#   - Column 1 = sample ID, format: Group_Age_Sex_Replicate
#                (e.g. Ctrl_young_M_1, Park_old_F_2)
#   - Columns 2–4 = the THREE housekeeping genes
#   - Columns 5+  = target genes
#
# Analysis design:
#   - Sex is ignored (M and F pooled within each Group_Age cell).
#   - Each Group_Age combination is one factor level; each is
#     contrasted against the reference "Ctrl_young".
#   - This gives clean pairwise contrasts but does NOT decompose
#     into main effects (group, age) and the group×age interaction.
#     If you ever want to ask "is the aging trajectory different in
#     Park vs Ctrl?", that's an interaction question — see the
#     commented factorial block at the end of section 5.
#
# Note on multi-HK normalization (Vandesompele 2002 / geNorm):
#   Normalization factor = geometric mean of HK expression values.
#   Because Ct is already on a log2 scale, the geometric mean of
#   linear-scale HK expression is equivalent to the arithmetic mean
#   of the HK Ct values. So we just average the 3 HK Cts per sample.
# ============================================================

# ---- 0. Packages ----
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggrepel)

# ---- 1. Settings ----
input_file <- file.choose()
ref_group  <- make.names("Ctrl_young")   # reference Group_Age combination
lod_ct     <- 30
alpha      <- 0.05
fc_cutoff  <- 1

out_dir <- dirname(input_file)
setwd(out_dir)

# ---- 2. Load & tidy ----
raw <- read_excel(input_file)
names(raw)[1] <- "sample"
# Case is preserved — group names like "Ctrl_young" are mixed-case.
# Make sure your Excel sample IDs use consistent casing.

# THREE housekeepers (columns 2:4), targets are columns 5+
housekeeping <- names(raw)[2:4]
target_genes <- names(raw)[-(1:4)]
cat("Housekeeping genes:", paste(housekeeping, collapse = " + "), "\n")
cat("Target genes (", length(target_genes), "):\n", sep = "")
cat(paste(target_genes, collapse = ", "), "\n\n")

# HK stability — Vandesompele geNorm M-value.
# For each HK j, compute pairwise log-ratios (= Ct differences, since
# Ct is log2) with the other HKs, take their SD across samples, then
# average. Rule of thumb: M < 0.5 = excellent (homogeneous samples);
# M < 1.0 = acceptable; M < 1.5 = acceptable for heterogeneous samples.
hk_ct <- as.data.frame(raw[, housekeeping])
M_values <- sapply(seq_along(housekeeping), function(j) {
  others <- setdiff(seq_along(housekeeping), j)
  pairwise_sd <- sapply(others, function(k) {
    sd(hk_ct[[j]] - hk_ct[[k]], na.rm = TRUE)
  })
  mean(pairwise_sd)
})
names(M_values) <- housekeeping
cat("HK stability (geNorm M-value, lower = more stable):\n")
print(round(M_values, 3))
if (any(M_values > 1.5)) {
  warning("One or more HKs have M > 1.5 — they may not be stably ",
          "expressed. Consider dropping the worst one and re-running, ",
          "or running NormFinder for a second opinion.")
}
cat("\n")

# Parse sample name → group, age, sex, replicate
meta <- raw %>%
  select(sample) %>%
  separate(sample, into = c("group", "age", "sex", "replicate"),
           sep = "_", remove = FALSE,
           extra = "merge", fill = "right")

bad_samples <- meta$sample[is.na(meta$replicate)]
if (length(bad_samples) > 0) {
  warning("These samples could not be parsed into Group_Age_Sex_Replicate ",
          "and will be dropped:\n  ", paste(bad_samples, collapse = "\n  "))
  meta <- meta %>% filter(!is.na(replicate))
  raw  <- raw  %>% filter(sample %in% meta$sample)
}

# Build combined Group_Age factor, ignoring sex
# make.names() ensures levels starting with a digit (e.g. "4MU") become
# valid R names (e.g. "X4MU"), which makeContrasts requires.
meta$group_age <- make.names(paste(meta$group, meta$age, sep = "_"))

all_lvls <- unique(meta$group_age)
if (!ref_group %in% all_lvls) {
  stop("ref_group '", ref_group, "' not found in sample names. ",
       "Available Group_Age combinations: ",
       paste(all_lvls, collapse = ", "))
}
other_lvls <- setdiff(all_lvls, ref_group)

# Sanity report — n per Group_Age cell, with sex breakdown for transparency
cat("Group_Age combinations (sex is shown for sanity but ignored in model):\n")
print(meta %>% count(group_age, sex) %>%
        pivot_wider(names_from = sex, values_from = n, values_fill = 0))
cat("\nReference level:", ref_group, "\n\n")

meta$group_age <- factor(meta$group_age,
                         levels = c(ref_group, other_lvls))

# Long format
long <- raw %>%
  pivot_longer(-sample, names_to = "gene", values_to = "Ct") %>%
  left_join(meta, by = "sample")

# ---- 3. ΔCt: subtract mean of the three housekeepers ----
hk <- long %>%
  filter(gene %in% housekeeping) %>%
  group_by(sample) %>%
  summarise(Ct_hk = mean(Ct, na.rm = TRUE), .groups = "drop")

targets <- long %>%
  filter(!gene %in% housekeeping) %>%
  left_join(hk, by = "sample") %>%
  mutate(dCt = Ct - Ct_hk,
         near_lod = Ct >= lod_ct)

# ---- 4. ΔΔCt vs reference (Ctrl_young), log2FC ----
ref_means <- targets %>%
  filter(group_age == ref_group) %>%
  group_by(gene) %>%
  summarise(dCt_ref = mean(dCt, na.rm = TRUE), .groups = "drop")

targets <- targets %>%
  left_join(ref_means, by = "gene") %>%
  mutate(ddCt   = dCt - dCt_ref,
         log2FC = -ddCt,
         FC     = 2^(-ddCt))

# ---- 5. Moderated stats with limma ----
# trend = FALSE: too few genes to reliably estimate a mean-variance trend.
expr_mat <- targets %>%
  mutate(neg_dCt = -dCt) %>%
  select(sample, gene, neg_dCt) %>%
  pivot_wider(names_from = sample, values_from = neg_dCt) %>%
  as.data.frame()
rownames(expr_mat) <- expr_mat$gene; expr_mat$gene <- NULL

meta_ord <- meta[match(colnames(expr_mat), meta$sample), ]
design <- model.matrix(~ 0 + group_age, data = meta_ord)
colnames(design) <- levels(meta_ord$group_age)

fit <- lmFit(as.matrix(expr_mat), design)

contrast_strs <- paste0(other_lvls, " - ", ref_group)
names(contrast_strs) <- paste0(other_lvls, "_vs_", ref_group)
contr <- makeContrasts(contrasts = contrast_strs, levels = design)
colnames(contr) <- names(contrast_strs)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = FALSE)

limma_res <- bind_rows(lapply(colnames(contr), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH",
                 sort.by = "none")
  tibble(gene = rownames(tt),
         contrast = cn,
         logFC = tt$logFC,
         t = tt$t,
         p = tt$P.Value,
         p.adj = tt$adj.P.Val)
}))

write.csv(limma_res, "qPCR_limma_results.csv", row.names = FALSE)

# --- Alternative designs (uncomment if you want them) ---
#
# (a) Block out sex without testing it — cheap insurance if sex shifts
#     baseline expression. Same contrasts as above but lower residual SD:
#       design <- model.matrix(~ 0 + group_age + sex, data = meta_ord)
#       colnames(design)[1:length(levels(meta_ord$group_age))] <-
#         levels(meta_ord$group_age)
#     (then build contrasts only over the group_age columns as above)
#
# (b) Factorial model — decompose into main effects + interaction:
#       design_fac <- model.matrix(~ group * age, data = meta_ord)
#       fit_fac    <- lmFit(as.matrix(expr_mat), design_fac)
#       fit_fac    <- eBayes(fit_fac, trend = FALSE)
#       topTable(fit_fac, coef = "groupPark:ageold", ...)   # interaction
#     Use this when you want to ask "does aging affect Park differently
#     than Ctrl?" rather than just "is Park_old different from Ctrl_young?"

# ---- 6. Wide log2FC table for heatmap + export ----
hm_mat <- limma_res %>%
  mutate(comp = sub(paste0("_vs_", ref_group, "$"), "", contrast)) %>%
  select(gene, comp, logFC) %>%
  pivot_wider(names_from = comp, values_from = logFC) %>%
  as.data.frame()
rownames(hm_mat) <- hm_mat$gene; hm_mat$gene <- NULL
hm_mat <- hm_mat[, other_lvls, drop = FALSE]

write.csv(hm_mat, "log2FC_matrix.csv")

# ---- 7. Heatmap ----
na_genes <- rownames(hm_mat)[apply(hm_mat, 1, function(x) any(is.na(x)))]
if (length(na_genes) > 0) {
  message("Heatmap: dropping ", length(na_genes),
          " gene(s) with NA logFC (missing data in some groups): ",
          paste(na_genes, collapse = ", "))
  hm_mat_plot <- hm_mat[!rownames(hm_mat) %in% na_genes, , drop = FALSE]
} else {
  hm_mat_plot <- hm_mat
}

max_abs <- max(abs(hm_mat_plot), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)

pheatmap(hm_mat_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         fontsize_row = 8,
         main = paste0("log2FC vs ", ref_group,
                       "  (HK: ", paste(housekeeping, collapse = "+"), ")"),
         filename = "heatmap_log2FC.pdf",
         width  = max(5, 1.2 * length(other_lvls) + 3),
         height = 9)

# ---- 8. Volcano plot ----
volcano_df <- limma_res %>%
  mutate(comp = sub(paste0("_vs_", ref_group, "$"), "", contrast),
         neglog10p = -log10(p.adj),
         sig = case_when(
           p.adj < alpha & logFC >  fc_cutoff ~ "Up",
           p.adj < alpha & logFC < -fc_cutoff ~ "Down",
           TRUE ~ "n.s."),
         sig = factor(sig, levels = c("Up", "Down", "n.s.")))

v <- ggplot(volcano_df, aes(logFC, neglog10p)) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
             linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = sig), size = 2, alpha = 0.8) +
  geom_text_repel(data = subset(volcano_df, sig != "n.s."),
                  aes(label = gene), size = 3,
                  max.overlaps = Inf, box.padding = 0.3,
                  segment.colour = "grey50", segment.size = 0.3) +
  scale_colour_manual(values = c(Up = "#C0392B", Down = "#2874A6", n.s. = "grey70")) +
  facet_wrap(~ comp, nrow = 1) +
  labs(x = bquote(log[2]~"fold change vs"~.(ref_group)),
       y = expression(-log[10]~adjusted~p~(limma)),
       colour = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

ggsave("volcano_vs_ref.pdf", v,
       width = max(3, 2.2 * length(other_lvls)), height = 4)

# ---- 9. Per-gene plots for significant hits ----
sig_genes <- limma_res %>%
  filter(p.adj < alpha) %>%
  pull(gene) %>% unique()

if (length(sig_genes) > 0) {
  p <- targets %>%
    filter(gene %in% sig_genes) %>%
    ggplot(aes(group_age, log2FC, fill = group_age)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.8, shape = 21) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(y = expression(log[2]~fold~change), x = NULL) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave("significant_genes.pdf", p,
         width  = min(12, 2.5 * ceiling(sqrt(length(sig_genes)))),
         height = min(10, 2.2 * ceiling(sqrt(length(sig_genes)))))
}

# ---- 10. Console report ----
cat("\n--- limma significant hits (p.adj <", alpha, ") ---\n")
sig_tbl <- limma_res %>% filter(p.adj < alpha) %>% arrange(p.adj)
if (nrow(sig_tbl) == 0) {
  cat("No genes pass p.adj <", alpha, ". Top 10 by adjusted p:\n")
  print(limma_res %>% arrange(p.adj) %>% head(10))
} else {
  print(sig_tbl)
}

cat("\nFiles written to:", out_dir, "\n")
cat("  qPCR_limma_results.csv  log2FC_matrix.csv\n")
cat("  heatmap_log2FC.pdf      volcano_vs_ref.pdf\n")
if (length(sig_genes) > 0) cat("  significant_genes.pdf\n")
