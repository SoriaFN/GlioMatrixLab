# ============================================================
# qPCR analysis — 2×2×2 factorial (Treatment × Age × Sex)
# Three housekeeping genes (ACTB, HPRT1, TBP)
# Sample ID format: Group_Age_Sex_Replicate
#   e.g. Ctrl_young_M_1, 4MU_aged_F_3
#
# Design philosophy
# -----------------
# We use cell-means parameterization: one factor `cell` with 8 levels
# (Ctrl_young_M, Ctrl_young_F, Ctrl_aged_M, Ctrl_aged_F,
#  4MU_young_M,  4MU_young_F,  4MU_aged_M,  4MU_aged_F).
# This is equivalent to ~ group*age*sex in terms of fit, but lets us
# build only the contrasts we actually care about — better power than
# blindly testing every coefficient of the full factorial.
#
# Contrasts computed (organized by question):
#   A. MAIN EFFECTS (averaged over the other two factors)
#        Main_Age          : aged   vs young
#        Main_4MU          : 4MU    vs Ctrl
#        Main_Sex          : F      vs M
#   B. 2-WAY INTERACTIONS (delta-of-delta)
#        Int_4MUxAge       : does 4MU's effect change with age?  *** key biological question
#        Int_4MUxSex       : does 4MU's effect differ between sexes?
#        Int_AgexSex       : does the aging effect differ between sexes?
#   C. STRATIFIED 4MU EFFECTS (drill-down for when interactions are significant)
#        4MU_vs_Ctrl_in_young_M / _young_F / _aged_M / _aged_F
#   D. STRATIFIED AGING EFFECTS
#        Aged_vs_Young_in_Ctrl_M / _Ctrl_F / _4MU_M / _4MU_F
#
# Note: with n = 4–6 per cell and 47 genes, main effects and 2-way
# interactions are reasonably powered; 3-way interactions would be
# underpowered and are deliberately omitted. eBayes shrinkage is
# weaker than with whole-transcriptome data but still helps.
# ============================================================

# ---- 0. Packages ----
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggrepel)

# ---- 1. Settings ----
input_file <- file.choose()
lod_ct      <- 30
alpha       <- 0.05
fc_cutoff   <- 1
hk_outlier_z <- 3   # flag samples whose mean HK Ct is > z SDs from the bulk
drop_hk_outliers <- FALSE  # set TRUE to auto-drop flagged samples

out_dir <- dirname(input_file)
setwd(out_dir)

# ---- 2. Load & tidy ----
raw <- read_excel(input_file)
names(raw)[1] <- "sample"
# Case-sensitive — make sure Excel sample IDs use consistent casing.

# Three housekeepers (cols 2:4), targets cols 5+
housekeeping <- names(raw)[2:4]
target_genes <- names(raw)[-(1:4)]
cat("Housekeeping genes:", paste(housekeeping, collapse = " + "), "\n")
cat("Target genes (", length(target_genes), "):\n", sep = "")
cat(paste(target_genes, collapse = ", "), "\n\n")

# ---- 2b. HK QC: stability (geNorm M-value) ----
# For each HK j, compute pairwise log-ratios (= Ct diffs, since Ct is log2)
# with each other HK, take their SD across samples, then average.
# M < 0.5 excellent (homogeneous), < 1.0 acceptable, < 1.5 acceptable (heterogeneous).
hk_ct <- as.data.frame(raw[, housekeeping])
M_values <- sapply(seq_along(housekeeping), function(j) {
  others <- setdiff(seq_along(housekeeping), j)
  mean(sapply(others, function(k) sd(hk_ct[[j]] - hk_ct[[k]], na.rm = TRUE)))
})
names(M_values) <- housekeeping
cat("HK stability (geNorm M-value, lower = more stable):\n")
print(round(M_values, 3))
if (any(M_values > 1.5)) {
  warning("One or more HKs have M > 1.5 — they may not be stably expressed. ",
          "Consider dropping the worst one and re-running.")
}
cat("\n")

# ---- 2c. HK QC: per-sample outlier detection ----
# Samples with uniformly high HK Cts across all 3 HKs are very likely
# poor preps (low input / degraded RNA). Geometric-mean normalization
# CANNOT rescue such samples if the inter-HK ratios are also distorted.
sample_hk_mean <- rowMeans(hk_ct, na.rm = TRUE)
hk_bulk_mean   <- mean(sample_hk_mean, na.rm = TRUE)
hk_bulk_sd     <- sd(sample_hk_mean,   na.rm = TRUE)
sample_hk_z    <- (sample_hk_mean - hk_bulk_mean) / hk_bulk_sd

hk_flag <- data.frame(
  sample   = raw$sample,
  mean_HK  = round(sample_hk_mean, 2),
  z        = round(sample_hk_z, 2)
)
flagged <- hk_flag[abs(hk_flag$z) > hk_outlier_z, ]
if (nrow(flagged) > 0) {
  cat("Samples flagged as HK outliers (|z| >", hk_outlier_z, "):\n")
  print(flagged, row.names = FALSE)
  if (drop_hk_outliers) {
    cat("→ Dropping flagged samples (drop_hk_outliers = TRUE)\n")
    raw <- raw[!raw$sample %in% flagged$sample, ]
  } else {
    cat("→ Keeping them in the analysis (set drop_hk_outliers = TRUE to drop)\n")
  }
} else {
  cat("No HK outlier samples flagged.\n")
}
cat("\n")

# ---- 2d. Parse sample IDs → group / age / sex / replicate ----
meta <- raw %>%
  select(sample) %>%
  separate(sample, into = c("group", "age", "sex", "replicate"),
           sep = "_", remove = FALSE, extra = "merge", fill = "right")

bad_samples <- meta$sample[is.na(meta$replicate)]
if (length(bad_samples) > 0) {
  warning("Unparseable samples (dropped):\n  ",
          paste(bad_samples, collapse = "\n  "))
  meta <- meta %>% filter(!is.na(replicate))
  raw  <- raw  %>% filter(sample %in% meta$sample)
}

# Build the 8-level cell factor (treatment × age × sex)
meta$cell <- paste(meta$group, meta$age, meta$sex, sep = "_")
meta$cell <- sub("^4MU", "X4MU", meta$cell)

cat("Cell counts (n per group × age × sex):\n")
print(meta %>% count(group, age, sex) %>%
        pivot_wider(names_from = sex, values_from = n, values_fill = 0))
cat("\n")

# Set factor level order with a sensible reference first (Ctrl_young_M)
ref_cell   <- "Ctrl_young_M"
all_cells  <- sort(unique(meta$cell))
if (!ref_cell %in% all_cells) stop("Reference cell '", ref_cell, "' not found.")
meta$cell  <- factor(meta$cell, levels = c(ref_cell, setdiff(all_cells, ref_cell)))

# ---- 2e. HK Ct per cell (does HK expression itself vary with the design?) ----
# If HKs differ systematically across experimental cells, normalization
# absorbs real biological signal. Print mean ± SD per cell for inspection.
cat("Mean HK Ct per cell (eyeball for systematic shifts):\n")
hk_per_cell <- raw %>%
  mutate(cell = meta$cell[match(sample, meta$sample)]) %>%
  group_by(cell) %>%
  summarise(across(all_of(housekeeping),
                   list(m = ~mean(.x, na.rm = TRUE), s = ~sd(.x, na.rm = TRUE))),
            .groups = "drop")
print(hk_per_cell)
cat("\n")

# ---- 3. Long format + ΔCt ----
long <- raw %>%
  pivot_longer(-sample, names_to = "gene", values_to = "Ct") %>%
  left_join(meta, by = "sample")

hk <- long %>%
  filter(gene %in% housekeeping) %>%
  group_by(sample) %>%
  summarise(Ct_hk = mean(Ct, na.rm = TRUE), .groups = "drop")

targets <- long %>%
  filter(!gene %in% housekeeping) %>%
  left_join(hk, by = "sample") %>%
  mutate(dCt = Ct - Ct_hk,
         near_lod = Ct >= lod_ct)

# ---- 4. ΔΔCt vs Ctrl_young (sex-pooled within cells of interest) ----
# Note: this ΔΔCt is for visualization only — limma below does the
# proper modeling with cell-means. Reference for ΔΔCt = mean of both
# Ctrl_young cells (M + F pooled).
ref_means <- targets %>%
  filter(group == "Ctrl", age == "young") %>%
  group_by(gene) %>%
  summarise(dCt_ref = mean(dCt, na.rm = TRUE), .groups = "drop")

targets <- targets %>%
  left_join(ref_means, by = "gene") %>%
  mutate(ddCt   = dCt - dCt_ref,
         log2FC = -ddCt,
         FC     = 2^(-ddCt))

# ---- 5. limma fit: cell-means model ----
expr_mat <- targets %>%
  mutate(neg_dCt = -dCt) %>%
  select(sample, gene, neg_dCt) %>%
  pivot_wider(names_from = sample, values_from = neg_dCt) %>%
  as.data.frame()
rownames(expr_mat) <- expr_mat$gene; expr_mat$gene <- NULL

meta_ord <- meta[match(colnames(expr_mat), meta$sample), ]
design <- model.matrix(~ 0 + cell, data = meta_ord)
colnames(design) <- levels(meta_ord$cell)

fit <- lmFit(as.matrix(expr_mat), design)

# ---- 6. Build curated contrasts ----
# Level names start with "4MU" (digit) → need backticks for makeContrasts to parse.
bq <- function(x) paste0("`", x, "`")
mean_str <- function(cells) {
  paste0("(", paste(bq(cells), collapse = " + "), ")/", length(cells))
}

cells   <- levels(meta_ord$cell)
aged    <- grep("_aged_",  cells, value = TRUE)
young   <- grep("_young_", cells, value = TRUE)
mu      <- grep("^X4MU_",   cells, value = TRUE)
ctrl    <- grep("^Ctrl_",  cells, value = TRUE)
female  <- grep("_F$",     cells, value = TRUE)
male    <- grep("_M$",     cells, value = TRUE)

# A) Main effects — equal weight per cell (tests the "balanced" main effect)
me_age <- paste(mean_str(aged),   "-", mean_str(young))
me_4mu <- paste(mean_str(mu),     "-", mean_str(ctrl))
me_sex <- paste(mean_str(female), "-", mean_str(male))

# B) 2-way interactions — delta-of-delta
int_4mu_age <- paste(
  "(", mean_str(intersect(mu,   aged)), "-", mean_str(intersect(mu,   young)), ")",
  "-",
  "(", mean_str(intersect(ctrl, aged)), "-", mean_str(intersect(ctrl, young)), ")"
)
int_4mu_sex <- paste(
  "(", mean_str(intersect(mu,   female)), "-", mean_str(intersect(mu,   male)), ")",
  "-",
  "(", mean_str(intersect(ctrl, female)), "-", mean_str(intersect(ctrl, male)), ")"
)
int_age_sex <- paste(
  "(", mean_str(intersect(aged,  female)), "-", mean_str(intersect(aged,  male)), ")",
  "-",
  "(", mean_str(intersect(young, female)), "-", mean_str(intersect(young, male)), ")"
)

# C) Stratified 4MU effects (one per Age × Sex stratum)
strat_4mu <- list(
  "X4MU_vs_Ctrl_in_young_M" = paste(bq("X4MU_young_M"), "-", bq("Ctrl_young_M")),
  "X4MU_vs_Ctrl_in_young_F" = paste(bq("X4MU_young_F"), "-", bq("Ctrl_young_F")),
  "X4MU_vs_Ctrl_in_aged_M"  = paste(bq("X4MU_aged_M"),  "-", bq("Ctrl_aged_M")),
  "X4MU_vs_Ctrl_in_aged_F"  = paste(bq("X4MU_aged_F"),  "-", bq("Ctrl_aged_F"))
)

# D) Stratified aging effects (one per Treatment × Sex stratum)
strat_age <- list(
  "Aged_vs_Young_in_Ctrl_M" = paste(bq("Ctrl_aged_M"), "-", bq("Ctrl_young_M")),
  "Aged_vs_Young_in_Ctrl_F" = paste(bq("Ctrl_aged_F"), "-", bq("Ctrl_young_F")),
  "Aged_vs_Young_in_X4MU_M"  = paste(bq("X4MU_aged_M"),  "-", bq("X4MU_young_M")),
  "Aged_vs_Young_in_X4MU_F"  = paste(bq("X4MU_aged_F"),  "-", bq("X4MU_young_F"))
)

all_contrasts <- c(
  list(Main_Age    = me_age,
       Main_4MU    = me_4mu,
       Main_Sex    = me_sex,
       Int_4MUxAge = int_4mu_age,
       Int_4MUxSex = int_4mu_sex,
       Int_AgexSex = int_age_sex),
  strat_4mu,
  strat_age
)

contrast_category <- c(
  rep("Main effect",  3),
  rep("Interaction",  3),
  rep("Stratified 4MU effect", 4),
  rep("Stratified aging effect", 4)
)
names(contrast_category) <- names(all_contrasts)

contr <- makeContrasts(contrasts = unlist(all_contrasts), levels = design)
colnames(contr) <- names(all_contrasts)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2, trend = FALSE)
# trend = FALSE: 47 genes is too few to fit a mean-variance trend reliably.
# If you ever scale up to a much larger panel, set trend = TRUE.

# ---- 7. Extract all results ----
limma_res <- bind_rows(lapply(colnames(contr), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, adjust.method = "BH", sort.by = "none")
  tibble(gene       = rownames(tt),
         contrast   = cn,
         category   = contrast_category[cn],
         logFC      = tt$logFC,
         t          = tt$t,
         p          = tt$P.Value,
         p.adj      = tt$adj.P.Val)
}))

# Multiple-testing note: limma's adj.P is BH-corrected within each contrast
# (across genes). If you want family-wise correction across contrasts too,
# apply a second BH on `p` grouped by category (uncomment below):
# limma_res <- limma_res %>%
#   group_by(category) %>%
#   mutate(p.adj_within_cat = p.adjust(p, method = "BH")) %>%
#   ungroup()

write.csv(limma_res, "qPCR_limma_results.csv", row.names = FALSE)

# ---- 8. Wide log2FC matrix for heatmaps ----
hm_mat <- limma_res %>%
  select(gene, contrast, logFC) %>%
  pivot_wider(names_from = contrast, values_from = logFC) %>%
  as.data.frame()
rownames(hm_mat) <- hm_mat$gene; hm_mat$gene <- NULL

# Order columns by category for visual grouping
col_order <- names(all_contrasts)
hm_mat <- hm_mat[, col_order, drop = FALSE]
write.csv(hm_mat, "log2FC_matrix.csv")

# Annotation bar for the heatmap
ann_col <- data.frame(Category = factor(contrast_category[col_order],
                                        levels = unique(contrast_category)))
rownames(ann_col) <- col_order
ann_colors <- list(Category = c(
  "Main effect"            = "#2C3E50",
  "Interaction"            = "#E67E22",
  "Stratified 4MU effect"  = "#C0392B",
  "Stratified aging effect" = "#27AE60"
))

# ---- 9. Heatmap (all contrasts side by side) ----
na_genes <- rownames(hm_mat)[apply(hm_mat, 1, function(x) any(is.na(x)))]
hm_plot <- if (length(na_genes) > 0) hm_mat[!rownames(hm_mat) %in% na_genes, ] else hm_mat

max_abs <- max(abs(hm_plot), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)

pheatmap(hm_plot,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = breaks,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         fontsize_row = 7,
         fontsize_col = 8,
         angle_col = 45,
         main = paste0("log2FC across contrasts  (HK: ",
                       paste(housekeeping, collapse = "+"), ")"),
         filename = "heatmap_all_contrasts.pdf",
         width  = max(8, 0.45 * ncol(hm_plot) + 4),
         height = max(7, 0.18 * nrow(hm_plot) + 2))

# ---- 10. Volcano plots, faceted by category ----
volcano_df <- limma_res %>%
  mutate(neglog10p = -log10(p.adj),
         sig = case_when(
           p.adj < alpha & logFC >  fc_cutoff ~ "Up",
           p.adj < alpha & logFC < -fc_cutoff ~ "Down",
           TRUE ~ "n.s."),
         sig = factor(sig, levels = c("Up", "Down", "n.s.")),
         contrast = factor(contrast, levels = names(all_contrasts)))

for (cat_i in unique(contrast_category)) {
  df_i <- volcano_df %>% filter(category == cat_i)
  n_facets <- length(unique(df_i$contrast))

  v <- ggplot(df_i, aes(logFC, neglog10p)) +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed", colour = "grey60") +
    geom_point(aes(colour = sig), size = 2, alpha = 0.8) +
    geom_text_repel(data = subset(df_i, sig != "n.s."),
                    aes(label = gene), size = 3,
                    max.overlaps = Inf, box.padding = 0.3,
                    segment.colour = "grey50", segment.size = 0.3) +
    scale_colour_manual(values = c(Up = "#C0392B", Down = "#2874A6", n.s. = "grey70")) +
    facet_wrap(~ contrast, nrow = 1) +
    labs(title = cat_i,
         x = expression(log[2]~"fold change"),
         y = expression(-log[10]~adjusted~p~(limma)),
         colour = NULL) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "top")

  fname <- paste0("volcano_", gsub("[^A-Za-z0-9]+", "_", tolower(cat_i)), ".pdf")
  ggsave(fname, v,
         width  = max(3, 2.2 * n_facets),
         height = 5.5,
         limitsize = FALSE)
}

# ---- 11. Per-gene boxplots for genes hit in any contrast ----
sig_genes <- limma_res %>% filter(p.adj < alpha) %>% pull(gene) %>% unique()

if (length(sig_genes) > 0) {
  p <- targets %>%
    filter(gene %in% sig_genes) %>%
    mutate(age = factor(age, levels = c("young", "aged")),
           group = factor(group, levels = c("Ctrl", "4MU")),
           sex = factor(sex, levels = c("M", "F"))) %>%
    ggplot(aes(interaction(group, age), log2FC, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_boxplot(alpha = 0.4, outlier.shape = NA) +
    geom_point(aes(shape = sex),
               position = position_jitterdodge(jitter.width = 0.15,
                                                dodge.width = 0.75),
               size = 1.8) +
    facet_wrap(~ gene, scales = "free_y") +
    scale_fill_manual(values = c(Ctrl = "#7F8C8D", `4MU` = "#E67E22")) +
    scale_shape_manual(values = c(M = 16, F = 17)) +
    labs(y = expression(log[2]~fold~change~vs~Ctrl[young]),
         x = NULL) +
    theme_classic(base_size = 10) +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top")

  ggsave("significant_genes.pdf", p,
         width  = min(14, 2.6 * ceiling(sqrt(length(sig_genes)))),
         height = min(12, 2.4 * ceiling(sqrt(length(sig_genes)))),
         limitsize = FALSE)
}

# ---- 12. Console summary ----
cat("\n=== Hits per contrast (p.adj <", alpha, ") ===\n")
hit_counts <- limma_res %>%
  group_by(category, contrast) %>%
  summarise(n_hits = sum(p.adj < alpha, na.rm = TRUE), .groups = "drop")
print(hit_counts)

cat("\n=== Top hits overall (lowest p.adj) ===\n")
print(limma_res %>% arrange(p.adj) %>% head(15))

cat("\nFiles written to:", out_dir, "\n")
cat("  qPCR_limma_results.csv     log2FC_matrix.csv\n")
cat("  heatmap_all_contrasts.pdf\n")
cat("  volcano_*.pdf              significant_genes.pdf\n")
