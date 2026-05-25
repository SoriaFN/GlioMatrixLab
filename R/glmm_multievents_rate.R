################################################################################
#
# GLMM analysis: Is a subset of events increased beyond what is expected
# from the global increase in total events?
#
# Federico Soria — Achucarro Basque Center for Neuroscience
# https://github.com/[your-repo]
#
# ──────────────────────────────────────────────────────────────────────────────
#
# PURPOSE
# -------
# In many experimental settings (calcium imaging, electrophysiology, behavioral
# scoring), a treatment increases total event counts. A natural follow-up 
# question is whether a specific *subtype* of events (e.g., synchronized,
# bursting, or long-duration events) is disproportionately increased — or 
# whether the increase merely tracks the overall rise in total activity.
#
# This script fits three Poisson Generalized Linear Mixed Models (GLMMs) to
# answer that question:
#
#   Model 1 — Offset model (primary)
#     SubEvents ~ Treatment + (1 | Experiment), offset = log(TotalEvents)
#
#     By placing log(TotalEvents) as an offset, we model the LOG-RATE:
#       log(E[SubEvents] / TotalEvents) = β0 + β_Treatment
#     A significant Treatment effect here means the RATE of sub-events per
#     total event differs across conditions — i.e., the increase in sub-events
#     is NOT simply a byproduct of more total events.
#
#   Model 2 — Covariate model
#     SubEvents ~ Treatment + TotalEvents + (1 | Experiment)
#
#     TotalEvents enters as a fixed-effect predictor rather than an offset.
#     This does not assume proportionality and is more flexible. A significant
#     Treatment term means treatment predicts sub-events even after adjusting
#     for total event count.
#
#   Model 3 — Interaction model (exploratory)
#     SubEvents ~ Treatment * TotalEvents + (1 | Experiment)
#
#     Tests whether the relationship between total and sub-events differs by
#     treatment. A significant interaction means the treatment changes the
#     coupling between overall activity and the sub-event type.
#
# The random intercept (1 | Experiment) handles paired/repeated-measures
# designs where the same experimental unit (slice, animal, session) is
# observed under multiple conditions.
#
# ──────────────────────────────────────────────────────────────────────────────
#
# INPUT
# -----
# An Excel (.xlsx) or CSV (.csv) file with at least 4 columns:
#   1. Experiment ID   — identifies the experimental unit (e.g., slice, animal)
#   2. Treatment/Group — the experimental condition
#   3. Sub-events      — count of the event subtype of interest
#   4. Total events    — count of all events
#
# The script will display the column names and prompt you to identify each one.
#
# OUTPUT (saved to the same directory as the input file)
# ------
#   - rate_by_treatment.pdf/.png  — dot plot of sub-event rates by treatment
#   - subevents_vs_total.pdf/.png — scatter plot of sub-events vs total events
#   - GLMM_summary.md             — Markdown report with formatted results tables
#   - Console output with full model summaries
#
# ──────────────────────────────────────────────────────────────────────────────
#
# REQUIREMENTS
# ------------
# install.packages(c("readxl", "lme4", "car", "emmeans", "ggplot2"))
#
################################################################################

# ============================================================================
# 0. LOAD LIBRARIES
# ============================================================================

library(readxl)
library(lme4)
library(car)
library(emmeans)
library(ggplot2)

# ============================================================================
# 1. DATA LOADING & COLUMN ASSIGNMENT
# ============================================================================

# --- Select file ---
filepath <- file.choose()
outdir   <- dirname(filepath)

# Read depending on extension
ext <- tolower(tools::file_ext(filepath))
if (ext %in% c("xlsx", "xls")) {
  dat_raw <- read_excel(filepath, sheet = 1)
} else if (ext == "csv") {
  dat_raw <- read.csv(filepath, stringsAsFactors = FALSE)
} else {
  stop("Unsupported file format. Please provide an .xlsx or .csv file.")
}

# --- Show columns and let user assign them ---
cat("====================================================================\n")
cat("Columns found in your file:\n")
cat("====================================================================\n\n")
for (i in seq_along(colnames(dat_raw))) {
  cat(sprintf("  [%d] %s\n", i, colnames(dat_raw)[i]))
}
cat("\n")

cat("Please identify which column corresponds to each variable.\n")
cat("Enter the column NUMBER for each prompt below.\n\n")

col_exp  <- as.integer(readline("  Experiment / Subject ID column: "))
col_trt  <- as.integer(readline("  Treatment / Group column:       "))
col_sub  <- as.integer(readline("  Sub-event count column:         "))
col_tot  <- as.integer(readline("  Total event count column:       "))

# Build a clean data frame with standardized names
dat <- data.frame(
  Experiment  = factor(dat_raw[[col_exp]]),
  Treatment   = trimws(as.character(dat_raw[[col_trt]])),
  SubEvents   = as.numeric(dat_raw[[col_sub]]),
  TotalEvents = as.numeric(dat_raw[[col_tot]])
)

# Store original column names for plot labels
orig_sub_label <- colnames(dat_raw)[col_sub]
orig_tot_label <- colnames(dat_raw)[col_tot]

# --- Set treatment factor levels ---
cat("\nTreatment levels detected:\n")
levs <- unique(dat$Treatment)
for (i in seq_along(levs)) cat(sprintf("  [%d] %s\n", i, levs[i]))

ref_idx <- as.integer(readline("\nWhich is the REFERENCE / control level? Enter number: "))
ref_level <- levs[ref_idx]
other_levels <- levs[levs != ref_level]
dat$Treatment <- factor(dat$Treatment, levels = c(ref_level, other_levels))

cat(sprintf("\nReference level set to: %s\n", ref_level))

# --- Summary ---
cat("\n====================================================================\n")
cat("Dataset loaded successfully\n")
cat("====================================================================\n\n")
print(as.data.frame(dat))
cat(sprintf("\n%d observations, %d experiments, %d treatments\n",
            nrow(dat), nlevels(dat$Experiment), nlevels(dat$Treatment)))

cat("\n=== Descriptives by Treatment ===\n")
print(aggregate(cbind(SubEvents, TotalEvents) ~ Treatment, data = dat,
                FUN = function(x) c(mean = mean(x), sd = sd(x))))

# ============================================================================
# 2. HANDLE ZERO TOTAL EVENTS
# ============================================================================
# Observations with zero total events cannot contribute to the offset model
# (log(0) is undefined). They are excluded from the primary analysis. If any
# exist, a sensitivity analysis using log(TotalEvents + 1) is run on the full
# dataset to verify that exclusion does not change the conclusions.

n_zero <- sum(dat$TotalEvents == 0)
cat(sprintf("\nObservations with zero total events: %d\n", n_zero))

if (n_zero > 0) {
  cat("These will be excluded from the offset model.\n")
  cat("A sensitivity analysis with log(TotalEvents + 1) is included below.\n")
}

dat_nozero <- dat[dat$TotalEvents > 0, ]

# ============================================================================
# 3. MODEL 1: OFFSET MODEL (primary analysis)
# ============================================================================
# By including log(TotalEvents) as an offset, we effectively model:
#   log(E[SubEvents] / TotalEvents) = β0 + β_Treatment
# This tests whether the RATE of sub-events per total event differs by
# treatment, which directly addresses whether the increase in sub-events is
# proportional to (and therefore explained by) the increase in total events.

cat("\n====================================================================\n")
cat("MODEL 1: GLMM with offset — rate of sub-events per total event\n")
cat("====================================================================\n\n")

m1 <- glmer(SubEvents ~ Treatment + (1 | Experiment),
            data = dat_nozero,
            family = poisson,
            offset = log(TotalEvents))

print(summary(m1))

cat("\n--- Type II Wald chi-square test for Treatment ---\n")
print(Anova(m1, type = "II"))

cat("\n--- Estimated marginal rates (sub-events per total event) ---\n")
emm1 <- emmeans(m1, ~ Treatment, type = "response", offset = 0)
print(emm1)

cat("\n--- Pairwise contrasts (ratio of rates) ---\n")
print(pairs(emm1, type = "response"))

# ============================================================================
# 4. MODEL 2: COVARIATE MODEL
# ============================================================================
# TotalEvents is included as a fixed-effect covariate rather than as an offset.
# This is more flexible because it does not assume a strict proportional
# relationship between sub-events and total events. The Treatment effect here
# means: "after adjusting for total event count, does treatment still predict
# sub-event count?"

cat("\n====================================================================\n")
cat("MODEL 2: GLMM with TotalEvents as covariate\n")
cat("====================================================================\n\n")

m2 <- glmer(SubEvents ~ Treatment + scale(TotalEvents) + (1 | Experiment),
            data = dat_nozero,
            family = poisson)

print(summary(m2))
cat("\n--- Type II Wald chi-square test ---\n")
print(Anova(m2, type = "II"))

# ============================================================================
# 5. MODEL 3: INTERACTION MODEL (exploratory)
# ============================================================================
# Adds a Treatment × TotalEvents interaction. A significant interaction would
# mean that the relationship between total events and sub-events itself differs
# across treatments — i.e., the treatment changes how strongly overall activity
# drives the sub-event type.

cat("\n====================================================================\n")
cat("MODEL 3: GLMM with Treatment x TotalEvents interaction\n")
cat("====================================================================\n\n")

m3 <- glmer(SubEvents ~ Treatment * scale(TotalEvents) + (1 | Experiment),
            data = dat_nozero,
            family = poisson)

print(summary(m3))
cat("\n--- Type II Wald chi-square test ---\n")
print(Anova(m3, type = "II"))

cat("\n--- LRT: Model 2 vs Model 3 (is the interaction needed?) ---\n")
print(anova(m2, m3))

# ============================================================================
# 6. SENSITIVITY: OFFSET MODEL WITH +1 (only if zero-total obs exist)
# ============================================================================

if (n_zero > 0) {
  cat("\n====================================================================\n")
  cat("SENSITIVITY: Offset model with log(TotalEvents + 1), all data\n")
  cat("====================================================================\n\n")

  m1s <- glmer(SubEvents ~ Treatment + (1 | Experiment),
               data = dat,
               family = poisson,
               offset = log(TotalEvents + 1))

  print(Anova(m1s, type = "II"))
}

# ============================================================================
# 7. OVERDISPERSION CHECK
# ============================================================================
# For a Poisson model, the Pearson dispersion ratio should be ~1. Values well
# above 1 indicate overdispersion (more variance than Poisson predicts), which
# inflates Type I error. If detected, we refit with an observation-level random
# effect (OLRE), which is equivalent to a log-normal-Poisson mixture and
# absorbs the extra variance.

cat("\n====================================================================\n")
cat("OVERDISPERSION DIAGNOSTICS\n")
cat("====================================================================\n\n")

resid_pearson <- sum(residuals(m1, type = "pearson")^2)
resid_df <- nrow(dat_nozero) - length(fixef(m1)) - 1
dispersion <- resid_pearson / resid_df
cat(sprintf("Pearson dispersion ratio (Model 1): %.3f\n", dispersion))
cat("Expected ~1 for Poisson. Values >> 1.5 suggest overdispersion.\n\n")

if (dispersion > 1.5) {
  cat("Overdispersion detected — fitting model with observation-level\n")
  cat("random effect (OLRE) to absorb extra-Poisson variance.\n\n")
  dat_nozero$obs <- factor(1:nrow(dat_nozero))
  m1_olre <- glmer(SubEvents ~ Treatment + (1 | Experiment) + (1 | obs),
                   data = dat_nozero,
                   family = poisson,
                   offset = log(TotalEvents))
  print(summary(m1_olre))
  cat("\n--- Type II Wald test (OLRE model) ---\n")
  print(Anova(m1_olre, type = "II"))
} else {
  cat("No substantial overdispersion — Poisson model is adequate.\n")
}

# ============================================================================
# 8. VISUALIZATION
# ============================================================================

dat_nozero$Rate <- dat_nozero$SubEvents / dat_nozero$TotalEvents

p1 <- ggplot(dat_nozero, aes(x = Treatment, y = Rate, color = Treatment)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.7) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.8) +
  labs(y = paste0(orig_sub_label, " / ", orig_tot_label),
       title = "Rate of sub-events by treatment") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(file.path(outdir, "rate_by_treatment.pdf"), p1, width = 5, height = 4)
ggsave(file.path(outdir, "rate_by_treatment.png"), p1, width = 5, height = 4, dpi = 300)

p2 <- ggplot(dat_nozero, aes(x = TotalEvents, y = SubEvents, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "glm", method.args = list(family = "poisson"),
              se = FALSE, linewidth = 1) +
  labs(x = orig_tot_label, y = orig_sub_label,
       title = paste0(orig_sub_label, " vs ", orig_tot_label, " by treatment")) +
  scale_x_continuous(trans = "log1p") +
  scale_y_continuous(trans = "log1p") +
  theme_classic(base_size = 14)

ggsave(file.path(outdir, "subevents_vs_total.pdf"), p2, width = 6, height = 4)
ggsave(file.path(outdir, "subevents_vs_total.png"), p2, width = 6, height = 4, dpi = 300)

cat("\nPlots saved to:", outdir, "\n")

# ============================================================================
# 9. SUMMARY TABLES → Markdown file
# ============================================================================

cat("\n====================================================================\n")
cat("GENERATING SUMMARY REPORT (Markdown)\n")
cat("====================================================================\n\n")

# --- Helper: format p-values ---
fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
  ifelse(p < 0.01,  sprintf("%.3f", p),
                     sprintf("%.3f", p)))
}

# --- Helper: data.frame to markdown table ---
md_table <- function(df) {
  header <- paste0("| ", paste(colnames(df), collapse = " | "), " |")
  sep    <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows   <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# -----------------------------------------------------------------------
# Table 1: Descriptives
# -----------------------------------------------------------------------
desc <- do.call(rbind, lapply(levels(dat$Treatment), function(trt) {
  sub <- dat_nozero[dat_nozero$Treatment == trt, ]
  data.frame(
    Treatment          = trt,
    n                  = nrow(sub),
    `Sub-events`       = sprintf("%.1f +/- %.1f", mean(sub$SubEvents), sd(sub$SubEvents)),
    `Total events`     = sprintf("%.1f +/- %.1f", mean(sub$TotalEvents), sd(sub$TotalEvents)),
    `Rate (sub/total)` = sprintf("%.3f +/- %.3f", mean(sub$Rate), sd(sub$Rate)),
    check.names = FALSE
  )
}))

# -----------------------------------------------------------------------
# Table 2: Omnibus Treatment effect from all models
# -----------------------------------------------------------------------
extract_wald <- function(model, label) {
  a <- Anova(model, type = "II")
  idx <- which(rownames(a) == "Treatment")
  data.frame(
    Model     = label,
    Chisq     = sprintf("%.2f", a$Chisq[idx]),
    df        = a$Df[idx],
    `p-value` = fmt_p(a$`Pr(>Chisq)`[idx]),
    check.names = FALSE
  )
}

wald_rows <- rbind(
  extract_wald(m1, "Model 1: Offset (rate)"),
  extract_wald(m2, "Model 2: Covariate"),
  extract_wald(m3, "Model 3: Interaction")
)

wald_rows$`What it tests` <- c(
  "Rate of sub-events per total event differs by treatment?",
  "Treatment predicts sub-events after adjusting for total events?",
  "Sub-event ~ total-event relationship changes by treatment?"
)

# -----------------------------------------------------------------------
# Table 3: Full fixed-effect coefficients for Models 2 and 3
# -----------------------------------------------------------------------
extract_fixef <- function(model, label) {
  cc  <- summary(model)$coefficients
  a   <- Anova(model, type = "II")
  out <- data.frame(
    Model          = label,
    Term           = rownames(cc),
    Estimate       = sprintf("%.3f", cc[, "Estimate"]),
    SE             = sprintf("%.3f", cc[, "Std. Error"]),
    z              = sprintf("%.2f", cc[, "z value"]),
    `p (coeff)`    = fmt_p(cc[, "Pr(>|z|)"]),
    check.names = FALSE
  )
  out$`p (Wald omnibus)` <- ""
  for (term in rownames(a)) {
    matching <- grepl(term, out$Term, fixed = TRUE)
    if (any(matching)) {
      first_match <- which(matching)[1]
      out$`p (Wald omnibus)`[first_match] <- fmt_p(a[term, "Pr(>Chisq)"])
    }
  }
  out
}

fixef_table <- rbind(
  extract_fixef(m2, "Model 2"),
  extract_fixef(m3, "Model 3")
)

# -----------------------------------------------------------------------
# Table 4: Pairwise comparisons from Model 1
# -----------------------------------------------------------------------
pw <- as.data.frame(pairs(emm1, type = "response"))

pw_clean <- data.frame(
  Contrast     = pw$contrast,
  `Rate ratio` = sprintf("%.3f", pw$ratio),
  SE           = sprintf("%.3f", pw$SE),
  z            = sprintf("%.2f", pw$z.ratio),
  `p-value`    = fmt_p(pw$p.value),
  check.names  = FALSE
)

# -----------------------------------------------------------------------
# Write markdown file
# -----------------------------------------------------------------------
md_path <- file.path(outdir, "GLMM_summary.md")

md <- c(
  "# GLMM Analysis: Sub-event Rate by Treatment",
  "",
  "## Question",
  "",
  "Is the increase in sub-events (e.g., synchronized events) explained by the",
  "global increase in total events, or does the treatment specifically enhance",
  "this event subtype beyond what is expected from higher overall activity?",
  "",
  "## Method",
  "",
  "We fitted three Poisson GLMMs with experiment as a random intercept.",
  "**Model 1** (offset) uses log(total events) as an offset, effectively modeling",
  "the rate of sub-events per total event. A significant Treatment effect in this",
  "model indicates that the increase in sub-events is disproportionate to the",
  "increase in total events.",
  "**Model 2** (covariate) includes total events as a fixed-effect covariate,",
  "testing whether treatment predicts sub-events after statistically adjusting",
  "for total event count.",
  "**Model 3** (interaction) adds a Treatment x TotalEvents interaction to test",
  "whether the relationship between total and sub-events itself differs across",
  "treatments.",
  "",
  "## Results",
  "",
  "### Table 1. Descriptive statistics by treatment (mean +/- SD)",
  "",
  md_table(desc),
  "",
  "### Table 2. Omnibus Treatment effect (Type II Wald chi-square test)",
  "",
  md_table(wald_rows),
  "",
  "### Table 3. Fixed-effect coefficients for Models 2 and 3",
  "",
  md_table(fixef_table),
  "",
  "### Table 4. Pairwise comparisons of sub-event rates (Model 1, Tukey-adjusted)",
  "",
  md_table(pw_clean),
  "",
  "## How to read these tables",
  "",
  "**Table 1** shows raw descriptive statistics.",
  "**Table 2** is the key result: a significant p-value in the Offset model",
  "(Model 1) means the *rate* of sub-events differs across treatments, ruling",
  "out that the increase is a trivial consequence of more total events.",
  "**Table 3** shows the full fixed-effect coefficients for Models 2 and 3,",
  "including the TotalEvents covariate and any interaction terms.",
  "**Table 4** shows1 pairwise rate-ratio comparisons from Model 1.",
  "",
  sprintf("**Overdispersion (Pearson ratio):** %.3f.", dispersion),
  ifelse(dispersion > 1.5,
    "Overdispersion was detected; an OLRE model was fitted (see console output).",
    "No substantial overdispersion; the Poisson model is adequate.")
)

writeLines(md, md_path)

cat(sprintf("Summary report saved: %s\n", md_path))
cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("All outputs saved to: %s\n", outdir))
