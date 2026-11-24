## ============================================================
## TCGA-UCEC APOBEC3B (A3B) Analysis -- SUBTYPE-STRATIFIED VERSION
##
## Produces, for each of the 4 TCGA molecular subtypes
## (POLE, MSI, CN-low, CN-high) SEPARATELY, plus overall:
##   0) Summary table: N, A3B high/low counts, events, median expr
##   (A) Boxplot of A3B expression by subtype (overall, unchanged)
##   (B1) KM curves, OS,  faceted by subtype (A3B high vs low)
##   (B2) KM curves, PFS, faceted by subtype (A3B high vs low)
##   (B3) Per-subtype log-rank p-values (table)
##   (C1) Multivariable Cox model: A3B adjusted for subtype
##   (C2) Cox model with A3B x SUBTYPE interaction term
##        (tests whether the A3B effect differs BY subtype)
##   (C3) Forest plot: subtype-specific HR for A3B (high vs low),
##        fit as separate Cox models within each subtype
##   (D) TMB vs A3B correlation, faceted + tested within each subtype
## ============================================================

## ---- 0. Packages ----
required_pkgs <- c("dplyr", "ggplot2", "survival", "survminer", "ggpubr", "broom")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

## ---- 1. Load + merge data ----
clinical <- read.csv("cbio_clinical_data.csv",   stringsAsFactors = FALSE)
expr     <- read.csv("cbio_clinical_data-2.csv", stringsAsFactors = FALSE)

# rename by position (raw header has spaces/colons that read.csv mangles)
names(expr)[5] <- "A3B_group"
names(expr)[6] <- "A3B_expr"
expr$A3B_expr <- as.numeric(expr$A3B_expr)

merged <- clinical %>%
  inner_join(expr, by = c("sampleId" = "sample")) %>%
  filter(SUBTYPE != "", !is.na(A3B_expr))

merged$SUBTYPE <- factor(
  merged$SUBTYPE,
  levels = c("UCEC_POLE", "UCEC_MSI", "UCEC_CN_LOW", "UCEC_CN_HIGH"),
  labels = c("POLE", "MSI", "CN-low", "CN-high")
)
merged$A3B_group <- factor(merged$A3B_group, levels = c("low", "high"))

merged$OS_time   <- merged$OS_MONTHS
merged$OS_event  <- ifelse(grepl("^1", merged$OS_STATUS), 1, 0)
merged$PFS_time  <- merged$PFS_MONTHS
merged$PFS_event <- ifelse(grepl("^1", merged$PFS_STATUS), 1, 0)
merged$TMB       <- as.numeric(merged$TMB_NONSYNONYMOUS)

cat("Merged dataset:", nrow(merged), "samples with clinical + A3B expression data\n\n")

## ============================================================
## 0. SUBTYPE-SPECIFIC SUMMARY TABLE
## ============================================================
summary_tbl <- merged %>%
  group_by(SUBTYPE) %>%
  summarise(
    N              = n(),
    N_A3B_high     = sum(A3B_group == "high"),
    N_A3B_low      = sum(A3B_group == "low"),
    Median_A3B     = round(median(A3B_expr, na.rm = TRUE), 2),
    OS_events      = sum(OS_event, na.rm = TRUE),
    PFS_events     = sum(PFS_event, na.rm = TRUE),
    Median_TMB     = round(median(TMB, na.rm = TRUE), 2),
    .groups = "drop"
  )
print(summary_tbl)
write.csv(summary_tbl, "subtype_summary_table.csv", row.names = FALSE)

## ============================================================
## (A) A3B expression by subtype -- overall boxplot (unchanged)
## ============================================================
pA <- ggboxplot(
  merged, x = "SUBTYPE", y = "A3B_expr",
  fill = "SUBTYPE", palette = "npg", add = "jitter",
  xlab = "TCGA Molecular Subtype", ylab = "APOBEC3B expression (log2 RSEM)"
) +
  stat_compare_means(method = "kruskal.test") +
  theme(legend.position = "none")
ggsave("A3B_by_subtype_boxplot.png", pA, width = 6, height = 5, dpi = 300)

## Pairwise subtype comparisons (Wilcoxon, BH-corrected) -- useful for
## the text/supplement even though the omnibus test is the headline
pairwise_res <- pairwise.wilcox.test(merged$A3B_expr, merged$SUBTYPE, p.adjust.method = "BH")
print(pairwise_res)
capture.output(pairwise_res, file = "A3B_subtype_pairwise_wilcox.txt")

## ============================================================
## (B1)/(B2) KM curves faceted by subtype: OS and PFS
## ============================================================
## NOTE: ggsurvplot_facet() re-evaluates the survfit call internally, so
## the formula must reference real column names directly (not variables
## built via paste0/as.formula inside a helper function -- that breaks
## survminer's internal call reconstruction). OS and PFS are written out
## explicitly below for that reason.

fit_os_facet <- survfit(Surv(OS_time, OS_event) ~ A3B_group, data = merged)
p_os_facet <- ggsurvplot_facet(
  fit_os_facet, data = merged, facet.by = "SUBTYPE",
  pval = TRUE, legend.title = "A3B expression",
  xlab = "Time (months)", ylab = "Overall Survival probability",
  palette = c("#2E9FDF", "#E7B800"), short.panel.labs = TRUE
)
ggsave("A3B_OS_KM_by_subtype.png", p_os_facet, width = 9, height = 7, dpi = 300)

fit_pfs_facet <- survfit(Surv(PFS_time, PFS_event) ~ A3B_group, data = merged)
p_pfs_facet <- ggsurvplot_facet(
  fit_pfs_facet, data = merged, facet.by = "SUBTYPE",
  pval = TRUE, legend.title = "A3B expression",
  xlab = "Time (months)", ylab = "Progression-Free Survival probability",
  palette = c("#2E9FDF", "#E7B800"), short.panel.labs = TRUE
)
ggsave("A3B_PFS_KM_by_subtype.png", p_pfs_facet, width = 9, height = 7, dpi = 300)

## (B3) Per-subtype log-rank p-values, collected into one table
logrank_by_subtype <- function(time_var, event_var) {
  merged %>%
    group_by(SUBTYPE) %>%
    group_modify(~{
      d <- .x
      if (length(unique(d$A3B_group)) < 2 || sum(d[[event_var]], na.rm = TRUE) < 2) {
        return(tibble(N = nrow(d), events = sum(d[[event_var]], na.rm = TRUE), p_logrank = NA))
      }
      sd <- survdiff(as.formula(paste0("Surv(", time_var, ",", event_var, ") ~ A3B_group")), data = d)
      p  <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
      tibble(N = nrow(d), events = sum(d[[event_var]], na.rm = TRUE), p_logrank = round(p, 3))
    }) %>% ungroup()
}
os_by_subtype  <- logrank_by_subtype("OS_time",  "OS_event")  %>% mutate(endpoint = "OS")
pfs_by_subtype <- logrank_by_subtype("PFS_time", "PFS_event") %>% mutate(endpoint = "PFS")
logrank_table  <- bind_rows(os_by_subtype, pfs_by_subtype)
print(logrank_table)
write.csv(logrank_table, "logrank_pvalues_by_subtype.csv", row.names = FALSE)

## ============================================================
## (C1) Multivariable Cox model: A3B adjusted for SUBTYPE
## ============================================================
cox_os_adj  <- coxph(Surv(OS_time, OS_event)   ~ A3B_expr + SUBTYPE, data = merged)
cox_pfs_adj <- coxph(Surv(PFS_time, PFS_event) ~ A3B_expr + SUBTYPE, data = merged)

cat("\n--- Cox model: OS ~ A3B_expr + SUBTYPE ---\n")
print(summary(cox_os_adj))
cat("\n--- Cox model: PFS ~ A3B_expr + SUBTYPE ---\n")
print(summary(cox_pfs_adj))

write.csv(tidy(cox_os_adj,  exponentiate = TRUE, conf.int = TRUE), "cox_OS_adjusted_for_subtype.csv",  row.names = FALSE)
write.csv(tidy(cox_pfs_adj, exponentiate = TRUE, conf.int = TRUE), "cox_PFS_adjusted_for_subtype.csv", row.names = FALSE)

## ============================================================
## (C2) Interaction test: does the A3B effect differ BY subtype?
## ============================================================
cox_os_int  <- coxph(Surv(OS_time, OS_event)   ~ A3B_expr * SUBTYPE, data = merged)
cox_pfs_int <- coxph(Surv(PFS_time, PFS_event) ~ A3B_expr * SUBTYPE, data = merged)

os_lrt  <- anova(cox_os_adj,  cox_os_int)
pfs_lrt <- anova(cox_pfs_adj, cox_pfs_int)
cat("\n--- LRT for A3B x SUBTYPE interaction (OS) ---\n");  print(os_lrt)
cat("\n--- LRT for A3B x SUBTYPE interaction (PFS) ---\n"); print(pfs_lrt)

capture.output(os_lrt, pfs_lrt, file = "interaction_test_A3B_x_subtype.txt")

## ============================================================
## (C3) Forest plot: subtype-specific HR (A3B high vs low)
## ============================================================
subtype_hr <- function(time_var, event_var) {
  merged %>%
    group_by(SUBTYPE) %>%
    group_modify(~{
      d <- .x
      if (length(unique(d$A3B_group)) < 2 || sum(d[[event_var]], na.rm = TRUE) < 3) {
        return(tibble(HR = NA, lower = NA, upper = NA, p = NA, N = nrow(d), events = sum(d[[event_var]])))
      }
      fit <- coxph(as.formula(paste0("Surv(", time_var, ",", event_var, ") ~ A3B_group")), data = d)
      s <- summary(fit)
      tibble(
        HR = s$coefficients[1, "exp(coef)"],
        lower = s$conf.int[1, "lower .95"],
        upper = s$conf.int[1, "upper .95"],
        p = s$coefficients[1, "Pr(>|z|)"],
        N = nrow(d), events = sum(d[[event_var]])
      )
    }) %>% ungroup()
}

os_hr  <- subtype_hr("OS_time",  "OS_event")  %>% mutate(endpoint = "OS")
pfs_hr <- subtype_hr("PFS_time", "PFS_event") %>% mutate(endpoint = "PFS")
forest_data <- bind_rows(os_hr, pfs_hr)
print(forest_data)
write.csv(forest_data, "subtype_specific_HR_A3B_high_vs_low.csv", row.names = FALSE)

## Some subtype/endpoint strata have too few events for a stable HR
## estimate (e.g., zero events in one arm -> "perfect separation",
## HR -> near-infinite). Flag and exclude those from the PLOT (they
## remain in the CSV table with their p-values for transparency).
forest_data <- forest_data %>%
  mutate(unstable = is.na(HR) | HR > 1000 | (upper / pmax(lower, 1e-6)) > 1000)
if (any(forest_data$unstable, na.rm = TRUE)) {
  cat("\nExcluded from forest PLOT due to separation/instability (see CSV for full values):\n")
  print(forest_data %>% filter(unstable) %>% select(SUBTYPE, endpoint, HR, lower, upper, p, N, events))
}
forest_plot_data <- forest_data %>% filter(!unstable)

pForest <- ggplot(forest_plot_data, aes(x = HR, y = SUBTYPE, color = endpoint)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2,
                  position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (A3B high vs low, 95% CI)", y = "TCGA Molecular Subtype",
       title = "Subtype-specific effect of A3B expression on survival") +
  theme_minimal(base_size = 13)
ggsave("A3B_forest_plot_by_subtype.png", pForest, width = 7, height = 5, dpi = 300)

## ============================================================
## (D) TMB vs A3B correlation, faceted by subtype + per-subtype stats
## ============================================================
pD <- ggscatter(
  merged, x = "A3B_expr", y = "TMB",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "spearman",
  xlab = "APOBEC3B expression (log2 RSEM)", ylab = "TMB (mutations/Mb)"
) + facet_wrap(~SUBTYPE, scales = "free_y")
ggsave("A3B_vs_TMB_by_subtype.png", pD, width = 8, height = 7, dpi = 300)

tmb_cor_by_subtype <- merged %>%
  group_by(SUBTYPE) %>%
  group_modify(~{
    ct <- cor.test(.x$A3B_expr, .x$TMB, method = "spearman", exact = FALSE)
    tibble(N = nrow(.x), rho = round(unname(ct$estimate), 3), p = round(ct$p.value, 3))
  }) %>% ungroup()
print(tmb_cor_by_subtype)
write.csv(tmb_cor_by_subtype, "TMB_correlation_by_subtype.csv", row.names = FALSE)

## ---- NOTE ----
## CD8+ T-cell infiltration was not present in either uploaded file.
## Once available (keyed by sample ID), merge it in the same way as
## `expr` above and repeat sections (D) with y = "<CD8 column>" to add
## a subtype-stratified immune-infiltration panel.

cat("\nDone. All figures + tables saved to working directory.\n")
