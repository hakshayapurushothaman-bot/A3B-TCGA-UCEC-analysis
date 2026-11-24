## Add-on to A3B_UCEC_analysis_subtype.R: run AFTER it (uses the `merged` object it creates).
## Not run by me (no R in my workspace) -- run and check for errors.
library(dplyr); library(survival); library(broom); library(rstatix)

## 1. Kruskal-Wallis: median (IQR) per subtype + Dunn post hoc, Holm-adjusted
merged %>% group_by(SUBTYPE) %>%
  summarise(N = n(), median = median(A3B_expr), Q1 = quantile(A3B_expr, .25), Q3 = quantile(A3B_expr, .75)) %>%
  print() %>% write.csv("A3B_median_IQR_by_subtype.csv", row.names = FALSE)
print(kruskal_test(A3B_expr ~ SUBTYPE, data = merged))
dunn <- dunn_test(A3B_expr ~ SUBTYPE, data = merged, p.adjust.method = "holm")
print(dunn); write.csv(dunn, "A3B_dunn_holm.csv", row.names = FALSE)
cat("Duplicate patients:", sum(duplicated(merged$patientId)), "\n")

## 2. Continuous A3B (per 1 SD) within each subtype, n and events
merged$A3B_z <- as.numeric(scale(merged$A3B_expr))
cont <- function(t, e) merged %>% filter(SUBTYPE != "POLE") %>% group_by(SUBTYPE) %>%
  group_modify(~{ f <- coxph(as.formula(paste0("Surv(", t, ",", e, ") ~ A3B_z")), data = .x)
    tidy(f, exponentiate = TRUE, conf.int = TRUE) %>% mutate(N = nrow(.x), events = sum(.x[[e]])) })
cont_res <- bind_rows(cont("OS_time","OS_event") %>% mutate(endpoint="OS"),
                      cont("PFS_time","PFS_event") %>% mutate(endpoint="PFS"))
print(cont_res); write.csv(cont_res, "A3B_continuous_HR_by_subtype.csv", row.names = FALSE)

## 3. Interaction (LRT), continuous A3B x subtype; POLE excluded for OS (0 OS events -> model cannot be fitted)
inter <- function(t, e, d) { a <- coxph(as.formula(paste0("Surv(",t,",",e,") ~ A3B_z + SUBTYPE")), data = d)
                             b <- coxph(as.formula(paste0("Surv(",t,",",e,") ~ A3B_z * SUBTYPE")), data = d); anova(a, b) }
d_os <- droplevels(filter(merged, SUBTYPE != "POLE"))
print(inter("OS_time","OS_event", d_os)); print(inter("PFS_time","PFS_event", merged))
merged$MSI <- as.integer(merged$SUBTYPE == "MSI")
print(anova(coxph(Surv(OS_time,OS_event) ~ A3B_z + MSI, merged),  coxph(Surv(OS_time,OS_event) ~ A3B_z * MSI, merged)))
print(anova(coxph(Surv(PFS_time,PFS_event) ~ A3B_z + MSI, merged), coxph(Surv(PFS_time,PFS_event) ~ A3B_z * MSI, merged)))

## 4. Pooled Cox adjusted for age, grade, histology, radiotherapy (complete cases)
merged$G3  <- as.integer(merged$GRADE %in% c("G3","High Grade"))
merged$SER <- as.integer(grepl("Serous", merged$CANCER_TYPE_DETAILED))
merged$RT  <- as.integer(merged$RADIATION_THERAPY == "Yes")
for (ep in c("OS","PFS")) {
  f <- coxph(as.formula(paste0("Surv(",ep,"_time,",ep,"_event) ~ A3B_group + AGE + G3 + SER + RT")), data = merged)
  cat(ep, "n =", f$n, "\n"); print(tidy(f, exponentiate = TRUE, conf.int = TRUE))
}
print(sessionInfo())   # gives R and package versions for the Methods
