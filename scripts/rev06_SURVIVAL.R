#!/usr/bin/env Rscript
# rev06 — Survival re-analysis (TCGA + CGGA-325 + CGGA-693) requested by review M6/DA8.
#
# Reviewer items addressed here:
#   * sample flow and cohort definition (patient-level, grade reported), and clarification
#     that the submitted "TCGA-GBM, n = 939" was a *sample*-level row count of a
#     GBM + LGG table (see rev05_cohort_flow.tsv);
#   * explicit event/time construction (rev05 prep: event = DECEASED, time = months * 30.4375);
#   * hazard ratios per SD of the z-scored MES score, with 95% CIs;
#   * proportional-hazards diagnostics (cox.zph, Schoenfeld);
#   * Kaplan-Meier with a *stated* cut-off rule (median), log-rank P and numbers at risk;
#   * random-effects meta-analysis as the primary pooled estimate, with the cohort range
#     and I^2, instead of foregrounding the fixed-effect estimate at I^2 = 89%;
#   * a purity / stromal-content adjustment (deconvolved malignant fraction + marker-based
#     stromal score) and a head-to-head benchmark against the published Neftel MES module.
#
# Inputs  : results/revision/rev05_{tcga,cgga325,cgga693}_cases.tsv
# Outputs : results/revision/rev06_*.tsv
#           figures/rev06_survival_forest.{pdf,jpeg}, rev06_km_*.{pdf,jpeg}
#
# Run: Rscript scripts/rev06_SURVIVAL.R      (from the project root)

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(metafor)
  library(ggplot2)
})

proj <- Sys.getenv("GLIOMA_PROJ", unset = "/run/media/meng/Seqdata/glioma/Paper_2")
res  <- file.path(proj, "results", "revision")
fig  <- file.path(proj, "figures")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

THEME <- theme_classic(base_size = 11) + theme(panel.grid = element_blank(),
                                               axis.line = element_line(linewidth = 0.4))

cohorts <- list(
  "TCGA"     = list(file = "rev05_tcga_cases.tsv",    time = "os_time_days", age = "age"),
  "CGGA-325" = list(file = "rev05_cgga325_cases.tsv", time = "os_days",      age = "age"),
  "CGGA-693" = list(file = "rev05_cgga693_cases.tsv", time = "os_days",      age = "age")
)

cox_one <- function(d, score, covars, label) {
  d <- d[is.finite(d[[score]]) & is.finite(d$time) & d$time > 0 & is.finite(d$event), ]
  d$z <- as.numeric(scale(d[[score]]))              # HR is *per SD of the z-scored score*
  cov <- c("z", covars)
  keep <- covars[vapply(covars, function(c) sum(is.finite(d[[c]])) > 50, logical(1))]
  cov <- c("z", keep)
  f <- as.formula(paste("Surv(time, event) ~", paste(cov, collapse = " + ")))
  m <- coxph(f, data = d)
  s <- summary(m)
  zph <- tryCatch(cox.zph(m), error = function(e) NULL)
  ph_p <- if (!is.null(zph)) zph$table["z", "p"] else NA_real_
  ci <- s$conf.int
  data.frame(model = label, score = score, n = m$n, events = m$nevent,
             HR_per_SD = ci["z", "exp(coef)"], lo = ci["z", "lower .95"],
             hi = ci["z", "upper .95"], p = s$coefficients["z", "Pr(>|z|)"],
             cox_zph_p = ph_p,
             covars = paste(cov, collapse = ","), row.names = NULL)
}

cohort_stats <- list()
cox_rows <- list()

for (cn in names(cohorts)) {
  spec <- cohorts[[cn]]
  d <- read.delim(file.path(res, spec$file), check.names = TRUE)
  d$time  <- as.numeric(d[[spec$time]])
  d$event <- as.numeric(d$event)
  d$age   <- as.numeric(d[[spec$age]])
  d$idh_mut <- suppressWarnings(as.numeric(d$idh_mut))
  d$grade_num <- suppressWarnings(as.numeric(d$grade_num))
  if (all(!is.finite(d$grade_num)) && "Grade" %in% names(d)) {
    rom <- c("II" = 2, "III" = 3, "IV" = 4)          # CGGA writes "WHO IV" etc.
    d$grade_num <- unname(rom[sub("^WHO ", "", as.character(d$Grade))])
    cat(sprintf("[%s] grade parsed from Roman numerals: %d/%d known\n",
                cn, sum(is.finite(d$grade_num)), nrow(d)))
  }
  d <- d[is.finite(d$time) & d$time > 0 & is.finite(d$event), ]

  # ---- cohort description (M6 item 1) ----
  cohort_stats[[cn]] <- data.frame(
    cohort = cn, n_patients = nrow(d), events = sum(d$event == 1),
    median_followup = median(d$time),
    age_median = median(d$age, na.rm = TRUE),
    grade_known = sum(is.finite(d$grade_num)),
    gbm_like = sum(d$grade_num == 4, na.rm = TRUE),
    idh_known = sum(is.finite(d$idh_mut)),
    mes_range = paste(round(range(d$MES_Neftel, na.rm = TRUE), 2), collapse = " to ")
  )

  # ---- Kaplan-Meier: median cut-off, log-rank, at-risk table ----
  d$grp <- ifelse(d$MES_Neftel > median(d$MES_Neftel, na.rm = TRUE), "High", "Low")
  sf  <- survfit(Surv(time, event) ~ grp, data = d)
  lr  <- survdiff(Surv(time, event) ~ grp, data = d)
  lr_p <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  cohort_stats[[cn]]$logrank_p <- lr_p

  km <- ggsurvplot(sf, data = d, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                   palette = c("#C0392B", "#2C3E50"), legend.title = "MES score",
                   xlab = "Time (days)", ylab = "Overall survival probability",
                   title = paste0(cn, ": MES (Neftel module) high vs low (median split)"))
  km$plot <- km$plot + THEME
  for (ext in c("pdf", "jpeg")) {
    if (ext == "pdf") {
      pdf(file.path(fig, paste0("rev06_km_", gsub("-", "", cn), ".pdf")), width = 7, height = 7)
    } else {
      jpeg(file.path(fig, paste0("rev06_km_", gsub("-", "", cn), ".jpeg")),
           width = 7, height = 7, units = "in", res = 300)
    }
    print(km)
    dev.off()
  }

  # ---- Cox models (primary score + benchmark + adjustments) ----
  cox_rows[[paste0(cn, "_primary")]]   <- cox_one(d, "MES_Neftel", c("age", "idh_mut"), paste0(cn, ": MES_Neftel + age + IDH"))
  cox_rows[[paste0(cn, "_grade")]]     <- cox_one(d, "MES_Neftel", c("age", "idh_mut", "grade_num"), paste0(cn, ": + WHO grade"))
  cox_rows[[paste0(cn, "_prog131")]]   <- cox_one(d, "MES_prog131", c("age", "idh_mut"), paste0(cn, ": MES_prog131 + age + IDH"))
  cox_rows[[paste0(cn, "_purity")]]    <- cox_one(d, "MES_Neftel", c("age", "idh_mut", "malignant_frac"), paste0(cn, ": + deconvolved malignant fraction"))
  cox_rows[[paste0(cn, "_stroma")]]    <- cox_one(d, "MES_Neftel", c("age", "idh_mut", "stromal_content", "myeloid_content"), paste0(cn, ": + stromal/myeloid content"))
  cox_rows[[paste0(cn, "_grade_only")]] <- {
    dd <- d[is.finite(d$grade_num), ]
    if (nrow(dd) > 60) cox_one(dd, "MES_Neftel", c("age", "grade_num"), paste0(cn, ": grade-stratified set")) else NULL
  }
}

cox_tab <- do.call(rbind, cox_rows)
write.table(cox_tab, file.path(res, "rev06_cox_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
desc <- do.call(rbind, cohort_stats)
write.table(desc, file.path(res, "rev06_cohort_description.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- meta-analysis: random-effects primary (DL), fixed-effects secondary ----
prim <- cox_tab[grepl("MES_Neftel \\+ age \\+ IDH$", cox_tab$model), ]
prim <- prim[!is.na(prim$HR_per_SD), ]
es <- escalc(measure = "RR", yi = log(prim$HR_per_SD),
             sei = (log(prim$hi) - log(prim$lo)) / (2 * qnorm(0.975)),
             slab = prim$model, data = prim)
re <- rma(yi, vi, data = es, method = "DL")
fe <- rma(yi, vi, data = es, method = "FE")
meta_tab <- data.frame(
  model = c("random-effects (primary)", "fixed-effects (secondary)"),
  HR = exp(c(as.numeric(re$b), as.numeric(fe$b))),
  lo = exp(c(re$ci.lb, fe$ci.lb)), hi = exp(c(re$ci.ub, fe$ci.ub)),
  p = c(re$pval, fe$pval), I2 = c(re$I2, fe$I2), tau2 = c(re$tau2, fe$tau2),
  cohort_HR_range = paste(round(range(prim$HR_per_SD), 2), collapse = "-")
)
write.table(meta_tab, file.path(res, "rev06_meta_analysis.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(meta_tab)

# ---- forest plot of all models (figure 5f) ----
ftab <- cox_tab[!is.na(cox_tab$HR_per_SD), ]
ftab$cohort <- sub(":.*", "", ftab$model)
ftab$grp <- ifelse(grepl("^CGGA", ftab$cohort), "CGGA", "TCGA")
ftab$model <- factor(ftab$model, levels = rev(ftab$model))
p <- ggplot(ftab, aes(x = HR_per_SD, y = model, colour = grp)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_colour_manual(values = c("CGGA" = "#C0392B", "TCGA" = "#2C3E50"), name = "Cohort") +
  labs(x = "Hazard ratio per SD of MES score (95% CI)", y = NULL,
       title = "Prognostic value of the MES program across cohorts and models") +
  THEME
ggsave(file.path(fig, "rev06_survival_forest.pdf"), p, width = 9, height = 5.5, device = cairo_pdf)
ggsave(file.path(fig, "rev06_survival_forest.jpeg"), p, width = 9, height = 5.5, dpi = 300)

# ---- grade trend (figure S2d) ----
gt <- read.delim(file.path(res, "rev05_tcga_cases.tsv"))
gt <- gt[is.finite(gt$grade_num) & is.finite(gt$MES_Neftel), ]
kt <- kruskal.test(MES_Neftel ~ grade_num, data = gt)
sp <- cor.test(gt$grade_num, gt$MES_Neftel, method = "spearman")
write.table(data.frame(test = c("Kruskal-Wallis MES by WHO grade (TCGA)", "Spearman grade vs MES (TCGA)"),
                       statistic = c(unname(kt$statistic), unname(sp$estimate)),
                       p = c(kt$p.value, sp$p.value)),
            file.path(res, "rev06_grade_trend.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
p2 <- ggplot(gt, aes(x = factor(grade_num), y = MES_Neftel)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey85") +
  geom_jitter(width = 0.14, size = 0.8, alpha = 0.5, colour = "#C0392B") +
  labs(x = "WHO grade", y = "MES score (Neftel module)",
       subtitle = sprintf("Kruskal-Wallis P = %.3g; Spearman rho = %.2f", kt$p.value, sp$estimate)) +
  THEME
ggsave(file.path(fig, "rev06_grade_trend.pdf"), p2, width = 5, height = 4, device = cairo_pdf)
ggsave(file.path(fig, "rev06_grade_trend.jpeg"), p2, width = 5, height = 4, dpi = 300)

writeLines("DONE rev06")
