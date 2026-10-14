#!/usr/bin/env Rscript
# rev09 — NicheNet ligand activity against the MES-transition programme (review M3 minor 6).
#
# Reviewer items addressed: the submitted text reported a "top ligand" list that did not
# match the plotted panel (LEFTY2, PSPN, LEAP2, SLITRK1, SPON1, VIP were shown but never
# discussed), the scoring formula was never given, and CSF1 was called the top ligand while
# its receptor was reported as "nearly absent" on malignant cells.
#
# Method (stated explicitly, base R on the NicheNet prior — no package install):
#   ligand activity = Pearson correlation between the ligand's regulatory-potential vector
#   (column of the ligand-target prior matrix) and a binary indicator of the MES-transition
#   programme, evaluated over the genes shared between the prior and the programme.
#   This is the NicheNet "ligand activity" definition, computed from the prior directly.
#
# Outputs: results/revision/rev09_nichenet_ligand_activity.tsv  (all ligands, ranked)
#          figures/rev09_nichenet_ligands.{pdf,jpeg}

suppressPackageStartupMessages({ library(ggplot2) })
proj <- Sys.getenv("GLIOMA_PROJ", unset = "/run/media/meng/Seqdata/glioma/Paper_2")
dat  <- Sys.getenv("GLIOMA_DATA", unset = "/run/media/meng/Seqdata/glioma/Datasets")
res  <- file.path(proj, "results", "revision"); fig <- file.path(proj, "figures")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

ltm <- readRDS(file.path(dat, "nichenet", "ligand_target_matrix_nsga2r_final.rds"))
lr  <- readRDS(file.path(dat, "nichenet", "lr_network_human_21122021.rds"))
cat(sprintf("ligand-target prior: %d genes x %d ligands\n", nrow(ltm), ncol(ltm)))

prog <- read.delim(file.path(proj, "results", "step04_traj", "mes_transition_program_core.tsv"))
prog_genes <- as.character(prog$gene[prog$rho > 0])
genes <- intersect(prog_genes, rownames(ltm))
cat(sprintf("MES-transition programme: %d genes, of which %d in the prior\n",
            length(prog_genes), length(genes)))

y <- as.numeric(rownames(ltm) %in% genes)
act <- apply(ltm[rownames(ltm) %in% rownames(ltm), , drop = FALSE], 2, function(x)
  suppressWarnings(cor(x, y, method = "pearson")))
act <- act[!is.na(act)]
d <- data.frame(ligand = names(act), ligand_activity = as.numeric(act))
d <- d[order(-d$ligand_activity), ]
write.table(d, file.path(res, "rev09_nichenet_ligand_activity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

top <- head(d, 20)
top$ligand <- factor(top$ligand, levels = rev(top$ligand))
p <- ggplot(top, aes(x = ligand_activity, y = ligand)) +
  geom_col(fill = "#2980B9") +
  labs(x = "NicheNet ligand activity vs the 128-gene MES-transition programme",
       y = NULL, title = "Ligands best predicting the MES programme (base-R NicheNet prior)") +
  theme_classic(base_size = 10) + theme(panel.grid = element_blank())
ggsave(file.path(fig, "rev09_nichenet_ligands.pdf"), p, width = 6.5, height = 5, device = cairo_pdf)
ggsave(file.path(fig, "rev09_nichenet_ligands.jpeg"), p, width = 6.5, height = 5, dpi = 300)
cat("top 15 ligands:\n"); print(head(d, 15))
writeLines("DONE rev09")
