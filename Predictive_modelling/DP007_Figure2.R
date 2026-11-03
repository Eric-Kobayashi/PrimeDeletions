# Core-feature subset of DP007_Figure2.R; manuscript point/fit styling follows
# ../Anchored screen/Anchored screen plots.R. Local prepared data replace old inputs.
# Called by run.R; unrelated historical exploratory plots are omitted.
library(ggplot2)
library(grid)

dipper <- rbind(read.csv(file.path(parameter_data_dir, "dipper_train.csv")),
                read.csv(file.path(parameter_data_dir, "dipper_test.csv")))
stopifnot(!anyDuplicated(dipper$unique))
dipper <- dipper[dipper$targetgroup != "control", ]
features <- c("pridict2_score_mS_HEK", "percGC_pbs", "DeepCas9",
              "length_mS_DPr04", "MCC_mS_coverage_bin1000")
labels <- c("PRIDICT2 score", "%GC (PBS)", "DeepCas9 score",
            "Deletion length [bp]", "Contact frequency")
groups <- list(lt10kb=dipper[dipper$length_mS_DPr04 < 10000, ],
               all=dipper,
               gt100kb=dipper[dipper$length_mS_DPr04 > 100000, ])
group_labels <- c(lt10kb="deletions <10 kb", all="all deletions", gt100kb="deletions >100 kb")
theme_pie <- theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        plot.background = element_blank())

outdir <- file.path(output, "core_features")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
plots <- list()
statistics <- list()
pdf(NULL) # measure the original 4-cm panels without an implicit Rplots.pdf
for (feature in features) {
  for (group in names(groups)) {
    dset <- groups[[group]]
    dset <- dset[is.finite(dset$L2FC) & is.finite(dset[[feature]]), ]
    correlation <- cor(dset$L2FC, dset[[feature]], method="pearson")
    log_y <- feature == "length_mS_DPr04" ||
      (feature == "MCC_mS_coverage_bin1000" && group != "gt100kb")
    n_complete <- nrow(dset)
    if (log_y) dset <- dset[dset[[feature]] > 0, ]
    # Original scatter/linear-fit grammar: measured LFC on x, parameter on y.
    p <- ggplot(dset, aes(x = L2FC, y = .data[[feature]])) +
      geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
      geom_point(shape = 16, alpha = 0.5, size = 0.6, color = "black") +
      labs(x = "LFC", y = labels[match(feature, features)], title = group_labels[group]) +
      annotate("text", x = Inf, y = if (log_y) min(dset[[feature]]) else -Inf,
               label = sprintf("R = %.2f", correlation), hjust = 1.05, vjust = -0.3,
               family = "Helvetica", size = 2) +
      theme_pie
    if (log_y) p <- p + scale_y_log10()
    # Fits use the same displayed scale as the original log-MCC panel.
    # The >100-kb MCC panel is linear and includes zero contact values.
    plots[[paste(feature, group, sep="_")]] <- p
    statistics[[length(statistics)+1]] <- data.frame(
      feature=feature, group=group, n_complete=n_complete, n_drawn=nrow(dset),
      log_y=log_y, pearson_r=correlation)
    panel <- egg::set_panel_size(p, width=unit(4,"cm"), height=unit(4,"cm"))
    panel <- gtable::gtable_add_padding(panel, unit(rep(2,4),"mm"))
    width <- convertWidth(sum(panel$widths), "in", valueOnly=TRUE)
    height <- convertHeight(sum(panel$heights), "in", valueOnly=TRUE)
    pdf(file.path(outdir,paste0(feature,"_",group,".pdf")),width=width,height=height,
        family="Helvetica",useDingbats=FALSE)
    grid.newpage(); grid.draw(panel); dev.off()
  }
}
# Five feature rows by three deletion-group columns.
overview <- egg::ggarrange(plots=plots, ncol=3, draw=FALSE)
ggsave(file.path(outdir,"core_features_vs_LFC.pdf"),overview,
       width=19,height=27,units="cm",device="pdf",useDingbats=FALSE)
ggsave(file.path(outdir,"core_features_vs_LFC.png"),overview,
       width=19,height=27,units="cm",dpi=200,bg="white")
write.csv(do.call(rbind,statistics),file.path(outdir,"plot_statistics.csv"),row.names=FALSE)
dev.off()
