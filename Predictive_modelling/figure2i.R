# Figure 2i: original upstream point-plot block, with local I/O and function scope.
# Source: Anchored screen/Anchored screen plots.R at e1a1a58711bf8237b0ed565d3155c8f70194fccc.
plot_figure2i <- function(predictions_file, output_file) {
  library(ggplot2)
  library(grid)
  test_data_corefeatures <- read.csv(predictions_file)
  theme_pie <-   theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(strip.background = element_blank(),
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black"),
          text = element_text(family = "Helvetica"),
          legend.text=element_text(size=7),
          legend.key.size = unit(0.3, "cm"),
          plot.background = element_blank())

  p <- ggplot(test_data_corefeatures, aes(x = L2FC, y = L2FC_pred)) +
    #geom_smooth(method = "lm", se = FALSE, color = "#D5CCAE", size = 0.2) +
    geom_abline(slope = 1, intercept = 0, color = "#D5CCAE", linetype = "dashed") +
    #geom_bin2d(bins = 75, show.legend = FALSE) +
    geom_point(shape = 16, alpha = 0.4, size = 0.8, color = "black") +
    scale_fill_continuous(type = "viridis") +
    coord_cartesian(xlim = c(-5,7.5), ylim =  c(-5, 7.5), expand = FALSE) +
    labs(x = "Observed L2FC", y = "Predicted L2FC", fill = "Density") +
    annotate("text", x = 6, y = -4.75, label = paste("R = ", round(cor(test_data_corefeatures$L2FC_pred, test_data_corefeatures$L2FC, method = "pearson", use = "complete.obs"), 2)), family = "Helvetica", size = 2) +
    theme_pie
  dir.create(dirname(output_file), recursive=TRUE, showWarnings=FALSE)
  # Keep the original 4-cm panel; size the PDF canvas to also include its labels.
  # The upstream 4-cm canvas clips labels around its 4-cm panel.
  pdf(NULL)
  panel <- egg::set_panel_size(p, width=unit(4, "cm"), height=unit(4, "cm"))
  panel <- gtable::gtable_add_padding(panel, unit(rep(2, 4), "mm"))
  width <- convertWidth(sum(panel$widths), "in", valueOnly=TRUE)
  height <- convertHeight(sum(panel$heights), "in", valueOnly=TRUE)
  dev.off()
  pdf(file=output_file, width=width, height=height)
  on.exit(dev.off(), add=TRUE)
  grid.draw(panel)
  invisible(p)
}
