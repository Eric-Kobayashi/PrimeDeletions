# Figure 2i-k: original upstream R scatter expression and original SHAP routines.
# Called by run.R with `output` and `context` set. No per-panel function wrappers.
# R source: ../Anchored screen/Anchored screen plots.R at e1a1a58711bf8237b0ed565d3155c8f70194fccc.


#### Figure 2i
library(ggplot2)
library(grid)
test_data_corefeatures <- read.csv(file.path(output, "test_predictions.csv"))
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
# Keep the original 4-cm panel; size the PDF canvas to also include its labels.
# The upstream 4-cm canvas clips labels around its 4-cm panel.
pdf(NULL)
panel <- egg::set_panel_size(p, width=unit(4, "cm"), height=unit(4, "cm"))
panel <- gtable::gtable_add_padding(panel, unit(rep(2, 4), "mm"))
width <- convertWidth(sum(panel$widths), "in", valueOnly=TRUE)
height <- convertHeight(sum(panel$heights), "in", valueOnly=TRUE)
dev.off()
pdf(file=file.path(output, "figure2i.pdf"), width=width, height=height)
grid.draw(panel)
dev.off()


#### Figure 2j
shap <- reticulate::import("shap", convert=FALSE)
plt <- reticulate::import("matplotlib.pyplot", convert=FALSE)
np <- reticulate::import("numpy", convert=FALSE)
np$random$seed(42L)
plt$figure()
shap$summary_plot(context$shap_values, context$features, context$labels,
                  show=FALSE, plot_size=reticulate::tuple(5.2, 4.0))
plt$savefig(file.path(output, "figure2j.pdf"), dpi=300L, transparent=TRUE, bbox_inches="tight")
plt$close()


#### Figure 2k
shap <- reticulate::import("shap", convert=FALSE)
plt <- reticulate::import("matplotlib.pyplot", convert=FALSE)
shap$summary_plot(context$shap_values, features=context$features,
                  feature_names=context$labels, plot_type="bar", show=FALSE,
                  color="#D5CCAE", plot_size=reticulate::tuple(4.5, 4.0))
plt$xlabel("Mean absolute SHAP value")
plt$setp(plt$gca()$patches, edgecolor="black", linewidth=0.5)
plt$savefig(file.path(output, "figure2k.pdf"), dpi=300L, transparent=TRUE, bbox_inches="tight")
plt$close()
