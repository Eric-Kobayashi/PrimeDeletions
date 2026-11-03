# Figure 2j: original DP9 cell 38 summary_plot beeswarm, called through R.
# Same values/algorithm; only labels, export size and layout seed are changed.
plot_figure2j <- function(context, output_file) {
  shap <- reticulate::import("shap", convert=FALSE)
  plt <- reticulate::import("matplotlib.pyplot", convert=FALSE)
  np <- reticulate::import("numpy", convert=FALSE)
  np$random$seed(42L)
  plt$figure()
  shap$summary_plot(context$shap_values, context$features, context$labels,
                    show=FALSE, plot_size=reticulate::tuple(5.2, 4.0))
  plt$savefig(output_file, dpi=300L, transparent=TRUE, bbox_inches="tight")
  plt$close()
  invisible(output_file)
}
