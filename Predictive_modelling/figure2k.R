# Figure 2k: original DP9 cell 39 mean-absolute-SHAP bar, called through R.
# Same SHAP summary_plot; short xlabel fixes original label clipping.
plot_figure2k <- function(context, output_file) {
  shap <- reticulate::import("shap", convert=FALSE)
  plt <- reticulate::import("matplotlib.pyplot", convert=FALSE)
  shap$summary_plot(context$shap_values, features=context$features,
                    feature_names=context$labels, plot_type="bar", show=FALSE,
                    color="#D5CCAE", plot_size=reticulate::tuple(4.5, 4.0))
  plt$xlabel("Mean absolute SHAP value")
  plt$setp(plt$gca()$patches, edgecolor="black", linewidth=0.5)
  plt$savefig(output_file, dpi=300L, transparent=TRUE, bbox_inches="tight")
  plt$close()
  invisible(output_file)
}
