# Predictive modelling and Figure 2i-k

Reproduce the five-feature core model and Figure 2i-k using the original analysis routines. The previous per-panel R wrappers have been replaced by one sequential plotting script, following the upstream script style.

| File | Purpose |
| --- | --- |
| `model.py` | Original core XGBoost constructor, saved-weight replay and exact training-row TreeSHAP export. |
| `Figure2_model_plots.R` | Combined Figure 2i-k plotting blocks; the upstream R scatter expression is unchanged and the original SHAP plotting algorithms are retained through reticulate. |
| `DP007_Figure2.R` | Byte-for-byte copy of the supplied original 696-line figure script; see the historical-input requirements below. |
| `run.R` | Local data/output adapter; loads committed weights by default and sources the combined plotting script. |
| `weights/core_model.json` | Verified five-feature fitted weights (approximately 1.2 MB), with original metadata alongside. |

Scripts, verified weights, documentation and locked dependencies are tracked. Prepared data, environments and generated outputs remain excluded.

## Environment

Install [Pixi](https://pixi.sh), enter this directory, then run:

```sh
pixi install --locked
```

The lock supports Apple Silicon macOS and Linux x86-64. It pins Python 3.8.20, XGBoost 2.0.3, SHAP 0.44.1, NumPy 1.24.4, pandas 2.0.3, and R 4.4 with ggplot2/reticulate. Micromamba provides the Conda-compatible Python runner required by reticulate inside Pixi; caches stay within the selected output directory. The environment was executed on Apple Silicon macOS; the Linux lock was resolved but not executed here. R calls the original Python SHAP library through `reticulate`, preserving its beeswarm and bar algorithms.

## Required inputs

Supply the archived, prepared `dipper_train.csv` (1,432 rows) and `dipper_test.csv` (614 rows) together in a data directory. The verified five-feature `core_model.json` is included under `weights/` and selected by default. The optional archived `dipper_test_with_results_core.csv` enables the historical prediction check. The prepared CSV inputs are **not included in this repository**; obtain them from the authors or the accompanying modelling package. The fitted weights and their metadata are committed here.

Expected SHA-256 checksums for the inputs used to verify this folder:

| File | SHA-256 |
| --- | --- |
| `dipper_train.csv` | `c66b665f7ddef5497418d1cb0677d03859f272a7a1e0133614ed41cde657f4b0` |
| `dipper_test.csv` | `60e8616cc5fa788a6714554381ffe70b7f63574899b1d22bcd92c29d00105dea` |
| Verified refit `core_model.json` | `91789bebb377419575303c036cc45381fb1dd5461b8fbb266baf566e5c579289` |

The five features, in original order, are `percGC_pbs`, `DeepCas9`, `pridict2_score_mS_HEK`, `length_mS_DPr04`, and `MCC_mS_coverage_bin1000`. The outcome is the archived centred `L2FC`. No new split, recentering or normalization is performed. This folder starts from the prepared archived data; it does not regenerate missing upstream annotations or raw-read processing.

## Replay saved weights and draw all three figures

Replace the example paths with your input locations:

```sh
pixi run --locked reproduce \
  --data-dir /path/to/prepared_data \
  --reference /path/to/dipper_test_with_results_core.csv \
  --output output/replay
```

The committed `weights/core_model.json` is loaded without fitting. Use `--model /path/to/another_model.json` to select other weights explicitly. Output includes `core_model.json`, `test_predictions.csv`, `training_features.csv`, `shap_values.npy`, `mean_absolute_shap.csv`, `metrics.json`, and separate `figure2i.pdf`, `figure2j.pdf`, `figure2k.pdf`. SHAP uses the 1,432 training rows, while panel i uses the 614 held-out rows. Reload equality, feature order, split separation, finite inputs, and SHAP additivity are checked. When a reference is supplied, the command fails if predictions differ by more than 1e-5.

To redraw existing analysis outputs without recomputing predictions/SHAP:

```sh
pixi run --locked reproduce --figures-only --output output/replay
```

To run just the model analysis without plotting:

```sh
pixi run --locked prepare --data-dir /path/to/prepared_data \
  --model weights/core_model.json --output output/replay
```

## Refit the original model

Use `--refit` to fit the original constructor on the archived split instead of loading the committed weights:

```sh
pixi run --locked reproduce --refit --data-dir /path/to/prepared_data --output output/refit
```

Or run `pixi run --locked train --data-dir /path/to/prepared_data --output output/refit` for model analysis only. The constructor retains 400 estimators, depth 6, learning rate 0.05, column/sample fractions 0.8, gamma 0.1 and squared-error regression, with the original unspecified defaults. No tuning or replacement R model is introduced.

**A fresh fit is not guaranteed to recover historical weights across platforms.** Previous same-version Mac fitting gave R 0.72066; the supplied verified refit gives R 0.72751735, matching archived test predictions within 1.364e-7. Its fitted object is a verified refit, not a recovered historical model object. For historical figure replay, use the default committed weights. Add `--reference` to a refit command to check rather than assume historical agreement; a mismatching candidate and its metrics are retained, and plotting stops.

## Source provenance and changes

Panel i copies the plotting expression and `theme_pie` directly from `../Anchored screen/Anchored screen plots.R`, upstream commit `e1a1a58711bf8237b0ed565d3155c8f70194fccc`. Changes are local input/output paths, qualified `egg::set_panel_size`, closing the PDF device reliably, and sizing the PDF canvas to include the labels around the original 4-cm panel with 2-mm edge padding (the upstream 4-cm canvas clips them). Point styling, scales, identity line, labels and R annotation retain the upstream code.

The model feature list, constructor, test evaluation and SHAP calls come from the manuscript version of `dels-main/scripts/submission1/DP9_XGBoost.ipynb` (SHA-256 `042d35361abef8688975468a01990f8e59c5c90b6e656a10551ccb6c0c8e1217`). The older standalone local `scripts/submission1` notebook has a different six-feature core list and would not reproduce these five manuscript features. The notebook itself is not duplicated here because only runnable scripts are needed.

Panels j and k retain the original `shap.summary_plot` algorithms, called sequentially from the combined R script through reticulate. Presentation edits use short manuscript labels, a smaller export size, seeded beeswarm layout, beige bars/thin borders, and a shorter bar-axis label to prevent clipping. No manual SHAP jitter or replacement importance metric is used. The final published figure assembly was not recovered; these are separate source-based panels rather than a claim of pixel-identical published page layout.

## Original DP007 figure script

`DP007_Figure2.R` is an **unchanged** copy of the supplied `scripts/DP007_Figure2.R`, not a replacement for the current manuscript model. Its SHA-256 is `f3c2648a46ee387f742545e46c480198bc5ac616f05a29e3516b31470cb3c31c`.

It preserves the original libraries, old feature names, plotting blocks, local paths, and statistical calculations. It requires historical inputs absent from the available prepared core-model data, including `51-14_counts_featurized.csv`, `51-14_counts_featurized_core.csv`, `DP102_mergeddata2_pbs.csv`, `DP102_mergeddata2_scaffoldsubset.csv`, and `DP102_mergeddata2_overhangs.csv`, plus the packages listed at its top. Its `/Users/jw38/repos/dipper/` paths must be adjusted to the location of those original files before running it. These dependencies are not claimed to be covered by the core-model Pixi environment.

The legacy script contains interactive exploratory sections and repeated export filenames; it has been preserved rather than silently repaired or remapped to incompatible prepared data. The `reproduce` task runs the consolidated manuscript-model figures only. Full execution of the historical DP007 script has not been verified because its original inputs are missing.
