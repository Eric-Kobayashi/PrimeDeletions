With [Pixi](https://pixi.sh) installed, place `dipper_train.csv` and `dipper_test.csv` in `data/` (inputs supplied separately), then run from this directory:

```sh
pixi run --locked reproduce --data-dir data
```

Uses the committed weights; writes analysis results and Figure 2i–k PDFs to `output/`. Add `--refit` to retrain; results may differ.

`DP007_Figure2.R` requires its original CSV inputs, libraries and local path adjustments.
