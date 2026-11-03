# Verified core fitted weights

`core_model.json` is the verified five-feature XGBoost refit previously supplied in `R1C4_predictive_modelling/outputs/core/`. It is copied unchanged. `core_metadata.json` is its original feature-order, parameter, version and input-checksum metadata. No original historical fitted object was recovered; this refit reproduces the archived 614 held-out predictions within 1.364e-7, with Pearson R = 0.7275173509756766.

SHA-256 checksums are in `SHA256SUMS`. The default `pixi run reproduce --data-dir PATH` loads these committed weights without retraining. The prepared input CSVs must be supplied separately. Fresh cross-platform fitting need not recover these same weights.
