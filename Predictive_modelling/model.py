"""Original DP9 core XGBoost fit, saved-weight replay, and exact TreeSHAP export."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import shap
import xgboost

# Original dels-main/scripts/submission1/DP9_XGBoost.ipynb, cell 31.
FEATURES = ['percGC_pbs', 'DeepCas9', 'pridict2_score_mS_HEK',
            'length_mS_DPr04', 'MCC_mS_coverage_bin1000']
LABELS = ['PBS %GC', 'DeepCas9', 'PRIDICT2', 'Deletion length', 'Contact frequency']


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# Original notebook cell 2: constructor and fit body retained.
def xgb(x, y, do_fit=False):
    model = xgboost.XGBRegressor(
        n_estimators=400,
        max_depth=6,
        learning_rate=0.05,
        colsample_bytree=0.8,
        subsample=0.8,
        gamma=0.1,
        objective='reg:squarederror',
        enable_categorical=True,
    )
    if do_fit:
        model.fit(x, y)
    return model


def reproduce(data_dir, output_dir, model_file=None, reference_file=None):
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    traindata = pd.read_csv(data_dir / 'dipper_train.csv')
    testdata = pd.read_csv(data_dir / 'dipper_test.csv')
    for data in (traindata, testdata):
        if not np.isfinite(data[FEATURES + ['L2FC']].to_numpy()).all():
            raise ValueError('Core inputs and outcomes must be finite.')
    for column in ('unique', 'spacer'):
        if column in traindata and column in testdata:
            if set(traindata[column]) & set(testdata[column]):
                raise ValueError('Train/test overlap in ' + column)
    X_train = traindata[FEATURES]
    y_train = traindata['L2FC']
    if model_file is None:
        model = xgb(X_train, y_train, do_fit=True)
    else:
        model = xgboost.XGBRegressor()
        model.load_model(str(model_file))
    if model.get_booster().feature_names != FEATURES:
        raise ValueError('Saved model feature order differs from the original core model.')

    # Original notebook cell 35: held-out evaluation, without a new split.
    testdata_model = testdata.copy().reset_index(drop=True)
    testdata_model['L2FC_pred'] = model.predict(testdata_model[FEATURES])
    corr, _ = stats.pearsonr(testdata_model['L2FC_pred'],
                           testdata_model['L2FC'].values.reshape(-1))
    testdata_model.to_csv(output_dir / 'test_predictions.csv', index=False)
    model.save_model(str(output_dir / 'core_model.json'))
    reloaded = xgboost.XGBRegressor()
    reloaded.load_model(str(output_dir / 'core_model.json'))
    if not np.array_equal(model.predict(testdata[FEATURES]), reloaded.predict(testdata[FEATURES])):
        raise AssertionError('Saved-weight reload changed predictions.')

    # Original notebook cell 38: exact TreeSHAP on the TRAINING rows.
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(traindata[FEATURES])
    additive_error = float(np.max(np.abs(
        shap_values.sum(axis=1) + explainer.expected_value - model.predict(X_train))))
    if not np.isfinite(shap_values).all() or additive_error >= 1e-5:
        raise AssertionError('SHAP additivity or finiteness check failed.')
    np.save(output_dir / 'shap_values.npy', shap_values)
    X_train.to_csv(output_dir / 'training_features.csv', index=False)
    importance = pd.DataFrame({'feature': FEATURES, 'label': LABELS,
                               'mean_absolute_shap': np.abs(shap_values).mean(axis=0)})
    importance.sort_values('mean_absolute_shap', ascending=False).to_csv(
        output_dir / 'mean_absolute_shap.csv', index=False)
    metadata = {
        'mode': 'original fit' if model_file is None else 'saved weights',
        'features': FEATURES, 'labels': LABELS,
        'train_rows': len(traindata), 'test_rows': len(testdata),
        'pearson_r': float(corr), 'shap_additivity_max_error': additive_error,
        'reload_predictions_identical': True,
        'input_sha256': {name: sha256(data_dir / name) for name in
                         ['dipper_train.csv', 'dipper_test.csv']},
        'source_model_sha256': None if model_file is None else sha256(model_file),
        'exported_model_sha256': sha256(output_dir / 'core_model.json'),
        'versions': {'python': platform.python_version(), 'xgboost': xgboost.__version__,
                     'shap': shap.__version__, 'numpy': np.__version__, 'pandas': pd.__version__},
    }
    if reference_file:
        reference = pd.read_csv(reference_file)
        if len(reference) != len(testdata):
            raise ValueError('Reference has a different number of test rows.')
        if 'unique' not in reference or not reference['unique'].equals(testdata['unique']):
            raise ValueError('Reference identifiers/order differ from held-out data.')
        if not np.allclose(reference.L2FC, testdata.L2FC, rtol=0, atol=1e-10):
            raise ValueError('Reference measured outcomes differ.')
        difference = float(np.max(np.abs(reference.L2FC_pred - testdata_model.L2FC_pred)))
        metadata['reference_max_prediction_difference'] = difference
        metadata['matches_reference_at_1e_minus_5'] = difference < 1e-5
    (output_dir / 'metrics.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))
    if reference_file and not metadata['matches_reference_at_1e_minus_5']:
        raise RuntimeError('Fit differs from archived predictions; candidate retained for inspection. '
                           'Use prepare with verified saved weights for historical replay.')


def load_plot_inputs(output_dir):
    output_dir = Path(output_dir)
    metadata = json.loads((output_dir / 'metrics.json').read_text())
    features = pd.read_csv(output_dir / 'training_features.csv')
    values = np.load(output_dir / 'shap_values.npy')
    if (metadata['features'] != FEATURES or list(features.columns) != FEATURES
            or values.shape != features.shape or len(features) != metadata['train_rows']):
        raise ValueError('SHAP values and training feature matrix do not align.')
    return {'features': features, 'shap_values': values, 'labels': LABELS}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['train', 'prepare'])
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('output'))
    parser.add_argument('--model', type=Path)
    parser.add_argument('--reference', type=Path)
    args = parser.parse_args()
    if args.mode == 'prepare' and args.model is None:
        parser.error('prepare requires --model')
    if args.mode == 'train' and args.model is not None:
        parser.error('train does not accept --model; use prepare to load weights')
    reproduce(args.data_dir, args.output, args.model, args.reference)
