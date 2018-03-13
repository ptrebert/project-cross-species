#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import time as ti

import pandas as pd
import numpy as np

from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer, r2_score, roc_auc_score
from sklearn.ensemble import RandomForestRegressor as RFreg
from sklearn.ensemble import RandomForestClassifier as RFcls
from scipy.stats import spearmanr, pearsonr

__ROOTDIR__ = '/TL/deep/fhgfs/projects/pebert/thesis/refdata/geneage'

# different feature sets are used - the below lists contain features
# that are NOT to be used in training

AGE_FEAT = ['ftorth_group', 'ftorth_pair', 'ftconf_bin_hgt',
            'ftage_bin_euk', 'ftage_bin_vert', 'ftage_bin_opist',
            'ftage_bin_eukbac', 'ftage_bin_mamm', 'ftage_bin_cellorg',
            'ftage_bin_eukarch', 'ftage_bin_eumeta',
            'ftcons_pct_body', 'ftcons_pct_reg5p']

BOTH_FEAT = ['ftorth_group', 'ftorth_pair', 'ftconf_bin_hgt',
             'ftage_bin_euk', 'ftage_bin_vert', 'ftage_bin_opist',
             'ftage_bin_eukbac', 'ftage_bin_mamm', 'ftage_bin_cellorg',
             'ftage_bin_eukarch', 'ftage_bin_eumeta']

DATA_FEAT = ['ftorth_group', 'ftorth_pair', 'ftconf_bin_hgt',
             'ftage_bin_euk', 'ftage_bin_vert', 'ftage_bin_opist',
             'ftage_bin_eukbac', 'ftage_bin_mamm', 'ftage_bin_cellorg',
             'ftage_bin_eukarch', 'ftage_bin_eumeta',
             'ftconf_abs_bimodality', 'ftconf_abs_entropy', 'ftconf_abs_nodeerror']

# joint predictions:
# chicken cow rat rhesus

# single predictions
# dog human  mouse opossum


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--train', action='store_true', dest='train', default=True)
    parser.add_argument('--jobs', '-j', type=int, default=4, dest='jobs')
    parser.add_argument('--root', type=str, default=__ROOTDIR__, dest='root')
    parser.add_argument('--data', type=str, default=os.path.join(__ROOTDIR__, 'norm'), dest='data')
    parser.add_argument('--model', type=str, default=os.path.join(__ROOTDIR__, 'model'), dest='model')
    parser.add_argument('--out', type=str, dest='outfile', required=True)
    parser.add_argument('--test', action='store_true', dest='test', default=False)
    parser.add_argument('--classify', action='store_true', default=False, dest='classify')
    args = parser.parse_args()
    return args


def spearman_score(a, b):
    """
    :param a:
    :param b:
    :return:
    """
    corr, _ = spearmanr(a, b)
    return np.round(corr, 3)


def pearson_score(a, b):
    """
    :param a:
    :param b:
    :return:
    """
    corr, _ = pearsonr(a, b)
    return np.round(corr, 3)


def coeff_det_score(a, b):
    """
    """
    r2 = r2_score(a, b)
    return np.round(r2, 3)


def train_model(dataset, args, unused_feat):
    """
    :param dataset:
    :return:
    """
    feat_cols = [c for c in dataset.columns if c.startswith('ft') and c not in unused_feat]

    rec_cols = ['n_estimators', 'min_samples_split',
                'oob_score', 'cv_score', 'test_score']
    rec_cols.extend(feat_cols)

    rec_run = pd.DataFrame(np.zeros((10, len(rec_cols)), dtype=np.float32),
                           index=range(10), columns=rec_cols)
    oob_est = pd.DataFrame(np.ones((dataset.shape[0], 10), dtype=np.float32),
                           columns=range(10), index=dataset['name'].values)
    oob_est *= -1

    if args.classify:
        cols = []
        for i in range(10):
            cols.append('{}_class'.format(i))
            cols.append('{}_pos_prob'.format(i))
        test_est = pd.DataFrame(np.ones((dataset.shape[0], 20), dtype=np.float32),
                                columns=cols, index=dataset['name'].values)
        test_est *= -1
    else:
        test_est = pd.DataFrame(np.ones((dataset.shape[0], 10), dtype=np.float32),
                                columns=range(10), index=dataset['name'].values)
        test_est *= -1
    for i in range(10):
        # shuffle dataset
        dataset = dataset.sample(frac=1).reset_index(drop=True)

        train = dataset.sample(frac=0.67, replace=False, axis=0)

        train_idx = dataset.index.isin(train.index)
        test_idx = ~train_idx

        test = dataset.loc[test_idx, :].copy()
        X_test = test[feat_cols]
        y_test = test['target']
        if args.classify:
            y_test = np.array(y_test > 0.5, dtype=np.int8)
        test_names = test['name'].values

        X_train = train[feat_cols]
        y_train = train['target']
        if args.classify:
            y_train = np.array(y_train > 0.5, dtype=np.int8)
        train_names = train['name'].values

        if args.classify:
            params = {'n_estimators': [250, 500, 750, 1000],
                      'min_samples_split': [8, 16, 32, 64, 128],
                      'criterion': ['gini'],
                      'oob_score': [True], 'bootstrap': [True],
                      'class_weight': ['balanced']}

            # params = {'n_estimators': [250],
            #           'min_samples_split': [8, 16],
            #           'criterion': ['gini'],
            #           'oob_score': [True], 'bootstrap': [True],
            #           'class_weight': ['balanced']}

            cv = GridSearchCV(RFcls(), params, scoring='roc_auc', n_jobs=args.jobs, cv=10, refit=True, verbose=1)
        else:
            #spr = make_scorer(spearman_score, greater_is_better=True)
            r2s = make_scorer(coeff_det_score, greater_is_better=True)
            params = {'n_estimators': [50, 75, 100, 150, 200],
                      'min_samples_split': [8, 16, 32, 64, 128],
                      #'criterion': ['mse', 'mae'],  # MAE criterion has memory issues
                      'criterion': ['mse'],
                      'oob_score': [True], 'bootstrap': [True]}

            cv = GridSearchCV(RFreg(), params, scoring=r2s, n_jobs=args.jobs, cv=10, refit=True, verbose=1)
        cv.fit(X_train, y_train)
        model = cv.best_estimator_
        if args.classify:
            # for classification, save probability of positive class
            # for OOB samples in analogy to what is done later for
            # the test set
            oob_pred = model.oob_decision_function_[:, 1]
        else:
            oob_pred = model.oob_prediction_
        oob_est.loc[train_names, i] = oob_pred
        rec_run.loc[i, 'cv_score'] = cv.best_score_
        for name, value in cv.best_params_.items():
            if name in ['bootstrap', 'class_weight']:
                continue
            elif name == 'oob_score':
                rec_run.loc[i, name] = model.oob_score_
            else:
                rec_run.loc[i, name] = value
        ftimp = model.feature_importances_
        ftimp = [(n, i) for n, i in zip(feat_cols, ftimp)]
        for name, imp in ftimp:
            rec_run.loc[i, name] = imp

        if args.classify:
            y_prob = pd.DataFrame(model.predict_proba(X_test),
                                  columns=[0, 1], index=test_names)
            y_hat = model.predict(X_test)
            pred_class_prob = y_prob.lookup(test_names, y_hat)
            pos_class_prob = pred_class_prob
            pos_class_prob[y_hat == 0] = (1 - pred_class_prob[y_hat == 0])
            perf = roc_auc_score(y_test, pos_class_prob)
            rec_run.loc[i, 'test_score'] = perf
            test_est.loc[test_names, '{}_pos_prob'.format(i)] = pos_class_prob
            test_est.loc[test_names, '{}_class'.format(i)] = y_hat
        else:
            y_pred = model.predict(X_test)
            perf = r2_score(y_test, y_pred)
            test_est.loc[test_names, i] = y_pred
            rec_run.loc[i, 'test_score'] = perf

    return rec_run, oob_est, test_est


def load_datasets(path, setting):
    """
    :param path:
    :return:
    """
    all_files = os.listdir(path)
    assert len(all_files) > 0, 'No feature files detected under path {}'.format(path)
    mrg_data = []
    for fn in all_files:
        date, spec, _ = fn.split('_')
        if setting == 'single':
            if spec not in ['dog', 'human', 'mouse', 'opossum']:
                continue
        elif setting == 'joint':
            if spec not in ['chicken', 'cow', 'rat', 'rhesus']:
                continue
        else:
            pass
        fp = os.path.join(path, fn)
        with pd.HDFStore(fp, 'r') as hdf:
            all_keys = [k for k in hdf.keys() if k.endswith('/final')]
            if any(['joint' in k for k in all_keys]):
                load_key = '/model/joint/final'
            else:
                assert len(all_keys) == 1, 'Too many keys: {}'.format(all_keys)
                load_key = all_keys[0]
            sub_data = hdf[load_key]
            assert not pd.isnull(sub_data).any(axis=0).any(), 'NULL data loaded: {}'.format(load_key)
            sub_data['species'] = spec
            mrg_data.append(sub_data)
    ml_data = pd.concat(mrg_data, axis=0, ignore_index=True)
    assert not pd.isnull(ml_data).any(axis=0).any(), 'NULL after concat'
    keep_cols = [c for c in ml_data.columns if c.startswith('ft')]
    keep_cols.extend(['name', 'target', 'species'])
    ml_data = ml_data.loc[:, keep_cols].copy()
    assert not pd.isnull(ml_data).any(axis=0).any(), 'NULL after drop'
    return ml_data


def main():
    """
    :return:
    """
    args = parse_command_line()
    filemode = 'w'
    if args.train:
        for grp_name, grp_feat in [('data', DATA_FEAT),
                                   ('age', AGE_FEAT),
                                   ('both', BOTH_FEAT)]:
            for setting in ['single', 'joint', 'merged']:
                ml_dataset = load_datasets(args.data, setting)
                perf, oob_est, test_est = train_model(ml_dataset, args, grp_feat)
                outpath = os.path.join(args.model, args.outfile)
                with pd.HDFStore(outpath, filemode, complevel=9) as hdf:
                    hdf.put('{}/{}/data'.format(grp_name, setting), ml_dataset, format='table')
                    hdf.put('{}/{}/perf'.format(grp_name, setting), perf, format='fixed')
                    hdf.put('{}/{}/oob'.format(grp_name, setting), oob_est, format='fixed')
                    hdf.put('{}/{}/test'.format(grp_name, setting), test_est, format='fixed')
                filemode = 'a'
    elif args.test:
        pass
    else:
        pass
    return 0


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
