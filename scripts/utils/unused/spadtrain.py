#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import re as re
import collections as col
import traceback as trb
import argparse as argp
import functools as fnt
import warnings as warn

import pandas as pd
import numpy as np

from sklearn.model_selection import GridSearchCV
from sklearn.metrics import f1_score, make_scorer, accuracy_score, r2_score
from scipy.stats import spearmanr as spr
from sklearn.ensemble import RandomForestClassifier as RFcls
from sklearn.ensemble import RandomForestRegressor as RFreg
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

"""
Some sandbox playing with species-adaptive co-training
This is not mature code!
"""

__ROOTDIR__ = '/TL/deep/fhgfs/projects/pebert/thesis/projects/cross_species/processing/playground/spadcotrain'
#__TRAINSET_CD4__ = 'G1240_EB16_TS30_mm9_ncd4.to.hg19.feat.seqout.h5'
#__TESTSET__ = 'G1240_EB16_TD29_hg19_ncd4.from.mm9.feat.seqout.h5'
#__TESTSET__ = 'G0151_EE12_TD21_hg19_liver.from.mm9.feat.seqout.h5'
__TRAINSET_LIVER__ = 'G0151_EE12_TE12_mm9_liver.to.hg19.feat.seqout.h5'
__TRAINSET_CH12__ = 'G0240_EE01_TE01_mm9_CH12.to.hg19.feat.seqout.h5'
#__TESTSET__ = 'G0240_EE01_TE04_hg19_CH12.from.mm9.feat.seqout.h5'
__TRAINSET__ = 'G0151_ED13_TD22_hg19_hepa.to.mm9.feat.seqout.h5'
__TESTSET__ = 'G0151_ED13_TS33_mm9_hepa.from.hg19.feat.seqout.h5'
#__STABLEGENES__ = 'mm9_stable_on_off.txt'
__STABLEGENES__ = 'mm9_stable_on_off.txt'
__ORTHOLOGS__ = 'odb9_6species.h5'

__MODEL__ = 'rfcls'

__TRAINMERGE__ = False

__CUTOFF__ = 1

#__CLUSTER_FEAT__ = 'ft(oecpg|cpg|gc)_\w+_reg5p$'
__CLUSTER_FEAT__ = 'ft((oecpg|cpg|gc)_\w+_reg5p$|msig_H3K(27ac|4me3)_abs_mean_reg5p$|msig_H3K36me3_abs_mean_body$)'


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--train', type=str, dest='train', nargs='+',
                        default=[os.path.join(__ROOTDIR__, __TRAINSET__),
                                 os.path.join(__ROOTDIR__, __TRAINSET_LIVER__),
                                 os.path.join(__ROOTDIR__, __TRAINSET_CH12__)])
    parser.add_argument('--test', type=str, dest='test', default=os.path.join(__ROOTDIR__, __TESTSET__))
    parser.add_argument('--init-genes', type=str, dest='initgenes', default=os.path.join(__ROOTDIR__, __STABLEGENES__))
    parser.add_argument('--orth', type=str, dest='orthologs', default=os.path.join(__ROOTDIR__, __ORTHOLOGS__))
    args = parser.parse_args()
    return args


def spearman(a, b):
    """
    :param a:
    :param b:
    :return:
    """
    r, _ = spr(a, b)
    return r


def pct(dataset, feat):
    """
    :param dataset:
    :param feat:
    :return:
    """
    vals = np.percentile(dataset[feat], [25, 50, 75, 90, 100])
    print('25pct >> ', vals, ' >> 100pct')
    return


def merge_traindata(datafiles, feat, subset):
    """
    :param datafiles:
    :param feat:
    :param subset:
    :return:
    """
    subsets = []
    for fp in datafiles:
        dsname = os.path.basename(fp).split('_')[1]
        data, labels, feat_cols = load_mldata(fp, feat, subset, 'train')
        data.index = data.index + '_' + dsname
        data['label'] = labels
        sub_pos = data.loc[data['label'] == 1, :].sample(n=700)
        sub_neg = data.loc[data['label'] == 0, :].sample(n=700)
        subsets.append(sub_pos)
        subsets.append(sub_neg)
    traindata = pd.concat(subsets, ignore_index=False, copy=True)
    print('Mrg traindata ', traindata.shape)
    trainlabels = traindata['label'].values
    traindata.drop('label', axis=1, inplace=True)
    return traindata[feat_cols], trainlabels, feat_cols


def load_mldata(fpath, feat, subset, mode):
    """
    :param fpath:
    :param feat:
    :return:
    """
    chromsets = []
    feat_cols = []
    with pd.HDFStore(fpath, 'r') as hdf:
        for k in hdf.keys():
            if k == '/metadata':
                continue
            chromdata = hdf[k]
            chromdata.index = chromdata['name']
            #print([c for c in chromdata.columns if c.startswith('ftmsig')])
            #raise
            if subset and mode == 'train':
                chromdata = chromdata.loc[chromdata.index.isin(subset), :]
            elif subset and mode == 'test':
                chromdata = chromdata.loc[~chromdata.index.isin(subset), :]
            else:
                pass
            if feat == 'all':
                feat_cols = sorted([c for c in chromdata.columns if c.startswith('ft')])
            else:
                match = re.compile(feat)
                feat_cols = sorted([c for c in chromdata.columns if match.match(c) if not None])
            chromdata = chromdata.loc[:, feat_cols + ['tpm_norm']]
            chromsets.append(chromdata)
    full_set = pd.concat(chromsets, axis=0, ignore_index=False)
    if mode == 'train':
        class_zero = (full_set['tpm_norm'] < __CUTOFF__).sum()
        pos_set = full_set.loc[full_set['tpm_norm'] >= __CUTOFF__, :].sample(min(class_zero, 2000))
        neg_set = full_set.loc[full_set['tpm_norm'] < __CUTOFF__, :].sample(min(class_zero, 2000))
        full_set = pd.concat([pos_set, neg_set], ignore_index=False, axis=0)
    if __MODEL__ == 'rfcls':
        labels = np.array((full_set['tpm_norm'] >= __CUTOFF__).values, dtype=np.int8)
    else:
        labels = np.log1p(np.array((full_set['tpm_norm']).values, dtype=np.float16))
    full_set.drop(['tpm_norm'], axis=1, inplace=True)
    assert labels.size == full_set.shape[0], 'Size mismatch'
    return full_set, labels, feat_cols


def load_orthologs(fpath):
    """
    :param fpath:
    :return:
    """
    with pd.HDFStore(fpath, 'r') as hdf:
        data = hdf['/auto/pairs/mouse/human']
    return data


def cluster_init_genes(dataset, labels):
    """
    :param dataset:
    :param labels:
    :return:
    """
    pca = PCA(n_components=3, whiten=True)
    red_space = pca.fit_transform(dataset)

    km = KMeans(n_clusters=2, max_iter=500, n_init=10, n_jobs=10,
                algorithm='elkan', precompute_distances=True)
    cdist = km.fit_transform(red_space)
    cdist = pd.DataFrame(cdist, dtype=np.float16, index=dataset.index, columns=['center0', 'center1'])
    cdist['cluster'] = km.labels_
    cdist['label'] = labels
    zero_zero = ((cdist['label'] == 0) & (cdist['cluster'] == 0)).sum()
    zero_one = ((cdist['label'] == 0) & (cdist['cluster'] == 1)).sum()
    if zero_zero > zero_one:
        cluster_off = 0
        cluster_on = 1
        sort_off, sort_on = 'center0', 'center1'
    else:
        cluster_off = 1
        cluster_on = 0
        sort_off, sort_on = 'center1', 'center0'
    cdist.sort_values(by=[sort_off], axis=0, ascending=True, inplace=True)
    genes_off = cdist.loc[(cdist['cluster'] == cluster_off), :].index.tolist()
    #genes_off = cdist.loc[(cdist['cluster'] == cluster_off) & (cdist['label'] == 0), :].index.tolist()

    cdist.sort_values(by=[sort_on], axis=0, ascending=True, inplace=True)
    genes_on = cdist.loc[(cdist['cluster'] == cluster_on), :].index.tolist()
    #genes_on = cdist.loc[(cdist['cluster'] == cluster_on) & cdist['label'] == 1, :].index.tolist()
    dataset['label'] = labels
    select_num = min(1000, len(genes_off), len(genes_on))
    select_genes = genes_on[:select_num] + genes_off[:select_num]
    seed_dataset = dataset.loc[select_genes, :]
    return seed_dataset


def cluster_comb_init_genes(traindata, trainlabels, testdata, testlabels, featcols):
    """
    :param traindata:
    :param trainlabels:
    :param testdata:
    :param testlabels:
    :param featcols:
    :return:
    """
    traindata['label'] = trainlabels
    testdata['label'] = -1

    joint = pd.concat([traindata, testdata], ignore_index=False, copy=True)

    pca = PCA(n_components=3, whiten=True)
    red_space = pd.DataFrame(pca.fit_transform(joint[featcols]), dtype=np.float16,
                             columns=['pc1', 'pc2', 'pc3'], index=joint.index)
    red_space['label'] = joint['label']
    mean_zero = red_space.loc[red_space['label'] == 0, ['pc1', 'pc2', 'pc3']].mean(axis=0)
    mean_one = red_space.loc[red_space['label'] == 1, ['pc1', 'pc2', 'pc3']].mean(axis=0)
    init_means = np.array([mean_zero, mean_one])

    subset_samples = set(red_space.index.tolist())
    for i, s in zip(range(4), [5000, 3500, 2500, 2000]):
        km = KMeans(n_clusters=2, max_iter=500, n_jobs=10, init=init_means,
                    algorithm='elkan', precompute_distances=True)
        subset = red_space.loc[red_space.index.isin(subset_samples), :].copy()
        cdist = km.fit_transform(subset[['pc1', 'pc2', 'pc3']])
        cdist = pd.DataFrame(cdist, dtype=np.float16, index=subset.index,
                             columns=['center0', 'center1'])
        cdist['cluster'] = km.labels_
        cdist['label'] = subset['label']

        cdist.sort_values(by=['center0'], axis=0, ascending=True, inplace=True)
        cluster0 = set(cdist.iloc[0:s, :].index.tolist())
        subset_samples = set(traindata.index.tolist()).union(cluster0)

        cdist.sort_values(by=['center1'], axis=0, ascending=True, inplace=True)
        cluster1 = set(cdist.iloc[0:s, :].index.tolist())
        subset_samples = subset_samples.union(cluster1)

    return subset_samples


def select_seed_samples(traindata, trainlabels, testdata, testlabels, initsamples,
                        featcols, remtrain, remlabels, orthos):
    """
    :param traindata:
    :return:
    """
    print('### Seed model')

    if __MODEL__ == 'rfcls':
        f1 = make_scorer(f1_score, **{'average': 'macro'})
        params = {'n_estimators': [500, 750, 1000, 1500], 'min_samples_split': [2, 4, 8],
                  'oob_score': [True], 'bootstrap': [True]}
        cv = GridSearchCV(RFcls(), params, scoring=f1, n_jobs=14, cv=5, refit=True)
    else:
        #r2 = make_scorer(r2_score)
        spear = make_scorer(spearman, greater_is_better=True)
        params = {'n_estimators': [200, 400, 600], 'min_samples_split': [2, 4, 8],
                  'oob_score': [True], 'bootstrap': [True]}
        cv = GridSearchCV(RFreg(), params, scoring=spear, n_jobs=14, cv=5, refit=True)
    cv.fit(traindata[featcols], trainlabels)
    print('Seed traindata shape ', traindata.shape)
    seed_model = cv.best_estimator_
    print(cv.best_params_)
    print('Seed CV perf ', cv.best_score_)
    print('Seed model remainder perf ', seed_model.score(remtrain[featcols], remlabels))
    print('Source testdata shape: ', testdata.shape)
    print('Seed model testdata perf ', seed_model.score(testdata[featcols], testlabels))
    print('=== against orthologs ')
    testdata['label'] = testlabels
    orthoset = testdata.loc[testdata.index.isin(orthos['mouse_name']), :]
    remtrain['label'] = remlabels
    traindata['label'] = trainlabels
    fulltrain = pd.concat([traindata, remtrain], ignore_index=False, copy=True)
    print('Model ortho perf ', seed_model.score(orthoset[featcols], orthoset['label']))

    if __MODEL__ in ['rfcls', 'rfreg']:
        seed_ftimp = seed_model.feature_importances_
        for n, i in sorted(zip(featcols, seed_ftimp), key=lambda x: x[1], reverse=True):
            print(n, i)
    # remain_test = testdata.copy()
    # satrain = None
    # for i in range(10):
    #     print('###==== ', i)
    #     with warn.catch_warnings():
    #         warn.simplefilter('ignore')
    #         test_fraction = cluster_comb_init_genes(traindata[featcols], trainlabels, remain_test, [], featcols)
    #     test_fraction = remain_test.loc[remain_test.index.isin(test_fraction), :].copy()
    #     print('Num test frac samples: ', test_fraction.shape)
    #     frac_labels = rf_seed_model.predict(test_fraction[featcols])
    #     frac_table = pd.DataFrame(rf_seed_model.predict_proba(test_fraction[featcols]),
    #                               dtype=np.float32, columns=['0', '1'],
    #                               index=test_fraction.index)
    #     print('Perf seed model test fraction: ', rf_seed_model.score(test_fraction[featcols], test_fraction['true_label']))
    #     cv = GridSearchCV(RFcls(), params, scoring=f1, n_jobs=14, cv=5, refit=True)
    #     test_fraction['label'] = frac_labels
    #     if satrain is None:
    #         traindata['label'] = trainlabels
    #         satrain = pd.concat([traindata, test_fraction], ignore_index=False, copy=True)
    #     else:
    #         satrain = pd.concat([satrain, test_fraction], ignore_index=False, copy=True)
    #     cv.fit(satrain[featcols], satrain['label'])
    #     sa_model = cv.best_estimator_
    #     print('SA model model remainder :', sa_model.score(remtrain[featcols], remlabels))
    #     remain_test = testdata.loc[~testdata.index.isin(satrain.index), :].copy()
    #     print('Remaining test data: ', remain_test.shape)
    #     if remain_test.empty:
    #         break
    return 0


def main():
    """
    :return:
    """
    args = parse_command_line()
    stable_genes = [gid.strip() for gid in open(args.initgenes, 'r').readlines()]
    orthos = load_orthologs(args.orthologs)
    if not __TRAINMERGE__:
        traindata, trainlabels, trainfeat = load_mldata(args.train[0], __CLUSTER_FEAT__, [], 'train')
    else:
        traindata, trainlabels, trainfeat = merge_traindata(args.train, __CLUSTER_FEAT__, stable_genes)
    #seed_train = cluster_init_genes(traindata, trainlabels)
    remtrain, remlabels, _ = load_mldata(args.train[0], __CLUSTER_FEAT__, traindata.index.tolist(), 'test')
    testdata, testlabels, testfeat = load_mldata(args.test, __CLUSTER_FEAT__, [], 'test')
    assert trainfeat == testfeat, 'Wrong feature sort order'
    with warn.catch_warnings():
        warn.simplefilter('ignore')
        #init_samples = cluster_comb_init_genes(traindata, trainlabels, testdata, testlabels, trainfeat)
    _ = select_seed_samples(traindata, trainlabels, testdata, testlabels, [], trainfeat,
                            remtrain, remlabels, orthos)
    #_ = select_seed_samples(seed_train, remtrain, remlabels, testdata, testlabels)
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
