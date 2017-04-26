#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import fnmatch as fnm
import json as js

import numpy as np
import numpy.random as rng

import pandas as pd
from sklearn.metrics import accuracy_score as accuracy
from scipy.stats import kendalltau as ktau
from sklearn.metrics import r2_score as r2s


ASSM_SPECIES_MAP = {'hg19': 'human', 'mm9': 'mouse', 'canFam3': 'dog',
                    'galGal3': 'chicken', 'bosTau7': 'cow', 'susScr2': 'pig'}


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--train-test', '-tt', nargs='+', type=str, dest='traintest')
    parser.add_argument('--orthologs', '-orth', nargs='+', type=str, dest='orthologs')

    parser.add_argument('--output', '-o', type=str, dest='outputfile')

    args = parser.parse_args()
    return args


def kendall_tau_scorer(x, y):
    """
    :param x:
    :param y:
    :return:
    """
    statistic, p_value = ktau(x, y, nan_policy='raise')
    return statistic


def get_model_subpath(modeltype):
    """
    :param modeltype:
    :return:
    """
    paths = {'rfcls.asig.all': 'status/asig', 'rfcls.seq.all': 'status/seq',
             'rfreg.rk.all': 'rank/all', 'rfreg.rk.act': 'rank/active',
             'rfreg.tpm.all': 'level/all', 'rfreg.tpm.act': 'level/active'}
    return paths[modeltype]


def aggregate_testing_data(mdfiles, trainperf, orthologs, outfile):
    """
    :param mdfiles:
    :param trainperf:
    :param orthologs:
    :param outfile:
    :return:
    """
    saved = 0
    for mdf in mdfiles:
        dump = js.load(open(mdf, 'r'))
        modelfile = dump['run_info']['model_file'].rsplit('.', 1)[0]
        if modelfile not in trainperf:
            # some skipped, see aggregate_training_data
            continue
        datafile = dump['run_info']['data_file'].rsplit('.', 1)[0]
        datainfo, _, qry, mtype = modelfile.split('.', 3)
        gid, _, _, trg, biosample = datainfo.split('_')
        trg_spec, qry_spec = ASSM_SPECIES_MAP[trg], ASSM_SPECIES_MAP[qry]
        if 'rfcls' in modelfile:
            class_0 = [t[0] for t in dump['testing_info']['probabilities']['all']]
            class_1 = [t[1] for t in dump['testing_info']['probabilities']['all']]
            data = pd.DataFrame({'prob_class_0': class_0,
                                 'prob_class_1': class_1,
                                 'true': dump['testing_info']['targets']['true'],
                                 'pred': dump['testing_info']['targets']['pred']},
                                index=dump['sample_info']['names'],
                                columns=['prob_class_0', 'prob_class_1', 'true', 'pred'])
        else:
            if dump['dataset_info']['target_type'].startswith('int'):
                true_labels = np.round(dump['testing_info']['targets']['true'], decimals=0)
                pred_labels = np.round(dump['testing_info']['targets']['pred'], decimals=0)
                datatype = np.int32
            elif dump['dataset_info']['target_type'].startswith('float'):
                true_labels = dump['testing_info']['targets']['true']
                pred_labels = dump['testing_info']['targets']['pred']
                datatype = np.float64
            else:
                raise ValueError('Unexpected data type for regressor: {}'.format(dump['dataset_info']['target_type']))
            data = pd.DataFrame({'true': true_labels,
                                 'pred': pred_labels},
                                index=dump['sample_info']['names'],
                                columns=['true', 'pred'],
                                dtype=datatype)
        qry_gene_name = '{}_name'.format(qry_spec)
        data.index.name = '{}_name'.format(qry_spec)
        try:
            metadata = {'modelfile': modelfile, 'target': trg_spec, 'query': qry_spec, 'datafile': datafile,
                        'biosample': biosample, 'estimators': str(dump['model_info']['params']['n_estimators']),
                        'perf_{}'.format(dump['testing_info']['scoring']): str(dump['testing_info']['performance'])}
        except KeyError as ke:
            raise KeyError('Could not find key {} in {}'.format(str(ke), modelfile))
        # add ortholog indicator columns
        for orth_name, annotations in orthologs.items():
            if orth_name == 'hcop' and trg_spec != 'human':
                continue
            pairs = annotations[(trg_spec, qry_spec)]
            pw_idx = np.array(data.index.isin(pairs[qry_gene_name]), dtype=np.int8)
            data['ortho_pair_{}'.format(orth_name)] = pw_idx
            groups = annotations['groups']
            grp_idx = np.array(data.index.isin(groups[qry_gene_name]), dtype=np.int8)
            data['ortho_group_{}'.format(orth_name)] = grp_idx
        if 'rfcls' in modelfile:
            perf_md = collect_test_class_perf(data, list(orthologs.keys()))
            metadata.update(perf_md)
        else:
            perf_md = collect_test_reg_perf(data, list(orthologs.keys()), dump['testing_info']['scoring'])
            metadata.update(perf_md)
        metadata = pd.DataFrame.from_dict(metadata, orient='index')
        modelfile = modelfile.replace('.', '_')
        datafile = datafile.replace('.', '_')
        basepath = os.path.join('testing', get_model_subpath(mtype), trg_spec, qry_spec, gid, datafile, modelfile)
        with pd.HDFStore(outfile, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put(os.path.join(basepath, 'data'), data, format='fixed')
            hdf.put(os.path.join(basepath, 'metadata'), metadata, format='table')
            hdf.flush()
            saved += 1
    return saved


def collect_test_reg_perf(dataset, ortho_names, scoring):
    """
    :param dataset:
    :param ortho_names:
    :param scoring:
    :return:
    """
    if scoring.startswith('kendall'):
        score_name = 'ktau'
        scorer = kendall_tau_scorer
    elif scoring.startswith('r2'):
        score_name = 'r2'
        scorer = r2s
    else:
        raise ValueError('Unknown score requested: {}'.format(scoring))
    num_samples = dataset.shape[0]
    raw_perf = scorer(dataset['true'], dataset['pred'])
    perf_md = {'raw_num_samples': str(num_samples), 'raw_perf_' + score_name: str(raw_perf)}
    for name in ortho_names:
        try:
            subset = dataset.loc[dataset['ortho_pair_' + name] == 1, :].copy()
            perf_md[name + '_pair_num_samples'] = str(subset.shape[0])
            perf_md[name + '_pair_perf_' + score_name] = str(scorer(subset['true'], subset['pred']))
            subset = dataset.loc[dataset['ortho_group_' + name] == 1, :].copy()
            perf_md[name + '_group_num_samples'] = str(subset.shape[0])
            perf_md[name + '_group_perf_' + score_name] = str(scorer(subset['true'], subset['pred']))
        except KeyError:
            # HCOP
            continue
    return perf_md


def collect_test_class_perf(dataset, ortho_names):
    """
    :param dataset:
    :param ortho_names:
    :return:
    """
    num_samples = dataset.shape[0]
    tp, fp, tn, fn, acc, _ = collect_confusion_metrics(dataset)
    perf_md = {'raw_num_samples': str(num_samples), 'raw_num_tp': str(tp), 'raw_num_fp': str(fp),
               'raw_num_tn': str(tn), 'raw_num_fn': str(fn), 'raw_wt_acc': str(acc)}
    for name in ortho_names:
        try:
            tp, fp, tn, fn, acc, nums = collect_confusion_metrics(dataset.loc[dataset['ortho_pair_' + name] == 1, :].copy())
            perf_md[name + '_pair_num_tp'] = str(tp)
            perf_md[name + '_pair_num_fp'] = str(fp)
            perf_md[name + '_pair_num_tn'] = str(tn)
            perf_md[name + '_pair_num_fn'] = str(fn)
            perf_md[name + '_pair_wt_acc'] = str(acc)
            perf_md[name + '_pair_num_samples'] = str(nums)
            tp, fp, tn, fn, acc, nums = collect_confusion_metrics(dataset.loc[dataset['ortho_group_' + name] == 1, :].copy())
            perf_md[name + '_group_num_tp'] = str(tp)
            perf_md[name + '_group_num_fp'] = str(fp)
            perf_md[name + '_group_num_tn'] = str(tn)
            perf_md[name + '_group_num_fn'] = str(fn)
            perf_md[name + '_group_wt_acc'] = str(acc)
            perf_md[name + '_group_num_samples'] = str(nums)
        except KeyError:
            # HCOP annotation only for target = human
            continue
    return perf_md


def collect_confusion_metrics(dataset):
    """
    :param dataset:
    :return:
    """
    num_samples = dataset.shape[0]
    true_labels = dataset['true'].tolist()
    wt_one = 0.5 / dataset['true'].sum()
    wt_zero = 0.5 / (num_samples - dataset['true'].sum())
    wt_vec = [wt_one if l == 1 else wt_zero for l in true_labels]
    wt_acc = accuracy(dataset['true'], dataset['pred'], sample_weight=wt_vec)
    tp = np.logical_and(dataset['true'] == 1, dataset['pred'] == 1).sum()
    tn = np.logical_and(dataset['true'] == 0, dataset['pred'] == 0).sum()
    fp = np.logical_and(dataset['true'] == 0, dataset['pred'] == 1).sum()
    fn = np.logical_and(dataset['true'] == 1, dataset['pred'] == 0).sum()
    return tp, fp, tn, fn, wt_acc, num_samples


def aggregate_training_data(mdfiles, orthofiles, outfile):
    """
    :param mdfiles:
    :param orthofiles:
    :param outfile:
    :return:
    """
    saved = 0
    collect_perf = dict()
    sample_weights = dict()
    orthos = dict()
    for mdf in mdfiles:
        dump = js.load(open(mdf, 'r'))
        modelfile = os.path.basename(mdf).rsplit('.', 1)[0]
        assert modelfile not in collect_perf, 'Duplicate model: {}'.format(modelfile)
        collect_perf[modelfile] = dump['training_info']['best_score']
        datainfo, _, qry, mtype = modelfile.split('.', 3)
        gid, _, _, trg, biosample = datainfo.split('_')
        trg_spec, qry_spec = ASSM_SPECIES_MAP[trg], ASSM_SPECIES_MAP[qry]
        if trg_spec == 'mouse' and qry_spec == 'human' and gid == 'G0951':
            # this skips over the hepa/kidney pairing
            continue
        if trg_spec == 'human' and qry_spec == 'mouse' and gid == 'G0951':
            # this skips over the redundant hepa/liver pairing (same as G0151)
            continue
        load_orthologs(orthofiles, (trg_spec, qry_spec), orthos)
        if 'rfcls' in mtype:
            class_0 = [t[0] for t in dump['attribute_info']['oob_decision_function_']]
            class_1 = [t[1] for t in dump['attribute_info']['oob_decision_function_']]
            pred_class = []
            for p0, p1 in zip(class_0, class_1):
                if p1 > 0.5:
                    pred_class.append(1)
                elif p1 == 0.5:
                    assert np.isclose(p0, p1), 'Numerical imprecisions: {} and {}'.format(p0, p1)
                    # turns out, this happens for a handful of data points per model
                    #sys.stderr.write('\nDatapoint on decision boundary: {}\n'.format(modelfile))
                    decide = rng.random()
                    if decide == 0.5:
                        raise ValueError('You got to be kidding me...')
                    elif decide < 0.5:
                        pred_class.append(0)
                    else:
                        pred_class.append(1)
                else:
                    pred_class.append(0)
            data = pd.DataFrame({'prob_class_0': class_0,
                                 'prob_class_1': class_1,
                                 'true': dump['sample_info']['targets'],
                                 'pred': pred_class},
                                index=dump['sample_info']['names'],
                                columns=['prob_class_0', 'prob_class_1', 'true', 'pred'])
        else:
            data = pd.DataFrame({'pred': dump['attribute_info']['oob_prediction_'],
                                 'true': dump['sample_info']['targets']},
                                index=dump['sample_info']['names'],
                                columns=['pred', 'true'])
        trg_gene_name = '{}_name'.format(trg_spec)
        data.index.name = trg_gene_name
        try:
            metadata = {'modelfile': modelfile, 'target': trg_spec, 'query': qry_spec,
                        'biosample': biosample, 'estimators': str(dump['model_info']['params']['n_estimators']),
                        'perf_{}'.format(dump['training_info']['scoring']): str(dump['training_info']['best_score'])}
        except KeyError as ke:
            raise KeyError('Could not find key {} in {}'.format(str(ke), modelfile))
        # if 'rfcls' in mtype:
        #     for ortho_name, annotations in orthos:
        #         ortho_md = calculate_training_weights(data, trg_gene_name, annotations[(trg_spec, qry_spec)],
        #                                               annotations['groups'], ortho_name)
        #         if modelfile not in sample_weights:
        #             sample_weights[modelfile] = {}
        #         sample_weights[modelfile]['raw_wt_zero'] = float(ortho_md['raw_wt_zero'])
        #         sample_weights[modelfile]['raw_wt_one'] = float(ortho_md['raw_wt_one'])
        #         sample_weights[modelfile]['pair_wt_zero_' + ortho_name] = float(ortho_md['pair_wt_zero_' + ortho_name])
        #         sample_weights[modelfile]['pair_wt_one_' + ortho_name] = float(ortho_md['pair_wt_one_' + ortho_name])
        #         sample_weights[modelfile]['group_wt_zero_' + ortho_name] = float(ortho_md['group_wt_zero_' + ortho_name])
        #         sample_weights[modelfile]['group_wt_one_' + ortho_name] = float(ortho_md['group_wt_one_' + ortho_name])
        #         metadata.update(ortho_md)

        feat = pd.DataFrame({'importances': dump['attribute_info']['feature_importances_']},
                            index=dump['feature_info']['order'])

        metadata = pd.DataFrame.from_dict(metadata, orient='index')
        modelfile = modelfile.replace('.', '_')
        basepath = os.path.join('training', get_model_subpath(mtype), trg_spec, qry_spec, gid, modelfile)
        with pd.HDFStore(outfile, 'a', complib='blosc', complevel=9) as hdf:
            hdf.put(os.path.join(basepath, 'data'), data, format='fixed')
            hdf.put(os.path.join(basepath, 'features'), feat, format='fixed')
            hdf.put(os.path.join(basepath, 'metadata'), metadata, format='table')
            hdf.flush()
            saved += 1
    return collect_perf, sample_weights, orthos, saved


def calculate_training_weights(dataset, gene_col, ortho_pair, ortho_group, ortho_name):
    """
    :param dataset:
    :param gene_col:
    :param ortho_pair:
    :param ortho_group:
    :param ortho_name:
    :return:
    """
    num_samples = dataset.shape[0]
    num_ones = dataset['true'].sum()
    num_zeros = num_samples - num_ones
    raw_zero_weight = 0.5 / num_zeros
    raw_one_weight = 0.5 / num_ones

    pair_data = dataset.loc[dataset.index.isin(ortho_pair[gene_col]), :].copy()
    num_pairs = pair_data.shape[0]
    pair_ones = pair_data['true'].sum()
    pair_zeros = num_pairs - pair_ones
    pair_zero_weight = 0.5 / pair_zeros
    pair_one_weight = 0.5 / pair_ones

    group_data = dataset.loc[dataset.index.isin(ortho_group[gene_col]), :].copy()
    num_groups = group_data.shape[0]
    group_ones = group_data['true'].sum()
    group_zeros = num_groups - group_ones
    group_zero_weight = 0.5 / group_zeros
    group_one_weight = 0.5 / group_ones
    md = {'num_samples': str(num_samples), 'num_class_one': str(num_ones), 'num_class_zero': str(num_zeros),
          'raw_wt_zero': str(raw_zero_weight), 'raw_wt_one': str(raw_one_weight),
          'num_pair_' + ortho_name: str(num_pairs), 'pair_class_one_' + ortho_name: str(pair_ones),
          'pair_class_zero_' + ortho_name: str(pair_zeros), 'pair_wt_one_' + ortho_name: str(pair_one_weight),
          'pair_wt_zero_' + ortho_name: str(pair_zero_weight),
          'num_group_' + ortho_name: str(num_groups), 'group_class_one_' + ortho_name: str(group_ones),
          'group_class_zero_' + ortho_name: str(group_zeros), 'group_wt_one_' + ortho_name: str(group_one_weight),
          'group_wt_zero_' + ortho_name: str(group_zero_weight)}
    return md


def load_orthologs(fpaths, pair, orthos):
    """
    :param fpaths:
    :param pair:
    :param orthos:
    :return:
    """
    pair_root = '/auto/pairs'
    group_path = '/auto/groups'
    for fp in fpaths:
        orth_name = os.path.basename(fp).split('_')[0]
        if orth_name == 'hcop' and pair[0] != 'human':
            # HCOP orthologs exist only for human
            continue
        if orth_name not in orthos:
            orthos[orth_name] = dict()
        with pd.HDFStore(fp, 'r') as hdf:
            load_path = os.path.join(pair_root, pair[0], pair[1])
            data = hdf[load_path]
            orthos[orth_name][pair] = data
            if 'groups' not in orthos[orth_name]:
                data = hdf[group_path]
                orthos[orth_name]['groups'] = data
    return


def main():
    """
    :return:
    """
    args = parse_command_line()
    # make sure output file is empty when running
    with pd.HDFStore(args.outputfile, 'w'):
        pass
    train_files = []
    test_files = []
    for i in args.traintest:
        if i.endswith('.pck'):
            train_md = i.rsplit('.', 1)[0] + '.json'
            assert os.path.isfile(train_md), 'Invalid path to train metadata file: {}'.format(train_md)
            train_files.append(train_md)
        else:
            test_files.append(i)
        # for debugging only
        # if '_trainmodel_' in i:
        #     train_files.append(i)
        # else:
        #     test_files.append(i)
    train_perf, _, orthos, saved = aggregate_training_data(sorted(train_files),
                                                           args.orthologs, args.outputfile)

    saved = aggregate_testing_data(test_files, train_perf, orthos, args.outputfile)

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
