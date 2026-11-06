#!/usr/bin/env python3
"""Independently reproduce published RIVER 133/140 paired-seed statistics.

Consumes the same 180 final per-seed TTA rows as Figure 3. Uses the historical
RIVER definition of exact two-sided signed-rank p by enumerating all 2^n signs.
Never trains models or loads patient records. Requires SciPy for t and rankdata.
"""
import argparse
import csv
import itertools
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev

from scipy.stats import rankdata, t

SEEDS = set(range(10))
CLASSES = ('normal', 'TypeA', 'TypeB')
METRICS = ('dice', 'iou', 'boundary_f1_tol3', 'hd95', 'assd')
N_TOL = 1e-10


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as file:
        return list(csv.DictReader(file))


def exact_wilcoxon(differences):
    nonzero = [float(x) for x in differences if math.isfinite(float(x)) and float(x) != 0]
    n = len(nonzero)
    if n == 0:
        return dict(wilcoxon_nonzero_n=0, wilcoxon_w_plus=0.,
                    wilcoxon_w_minus=0., wilcoxon_exact_p=1.)
    ranks = [float(x) for x in rankdata([abs(x) for x in nonzero], method='average')]
    wplus = sum(rank for x, rank in zip(nonzero, ranks) if x > 0)
    total = sum(ranks)
    wminus = total - wplus
    min_observed = min(wplus, wminus)
    # The exact conditional distribution uses the averaged ranks even with ties.
    count = 0
    for signs in itertools.product((0, 1), repeat=n):
        value = sum(rank for rank, sign in zip(ranks, signs) if sign)
        if min(value, total - value) <= min_observed + 1e-12:
            count += 1
    return dict(wilcoxon_nonzero_n=n, wilcoxon_w_plus=wplus,
                wilcoxon_w_minus=wminus, wilcoxon_exact_p=count / (2 ** n))


def holm(pvalues):
    order = sorted(range(len(pvalues)), key=lambda i: (pvalues[i], i))
    result = [None] * len(pvalues)
    running = 0.
    for rank, i in enumerate(order):
        running = max(running, min(1., float(pvalues[i]) * (len(pvalues)-rank)))
        result[i] = running
    return result


def summarize(values):
    if len(values) != 10 or not all(math.isfinite(x) for x in values):
        raise ValueError('Expected 10 finite paired-seed values')
    mu = mean(values)
    sd = stdev(values)
    margin = float(t.ppf(0.975, 9)) * sd / math.sqrt(10)
    return dict(n_seeds=10, delta_mean=mu, delta_sd=sd,
                ci95_low_t=mu-margin, ci95_high_t=mu+margin)


def input_metrics(path):
    rows = read(path)
    if len(rows) != 180:
        raise ValueError('Expected 180 rows, got %d' % len(rows))
    values = {}
    for row in rows:
        seed = int(row['seed'])
        key = (seed, row['condition'], row['eval_class'])
        if seed not in SEEDS or row['eval_class'] not in CLASSES or key in values:
            raise ValueError('Invalid or duplicate seed/condition/class: %s' % (key,))
        if row['mode'] != 'tta':
            raise ValueError('Only final TTA rows supported: %s' % (key,))
        values[key] = {name:float(row[name]) for name in METRICS}
    for condition in ('baseline_2000','robust_05','robust_15','robust_30','robust_50','robust_max643'):
        for cls in CLASSES:
            if {s for s,cond,c in values if cond==condition and c==cls} != SEEDS:
                raise ValueError('Incomplete seed group: %s / %s' % (condition,cls))
    return values


def deltas(metrics, left_cond, left_class, right_cond, right_class, metric):
    return [metrics[(seed,right_cond,right_class)][metric] -
            metrics[(seed,left_cond,left_class)][metric] for seed in range(10)]


def compare_numeric(actual, expected, columns, key):
    for col in columns:
        va, vb = float(actual[col]), float(expected[col])
        if not math.isclose(va,vb,abs_tol=N_TOL,rel_tol=0):
            raise ValueError('%s %s mismatch: %s != %s' % (key,col,va,vb))


def compare_133(metrics, reference_path):
    reference = read(reference_path)
    if len(reference) != 65:
        raise ValueError('Expected 65 primary statistic rows')
    by_family = defaultdict(list)
    for row in reference:
        by_family[row['analysis_family']].append(row)
    if len(by_family) != 13:
        raise ValueError('Expected 13 primary comparison families')
    checked = 0
    for family, rows in by_family.items():
        if len(rows)!=5 or {r['metric'] for r in rows}!=set(METRICS):
            raise ValueError('Invalid primary Holm family '+family)
        pvalues=[]
        for row in rows:
            if row['comparison'] in ('Baseline2000_to_Robust15','Baseline2000_to_Robust50','Robust15_to_Robust50'):
                pair={'Baseline2000_to_Robust15':('baseline_2000','robust_15'),
                      'Baseline2000_to_Robust50':('baseline_2000','robust_50'),
                      'Robust15_to_Robust50':('robust_15','robust_50')}[row['comparison']]
                cls=row['condition_context']
                diff=deltas(metrics,pair[0],cls,pair[1],cls,row['metric'])
            elif row['comparison'] in ('TypeA_minus_normal','TypeB_minus_normal'):
                cond=row['condition_context']; cls=row['comparison'].split('_')[0]
                diff=deltas(metrics,cond,'normal',cond,cls,row['metric'])
            else:
                raise ValueError('Unknown comparison '+row['comparison'])
            s=summarize(diff);w=exact_wilcoxon(diff)
            actual={**s, **w, 'wilcoxon_w_statistic_min':min(w['wilcoxon_w_plus'],w['wilcoxon_w_minus'])}
            compare_numeric(actual,row,('n_seeds','delta_mean','ci95_low_t','ci95_high_t',
                    'wilcoxon_nonzero_n','wilcoxon_w_plus','wilcoxon_w_minus',
                    'wilcoxon_w_statistic_min','wilcoxon_exact_p'),family+'/'+row['metric'])
            pvalues.append(w['wilcoxon_exact_p'])
            checked+=1
        for row,adjusted in zip(rows,holm(pvalues)):
            compare_numeric({'holm_p_within_5_main_metrics':adjusted},row,
                            ('holm_p_within_5_main_metrics',),family+'/'+row['metric'])
            if str(adjusted<0.05).lower()!=row['significant_holm_0p05'].lower():
                raise ValueError('Holm significance flag mismatch '+family+'/'+row['metric'])
    return checked


def compare_140(metrics, reference_path):
    reference=read(reference_path)
    if len(reference)!=15 or {r['eval_class'] for r in reference}!=set(CLASSES):
        raise ValueError('Expected 15 upper-bound statistic rows')
    pvalues=[]
    for row in reference:
        metric=row['metric']; cls=row['eval_class']
        if metric not in METRICS:
            raise ValueError('Unexpected metric '+metric)
        diff=deltas(metrics,'robust_50',cls,'robust_max643',cls,metric)
        s=summarize(diff); w=exact_wilcoxon(diff)
        actual={**s,**w}
        compare_numeric(actual,row,('n_seeds','delta_mean','delta_sd','ci95_low_t','ci95_high_t',
                                    'wilcoxon_nonzero_n','wilcoxon_w_plus','wilcoxon_w_minus',
                                    'wilcoxon_exact_p'),cls+'/'+metric)
        pvalues.append(w['wilcoxon_exact_p'])
    adjusted15=holm(pvalues)
    for cls in CLASSES:
        indexes=[i for i,r in enumerate(reference) if r['eval_class']==cls]
        adjusted5=holm([pvalues[i] for i in indexes])
        for j,i in enumerate(indexes):
            r=reference[i]
            compare_numeric({'holm_p_within_5_metrics':adjusted5[j],
                             'holm_p_across_all_15_tests':adjusted15[i]},r,
                            ('holm_p_within_5_metrics','holm_p_across_all_15_tests'),cls+'/'+r['metric'])
            for flag,score in (('significant_holm5_0p05',adjusted5[j]),
                               ('significant_holm15_0p05',adjusted15[i])):
                if str(score<0.05).lower()!=r[flag].lower():
                    raise ValueError('Significance mismatch '+cls+'/'+r['metric']+'/'+flag)
    return len(reference)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args(); root=args.repo/'release_data'
    metrics=input_metrics(root/'139_seed_level_6level_primary5_metrics.csv')
    n133=compare_133(metrics,root/'133_all_primary_statistical_comparisons.csv')
    n140=compare_140(metrics,root/'140_upper_bound_50_to_64p3_comparison.csv')
    print('PASS: 133 primary statistics (%d rows) and 140 upper bound (%d rows) rederived' % (n133,n140))

if __name__=='__main__': main()
