#!/usr/bin/env python3
"""Regenerate Figure 3's 90 numeric rows from the final 180 seed-level rows."""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev

METRICS=('dice','iou','boundary_f1_tol3','hd95','assd')


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))


def recreate(input_csv,output_csv):
    try: from scipy.stats import t
    except ImportError as e: raise SystemExit('scipy required for Student t 95% CI') from e
    a=read(input_csv)
    if len(a)!=180: raise ValueError('requires 180 seed-level rows')
    groups=defaultdict(dict)
    for r in a:
        k=(r['condition'],r['eval_class'])
        seed=int(r['seed'])
        if seed in groups[k]: raise ValueError('duplicate seed in group')
        groups[k][seed]=r
    if len(groups)!=18: raise ValueError('requires 6 conditions x 3 subgroups')
    result=[]
    for (cond,cls),data in sorted(groups.items()):
        if set(data)!=set(range(10)): raise ValueError('seed set is not 0..9')
        for m in METRICS:
            values=[float(data[s][m]) for s in range(10)]
            mu=mean(values);sd=stdev(values);margin=t.ppf(0.975,9)*sd/math.sqrt(10)
            result.append({'condition':cond,'exposure_percent':data[0]['exposure_percent'],
                           'eval_class':cls,'metric':m,'n_seeds':'10','mean':repr(mu),
                           'sd_across_seeds':repr(float(sd)),'ci95_low_t':repr(float(mu-margin)),'ci95_high_t':repr(float(mu+margin))})
    output_csv.parent.mkdir(parents=True,exist_ok=True)
    if output_csv.exists(): raise FileExistsError(str(output_csv))
    with output_csv.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    return result


def check_expected(rows,expected_csv,tol=1e-10):
    reference=read(expected_csv)
    if len(reference)!=90: raise ValueError('invalid reference Figure 3 source row count')
    get_key=lambda r:(r['condition'],r['eval_class'],r['metric'])
    a={get_key(r):r for r in rows};b={get_key(r):r for r in reference}
    if set(a)!=set(b):raise ValueError('Figure 3 row identity mismatch')
    for k in a:
        for v in ('mean','sd_across_seeds','ci95_low_t','ci95_high_t'):
            if abs(float(a[k][v])-float(b[k][v]))>tol:
                raise ValueError('Figure 3 discrepancy at %s field=%s'%(k,v))


def main():
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--output',type=Path,default=Path('figure3_recalculated.csv'))
    a=p.parse_args();data=a.repo/'release_data'
    r=recreate(data/'139_seed_level_6level_primary5_metrics.csv',a.output)
    check_expected(r,data/'139_exposure_curve_6level_primary5_summary.csv')
    print('PASS: 90 Figure 3 rows reproduced and compared with released source')

if __name__=='__main__':main()
