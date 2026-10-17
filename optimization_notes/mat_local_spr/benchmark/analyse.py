#!/usr/bin/env python3
"""
Summarises the big-3 comparison runs (see run_benchmark.sh): for each run, method and parameter of the
basket, the posterior mean, 95% HPD interval, effective sample size (ESS) and ESS per minute, and whether
the replicates of each method, and each method and the baseline, agree.

ESS is computed as Tracer does (autocovariances summed in pairs until a pair is negative). ESS per minute
is per minute of sampling after the burn-in (the run's wall time times the proportion of states kept).
Agreement is the difference between means in units of its Monte Carlo standard error
(sqrt(se1^2 + se2^2), se = sd / sqrt(ESS)): |z| > 3 is flagged.

usage: analyse.py RUN_DIR [--burnin 0.25] [--baseline baseline] [--markdown out.md]
"""
import argparse
import glob
import math
import os

import numpy as np

BASKET = ['rootHeight', 'treeLength', 'constant.popSize', 'kappa', 'clock.rate',
          'frequencies1', 'frequencies2', 'frequencies3', 'frequencies4', 'coalescent', 'prior']


def read_log(path):
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if header is None:
                header = fields
            else:
                rows.append([float(x) for x in fields])
    data = np.array(rows)
    return {name: data[:, i] for i, name in enumerate(header)}


def ess(x):
    """Effective sample size as in Tracer (initial positive sequence of paired autocovariances)."""
    n = len(x)
    if n < 4 or np.var(x) == 0.0:
        return float('nan')
    y = x - np.mean(x)
    size = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(y, size)
    gamma = np.fft.irfft(f * np.conj(f), size)[:n] / n
    variance = gamma[0]
    lag = 1
    while lag + 1 < n:
        pair = gamma[lag] + gamma[lag + 1]
        if pair <= 0.0:
            break
        variance += 2.0 * pair
        lag += 2
    act = variance / gamma[0]
    return n / act


def hpd(x, mass=0.95):
    s = np.sort(x)
    k = int(math.floor(mass * len(s)))
    widths = s[k:] - s[:len(s) - k]
    i = int(np.argmin(widths))
    return s[i], s[i + k]


def summarise(run_dir, burnin):
    runs = {}
    for path in sorted(glob.glob(os.path.join(run_dir, '*', 'time.txt'))):
        directory = os.path.dirname(path)
        name = os.path.basename(directory)
        method, seed = name.rsplit('-', 1)
        logs = [p for p in glob.glob(os.path.join(directory, '*.log'))]
        if not logs:
            continue
        log = read_log(logs[0])
        with open(path) as f:
            seconds = float(f.read().strip())
        states = log['state']
        keep = states >= burnin * states[-1]
        minutes = seconds / 60.0 * (1.0 - burnin)
        stats = {}
        for p in BASKET:
            if p not in log:
                continue
            x = log[p][keep]
            e = ess(x)
            stats[p] = {'mean': float(np.mean(x)), 'sd': float(np.std(x, ddof=1)), 'ess': e,
                        'ess_per_minute': e / minutes, 'hpd': hpd(x),
                        # first and second halves of the kept samples, for stationarity
                        'halves': (float(np.mean(x[:len(x) // 2])), float(np.mean(x[len(x) // 2:])))}
        runs[(method, seed)] = {'seconds': seconds, 'states': int(states[-1]), 'samples': int(keep.sum()),
                                'stats': stats}
    return runs


def z_score(a, b):
    se = math.sqrt(a['sd'] ** 2 / a['ess'] + b['sd'] ** 2 / b['ess'])
    return (a['mean'] - b['mean']) / se if se > 0 else float('nan')


def pooled(stats_list):
    """Pools replicates: mean weighted by ESS, and the summed ESS."""
    ess_total = sum(s['ess'] for s in stats_list)
    mean = sum(s['mean'] * s['ess'] for s in stats_list) / ess_total
    sd = math.sqrt(sum(s['sd'] ** 2 * s['ess'] for s in stats_list) / ess_total)
    return {'mean': mean, 'sd': sd, 'ess': ess_total}


def fmt(x):
    if x != x:
        return 'nan'
    if abs(x) >= 1000 or (abs(x) < 0.01 and x != 0):
        return '%.4g' % x
    return '%.4f' % x


def report(runs, baseline, burnin):
    lines = []
    out = lines.append
    methods = sorted({m for m, _ in runs}, key=lambda m: (m != baseline, m))
    out('Burn-in: %d%% of each run.\n' % round(100 * burnin))
    out('| Run | States | Wall time (min) | Time per million states (min) | Samples kept |')
    out('|---|---|---|---|---|')
    for (m, s), r in sorted(runs.items(), key=lambda kv: (methods.index(kv[0][0]), kv[0][1])):
        out('| %s %s | %d | %.1f | %.2f | %d |' % (m, s, r['states'], r['seconds'] / 60,
                                                   r['seconds'] / 60 / (r['states'] / 1e6), r['samples']))
    out('')

    out('## Posterior means (95% HPD) by run\n')
    out('| Parameter | ' + ' | '.join('%s %s' % k for k in sorted(runs, key=lambda k: (methods.index(k[0]), k[1]))) + ' |')
    out('|---|' + '---|' * len(runs))
    keys = sorted(runs, key=lambda k: (methods.index(k[0]), k[1]))
    for p in BASKET:
        cells = []
        for k in keys:
            st = runs[k]['stats'].get(p)
            cells.append('%s (%s, %s)' % (fmt(st['mean']), fmt(st['hpd'][0]), fmt(st['hpd'][1])) if st else '')
        out('| %s | %s |' % (p, ' | '.join(cells)))
    out('')

    out('## ESS (ESS per minute) by run\n')
    out('| Parameter | ' + ' | '.join('%s %s' % k for k in keys) + ' |')
    out('|---|' + '---|' * len(keys))
    for p in BASKET:
        cells = []
        for k in keys:
            st = runs[k]['stats'].get(p)
            cells.append('%.0f (%.2f)' % (st['ess'], st['ess_per_minute']) if st else '')
        out('| %s | %s |' % (p, ' | '.join(cells)))
    out('')

    out('## ESS per minute by method (mean over replicates; minimum over the basket)\n')
    out('| Method | ' + ' | '.join(BASKET) + ' | minimum |')
    out('|---|' + '---|' * (len(BASKET) + 1))
    for m in methods:
        reps = [runs[k] for k in keys if k[0] == m]
        cells = []
        values = []
        for p in BASKET:
            v = [r['stats'][p]['ess_per_minute'] for r in reps if p in r['stats']]
            value = sum(v) / len(v) if v else float('nan')
            values.append(value)
            cells.append('%.2f' % value)
        out('| %s | %s | %.2f |' % (m, ' | '.join(cells), min(values)))
    out('')

    out('## Agreement (z = difference / Monte Carlo SE; |z| > 3 flagged)\n')
    out('Within runs: second half of the kept samples against the first (stationarity). Between replicates of a '
        'method. Between each method (replicates pooled) and %s.\n' % baseline)
    out('| Parameter | ' + ' | '.join(['halves %s %s' % k for k in keys] +
                                       ['reps %s' % m for m in methods] +
                                       ['%s vs %s' % (m, baseline) for m in methods if m != baseline]) + ' |')
    ncols = len(keys) + len(methods) + len(methods) - 1
    out('|---|' + '---|' * ncols)
    for p in BASKET:
        cells = []
        for k in keys:
            st = runs[k]['stats'].get(p)
            if st:
                half = {'mean': st['halves'][1], 'sd': st['sd'], 'ess': st['ess'] / 2}
                first = {'mean': st['halves'][0], 'sd': st['sd'], 'ess': st['ess'] / 2}
                cells.append(flag(z_score(half, first)))
            else:
                cells.append('')
        for m in methods:
            reps = [runs[k]['stats'][p] for k in keys if k[0] == m and p in runs[k]['stats']]
            cells.append(flag(z_score(reps[0], reps[1])) if len(reps) >= 2 else '')
        base = [runs[k]['stats'][p] for k in keys if k[0] == baseline and p in runs[k]['stats']]
        for m in methods:
            if m == baseline:
                continue
            other = [runs[k]['stats'][p] for k in keys if k[0] == m and p in runs[k]['stats']]
            cells.append(flag(z_score(pooled(other), pooled(base))) if other and base else '')
        out('| %s | %s |' % (p, ' | '.join(cells)))
    return '\n'.join(lines)


def flag(z):
    if z != z:
        return 'nan'
    return ('**%.1f**' if abs(z) > 3 else '%.1f') % z


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('run_dir')
    parser.add_argument('--burnin', type=float, default=0.25)
    parser.add_argument('--baseline', default='baseline')
    parser.add_argument('--markdown')
    args = parser.parse_args()
    runs = summarise(args.run_dir, args.burnin)
    text = report(runs, args.baseline, args.burnin)
    print(text)
    if args.markdown:
        with open(args.markdown, 'w') as f:
            f.write(text + '\n')


if __name__ == '__main__':
    main()
