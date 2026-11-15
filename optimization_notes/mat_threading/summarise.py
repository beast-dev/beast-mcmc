#!/usr/bin/env python3
"""
Summarises the profiling study's runs as markdown tables: per run the chain's speed, the step profile of
the subtreeLeap operator (and of all operators), the MAT phases per subtreeLeap step and per evaluation,
and the counters per move; plus the overhead of the profiles (profiled against unprofiled runs).

usage: summarise.py RUNS_DIR [--markdown OUT.md]
(RUNS_DIR holds profiled/<name>/out.txt, and optionally unprofiled/<name>/out.txt and jfr/<name>/out.txt)
"""
import argparse
import os
import re
import sys


def block_after(lines, prefix, skip):
    """@return the lines after the first line starting with prefix (skipping `skip` header lines) up to a blank line"""
    for i, line in enumerate(lines):
        if line.startswith(prefix):
            block = []
            for l in lines[i + 1 + skip:]:
                if not l.strip():
                    break
                block.append(l)
            return block
    return []


def parse_out(path):
    text = open(path).read()
    lines = text.splitlines()
    run = {'name': os.path.basename(os.path.dirname(path))}
    m = re.search(r'MAT likelihood: (\d+) variable sites, (\d+) constant sites', text)
    if m:
        run['variable'] = int(m.group(1))
        run['constant'] = int(m.group(2))
    rates = re.findall(r'([0-9.]+) minutes/million states', text)
    run['min_per_million'] = float(rates[-1]) if rates else None
    m = re.search(r'^([0-9.]+) seconds$', text, re.M)
    run['seconds'] = float(m.group(1)) if m else None
    ops = {}
    for line in block_after(lines, 'Operator analysis', 1):
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            ops[parts[0]] = {'count': int(parts[-4]), 'time_s': float(parts[-3]), 'us_per_op': float(parts[-2]),
                             'accept': float(parts[-1])}
        except ValueError:
            pass
    run['operators'] = ops
    steps = {}
    for line in block_after(lines, 'Step profile', 1):
        parts = line.split()
        if len(parts) < 8:
            continue
        try:
            steps[' '.join(parts[:-7])] = {'steps': int(parts[-7]), 'store': float(parts[-6]), 'propose': float(parts[-5]),
                                           'evaluate': float(parts[-4]), 'acceptRestore': float(parts[-3]),
                                           'adapt': float(parts[-2]), 'total': float(parts[-1])}
        except ValueError:
            pass
    run['steps'] = steps
    phases = {}
    for line in block_after(lines, 'MAT profile', 1):
        parts = line.split()
        try:
            phases[parts[0]] = {'calls': int(parts[1]), 'total_s': float(parts[2]), 'us_per_call': float(parts[3])}
        except (ValueError, IndexError):
            pass
    counters = {}
    for line in block_after(lines, 'Counter ', 0):
        parts = line.split()
        try:
            counters[parts[0]] = int(parts[1])
        except (ValueError, IndexError):
            pass
    run['phases'] = phases
    run['counters'] = counters
    # per-operator MAT profiles
    per_op = {}
    for i, line in enumerate(lines):
        m = re.match(r'MAT profile of (.+) \((\d+) steps;', line)
        if not m:
            continue
        op, steps_n = m.group(1), int(m.group(2))
        op_phases, op_counters = {}, {}
        in_counters = False
        for l in lines[i + 2:]:
            if not l.strip():
                break
            parts = l.split()
            if l.startswith('Counter'):
                in_counters = True
                continue
            try:
                if in_counters:
                    op_counters[parts[0]] = (int(parts[1]), float(parts[2]))
                else:
                    op_phases[parts[0]] = {'calls': int(parts[1]), 'total_s': float(parts[2]),
                                           'us_per_call': float(parts[3]), 'us_per_step': float(parts[4])}
            except (ValueError, IndexError):
                pass
        per_op[op] = {'steps': steps_n, 'phases': op_phases, 'counters': op_counters}
    run['per_operator'] = per_op
    return run


def load(runs_dir, group):
    runs = []
    d = os.path.join(runs_dir, group)
    if not os.path.isdir(d):
        return runs
    for name in sorted(os.listdir(d)):
        path = os.path.join(d, name, 'out.txt')
        if os.path.exists(path) and os.path.exists(os.path.join(d, name, 'time.txt')):
            runs.append(parse_out(path))
    return runs


def stl(run):
    for name, s in run['steps'].items():
        if name.startswith('subtreeLeap'):
            return s
    return None


def table(header, rows):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '|'.join(['---'] * len(header)) + '|']
    for row in rows:
        out.append('| ' + ' | '.join(str(x) for x in row) + ' |')
    return '\n'.join(out) + '\n'


def fmt(x, digits=2):
    return '' if x is None else ('%.*f' % (digits, x))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('runs_dir')
    parser.add_argument('--markdown', default=None)
    args = parser.parse_args()
    profiled = load(args.runs_dir, 'profiled')
    unprofiled = {r['name']: r for r in load(args.runs_dir, 'unprofiled')}
    out = []

    out.append('## Runs\n')
    rows = []
    for r in profiled:
        s = stl(r)
        rows.append([r['name'], r.get('variable', ''), fmt(r['min_per_million']), fmt(r['seconds'], 0),
                     fmt(s['total']) if s else '', fmt(s['propose']) if s else '', fmt(s['evaluate']) if s else '',
                     fmt(s['store'] + s['acceptRestore'] + s['adapt']) if s else ''])
    out.append(table(['Run', 'Variable sites', 'min/M states', 'Wall (s)', 'STL step (us)', 'propose', 'evaluate',
                      'store+accept/restore+adapt'], rows))

    out.append('## Step profile of every operator (microseconds per step)\n')
    for r in profiled:
        out.append('### %s\n' % r['name'])
        rows = []
        for name, s in r['steps'].items():
            rows.append([name, s['steps'], fmt(s['store']), fmt(s['propose']), fmt(s['evaluate']),
                         fmt(s['acceptRestore']), fmt(s['adapt']), fmt(s['total'])])
        out.append(table(['Operator', 'Steps', 'store', 'propose', 'evaluate', 'accept/restore', 'adapt', 'total'], rows))

    out.append('## MAT phases\n')
    out.append('Microseconds per call, and per subtree move (local.before + local.after, or path.*) or per '
               'evaluation (lik.*), from the whole run.\n')
    for r in profiled:
        out.append('### %s\n' % r['name'])
        rows = []
        moves = r['counters'].get('local_moves', 0) + r['counters'].get('path_moves', 0)
        evaluations = r['counters'].get('evaluations', 0)
        for name, p in r['phases'].items():
            per = ''
            if name.startswith('lik.') and evaluations:
                per = fmt(p['total_s'] * 1e6 / evaluations)
            elif moves:
                per = fmt(p['total_s'] * 1e6 / moves)
            rows.append([name, p['calls'], fmt(p['total_s'], 3), fmt(p['us_per_call']), per])
        out.append(table(['Phase', 'Calls', 'Total (s)', 'us/call', 'us per move or evaluation'], rows))
        rows = []
        for name, c in r['counters'].items():
            per = fmt(c / moves) if moves and name not in ('local_moves', 'path_moves', 'gibbs_updates', 'evaluations') else ''
            rows.append([name, c, per])
        out.append(table(['Counter', 'Total', 'per move'], rows))

    out.append('## MAT phases per step of each tree operator (microseconds per step)\n')
    for r in profiled:
        if not r.get('per_operator'):
            continue
        out.append('### %s\n' % r['name'])
        ops = [op for op in r['per_operator'] if op.startswith(('subtreeLeap', 'subtreeJump', 'fixedHeight', 'mutationHistoryGibbs'))]
        names = []
        for op in ops:
            for ph in r['per_operator'][op]['phases']:
                if ph not in names:
                    names.append(ph)
        rows = []
        for ph in names:
            rows.append([ph] + [fmt(r['per_operator'][op]['phases'].get(ph, {}).get('us_per_step')) for op in ops])
        header = ['Phase'] + ['%s (%d steps)' % (op.split('(')[0], r['per_operator'][op]['steps']) for op in ops]
        out.append(table(header, rows))
        cnames = []
        for op in ops:
            for c in r['per_operator'][op]['counters']:
                if c not in cnames:
                    cnames.append(c)
        rows = [[c] + [fmt(r['per_operator'][op]['counters'].get(c, (0, None))[1]) for op in ops] for c in cnames]
        out.append(table(['Counter per step'] + [op.split('(')[0] for op in ops], rows))

    if unprofiled:
        out.append('## Overhead of the profiles\n')
        rows = []
        for r in profiled:
            u = unprofiled.get(r['name'])
            if u is None:
                continue
            s = stl(r)
            su = u['operators'].get('subtreeLeap(treeModel)')
            rows.append([r['name'], fmt(r['min_per_million']), fmt(u['min_per_million']),
                         fmt(100.0 * (r['min_per_million'] / u['min_per_million'] - 1.0), 1)
                         if r['min_per_million'] and u['min_per_million'] else '',
                         fmt(s['total']) if s else '', fmt(su['us_per_op']) if su else ''])
        out.append(table(['Run', 'min/M (profiled)', 'min/M (unprofiled)', 'Overhead %', 'STL us (profiled)',
                          'STL us (unprofiled)'], rows))

    text = '\n'.join(out)
    if args.markdown:
        open(args.markdown, 'w').write(text)
        print('wrote', args.markdown)
    else:
        print(text)


if __name__ == '__main__':
    main()
