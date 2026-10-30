#!/usr/bin/env python3
"""
Makes the XMLs of the MAT profiling study from the beast-benchmark MAT XMLs: for each dataset a run with
the local block and one with the path block of the mutation history resampler (the gamma variant with the
local block only), with the chain length and logging set for timing. The data, model, priors and operators
are those of the benchmark XMLs.

usage: make_profile_xmls.py BENCHMARKS_DIR OUT_DIR [--chain-length N] [--path-chain-length N] [--log-every N]
"""
import argparse
import os
import re

SOURCES = {
    'small-2': 'Benchmark-small/benchmark-small-2-mat.xml',
    'big-3': 'Benchmark-big/benchmark-big-3-mat.xml',
    'big-3-gamma': 'Benchmark-big/benchmark-big-3-gamma-mat.xml',
    'big-2': 'Benchmark-big/benchmark-big-2-mat.xml',
    'wide-1': 'Benchmark-wide/benchmark-wide-1-mat.xml',
}
BLOCKS = {'local': ('local', 'path'), 'gamma': ('local',)}


def make(name, xml, block, chain_length, log_every):
    xml = re.sub(r'<mutationHistoryResampler id="matResampler"[^>]*>',
                 '<mutationHistoryResampler id="matResampler" block="%s">' % block, xml, count=1)
    xml = re.sub(r'chainLength="\d+"', 'chainLength="%d"' % chain_length, xml)
    xml = re.sub(r' operatorAnalysis="[^"]*"', '', xml)
    xml = re.sub(r'(<log id="(?:screenLog|fileLog)" logEvery=)"\d+"', r'\1"%d"' % log_every, xml)
    xml = re.sub(r'(<logTree id="treeFileLog" logEvery=)"\d+"', r'\1"%d"' % max(log_every, chain_length // 10), xml)
    xml = re.sub(r'fileName="[\w-]+\.log"', 'fileName="%s.log"' % name, xml)
    xml = re.sub(r'fileName="[\w-]+\.trees"', 'fileName="%s.trees"' % name, xml)
    return xml


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('benchmarks_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--chain-length', type=int, default=1000000)
    parser.add_argument('--path-chain-length', type=int, default=500000)
    parser.add_argument('--log-every', type=int, default=10000)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    for dataset, source in SOURCES.items():
        with open(os.path.join(args.benchmarks_dir, source)) as f:
            xml = f.read()
        for block in BLOCKS['gamma' if dataset.endswith('gamma') else 'local']:
            name = '%s-%s' % (dataset, block)
            length = args.path_chain_length if block == 'path' else args.chain_length
            path = os.path.join(args.out_dir, name + '.xml')
            with open(path, 'w') as f:
                f.write(make(name, xml, block, length, args.log_every))
            print('wrote', path)


if __name__ == '__main__':
    main()
