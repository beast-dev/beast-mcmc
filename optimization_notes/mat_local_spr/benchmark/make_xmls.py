#!/usr/bin/env python3
"""
Makes the XMLs of the big-3 comparison of the Felsenstein likelihood (BEAGLE) with the MAT likelihood
(integrated hits) using the path and the local resampling blocks, with and without the Gibbs operator on
the mutational history (mutationHistoryGibbsOperator, weight GIBBS_WEIGHT, blocks of up to 3 nodes).

The sources are the big-3 XMLs of the beast-benchmark repository (899 BDBV genomes, HKY, strict clock,
constant-size coalescent). All three use the same data, model, priors and operators; they differ only in
the likelihood and, for the MAT, the resampling block. The MAT patterns use constantOnly="true" WITHOUT
ambiguousConstant, so that ambiguous states are treated exactly as BEAGLE treats them: with the ancestral
states summed out the MAT (integrated hits) and Felsenstein posteriors are then the same.

usage: make_xmls.py SOURCE_DIR OUT_DIR [--chain-length N] [--log-every N] [--tree-log-every N]
"""
import argparse
import os
import re

SOURCES = {
    'baseline': 'benchmark-big-3-min.xml',
    'mat-path': 'benchmark-big-3-mat.xml',
    'mat-local': 'benchmark-big-3-mat.xml',
    'mat-path-gibbs': 'benchmark-big-3-mat.xml',
    'mat-local-gibbs': 'benchmark-big-3-mat.xml',
}

GIBBS_WEIGHT = 90


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError('expected exactly one occurrence of: ' + old)
    return text.replace(old, new)


def make(name, source, chain_length, log_every, tree_log_every):
    xml = source
    xml = re.sub(r'chainLength="\d+"', 'chainLength="%d"' % chain_length, xml)
    xml = re.sub(r'(<log id="(?:screenLog|fileLog)" logEvery=)"\d+"', r'\1"%d"' % log_every, xml)
    xml = re.sub(r'(<logTree id="treeFileLog" logEvery=)"\d+"', r'\1"%d"' % tree_log_every, xml)
    xml = re.sub(r'fileName="benchmark-big-3-[\w-]+\.log"', 'fileName="%s.log"' % name, xml)
    xml = re.sub(r'fileName="benchmark-big-3-[\w-]+\.trees"', 'fileName="%s.trees"' % name, xml)
    if name.startswith('mat'):
        # treat ambiguities exactly as BEAGLE does
        xml = replace_once(xml, ' ambiguousConstant="true" ambiguityThreshold="0.25"', '')
        block = 'local' if name.startswith('mat-local') else 'path'
        xml = replace_once(xml, '<mutationHistoryResampler id="matResampler" siteSelection="all">',
                           '<mutationHistoryResampler id="matResampler" block="%s">' % block)
        if name.endswith('-gibbs'):
            xml = replace_once(xml, '\t</operators>',
                               '\t\t<mutationHistoryGibbsOperator weight="%d" size="3">\n'
                               '\t\t\t<mutationHistoryResampler idref="matResampler"/>\n'
                               '\t\t</mutationHistoryGibbsOperator>\n'
                               '\t</operators>' % GIBBS_WEIGHT)
    return xml


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('source_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--chain-length', type=int, default=20000000)
    parser.add_argument('--log-every', type=int, default=5000)
    parser.add_argument('--tree-log-every', type=int, default=100000)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    for name, source in SOURCES.items():
        with open(os.path.join(args.source_dir, source)) as f:
            xml = make(name, f.read(), args.chain_length, args.log_every, args.tree_log_every)
        path = os.path.join(args.out_dir, name + '.xml')
        with open(path, 'w') as f:
            f.write(xml)
        print('wrote', path)


if __name__ == '__main__':
    main()
