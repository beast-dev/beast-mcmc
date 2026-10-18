#!/usr/bin/env python3
"""
Makes the XMLs of a validation of the MAT local block, with and without the Gibbs operator on the mutational
history, against the Felsenstein likelihood (BEAGLE) on a dataset small enough for the tree parameters to
mix well (by default small-2 of the beast-benchmark repository: 165 influenza A HA sequences, JC69, strict
clock, constant-size coalescent), so that a small bias in either MAT method would show.

The MAT patterns of the source use constantOnly="true" without ambiguousConstant, so the MAT (integrated
hits) and Felsenstein posteriors are the same.

usage: make_validation_xmls.py SOURCE_DIR OUT_DIR [--name benchmark-small-2] [--chain-length N]
                               [--log-every N] [--tree-log-every N] [--gibbs-weight W]
"""
import argparse
import os
import re

from make_xmls import replace_once


def make(method, source, name, chain_length, log_every, tree_log_every, gibbs_weight):
    xml = re.sub(r'chainLength="\d+"', 'chainLength="%d"' % chain_length, source)
    xml = re.sub(r'(<log id="(?:screenLog|fileLog)" logEvery=)"\d+"', r'\1"%d"' % log_every, xml)
    xml = re.sub(r'(<logTree id="treeFileLog" logEvery=)"\d+"', r'\1"%d"' % tree_log_every, xml)
    xml = re.sub(r'fileName="%s(?:-mat)?\.log"' % name, 'fileName="%s.log"' % method, xml)
    xml = re.sub(r'fileName="%s(?:-mat)?\.trees"' % name, 'fileName="%s.trees"' % method, xml)
    if method.startswith('mat'):
        xml = replace_once(xml, '<mutationHistoryResampler id="matResampler">',
                           '<mutationHistoryResampler id="matResampler" block="local">')
        if method.endswith('-gibbs'):
            xml = replace_once(xml, '\t</operators>',
                               '\t\t<mutationHistoryGibbsOperator weight="%g" size="3">\n'
                               '\t\t\t<mutationHistoryResampler idref="matResampler"/>\n'
                               '\t\t</mutationHistoryGibbsOperator>\n'
                               '\t</operators>' % gibbs_weight)
    return xml


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('source_dir')
    parser.add_argument('out_dir')
    parser.add_argument('--name', default='benchmark-small-2')
    parser.add_argument('--chain-length', type=int, default=20000000)
    parser.add_argument('--log-every', type=int, default=10000)
    parser.add_argument('--tree-log-every', type=int, default=50000)
    parser.add_argument('--gibbs-weight', type=float, default=16.5)
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sources = {'baseline': args.name + '-min.xml', 'mat-local': args.name + '-mat.xml',
               'mat-local-gibbs': args.name + '-mat.xml'}
    for method, source in sources.items():
        with open(os.path.join(args.source_dir, source)) as f:
            xml = make(method, f.read(), args.name, args.chain_length, args.log_every, args.tree_log_every,
                       args.gibbs_weight)
        path = os.path.join(args.out_dir, method + '.xml')
        with open(path, 'w') as f:
            f.write(xml)
        print('wrote', path)


if __name__ == '__main__':
    main()
