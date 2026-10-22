#!/usr/bin/env python3
"""
Makes the XMLs of a validation of the MAT local block, with and without the Gibbs operator on the mutational
history, against the Felsenstein likelihood (BEAGLE) on a dataset small enough for the tree parameters to
mix well (by default small-2 of the beast-benchmark repository: 165 influenza A HA sequences, JC69, strict
clock, constant-size coalescent), so that a small bias in either MAT method would show.

The MAT patterns of the source use constantOnly="true" without ambiguousConstant, so the MAT (integrated
hits) and Felsenstein posteriors are the same.

With --gamma, every method also has a discrete gamma site model (4 categories, shape alpha with an
exponential prior of mean 0.5 and a scale operator); the MAT then samples each variable site's category,
with a siteRateCategoryGibbsOperator (weight --category-weight, by default 1% of the SubtreeLeap weight: each
update draws every site's category, which mix quickly), and sums the constant sites over them.

usage: make_validation_xmls.py SOURCE_DIR OUT_DIR [--name benchmark-small-2] [--chain-length N]
                               [--log-every N] [--tree-log-every N] [--gibbs-weight W] [--gamma]
                               [--category-weight W] [--methods M,...]
"""
import argparse
import os
import re

from make_xmls import replace_once


GAMMA_SITE_MODEL = '''		<gammaShape gammaCategories="4">
			<parameter id="alpha" value="0.5" lower="0.0"/>
		</gammaShape>
	</siteModel>'''


def add_gamma(xml):
    """A discrete gamma site model, with a prior and an operator on its shape, logged."""
    xml = replace_once(xml, '''		</substitutionModel>
	</siteModel>''', '''		</substitutionModel>
''' + GAMMA_SITE_MODEL)
    xml = replace_once(xml, '''				<coalescentLikelihood idref="coalescent"/>''', '''				<exponentialPrior mean="0.5" offset="0.0">
					<parameter idref="alpha"/>
				</exponentialPrior>
				<coalescentLikelihood idref="coalescent"/>''')
    xml = replace_once(xml, '''	</operators>''', '''		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="alpha"/>
		</scaleOperator>
	</operators>''')
    xml = replace_once(xml, '''			<parameter idref="clock.rate"/>
			<rateStatistic idref="meanRate"/>''', '''			<parameter idref="clock.rate"/>
			<parameter idref="alpha"/>
			<rateStatistic idref="meanRate"/>''')
    return xml


def make(method, source, name, chain_length, log_every, tree_log_every, gibbs_weight, gamma=False,
         category_weight=None):
    xml = re.sub(r'chainLength="\d+"', 'chainLength="%d"' % chain_length, source)
    xml = re.sub(r'(<log id="(?:screenLog|fileLog)" logEvery=)"\d+"', r'\1"%d"' % log_every, xml)
    xml = re.sub(r'(<logTree id="treeFileLog" logEvery=)"\d+"', r'\1"%d"' % tree_log_every, xml)
    xml = re.sub(r'fileName="%s(?:-mat)?\.log"' % name, 'fileName="%s.log"' % method, xml)
    xml = re.sub(r'fileName="%s(?:-mat)?\.trees"' % name, 'fileName="%s.trees"' % method, xml)
    if method.startswith('mat'):
        xml = replace_once(xml, '<mutationHistoryResampler id="matResampler">',
                           '<mutationHistoryResampler id="matResampler" block="local">')
        if method.endswith('-gibbs') and 'mutationHistoryGibbsOperator' not in xml:
            xml = replace_once(xml, '\t</operators>',
                               '\t\t<mutationHistoryGibbsOperator weight="%g" size="3">\n'
                               '\t\t\t<mutationHistoryResampler idref="matResampler"/>\n'
                               '\t\t</mutationHistoryGibbsOperator>\n'
                               '\t</operators>' % gibbs_weight)
        elif not method.endswith('-gibbs'):
            # (the source may already have the history Gibbs operator)
            xml = re.sub(r'\t\t<mutationHistoryGibbsOperator[^>]*>\n[^\n]*\n\t\t</mutationHistoryGibbsOperator>\n', '', xml)
            assert 'mutationHistoryGibbsOperator' not in xml
    if gamma:
        xml = add_gamma(xml)
        if method.startswith('mat'):
            xml = replace_once(xml, '\t</operators>',
                               '\t\t<siteRateCategoryGibbsOperator weight="%g">\n'
                               '\t\t\t<matDataLikelihood idref="treeLikelihood"/>\n'
                               '\t\t</siteRateCategoryGibbsOperator>\n'
                               '\t</operators>' % category_weight)
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
    parser.add_argument('--gamma', action='store_true')
    parser.add_argument('--category-weight', type=float, default=None,
                        help='weight of the siteRateCategoryGibbsOperator (default 1%% of the SubtreeLeap weight)')
    parser.add_argument('--methods', default='baseline,mat-local,mat-local-gibbs')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sources = {'baseline': args.name + '-min.xml', 'mat-local': args.name + '-mat.xml',
               'mat-local-gibbs': args.name + '-mat.xml'}
    for method, source in sources.items():
        if method not in args.methods.split(','):
            continue
        with open(os.path.join(args.source_dir, source)) as f:
            text = f.read()
        category_weight = args.category_weight
        if category_weight is None:
            category_weight = round(0.01 * float(re.search(r'<subtreeLeap [^>]*weight="([^"]+)"', text).group(1)), 3)
        xml = make(method, text, args.name, args.chain_length, args.log_every, args.tree_log_every,
                   args.gibbs_weight, args.gamma, category_weight)
        path = os.path.join(args.out_dir, method + '.xml')
        with open(path, 'w') as f:
            f.write(xml)
        print('wrote', path)


if __name__ == '__main__':
    main()
