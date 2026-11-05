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

With --clock ucld, every method has an uncorrelated lognormal relaxed clock (discretizedBranchRates, as BEAUti
writes it: ucld.mean with the CTMC scale prior, ucld.stdev with an exponential prior of mean 1/3, scale
operators on both, swap and uniform-integer operators on the rate categories; the coefficient of variation and
the rate covariance are logged) in place of the strict clock.

usage: make_validation_xmls.py SOURCE_DIR OUT_DIR [--name benchmark-small-2] [--chain-length N]
                               [--log-every N] [--tree-log-every N] [--gibbs-weight W] [--gamma]
                               [--category-weight W] [--clock strict|ucld] [--methods M,...]
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


STRICT_CLOCK = '''	<strictClockBranchRates id="branchRates">
		<rate>
			<parameter id="clock.rate" value="0.001" lower="0.0"/>
		</rate>
	</strictClockBranchRates>'''

UCLD_CLOCK = '''	<discretizedBranchRates id="branchRates">
		<treeModel idref="treeModel"/>
		<distribution>
			<logNormalDistributionModel meanInRealSpace="true">
				<mean>
					<parameter id="ucld.mean" value="0.001" lower="0.0"/>
				</mean>
				<stdev>
					<parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0"/>
				</stdev>
			</logNormalDistributionModel>
		</distribution>
		<rateCategories>
			<parameter id="branchRates.categories"/>
		</rateCategories>
	</discretizedBranchRates>'''

UCLD_STATISTICS = '''	<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
	<rateCovarianceStatistic id="covariance" name="covariance">
		<treeModel idref="treeModel"/>
		<discretizedBranchRates idref="branchRates"/>
	</rateCovarianceStatistic>
'''

UCLD_OPERATORS = '''		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="ucld.mean"/>
		</scaleOperator>
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="ucld.stdev"/>
		</scaleOperator>
		<swapOperator size="1" weight="30" autoOptimize="false">
			<parameter idref="branchRates.categories"/>
		</swapOperator>
		<uniformIntegerOperator weight="30">
			<parameter idref="branchRates.categories"/>
		</uniformIntegerOperator>'''


def add_ucld(xml):
    """The uncorrelated lognormal relaxed clock in place of the strict clock, with its priors, operators and
    statistics; clock.rate becomes ucld.mean wherever it is referred to."""
    xml = replace_once(xml, STRICT_CLOCK, UCLD_CLOCK)
    xml = xml.replace('<strictClockBranchRates idref="branchRates"/>', '<discretizedBranchRates idref="branchRates"/>')
    xml = replace_once(xml, '''		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
''', '''		<discretizedBranchRates idref="branchRates"/>
	</rateStatistic>
''' + UCLD_STATISTICS)
    xml = replace_once(xml, '''		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="clock.rate"/>
		</scaleOperator>''', UCLD_OPERATORS)
    xml = xml.replace('idref="clock.rate"', 'idref="ucld.mean"').replace('label="clock.rate"', 'label="ucld.mean"')
    xml = replace_once(xml, '''				</ctmcScalePrior>''', '''				</ctmcScalePrior>
				<exponentialPrior mean="0.3333333333333333" offset="0.0">
					<parameter idref="ucld.stdev"/>
				</exponentialPrior>''')
    xml = replace_once(xml, '''			<rateStatistic idref="meanRate"/>''', '''			<parameter idref="ucld.stdev"/>
			<rateStatistic idref="meanRate"/>
			<rateStatistic idref="coefficientOfVariation"/>
			<rateCovarianceStatistic idref="covariance"/>''')
    assert 'clock.rate' not in xml
    return xml


def make(method, source, name, chain_length, log_every, tree_log_every, gibbs_weight, gamma=False,
         category_weight=None, clock='strict'):
    xml = re.sub(r'chainLength="\d+"', 'chainLength="%d"' % chain_length, source)
    if clock == 'ucld':
        xml = add_ucld(xml)
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
    parser.add_argument('--clock', choices=['strict', 'ucld'], default='strict',
                        help='the strict clock of the source, or an uncorrelated lognormal relaxed clock')
    parser.add_argument('--methods', default='baseline,mat-local,mat-local-gibbs')
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sources = {'baseline': args.name + '-min.xml', 'mat-local': args.name + '-mat.xml',
               'mat-local-gibbs': args.name + '-mat.xml',
               # the same as mat-local-gibbs, named for runs made since the resampler moved into the likelihood
               'mat-lazy-gibbs': args.name + '-mat.xml'}
    for method, source in sources.items():
        if method not in args.methods.split(','):
            continue
        with open(os.path.join(args.source_dir, source)) as f:
            text = f.read()
        category_weight = args.category_weight
        if category_weight is None:
            category_weight = round(0.01 * float(re.search(r'<subtreeLeap [^>]*weight="([^"]+)"', text).group(1)), 3)
        xml = make(method, text, args.name, args.chain_length, args.log_every, args.tree_log_every,
                   args.gibbs_weight, args.gamma, category_weight, args.clock)
        path = os.path.join(args.out_dir, method + '.xml')
        with open(path, 'w') as f:
            f.write(xml)
        print('wrote', path)


if __name__ == '__main__':
    main()
