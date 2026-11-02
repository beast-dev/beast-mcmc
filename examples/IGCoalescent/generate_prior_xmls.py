#!/usr/bin/env python3
"""
Generates paired "full" (sample Ne with a ScaleOperator under an inverse-Gamma
prior) vs "marginalised" (integrate Ne out analytically via coalescentLikelihoodIG)
prior-only BEAST XMLs, for the scenarios described in PR #1075:

  - 5 contemporaneous taxa
  - 50 contemporaneous taxa
  - 17 heterochronous taxa (the Dengue-4 tip dates from
    examples/TestXML/TreePriors/testConstantSize.xml)

Alpha/beta are fixed at 2.0/10.0 throughout (E[Ne] = 10, as in the PR).
No sequence data: the MCMC target is the coalescent prior alone (plus the
population-size prior in the "full" variant).
"""

import os

ALPHA = 2.0
BETA = 10.0
CHAIN_LENGTH = 10_000_000
LOG_EVERY = 1000

HETEROCHRONOUS_TAXA = [
    ("D4Brazi82", "1982.0"), ("D4ElSal83", "1983.0"), ("D4ElSal94", "1994.0"),
    ("D4Indon76", "1976.0"), ("D4Indon77", "1977.0"), ("D4Mexico84", "1984.0"),
    ("D4NewCal81", "1981.0"), ("D4Philip64", "1964.0"), ("D4Philip56", "1956.0"),
    ("D4Philip84", "1984.0"), ("D4PRico86", "1986.0"), ("D4SLanka78", "1978.0"),
    ("D4Tahiti79", "1979.0"), ("D4Tahiti85", "1985.0"), ("D4Thai63", "1963.0"),
    ("D4Thai78", "1978.0"), ("D4Thai84", "1984.0"),
]


def taxa_block(n, heterochronous):
    lines = ['\t<taxa id="taxa">']
    if heterochronous:
        for taxon_id, date in HETEROCHRONOUS_TAXA:
            lines.append(f'\t\t<taxon id="{taxon_id}">')
            lines.append(f'\t\t\t<date value="{date}" direction="forwards" units="years"/>')
            lines.append('\t\t</taxon>')
    else:
        for i in range(1, n + 1):
            lines.append(f'\t\t<taxon id="taxon{i}"/>')
    lines.append('\t</taxa>')
    return '\n'.join(lines)


def shared_prefix(n, heterochronous):
    return f"""<?xml version="1.0" standalone="yes"?>

<!-- Prior-only run for PR #1075 (integrated inverse-Gamma coalescent prior). -->
<!-- {n} {'heterochronous' if heterochronous else 'contemporaneous'} taxa. -->
<!-- No sequence data: the MCMC samples the coalescent (and, in the "full"  -->
<!-- variant, the population size) prior alone.                             -->
<beast>

{taxa_block(n, heterochronous)}

\t<constantSize id="demo" units="years">
\t\t<populationSize>
\t\t\t<parameter id="demo.popSize" value="{BETA}" lower="0.0"/>
\t\t</populationSize>
\t</constantSize>

\t<coalescentSimulator id="startingTree">
\t\t<taxa idref="taxa"/>
\t\t<constantSize idref="demo"/>
\t</coalescentSimulator>

\t<treeModel id="treeModel">
\t\t<coalescentTree idref="startingTree"/>
\t\t<rootHeight>
\t\t\t<parameter id="treeModel.rootHeight"/>
\t\t</rootHeight>
\t\t<nodeHeights internalNodes="true">
\t\t\t<parameter id="treeModel.internalNodeHeights"/>
\t\t</nodeHeights>
\t\t<nodeHeights internalNodes="true" rootNode="true">
\t\t\t<parameter id="treeModel.allInternalNodeHeights"/>
\t\t</nodeHeights>
\t</treeModel>

\t<treeLengthStatistic id="treeLength">
\t\t<treeModel idref="treeModel"/>
\t</treeLengthStatistic>

\t<tmrcaStatistic id="age(root)" absolute="true">
\t\t<treeModel idref="treeModel"/>
\t</tmrcaStatistic>
"""


TREE_OPERATORS = """\t\t<subtreeSlide weight="15" gaussian="true" size="1.0">
\t\t\t<treeModel idref="treeModel"/>
\t\t</subtreeSlide>
\t\t<narrowExchange weight="15">
\t\t\t<treeModel idref="treeModel"/>
\t\t</narrowExchange>
\t\t<wideExchange weight="3">
\t\t\t<treeModel idref="treeModel"/>
\t\t</wideExchange>
\t\t<wilsonBalding weight="3">
\t\t\t<treeModel idref="treeModel"/>
\t\t</wilsonBalding>
\t\t<uniformOperator weight="30">
\t\t\t<parameter idref="treeModel.internalNodeHeights"/>
\t\t</uniformOperator>
\t\t<scaleOperator scaleFactor="0.75" weight="3">
\t\t\t<parameter idref="treeModel.rootHeight"/>
\t\t</scaleOperator>"""


def full_xml(n, heterochronous, filename_stem):
    return shared_prefix(n, heterochronous) + f"""
\t<coalescentLikelihood id="coalescent">
\t\t<model>
\t\t\t<constantSize idref="demo"/>
\t\t</model>
\t\t<populationTree>
\t\t\t<treeModel idref="treeModel"/>
\t\t</populationTree>
\t</coalescentLikelihood>

\t<mcmc id="mcmc" chainLength="{CHAIN_LENGTH}" autoOptimize="true">
\t\t<prior id="prior">
\t\t\t<invgammaPrior shape="{ALPHA}" scale="{BETA}">
\t\t\t\t<parameter idref="demo.popSize"/>
\t\t\t</invgammaPrior>
\t\t\t<coalescentLikelihood idref="coalescent"/>
\t\t</prior>

\t\t<log id="screenLog" logEvery="{LOG_EVERY * 10}">
\t\t\t<column label="prior" dp="4" width="12">
\t\t\t\t<prior idref="prior"/>
\t\t\t</column>
\t\t\t<column label="age(root)" sf="6" width="12">
\t\t\t\t<tmrcaStatistic idref="age(root)"/>
\t\t\t</column>
\t\t\t<column label="popSize" sf="6" width="12">
\t\t\t\t<parameter idref="demo.popSize"/>
\t\t\t</column>
\t\t</log>

\t\t<log id="fileLog" logEvery="{LOG_EVERY}" fileName="{filename_stem}.log">
\t\t\t<prior idref="prior"/>
\t\t\t<parameter idref="treeModel.rootHeight"/>
\t\t\t<tmrcaStatistic idref="age(root)"/>
\t\t\t<treeLengthStatistic idref="treeLength"/>
\t\t\t<parameter idref="demo.popSize"/>
\t\t</log>

\t\t<operators id="operators" optimizationSchedule="default">
\t\t\t<scaleOperator scaleFactor="0.75" weight="3">
\t\t\t\t<parameter idref="demo.popSize"/>
\t\t\t</scaleOperator>
{TREE_OPERATORS}
\t\t</operators>
\t</mcmc>

\t<report>
\t\t<property name="timer">
\t\t\t<mcmc idref="mcmc"/>
\t\t</property>
\t</report>

</beast>
"""


def marginalised_xml(n, heterochronous, filename_stem):
    return shared_prefix(n, heterochronous) + f"""
\t<coalescentLikelihoodIG id="coalescentIG" alpha="{ALPHA}" beta="{BETA}">
\t\t<populationTree>
\t\t\t<treeModel idref="treeModel"/>
\t\t</populationTree>
\t</coalescentLikelihoodIG>

\t<mcmc id="mcmc" chainLength="{CHAIN_LENGTH}" autoOptimize="true">
\t\t<coalescentLikelihoodIG idref="coalescentIG"/>

\t\t<log id="screenLog" logEvery="{LOG_EVERY * 10}">
\t\t\t<column label="coalescentIG" dp="4" width="12">
\t\t\t\t<coalescentLikelihoodIG idref="coalescentIG"/>
\t\t\t</column>
\t\t\t<column label="age(root)" sf="6" width="12">
\t\t\t\t<tmrcaStatistic idref="age(root)"/>
\t\t\t</column>
\t\t</log>

\t\t<log id="fileLog" logEvery="{LOG_EVERY}" fileName="{filename_stem}.log">
\t\t\t<coalescentLikelihoodIG idref="coalescentIG"/>
\t\t\t<parameter idref="treeModel.rootHeight"/>
\t\t\t<tmrcaStatistic idref="age(root)"/>
\t\t\t<treeLengthStatistic idref="treeLength"/>
\t\t</log>

\t\t<operators id="operators" optimizationSchedule="default">
{TREE_OPERATORS}
\t\t</operators>
\t</mcmc>

\t<report>
\t\t<property name="timer">
\t\t\t<mcmc idref="mcmc"/>
\t\t</property>
\t</report>

</beast>
"""


SCENARIOS = [
    (5, False, "prior_n5"),
    (50, False, "prior_n50"),
    (17, True, "prior_n17_hetero"),
]

if __name__ == "__main__":
    out_dir = os.path.dirname(os.path.abspath(__file__))
    for n, hetero, stem in SCENARIOS:
        full_path = os.path.join(out_dir, f"{stem}_full.xml")
        marg_path = os.path.join(out_dir, f"{stem}_marginalised.xml")
        with open(full_path, "w") as f:
            f.write(full_xml(n, hetero, f"{stem}_full"))
        with open(marg_path, "w") as f:
            f.write(marginalised_xml(n, hetero, f"{stem}_marginalised"))
        print(f"wrote {full_path}\nwrote {marg_path}")
