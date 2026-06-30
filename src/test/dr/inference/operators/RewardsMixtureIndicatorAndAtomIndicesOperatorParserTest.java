/*
 * RewardsMixtureIndicatorAndAtomIndicesOperatorParserTest.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.operators;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RewardsAwarePartialsRepresentation;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.RewardsMixtureIndicatorAndAtomIndicesOperator;
import dr.inferencexml.operators.RewardsMixtureIndicatorAndAtomIndicesOperatorParser;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import test.dr.math.MathTestCase;

import javax.xml.parsers.DocumentBuilderFactory;
import java.lang.reflect.Method;

/**
 * XML-level smoke test for parsing the reward-mixture operator with dependent
 * CTMC likelihoods and executing one proposal.
 *
 * @author Filippo Monti
 */
public class RewardsMixtureIndicatorAndAtomIndicesOperatorParserTest extends MathTestCase {

    public void testParserAcceptsDependentCtmcLikelihoodsAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260610);

        final Fixture fixture = createFixture();
        final RewardsMixtureIndicatorAndAtomIndicesOperatorParser parser =
                new RewardsMixtureIndicatorAndAtomIndicesOperatorParser();

        final XMLObject xo = operatorXmlObject(fixture);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof RewardsMixtureIndicatorAndAtomIndicesOperator);

        final RewardsMixtureIndicatorAndAtomIndicesOperator operator =
                (RewardsMixtureIndicatorAndAtomIndicesOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(fixture.dependentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());
        assertIndicatorValuesAreBinary(fixture.indicator);
        assertAtomIndicesAreInRange(fixture.atomIndices, 4);
    }

    private static XMLObject operatorXmlObject(final Fixture fixture) throws Exception {
        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                RewardsMixtureIndicatorAndAtomIndicesOperatorParser.OPERATOR_NAME,
                MCMCOperator.WEIGHT, "1.0",
                "updateProportion", "1.0",
                "autoOptimize", "false",
                "useClusterMoves", "false");

        addChild(xo, wrapper(document, "indicatorZ", nativeObject(document, "parameter", fixture.indicator)));
        addChild(xo, wrapper(document, "atomIndex", nativeObject(document, "parameter", fixture.atomIndices)));
        addChild(xo, nativeObject(document, "rewardsAwareBranchModel", fixture.rewardsAwareBranchModel));
        addChild(xo, nativeObject(document, "treeDataLikelihood", fixture.independentLikelihood));
        addChild(xo, wrapper(document, "dependentCtmcLikelihoods",
                nativeObject(document, "treeDataLikelihood", fixture.dependentLikelihood)));

        return xo;
    }

    private static Fixture createFixture() {
        final TreeModel tree = createThreeTipTree();
        final SitePatterns independentPatterns = createSitePatterns("A", "C", "G");
        final SitePatterns dependentPatterns = createSitePatterns("ACGT", "TGCA", "CAGT");
        final SubstitutionModel substitutionModel = createNucleotideSubstitutionModel();
        final GammaSiteRateModel independentSiteRateModel = new GammaSiteRateModel("independentRewardRateModel");
        final GammaSiteRateModel dependentSiteRateModel = new GammaSiteRateModel("dependentRewardRateModel");

        final Parameter ctsRewards = new Parameter.Default("rewardCts", new double[]{0.50, 0.50, 0.50, 0.50});
        final Parameter indicator = new Parameter.Default("rewardIndicator", new double[]{0.0, 0.0, 0.0, 0.0});
        final Parameter atomIndices = new Parameter.Default("rewardAtomIndices", new double[]{0.0, 1.0, 2.0, 3.0});
        final RewardRates rewardRates = new RewardRates(
                new Parameter.Default("rewardRates", new double[]{0.20, 0.40, 0.60, 0.80}),
                null,
                new Parameter.Default("rewardRatesInternal", new double[0]),
                new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );
        final RewardsAwareMixtureBranchRates rewardBranchRates = new RewardsAwareMixtureBranchRates(
                tree,
                ctsRewards,
                indicator,
                atomIndices,
                rewardRates,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false,
                TreeParameterModel.Type.WITHOUT_ROOT
        );
        final RewardsAwareBranchModel rewardsAwareBranchModel = new RewardsAwareBranchModel(
                tree,
                substitutionModel,
                rewardRates,
                indicator,
                rewardBranchRates,
                atomIndices,
                false
        );

        final DiscreteDataLikelihoodDelegate independentDelegate = new DiscreteDataLikelihoodDelegate(
                tree,
                independentPatterns,
                rewardsAwareBranchModel,
                independentSiteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, false, false, false),
                new RewardsAwarePartialsRepresentation(rewardsAwareBranchModel),
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                null,
                null
        );
        final TreeDataLikelihood independentLikelihood =
                new TreeDataLikelihood(independentDelegate, tree, new DefaultBranchRateModel());

        final BeagleDataLikelihoodDelegate dependentDelegate = new BeagleDataLikelihoodDelegate(
                tree,
                dependentPatterns,
                new HomogeneousBranchModel(substitutionModel),
                dependentSiteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, true, false, false, false)
        );
        final TreeDataLikelihood dependentLikelihood =
                new TreeDataLikelihood(dependentDelegate, tree, rewardBranchRates);

        return new Fixture(rewardsAwareBranchModel, independentLikelihood, dependentLikelihood, indicator, atomIndices);
    }

    private static SitePatterns createSitePatterns(final String a, final String b, final String c) {
        final DataType dataType = Nucleotides.INSTANCE;
        final SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);
        alignment.addSequence(sequence("a", a, dataType));
        alignment.addSequence(sequence("b", b, dataType));
        alignment.addSequence(sequence("c", c, dataType));
        return new SitePatterns(alignment, null, 0, -1, 1, true);
    }

    private static Sequence sequence(final String taxonName, final String sequenceString, final DataType dataType) {
        final Sequence sequence = new Sequence(new Taxon(taxonName), sequenceString);
        sequence.setDataType(dataType);
        return sequence;
    }

    private static TreeModel createThreeTipTree() {
        final SimpleNode left = new SimpleNode();
        left.setTaxon(new Taxon("a"));
        left.setHeight(0.0);

        final SimpleNode right = new SimpleNode();
        right.setTaxon(new Taxon("b"));
        right.setHeight(0.0);

        final SimpleNode internal = new SimpleNode();
        internal.setHeight(0.4);
        internal.addChild(left);
        internal.addChild(right);

        final SimpleNode outgroup = new SimpleNode();
        outgroup.setTaxon(new Taxon("c"));
        outgroup.setHeight(0.0);

        final SimpleNode root = new SimpleNode();
        root.setHeight(1.0);
        root.addChild(internal);
        root.addChild(outgroup);

        final Tree tree = new SimpleTree(root);
        return new DefaultTreeModel("rewardMixtureOperatorParserTree", tree);
    }

    private static SubstitutionModel createNucleotideSubstitutionModel() {
        final DataType dataType = Nucleotides.INSTANCE;
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.25, 0.25, 0.25, 0.25});
        return new HKY(new Parameter.Default("rewardMixtureParserKappa", 2.0), frequencyModel);
    }

    private static XMLObject xmlObject(final Document document, final String name, final String... attributes) {
        final Element element = document.createElement(name);
        for (int i = 0; i < attributes.length; i += 2) {
            element.setAttribute(attributes[i], attributes[i + 1]);
        }
        return new XMLObject(element, null);
    }

    private static XMLObject nativeObject(final Document document, final String name, final Object object) {
        final XMLObject xo = xmlObject(document, name);
        xo.setNativeObject(object);
        return xo;
    }

    private static XMLObject wrapper(final Document document, final String name, final XMLObject child) throws Exception {
        final XMLObject xo = xmlObject(document, name);
        addChild(xo, child);
        return xo;
    }

    private static void addChild(final XMLObject parent, final XMLObject child) throws Exception {
        final Method addChild = XMLObject.class.getDeclaredMethod("addChild", Object.class);
        addChild.setAccessible(true);
        addChild.invoke(parent, child);
    }

    private static void assertFinite(final double value) {
        assertTrue("Expected finite value but found " + value,
                !Double.isNaN(value) && !Double.isInfinite(value));
    }

    private static void assertIndicatorValuesAreBinary(final Parameter indicator) {
        for (int i = 0; i < indicator.getDimension(); i++) {
            final double value = indicator.getParameterValue(i);
            assertTrue("Indicator value must be 0 or 1, found " + value,
                    value == 0.0 || value == 1.0);
        }
    }

    private static void assertAtomIndicesAreInRange(final Parameter atomIndices, final int stateCount) {
        for (int i = 0; i < atomIndices.getDimension(); i++) {
            final int value = (int) Math.round(atomIndices.getParameterValue(i));
            assertTrue("Atom index out of range: " + value, value >= 0 && value < stateCount);
        }
    }

    private static final class Fixture {
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        final TreeDataLikelihood independentLikelihood;
        final TreeDataLikelihood dependentLikelihood;
        final Parameter indicator;
        final Parameter atomIndices;

        private Fixture(final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final TreeDataLikelihood independentLikelihood,
                        final TreeDataLikelihood dependentLikelihood,
                        final Parameter indicator,
                        final Parameter atomIndices) {
            this.rewardsAwareBranchModel = rewardsAwareBranchModel;
            this.independentLikelihood = independentLikelihood;
            this.dependentLikelihood = dependentLikelihood;
            this.indicator = indicator;
            this.atomIndices = atomIndices;
        }
    }
}
