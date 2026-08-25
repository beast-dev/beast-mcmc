/*
 * OneZeroOneShuffleGibbsOperatorParserTest.java
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
import dr.evomodel.branchmodel.RewardMixtureAtomicPseudoPrior;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RewardsAwarePartialsRepresentation;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.OneZeroOneShuffleGibbsOperator;
import dr.inferencexml.operators.OneZeroOneShuffleGibbsOperatorParser;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import test.dr.math.MathTestCase;

import javax.xml.parsers.DocumentBuilderFactory;
import java.lang.reflect.Method;

/**
 * XML-level smoke test for the full-conditional reward-rate mapping Gibbs
 * refresh operator.
 *
 * @author Filippo Monti
 */
public class OneZeroOneShuffleGibbsOperatorParserTest extends MathTestCase {

    public void testParserAcceptsRewardRatesAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260818);

        final Fixture fixture = createFixture();
        final OneZeroOneShuffleGibbsOperatorParser parser = new OneZeroOneShuffleGibbsOperatorParser();

        final XMLObject xo = operatorXmlObject(fixture);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof OneZeroOneShuffleGibbsOperator);
        final OneZeroOneShuffleGibbsOperator operator = (OneZeroOneShuffleGibbsOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());

        assertIsValidPermutation(fixture.rewardRatesMapping);
        assertEquals(0.0, fixture.rewardRatesValues.getParameterValue(0), 1e-12);
        assertEquals(1.0, fixture.rewardRatesValues.getParameterValue(1), 1e-12);
    }

    public void testParserAcceptsAtomicPseudoPriorAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260818);

        final Fixture fixture = createFixture();
        fixture.ctsRewards.setParameterValue(0, 0.05);
        fixture.categoryState.setParameterValue(0, 1.50);

        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);

        final OneZeroOneShuffleGibbsOperatorParser parser = new OneZeroOneShuffleGibbsOperatorParser();
        final XMLObject xo = operatorXmlObject(fixture, pseudoPrior);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof OneZeroOneShuffleGibbsOperator);
        final OneZeroOneShuffleGibbsOperator operator = (OneZeroOneShuffleGibbsOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(pseudoPrior.getLogLikelihood());
        assertFinite(operator.doOperation());

        assertIsValidPermutation(fixture.rewardRatesMapping);
    }

    private static XMLObject operatorXmlObject(final Fixture fixture) throws Exception {
        return operatorXmlObject(fixture, null);
    }

    private static XMLObject operatorXmlObject(final Fixture fixture,
                                              final RewardMixtureAtomicPseudoPrior pseudoPrior) throws Exception {
        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                OneZeroOneShuffleGibbsOperatorParser.OPERATOR,
                MCMCOperator.WEIGHT, "1.0");

        addChild(xo, nativeObject(document, "rewardRates", fixture.rewardRates));
        addChild(xo, nativeObject(document, "treeDataLikelihood", fixture.independentLikelihood));
        if (pseudoPrior != null) {
            addChild(xo, nativeObject(document, "rewardMixtureAtomicPseudoPrior", pseudoPrior));
        }

        return xo;
    }

    private static Fixture createFixture() {
        final TreeModel tree = createThreeTipTree();
        final SitePatterns patterns = createSitePatterns("A", "C", "G");
        final SubstitutionModel substitutionModel = createNucleotideSubstitutionModel();
        final GammaSiteRateModel siteRateModel = new GammaSiteRateModel("shuffleGibbsSiteRateModel");

        final Parameter ctsRewards =
                new Parameter.Default("rewardCts", new double[]{0.50, 0.50, 0.50, 0.50});
        ctsRewards.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ctsRewards.getDimension()));
        final Parameter categoryState =
                new Parameter.Default("rewardCategory", new double[]{0.50, 0.50, 0.50, 0.50});
        final Parameter categoryCuts =
                new Parameter.Default("rewardCategoryCuts", new double[]{0, 1, 2, 3, 4, 5});

        final Parameter rewardRatesValues =
                new Parameter.Default("rewardRatesValues", new double[]{0.0, 1.0, 0.4, 0.6});
        final Parameter rewardRatesMapping =
                new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0});
        final RewardRates rewardRates = new RewardRates(
                rewardRatesValues,
                null,
                new Parameter.Default("rewardRatesVarying", new double[]{0.4, 0.6}),
                rewardRatesMapping
        );

        final RewardsAwareCategoricalMixtureBranchRates rewardBranchRates =
                new RewardsAwareCategoricalMixtureBranchRates(
                        tree,
                        ctsRewards,
                        categoryState,
                        categoryCuts,
                        rewardRates,
                        new ArbitraryBranchRates.BranchRateTransform.None(),
                        false,
                        TreeParameterModel.Type.WITHOUT_ROOT);
        final RewardsAwareBranchModel rewardsAwareBranchModel = new RewardsAwareBranchModel(
                tree,
                substitutionModel,
                rewardRates,
                categoryState,
                categoryCuts,
                rewardBranchRates,
                false
        );

        final DiscreteDataLikelihoodDelegate independentDelegate = new DiscreteDataLikelihoodDelegate(
                tree,
                patterns,
                rewardsAwareBranchModel,
                siteRateModel,
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
                new TreeDataLikelihood(independentDelegate, tree, rewardBranchRates);

        return new Fixture(independentLikelihood, rewardsAwareBranchModel, rewardBranchRates,
                ctsRewards, categoryState, rewardRates, rewardRatesValues, rewardRatesMapping);
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

        final Tree simpleTree = new SimpleTree(root);
        return new DefaultTreeModel("shuffleGibbsParserTree", simpleTree);
    }

    private static SubstitutionModel createNucleotideSubstitutionModel() {
        final DataType dataType = Nucleotides.INSTANCE;
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.25, 0.25, 0.25, 0.25});
        return new HKY(new Parameter.Default("shuffleGibbsKappa", 2.0), frequencyModel);
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

    private static void addChild(final XMLObject parent, final XMLObject child) throws Exception {
        final Method addChild = XMLObject.class.getDeclaredMethod("addChild", Object.class);
        addChild.setAccessible(true);
        addChild.invoke(parent, child);
    }

    private static void assertFinite(final double value) {
        assertTrue("Expected finite value but found " + value,
                !Double.isNaN(value) && !Double.isInfinite(value));
    }

    private static void assertIsValidPermutation(final Parameter mapping) {
        final int n = mapping.getDimension();
        final boolean[] seen = new boolean[n];
        for (int i = 0; i < n; i++) {
            final int value = (int) Math.round(mapping.getParameterValue(i));
            assertTrue("Mapping value out of range: " + value, value >= 0 && value < n);
            assertTrue("Mapping value repeated: " + value, !seen[value]);
            seen[value] = true;
        }
    }

    private static final class Fixture {
        final TreeDataLikelihood independentLikelihood;
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        final ArbitraryBranchRates rewardBranchRates;
        final Parameter ctsRewards;
        final Parameter categoryState;
        final RewardRates rewardRates;
        final Parameter rewardRatesValues;
        final Parameter rewardRatesMapping;

        private Fixture(final TreeDataLikelihood independentLikelihood,
                        final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final ArbitraryBranchRates rewardBranchRates,
                        final Parameter ctsRewards,
                        final Parameter categoryState,
                        final RewardRates rewardRates,
                        final Parameter rewardRatesValues,
                        final Parameter rewardRatesMapping) {
            this.independentLikelihood = independentLikelihood;
            this.rewardsAwareBranchModel = rewardsAwareBranchModel;
            this.rewardBranchRates = rewardBranchRates;
            this.ctsRewards = ctsRewards;
            this.categoryState = categoryState;
            this.rewardRates = rewardRates;
            this.rewardRatesValues = rewardRatesValues;
            this.rewardRatesMapping = rewardRatesMapping;
        }
    }
}
