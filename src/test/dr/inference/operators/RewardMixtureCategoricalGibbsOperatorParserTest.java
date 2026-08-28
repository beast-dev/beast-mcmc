/*
 * RewardMixtureCategoricalGibbsOperatorParserTest.java
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
import dr.evomodel.branchratemodel.PerBranchRewardMixtureCategoryDecoder;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRates;
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesDynamic;
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
import dr.inference.operators.RewardMixtureBranchJointCtsCategoryOperator;
import dr.inference.operators.RewardMixtureCategoricalGibbsOperator;
import dr.inference.operators.RewardMixtureContinuousBranchSliceOperator;
import dr.inference.operators.RewardsMixtureBranchWeightProvider;
import dr.inferencexml.operators.RewardMixtureBranchJointCtsCategoryOperatorParser;
import dr.inferencexml.operators.RewardMixtureCategoricalGibbsOperatorParser;
import dr.inferencexml.operators.RewardMixtureContinuousBranchSliceOperatorParser;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.xml.XMLObject;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import test.dr.math.MathTestCase;

import javax.xml.parsers.DocumentBuilderFactory;
import java.lang.reflect.Method;

/**
 * XML-level smoke test for the embedded categorical reward-mixture Gibbs
 * refresh operator.
 *
 * @author Filippo Monti
 */
public class RewardMixtureCategoricalGibbsOperatorParserTest extends MathTestCase {

    public void testContinuousBranchSliceParserAcceptsStateAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260826);

        final Fixture fixture = createFixture();
        final RewardMixtureContinuousBranchSliceOperatorParser parser =
                new RewardMixtureContinuousBranchSliceOperatorParser();

        final XMLObject xo = continuousSliceOperatorXmlObject(fixture);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof RewardMixtureContinuousBranchSliceOperator);

        final RewardMixtureContinuousBranchSliceOperator operator =
                (RewardMixtureContinuousBranchSliceOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());
        assertCategoryValuesAreInsideCuts(fixture.categoryState, fixture.categoryCuts);
        assertCtsValuesAreInsideBounds(fixture.ctsRewards);
    }

    public void testJointBranchCtsCategoryParserAcceptsStateAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260828);

        final Fixture fixture = createFixture();
        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);
        final RewardMixtureBranchJointCtsCategoryOperatorParser parser =
                new RewardMixtureBranchJointCtsCategoryOperatorParser();

        final XMLObject xo = jointOperatorXmlObject(fixture, pseudoPrior);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof RewardMixtureBranchJointCtsCategoryOperator);

        final RewardMixtureBranchJointCtsCategoryOperator operator =
                (RewardMixtureBranchJointCtsCategoryOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());
        assertCategoryValuesAreInsideCuts(fixture.categoryState, fixture.categoryCuts);
        assertCtsValuesAreInsideBounds(fixture.ctsRewards);
    }

    public void testJointBranchCtsCategoryOperatorRunsWithDynamicDecoder() {
        MathUtils.setSeed(20260830);

        final Fixture fixture = createFixture(true);
        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);
        final RewardMixtureBranchJointCtsCategoryOperator operator =
                new RewardMixtureBranchJointCtsCategoryOperator(
                        fixture.categoryState,
                        fixture.categoryCuts,
                        fixture.rewardsAwareBranchModel,
                        fixture.independentLikelihood,
                        new TreeDataLikelihood[0],
                        new TreeDataLikelihood[0],
                        pseudoPrior,
                        0.25,
                        10,
                        100,
                        1.0);

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());
        fixture.rewardsAwareBranchModel.refreshCategoryDecoderEmbedding();
        assertCategoryValuesAreInsideCuts(fixture.categoryState, fixture.categoryCuts);
        assertCtsValuesAreInsideBounds(fixture.ctsRewards);
        for (int i = 0; i < fixture.categoryState.getDimension(); i++) {
            final int category = fixture.rewardsAwareBranchModel.getCategoryDecoder()
                    .getCategoryForParameterIndex(i);
            assertTrue("Decoded category out of range: " + category,
                    category >= 0 && category < 5);
        }
    }

    public void testParserAcceptsCategoricalRewardStateAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260701);

        final Fixture fixture = createFixture();
        final RewardMixtureCategoricalGibbsOperatorParser parser =
                new RewardMixtureCategoricalGibbsOperatorParser();

        final XMLObject xo = operatorXmlObject(fixture);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof RewardMixtureCategoricalGibbsOperator);

        final RewardMixtureCategoricalGibbsOperator operator =
                (RewardMixtureCategoricalGibbsOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());
        assertCategoryValuesAreInsideCuts(fixture.categoryState, fixture.categoryCuts);
    }

    public void testAtomicPseudoPriorPullsOnlyAtomicBranches() {
        final Fixture fixture = createFixture();
        fixture.ctsRewards.setParameterValue(0, 0.50);
        fixture.categoryState.setParameterValue(0, 2.50);

        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);

        final double expectedLogDensity = NormalDistribution.logPdf(0.50, 0.40, 0.10) -
                Math.log(NormalDistribution.cdf(1.0, 0.40, 0.10) -
                        NormalDistribution.cdf(0.0, 0.40, 0.10));

        assertEquals(expectedLogDensity, pseudoPrior.getLogDensityForCurrentCategory(0), 1E-10);
        assertEquals(0.0, pseudoPrior.getLogDensityForCurrentCategory(1), 1E-10);

        final double[] gradient = pseudoPrior.getGradientLogDensity();
        assertEquals(-10.0, gradient[0], 1E-10);
        for (int i = 1; i < gradient.length; i++) {
            assertEquals(0.0, gradient[i], 1E-10);
        }
    }

    public void testAtomicPseudoPriorCandidateValueMatchesTemporarySetValue() {
        final Fixture fixture = createFixture();
        fixture.ctsRewards.setParameterValue(0, 0.50);
        fixture.categoryState.setParameterValue(0, 2.50);

        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);

        final double candidate = 0.55;
        final int category = 2;
        final double candidateLogDensity =
                pseudoPrior.getLogDensityForCategoryAtValue(0, category, candidate);
        final double currentLogDensity =
                pseudoPrior.getLogDensityForCategory(0, category);
        assertTrue("Candidate density should not read the stored CTS value",
                Math.abs(candidateLogDensity - currentLogDensity) > 1E-8);

        fixture.ctsRewards.setParameterValue(0, candidate);
        assertEquals(candidateLogDensity, pseudoPrior.getLogDensityForCategory(0, category), 1E-10);
        assertEquals(-15.0, pseudoPrior.getGradientForCategoryAtValue(0, category, candidate), 1E-10);
    }

    public void testDynamicDecoderCandidateCtsCutsMatchRefreshedEmbedding() {
        final Parameter ctsRewards =
                new Parameter.Default("dynamicDecoderCts", new double[]{0.30});
        final Parameter categoryState =
                new Parameter.Default("dynamicDecoderCategory", new double[]{0.50});
        final Parameter categoryCuts =
                new Parameter.Default("dynamicDecoderCuts", new double[]{0, 1, 2, 3, 4, 5});
        final RewardRates rewardRates = createRewardRates();
        final PerBranchRewardMixtureCategoryDecoder decoder =
                new PerBranchRewardMixtureCategoryDecoder(
                        categoryState,
                        categoryCuts,
                        ctsRewards,
                        rewardRates,
                        4,
                        1);

        assertEquals(1.0, decoder.getLowerCut(0, 0), 0.0);
        assertEquals(2.0, decoder.getUpperCut(0, 0), 0.0);

        final double candidate = 0.70;
        final double candidateLower = decoder.getLowerCutForCategoryAtCtsValue(0, 0, candidate);
        final double candidateUpper = decoder.getUpperCutForCategoryAtCtsValue(0, 0, candidate);
        assertEquals(3.0, candidateLower, 0.0);
        assertEquals(4.0, candidateUpper, 0.0);

        ctsRewards.setParameterValue(0, candidate);
        decoder.refreshEmbedding();
        assertEquals(candidateLower, decoder.getLowerCut(0, 0), 0.0);
        assertEquals(candidateUpper, decoder.getUpperCut(0, 0), 0.0);
    }

    public void testJointOperatorDoesNotMutateWhenAllCategoryMassesAreInvalid() {
        MathUtils.setSeed(20260829);

        final Fixture fixture = createFixture();
        for (int i = 0; i < fixture.ctsRewards.getDimension(); i++) {
            fixture.ctsRewards.setParameterValue(i, Double.NaN);
            fixture.categoryState.setParameterValue(i, 2.50);
        }
        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);
        final RewardMixtureBranchJointCtsCategoryOperator operator =
                new RewardMixtureBranchJointCtsCategoryOperator(
                        fixture.categoryState,
                        fixture.categoryCuts,
                        fixture.rewardsAwareBranchModel,
                        fixture.independentLikelihood,
                        new TreeDataLikelihood[0],
                        new TreeDataLikelihood[0],
                        pseudoPrior,
                        0.05,
                        10,
                        100,
                        1.0);

        final double result = operator.doOperation();
        assertTrue("Expected invalid operation but found " + result,
                Double.isInfinite(result) && result < 0.0);
        for (int i = 0; i < fixture.ctsRewards.getDimension(); i++) {
            assertTrue(Double.isNaN(fixture.ctsRewards.getParameterValue(i)));
            assertEquals(2.50, fixture.categoryState.getParameterValue(i), 0.0);
        }
    }

    public void testContinuousBranchLocalWeightDeltaMatchesFullLikelihoodDelta() {
        final Fixture fixture = createFixture();
        final RewardsMixtureBranchWeightProvider provider =
                new RewardsMixtureBranchWeightProvider(
                        fixture.rewardsAwareBranchModel,
                        fixture.independentLikelihood,
                        new TreeDataLikelihood[0],
                        new TreeDataLikelihood[0]);

        final int parameterIndex = 0;
        final double current = fixture.ctsRewards.getParameterValue(parameterIndex);
        final double proposed = 0.55;

        final double fullCurrent = fixture.independentLikelihood.getLogLikelihood();
        provider.beginSingleBranchOperationCache(parameterIndex);
        final double localCurrent =
                provider.computeContinuousLogWeightForParameterIndex(parameterIndex, current);
        final double localProposed =
                provider.computeContinuousLogWeightForParameterIndex(parameterIndex, proposed);

        fixture.ctsRewards.setParameterValue(parameterIndex, proposed);
        fixture.independentLikelihood.makeDirty();
        final double fullProposed = fixture.independentLikelihood.getLogLikelihood();

        assertFinite(fullCurrent);
        assertFinite(fullProposed);
        assertFinite(localCurrent);
        assertFinite(localProposed);
        assertEquals(fullProposed - fullCurrent, localProposed - localCurrent, 1E-8);
    }

    public void testParserAcceptsAtomicPseudoPriorAndProposalRuns() throws Exception {
        MathUtils.setSeed(20260702);

        final Fixture fixture = createFixture();
        final RewardMixtureAtomicPseudoPrior pseudoPrior =
                new RewardMixtureAtomicPseudoPrior(
                        fixture.rewardsAwareBranchModel,
                        fixture.rewardBranchRates,
                        0.10);
        final RewardMixtureCategoricalGibbsOperatorParser parser =
                new RewardMixtureCategoricalGibbsOperatorParser();

        final XMLObject xo = operatorXmlObject(fixture, pseudoPrior);
        final Object parsed = parser.parseXMLObject(xo);

        assertTrue(parsed instanceof RewardMixtureCategoricalGibbsOperator);

        final RewardMixtureCategoricalGibbsOperator operator =
                (RewardMixtureCategoricalGibbsOperator) parsed;

        assertFinite(fixture.independentLikelihood.getLogLikelihood());
        assertFinite(operator.doOperation());
        assertCategoryValuesAreInsideCuts(fixture.categoryState, fixture.categoryCuts);
    }

    private static XMLObject operatorXmlObject(final Fixture fixture) throws Exception {
        return operatorXmlObject(fixture, null);
    }

    private static XMLObject operatorXmlObject(final Fixture fixture,
                                               final RewardMixtureAtomicPseudoPrior pseudoPrior) throws Exception {
        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                RewardMixtureCategoricalGibbsOperatorParser.OPERATOR_NAME,
                MCMCOperator.WEIGHT, "1.0",
                RewardMixtureCategoricalGibbsOperatorParser.UPDATE_PROPORTION, "1.0");

        addChild(xo, wrapper(document, "categoryState", nativeObject(document, "parameter", fixture.categoryState)));
        addChild(xo, wrapper(document, "categoryCuts", nativeObject(document, "parameter", fixture.categoryCuts)));
        addChild(xo, nativeObject(document, "rewardsAwareBranchModel", fixture.rewardsAwareBranchModel));
        addChild(xo, nativeObject(document, "treeDataLikelihood", fixture.independentLikelihood));
        if (pseudoPrior != null) {
            addChild(xo, nativeObject(document, "rewardMixtureAtomicPseudoPrior", pseudoPrior));
        }

        return xo;
    }

    private static XMLObject jointOperatorXmlObject(final Fixture fixture,
                                                    final RewardMixtureAtomicPseudoPrior pseudoPrior) throws Exception {
        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                RewardMixtureBranchJointCtsCategoryOperatorParser.OPERATOR_NAME,
                MCMCOperator.WEIGHT, "1.0",
                RewardMixtureBranchJointCtsCategoryOperatorParser.WINDOW_SIZE, "0.05");

        addChild(xo, wrapper(document, "categoryState", nativeObject(document, "parameter", fixture.categoryState)));
        addChild(xo, wrapper(document, "categoryCuts", nativeObject(document, "parameter", fixture.categoryCuts)));
        addChild(xo, nativeObject(document, "rewardsAwareBranchModel", fixture.rewardsAwareBranchModel));
        addChild(xo, nativeObject(document, "treeDataLikelihood", fixture.independentLikelihood));
        if (pseudoPrior != null) {
            addChild(xo, nativeObject(document, "rewardMixtureAtomicPseudoPrior", pseudoPrior));
        }

        return xo;
    }

    private static XMLObject continuousSliceOperatorXmlObject(final Fixture fixture) throws Exception {
        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                RewardMixtureContinuousBranchSliceOperatorParser.OPERATOR_NAME,
                MCMCOperator.WEIGHT, "1.0",
                RewardMixtureContinuousBranchSliceOperatorParser.WINDOW_SIZE, "0.05");

        addChild(xo, wrapper(document, "categoryState", nativeObject(document, "parameter", fixture.categoryState)));
        addChild(xo, wrapper(document, "categoryCuts", nativeObject(document, "parameter", fixture.categoryCuts)));
        addChild(xo, nativeObject(document, "rewardsAwareBranchModel", fixture.rewardsAwareBranchModel));
        addChild(xo, nativeObject(document, "treeDataLikelihood", fixture.independentLikelihood));

        return xo;
    }

    private static Fixture createFixture() {
        return createFixture(false);
    }

    private static Fixture createFixture(final boolean dynamicCategoryOrdering) {
        final TreeModel tree = createThreeTipTree();
        final SitePatterns patterns = createSitePatterns("A", "C", "G");
        final SubstitutionModel substitutionModel = createNucleotideSubstitutionModel();
        final GammaSiteRateModel siteRateModel = new GammaSiteRateModel("categoricalGibbsSiteRateModel");

        final Parameter ctsRewards =
                new Parameter.Default("rewardCts", new double[]{0.50, 0.50, 0.50, 0.50});
        ctsRewards.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ctsRewards.getDimension()));
        final double initialCategoryValue = dynamicCategoryOrdering ? 2.50 : 0.50;
        final Parameter categoryState =
                new Parameter.Default("rewardCategory", new double[]{
                        initialCategoryValue,
                        initialCategoryValue,
                        initialCategoryValue,
                        initialCategoryValue});
        final Parameter categoryCuts =
                new Parameter.Default("rewardCategoryCuts", new double[]{0, 1, 2, 3, 4, 5});
        final RewardRates rewardRates = createRewardRates();

        final ArbitraryBranchRates rewardBranchRates;
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        if (dynamicCategoryOrdering) {
            final RewardsAwareCategoricalMixtureBranchRatesDynamic dynamicRewardBranchRates =
                    new RewardsAwareCategoricalMixtureBranchRatesDynamic(
                            tree,
                            ctsRewards,
                            categoryState,
                            categoryCuts,
                            rewardRates,
                            new ArbitraryBranchRates.BranchRateTransform.None(),
                            false,
                            TreeParameterModel.Type.WITHOUT_ROOT);
            rewardBranchRates = dynamicRewardBranchRates;
            rewardsAwareBranchModel = new RewardsAwareBranchModel(
                    tree,
                    substitutionModel,
                    rewardRates,
                    dynamicRewardBranchRates.getCategoryDecoder(),
                    rewardBranchRates,
                    false,
                    RewardsAwareBranchModel.DEFAULT_SERICOLA_SERIES_RESCALING);
        } else {
            final RewardsAwareCategoricalMixtureBranchRates staticRewardBranchRates =
                    new RewardsAwareCategoricalMixtureBranchRates(
                        tree,
                        ctsRewards,
                        categoryState,
                        categoryCuts,
                        rewardRates,
                        new ArbitraryBranchRates.BranchRateTransform.None(),
                        false,
                        TreeParameterModel.Type.WITHOUT_ROOT);
            rewardBranchRates = staticRewardBranchRates;
            rewardsAwareBranchModel = new RewardsAwareBranchModel(
                    tree,
                    substitutionModel,
                    rewardRates,
                    categoryState,
                    categoryCuts,
                    rewardBranchRates,
                    false
            );
        }

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

        return new Fixture(
                rewardsAwareBranchModel,
                rewardBranchRates,
                independentLikelihood,
                ctsRewards,
                categoryState,
                categoryCuts);
    }

    private static RewardRates createRewardRates() {
        return new RewardRates(
                new Parameter.Default("rewardRates", new double[]{0.20, 0.40, 0.60, 0.80}),
                null,
                new Parameter.Default("rewardRatesInternal", new double[0]),
                new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );
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
        return new DefaultTreeModel("categoricalGibbsParserTree", simpleTree);
    }

    private static SubstitutionModel createNucleotideSubstitutionModel() {
        final DataType dataType = Nucleotides.INSTANCE;
        final FrequencyModel frequencyModel = new FrequencyModel(dataType, new double[]{0.25, 0.25, 0.25, 0.25});
        return new HKY(new Parameter.Default("categoricalGibbsKappa", 2.0), frequencyModel);
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

    private static XMLObject wrapper(final Document document, final String name, final XMLObject child)
            throws Exception {
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

    private static void assertCategoryValuesAreInsideCuts(final Parameter categoryState,
                                                          final Parameter categoryCuts) {
        final double lower = categoryCuts.getParameterValue(0);
        final double upper = categoryCuts.getParameterValue(categoryCuts.getDimension() - 1);
        for (int i = 0; i < categoryState.getDimension(); i++) {
            final double value = categoryState.getParameterValue(i);
            assertTrue("Category value out of cut support: " + value,
                    value >= lower && value <= upper);
        }
    }

    private static void assertCtsValuesAreInsideBounds(final Parameter ctsRewards) {
        for (int i = 0; i < ctsRewards.getDimension(); i++) {
            final double value = ctsRewards.getParameterValue(i);
            assertTrue("CTS value below lower bound: " + value,
                    value >= ctsRewards.getBounds().getLowerLimit(i));
            assertTrue("CTS value above upper bound: " + value,
                    value <= ctsRewards.getBounds().getUpperLimit(i));
        }
    }

    private static final class Fixture {
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        final ArbitraryBranchRates rewardBranchRates;
        final TreeDataLikelihood independentLikelihood;
        final Parameter ctsRewards;
        final Parameter categoryState;
        final Parameter categoryCuts;

        private Fixture(final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final ArbitraryBranchRates rewardBranchRates,
                        final TreeDataLikelihood independentLikelihood,
                        final Parameter ctsRewards,
                        final Parameter categoryState,
                        final Parameter categoryCuts) {
            this.rewardsAwareBranchModel = rewardsAwareBranchModel;
            this.rewardBranchRates = rewardBranchRates;
            this.independentLikelihood = independentLikelihood;
            this.ctsRewards = ctsRewards;
            this.categoryState = categoryState;
            this.categoryCuts = categoryCuts;
        }
    }
}
