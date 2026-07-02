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
import dr.inference.operators.RewardMixtureCategoricalGibbsOperator;
import dr.inferencexml.operators.RewardMixtureCategoricalGibbsOperatorParser;
import dr.math.MathUtils;
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

    private static XMLObject operatorXmlObject(final Fixture fixture) throws Exception {
        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                RewardMixtureCategoricalGibbsOperatorParser.OPERATOR_NAME,
                MCMCOperator.WEIGHT, "1.0",
                RewardMixtureCategoricalGibbsOperatorParser.UPDATE_PROPORTION, "1.0");

        addChild(xo, wrapper(document, "categoryState", nativeObject(document, "parameter", fixture.categoryState)));
        addChild(xo, wrapper(document, "categoryCuts", nativeObject(document, "parameter", fixture.categoryCuts)));
        addChild(xo, nativeObject(document, "rewardsAwareBranchModel", fixture.rewardsAwareBranchModel));
        addChild(xo, nativeObject(document, "treeDataLikelihood", fixture.independentLikelihood));

        return xo;
    }

    private static Fixture createFixture() {
        final TreeModel tree = createThreeTipTree();
        final SitePatterns patterns = createSitePatterns("A", "C", "G");
        final SubstitutionModel substitutionModel = createNucleotideSubstitutionModel();
        final GammaSiteRateModel siteRateModel = new GammaSiteRateModel("categoricalGibbsSiteRateModel");

        final Parameter ctsRewards =
                new Parameter.Default("rewardCts", new double[]{0.50, 0.50, 0.50, 0.50});
        final Parameter categoryState =
                new Parameter.Default("rewardCategory", new double[]{0.50, 0.50, 0.50, 0.50});
        final Parameter categoryCuts =
                new Parameter.Default("rewardCategoryCuts", new double[]{0, 1, 2, 3, 4, 5});
        final RewardRates rewardRates = new RewardRates(
                new Parameter.Default("rewardRates", new double[]{0.20, 0.40, 0.60, 0.80}),
                null,
                new Parameter.Default("rewardRatesInternal", new double[0]),
                new Parameter.Default("rewardRatesMapping", new double[]{0.0, 1.0, 2.0, 3.0})
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

        return new Fixture(rewardsAwareBranchModel, independentLikelihood, categoryState, categoryCuts);
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

    private static final class Fixture {
        final RewardsAwareBranchModel rewardsAwareBranchModel;
        final TreeDataLikelihood independentLikelihood;
        final Parameter categoryState;
        final Parameter categoryCuts;

        private Fixture(final RewardsAwareBranchModel rewardsAwareBranchModel,
                        final TreeDataLikelihood independentLikelihood,
                        final Parameter categoryState,
                        final Parameter categoryCuts) {
            this.rewardsAwareBranchModel = rewardsAwareBranchModel;
            this.independentLikelihood = independentLikelihood;
            this.categoryState = categoryState;
            this.categoryCuts = categoryCuts;
        }
    }
}
