/*
 * StaleScaleFactorGradientTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
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
 *
 */

package test.dr.evomodel.treedatalikelihood;

import beagle.BeagleInfo;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.BeagleFunctionality;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode;
import dr.evomodel.treedatalikelihood.discrete.LogCtmcRateGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * The dynamic rescaling scheme reuses the scale factors of an earlier evaluation. After a large change, partials
 * scaled with those stale factors can be so large that the pre-order traversal overflows while the likelihood is
 * still finite. The gradient must then be computed again with fresh scale factors.
 * <p>
 * The tree is a spine of mismatched cherries above a ladder of mismatched tips. In state A the cherries' branches
 * are nearly zero, so their scale factors are tiny; in state B the ladder's are. With A's factors in B, the cherries'
 * partials are large and the ladder's small: the likelihood is finite, but the pre-order partial below the cherries,
 * a product of theirs, overflows.
 */
public class StaleScaleFactorGradientTest extends TestCase {

    private static final int CHERRIES = 60;
    private static final int LADDER_TIPS = 21;
    private static final double SHORT = 1E-6;

    public StaleScaleFactorGradientTest(String name) {
        super(name);
    }

    @Override
    public void runBare() throws Throwable {
        if (!BeagleFunctionality.IS_MATRIX_VECTOR_SUPPORTED()) {
            System.err.println("SKIPPED " + getClass().getSimpleName() + "." + getName() +
                    ": needs native BEAGLE 4.2.0 or later, found " + BeagleInfo.getVersion());
            return;
        }
        super.runBare();
    }

    public void testGradientWithStaleScaleFactors() {

        Fixture dynamic = new Fixture(PartialsRescalingScheme.DYNAMIC);
        Fixture always = new Fixture(PartialsRescalingScheme.ALWAYS);

        for (Fixture f : new Fixture[]{dynamic, always}) {
            f.setBranchRates(false);
            for (int i = 0; i < 4; ++i) { // past the dynamic scheme's initial rescalings, which are in state A
                f.likelihood.makeDirty();
                f.likelihood.getLogLikelihood();
            }
            f.setBranchRates(true);
        }

        final double expected = always.likelihood.getLogLikelihood();
        assertEquals(expected, dynamic.likelihood.getLogLikelihood(), 1E-8 * Math.abs(expected));

        final double[] reference = always.gradient.getGradientLogDensity();
        final double[] gradient = dynamic.gradient.getGradientLogDensity();
        for (int k = 0; k < reference.length; ++k) {
            assertTrue("reference gradient " + k + " is " + reference[k], Double.isFinite(reference[k]));
            assertEquals("gradient " + k, reference[k], gradient[k], 1E-8 * Math.max(1.0, Math.abs(reference[k])));
        }
    }

    private static class Fixture {

        final TreeDataLikelihood likelihood;
        final LogCtmcRateGradient gradient;
        final DefaultTreeModel tree;
        final ArbitraryBranchRates branchRates;
        final Parameter rates;

        Fixture(PartialsRescalingScheme scheme) {

            SimpleAlignment alignment = new SimpleAlignment();
            alignment.setDataType(Nucleotides.INSTANCE);

            final String states = "ACGT";
            SimpleNode node = tip(alignment, "ladder0", "A");
            double height = 0.0;
            for (int i = 1; i < LADDER_TIPS; ++i) {
                height += 0.1;
                node = parent(node, tip(alignment, "ladder" + i, states.substring(i % 4, i % 4 + 1)), height);
            }
            for (int i = 0; i < CHERRIES; ++i) {
                SimpleNode cherry = parent(tip(alignment, "cherry" + i + "a", "A"),
                        tip(alignment, "cherry" + i + "c", "C"), 0.1);
                height += 0.1;
                node = parent(node, cherry, height);
            }
            tree = new DefaultTreeModel(new SimpleTree(node));

            double[] logRates = new double[12];
            for (int k = 0; k < logRates.length; ++k) {
                logRates[k] = 0.1 * (k % 5) - 0.2;
            }
            LogRateSubstitutionModel model = new LogRateSubstitutionModel("model", Nucleotides.INSTANCE,
                    new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.3, 0.2, 0.2, 0.3}),
                    new LogAdditiveCtmcRateProvider.DataAugmented.Basic("logRates",
                            new Parameter.Default("logRates", logRates)));
            model.setNormalization(true);
            model.setScaleRatesByFrequencies(false);

            rates = new Parameter.Default(tree.getNodeCount() - 1, 1.0);
            branchRates = new ArbitraryBranchRates(tree, rates, ArbitraryBranchRates.make(false, false, false),
                    false);

            BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree,
                    new SitePatterns(alignment), new HomogeneousBranchModel(model),
                    new GammaSiteRateModel("siteModel", 0.5, 1), false, false, scheme, true,
                    new PreOrderSettings(true, false, false, false, true, false, false));

            likelihood = new TreeDataLikelihood(delegate, tree, branchRates);
            gradient = new LogCtmcRateGradient("test", likelihood, delegate, model, ApproximationMode.EXACT_SPECTRAL,
                    false);
        }

        /**
         * State A (stateB false): the cherries' branches are short. State B: the ladder's are.
         */
        void setBranchRates(boolean stateB) {
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                NodeRef node = tree.getNode(i);
                if (tree.isRoot(node)) {
                    continue;
                }
                final boolean ladder = inLadder(node);
                final boolean cherryTip = tree.isExternal(node) && !ladder;
                final double rate = (stateB ? ladder : cherryTip) ? SHORT : 1.0;
                rates.setParameterValue(branchRates.getParameterIndexFromNode(node), rate);
            }
        }

        private boolean inLadder(NodeRef node) {
            if (tree.isExternal(node)) {
                return tree.getNodeTaxon(node).getId().startsWith("ladder");
            }
            for (int c = 0; c < tree.getChildCount(node); ++c) {
                if (!inLadder(tree.getChild(node, c))) {
                    return false;
                }
            }
            return true;
        }

        private static SimpleNode tip(SimpleAlignment alignment, String name, String state) {
            Taxon taxon = new Taxon(name);
            alignment.addSequence(new Sequence(taxon, state));
            SimpleNode tip = new SimpleNode();
            tip.setTaxon(taxon);
            tip.setHeight(0.0);
            return tip;
        }

        private static SimpleNode parent(SimpleNode left, SimpleNode right, double height) {
            SimpleNode parent = new SimpleNode();
            parent.setHeight(height);
            parent.addChild(left);
            parent.addChild(right);
            return parent;
        }
    }
}
