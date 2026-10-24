/*
 * LocalClockGradientTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package test.dr.evomodel.branchratemodel;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.LocalClockModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Gradients of the fixed local clock wrt its global and local rates (plan 3, E4): the chain rule from branch rates to
 * clock parameters against finite differences of getBranchRate, and the full gradient through a tree likelihood.
 */
public class LocalClockGradientTest extends TestCase {

    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";
    private static final String[][] SEQUENCES = {
            {"A", "ACGTACGTTACGATCGATCGTAGCTAGCTAGCATCG"},
            {"B", "ACGTACGTAACGATCGTTCGTAGCTAGCTAGCATCA"},
            {"C", "ACGAACGTTACGTTCGATCGTAGGTAGCTAGGATCG"},
            {"D", "TCGTACGATACGATCGATCCTAGCTAGCAAGCATCG"},
            {"E", "TCGTACGATACGATCGATCCTAGCTTGCAAGCATGG"},
    };
    private static final double H = 1E-6;

    public LocalClockGradientTest(String name) {
        super(name);
    }

    public void testRelativeCladeClock() throws Exception {
        TreeModel tree = tree();
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "C", "D", "E"), new Parameter.Default(2.0), true, null, false, false);
        checkChainRule(tree, clock, "relative clade");
    }

    public void testAbsoluteCladeClockWithStem() throws Exception {
        TreeModel tree = tree();
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "C", "D", "E"), new Parameter.Default(1.6), false, new Parameter.Default(0.3),
                false, false);
        checkChainRule(tree, clock, "absolute clade with stem");
    }

    public void testStemAsTime() throws Exception {
        TreeModel tree = tree();
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "C", "D", "E"), new Parameter.Default(2.0), true, new Parameter.Default(0.05),
                true, false);
        checkChainRule(tree, clock, "stem as time");
    }

    public void testNestedCladesAndTips() throws Exception {
        TreeModel tree = tree();
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "C", "D", "E"), new Parameter.Default(2.0), true, new Parameter.Default(0.4),
                false, false);
        clock.addCladeClock(taxa(tree, "D", "E"), new Parameter.Default(0.5), false, new Parameter.Default(0.2),
                false, false);
        clock.addExternalBranchClock(taxa(tree, "A"), new Parameter.Default(3.0), true);
        checkChainRule(tree, clock, "nested clades and a tip");
    }

    public void testExcludedCladeAndTrunk() throws Exception {
        TreeModel tree = tree();
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "A", "B"), new Parameter.Default(1.8), true, null, false, true);
        clock.addTrunkClock(taxa(tree, "D"), new Parameter.Default(0.4), null, true);
        checkChainRule(tree, clock, "excluded clade and trunk");
    }

    public void testSharedRateParameter() throws Exception {
        TreeModel tree = tree();
        Parameter shared = new Parameter.Default(2.0);
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "A", "B"), shared, true, null, false, false);
        clock.addCladeClock(taxa(tree, "D", "E"), shared, true, null, false, false);
        assertEquals(2, clock.getRateParameter().getDimension());
        checkChainRule(tree, clock, "shared rate");
    }

    public void testThroughTreeLikelihood() throws Exception {
        TreeModel tree = tree();
        LocalClockModel clock = new LocalClockModel(tree, new Parameter.Default(0.7));
        clock.addCladeClock(taxa(tree, "C", "D", "E"), new Parameter.Default(2.0), true, new Parameter.Default(0.4),
                false, false);
        clock.addExternalBranchClock(taxa(tree, "A"), new Parameter.Default(3.0), false);

        TreeDataLikelihood likelihood = likelihood(tree, clock);
        Parameter rates = clock.getRateParameter();
        double[] analytic = new BranchRateGradientForDiscreteTrait(null, likelihood,
                (BeagleDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate(), rates, false)
                .getGradientLogDensity();
        assertEquals(rates.getDimension(), analytic.length);
        for (int k = 0; k < rates.getDimension(); k++) {
            final double value = rates.getParameterValue(k);
            rates.setParameterValue(k, value + H);
            final double plus = likelihood.getLogLikelihood();
            rates.setParameterValue(k, value - H);
            final double minus = likelihood.getLogLikelihood();
            rates.setParameterValue(k, value);
            final double numeric = (plus - minus) / (2 * H);
            assertEquals("rate " + k, numeric, analytic[k], 1E-6 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    // sum_e w_e d rate_e / d theta_k against finite differences, with arbitrary branch weights w_e
    private static void checkChainRule(TreeModel tree, LocalClockModel clock, String label) {
        Parameter rates = clock.getRateParameter();
        double[] weights = new double[tree.getNodeCount() - 1];
        for (int v = 0; v < weights.length; v++) {
            weights[v] = 1.0 + 0.37 * v;
        }
        double[] analytic = clock.updateGradientLogDensity(weights, null, 0, weights.length);
        assertEquals(rates.getDimension(), analytic.length);
        for (int k = 0; k < rates.getDimension(); k++) {
            final double value = rates.getParameterValue(k);
            rates.setParameterValue(k, value + H);
            final double plus = weightedRates(tree, clock, weights);
            rates.setParameterValue(k, value - H);
            final double minus = weightedRates(tree, clock, weights);
            rates.setParameterValue(k, value);
            final double numeric = (plus - minus) / (2 * H);
            assertEquals(label + " rate " + k, numeric, analytic[k], 1E-7 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static double weightedRates(TreeModel tree, LocalClockModel clock, double[] weights) {
        double sum = 0.0;
        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                sum += weights[v++] * clock.getBranchRate(tree, node);
            }
        }
        return sum;
    }

    private static TreeDataLikelihood likelihood(TreeModel tree, LocalClockModel clock) {
        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);
        for (String[] row : SEQUENCES) {
            alignment.addSequence(new Sequence(tree.getTaxon(tree.getTaxonIndex(row[0])), row[1]));
        }
        HKY hky = new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.1, 0.2, 0.3, 0.4}));
        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree,
                new SitePatterns(alignment, null, 0, -1, 1, true), new HomogeneousBranchModel(hky),
                new GammaSiteRateModel("siteModel"), true, false, PartialsRescalingScheme.ALWAYS, false,
                new PreOrderSettings(true, false, true, true));
        return new TreeDataLikelihood(delegate, tree, clock);
    }

    private static Taxa taxa(TreeModel tree, String... names) {
        Taxa taxa = new Taxa();
        for (String name : names) {
            taxa.addTaxon(tree.getTaxon(tree.getTaxonIndex(name)));
        }
        return taxa;
    }

    private static TreeModel tree() throws Exception {
        return new DefaultTreeModel("treeModel", new NewickImporter(NEWICK).importTree(null));
    }
}
