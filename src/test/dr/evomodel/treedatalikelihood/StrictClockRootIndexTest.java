/*
 * StrictClockRootIndexTest.java
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
package test.dr.evomodel.treedatalikelihood;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.BranchRateGradientForDiscreteTrait;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Branch-rate gradients when the root is not node 2n-2 (plan 3, B7): the gradient array is indexed by branch with
 * the root skipped, so the strict clock (whose index used to be the raw node number) neither overflows nor
 * double-counts, and per-branch rates keep their TreeParameterModel indices after the root moves.
 */
public class StrictClockRootIndexTest extends TestCase {

    private static final String[][] SEQUENCES = {
            {"A", "AAAAAACCCCGGGGGTTTTTTTAACAAAGTTACGARACCCGGTA"},
            {"B", "AAAAAACCCCGGGGGTTTTTTTAACGGCGTTATATAACCCGGCC"},
            {"C", "AAAAAACCCCGGGGGTTTTTTTGGTAAGGCCACGAA-CCCAATC"},
            {"D", "AAAAAACCCCGGGGGTTTTTTTGGTAATATTATATGGCCCGGCC"},
            {"E", "AAAAAACCCCGGGGGTTTTTTTGGTAAAATTGCGAGGTTTGGCA"},
    };
    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";

    public StrictClockRootIndexTest(String name) {
        super(name);
    }

    public void testStrictClockWithRootNotLast() throws Exception {
        for (boolean moveRoot : new boolean[]{false, true}) {
            Parameter clockRate = new Parameter.Default("clock.rate", 1.0);
            DefaultTreeModel tree = tree(moveRoot);
            assertEquals(moveRoot ? 7 : 8, tree.getRoot().getNumber());
            TreeDataLikelihood likelihood = likelihood(tree, new StrictClockBranchRates(clockRate));

            BranchRateGradientForDiscreteTrait gradient = new BranchRateGradientForDiscreteTrait("strict",
                    likelihood, (BeagleDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate(), clockRate, false);
            double[] analytic = gradient.getGradientLogDensity();
            assertEquals(1, analytic.length);
            assertEquals("root " + tree.getRoot().getNumber(), numeric(likelihood, clockRate, 0), analytic[0],
                    1E-6 * Math.max(1.0, Math.abs(analytic[0])));
        }
    }

    public void testArbitraryRatesAfterRootMove() throws Exception {
        DefaultTreeModel tree = tree(false);
        Parameter rates = new Parameter.Default("rates", tree.getNodeCount() - 1, 1.0);
        for (int i = 0; i < rates.getDimension(); i++) {
            rates.setParameterValue(i, 0.7 + 0.15 * i);
        }
        ArbitraryBranchRates branchRates = new ArbitraryBranchRates(tree, rates,
                new ArbitraryBranchRates.BranchRateTransform.None(), false);
        TreeDataLikelihood likelihood = likelihood(tree, branchRates);
        moveRoot(tree);

        BranchRateGradientForDiscreteTrait gradient = new BranchRateGradientForDiscreteTrait("arbitrary",
                likelihood, (BeagleDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate(), rates, false);
        double[] analytic = gradient.getGradientLogDensity();
        for (int i = 0; i < rates.getDimension(); i++) {
            assertEquals("rate " + i, numeric(likelihood, rates, i), analytic[i],
                    1E-6 * Math.max(1.0, Math.abs(analytic[i])));
        }
    }

    private static double numeric(TreeDataLikelihood likelihood, Parameter parameter, int i) {
        final double value = parameter.getParameterValue(i);
        final double h = 1E-6 * Math.max(1.0, value);
        parameter.setParameterValue(i, value + h);
        final double plus = likelihood.getLogLikelihood();
        parameter.setParameterValue(i, value - h);
        final double minus = likelihood.getLogLikelihood();
        parameter.setParameterValue(i, value);
        return (plus - minus) / (2 * h);
    }

    private static DefaultTreeModel tree(boolean moveRoot) throws Exception {
        Tree newick = new NewickImporter(NEWICK).importTree(null);
        DefaultTreeModel tree = new DefaultTreeModel("treeModel", newick);
        if (moveRoot) {
            moveRoot(tree);
        }
        return tree;
    }

    // ((A,B)5,(C,(D,E)6)7)8 becomes (((A,B)5,C)8,(D,E)6)7, so the root is node 7 rather than 2n-2 = 8
    private static void moveRoot(DefaultTreeModel tree) {
        NodeRef oldRoot = tree.getNode(8);
        NodeRef newRoot = tree.getNode(7);
        NodeRef c = tree.getNode(2);
        tree.beginTreeEdit();
        tree.removeChild(oldRoot, newRoot);
        tree.removeChild(newRoot, c);
        tree.addChild(oldRoot, c);
        tree.addChild(newRoot, oldRoot);
        tree.setRoot(newRoot);
        tree.setNodeHeight(oldRoot, 0.3);
        tree.setNodeHeight(newRoot, 0.43);
        tree.endTreeEdit();
    }

    private static TreeDataLikelihood likelihood(DefaultTreeModel tree, BranchRateModel branchRates) throws Exception {
        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);
        for (String[] row : SEQUENCES) {
            alignment.addSequence(new Sequence(new Taxon(row[0]), row[1]));
        }
        HKY hky = new HKY(new Parameter.Default(2.5), new FrequencyModel(Nucleotides.INSTANCE,
                new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4})));
        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree,
                new SitePatterns(alignment, null, 0, -1, 1, true), new HomogeneousBranchModel(hky),
                new GammaSiteRateModel("siteModel", 0.5, 4), true, false, PartialsRescalingScheme.NONE, false,
                new PreOrderSettings(true, true, false, true));
        return new TreeDataLikelihood(delegate, tree, branchRates);
    }
}
