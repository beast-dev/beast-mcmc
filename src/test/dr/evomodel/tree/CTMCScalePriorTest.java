/*
 * CTMCScalePriorTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.evomodel.tree;

import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.operators.FixedHeightSubtreePruneRegraftOperator;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.tree.CTMCScalePrior;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.GammaFunction;
import dr.math.MathUtils;
import junit.framework.TestCase;

/**
 * Checks that CTMCScalePrior, which keeps the tree length up to date from tree changed events, gives
 * exactly the same log density as calculating the tree length in full, through store/accept/restore
 * cycles of tree moves (including root changes and moves of all heights) and moves of the scale.
 */
public class CTMCScalePriorTest extends TestCase {

    public CTMCScalePriorTest(String name) {
        super(name);
    }

    public void testIncrementalTreeLength() {
        MathUtils.setSeed(1234);

        Taxa taxa = new Taxa();
        for (int i = 0; i < 40; i++) {
            Taxon taxon = new Taxon("t" + i);
            // a few tips share a date
            taxon.setDate(new Date(MathUtils.nextInt(10) * 0.1, Units.Type.YEARS, true));
            taxa.addTaxon(taxon);
        }
        ConstantPopulation population = new ConstantPopulation(Units.Type.YEARS);
        population.setN0(1.0);
        Tree start = new CoalescentSimulator().simulateTree(taxa, population);
        DefaultTreeModel tree = new DefaultTreeModel("treeModel", start);

        Parameter scale = new Parameter.Default(0.001);
        scale.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        CTMCScalePrior prior = new CTMCScalePrior("ctmcScalePrior", scale, tree, null, false, null, false);

        Parameter internalHeights = tree.createNodeHeightsParameter(false, true, false);
        Parameter allInternalHeights = tree.createNodeHeightsParameter(true, true, false);
        SimpleMCMCOperator[] operators = {
                new SubtreeLeapOperator(tree, 1.0, 0.5, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                        AdaptationMode.ADAPTATION_OFF, 0.234),
                new FixedHeightSubtreePruneRegraftOperator(tree, 1.0),
                new UniformOperator(internalHeights, 1.0),
                new ScaleOperator(tree.getRootHeightParameter(), 0.75, AdaptationMode.ADAPTATION_OFF, 1.0),
                new UpDownOperator(new Scalable[]{new Scalable.Default(allInternalHeights)},
                        new Scalable[]{new Scalable.Default(scale)}, 0.75, 1.0, AdaptationMode.ADAPTATION_OFF),
                new ScaleOperator(scale, 0.75, AdaptationMode.ADAPTATION_OFF, 1.0)
        };

        assertEquals(expected(tree, scale), prior.getLogLikelihood(), 0.0);

        for (int iteration = 0; iteration < 20000; iteration++) {
            tree.storeModelState();
            scale.storeParameterValues();
            prior.storeModelState();

            // usually one move, sometimes several; occasionally the prior is made dirty
            int moves = MathUtils.nextInt(4) == 0 ? 2 + MathUtils.nextInt(3) : 1;
            for (int m = 0; m < moves; m++) {
                try {
                    operators[MathUtils.nextInt(operators.length)].doOperation();
                } catch (RuntimeException e) {
                    // a move that cannot be made: the state is restored below
                }
                if (MathUtils.nextInt(50) == 0) {
                    prior.makeDirty();
                }
                // sometimes the prior is not evaluated until after several moves
                if (MathUtils.nextBoolean()) {
                    assertEquals("proposed " + iteration, expected(tree, scale), prior.getLogLikelihood(), 0.0);
                }
            }
            assertEquals("proposed " + iteration, expected(tree, scale), prior.getLogLikelihood(), 0.0);

            if (MathUtils.nextBoolean()) {
                tree.acceptModelState();
                scale.acceptParameterValues();
                prior.acceptModelState();
            } else {
                tree.restoreModelState();
                scale.restoreParameterValues();
                prior.restoreModelState();
                assertEquals("restored " + iteration, expected(tree, scale), prior.getLogLikelihood(), 0.0);
            }
        }
    }

    /**
     * The CTMC scale reference prior, calculated from the full tree length.
     */
    private static double expected(Tree tree, Parameter scale) {
        double treeLength = TreeUtils.getTreeLength(tree);
        // as in CTMCScalePrior (one tree, proportion 1)
        treeLength /= 1;
        treeLength *= 1.0;
        double logNormalization = 0.5 * Math.log(treeLength) - GammaFunction.lnGamma(0.5);
        double logLike = 0;
        for (int i = 0; i < scale.getDimension(); ++i) {
            double ab = scale.getParameterValue(i);
            logLike += logNormalization - 0.5 * Math.log(ab) - ab * treeLength;
        }
        return logLike;
    }
}
