/*
 * CharacterTreePriorTest.java
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

package test.dr.evomodel.substmodel;

import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.operators.CharacterTreeKeepFlip;
import dr.evomodel.operators.CharacterTreeNNI;
import dr.evomodel.operators.CharacterTreeSPR;
import dr.evomodel.substmodel.CharacterTreeModel;
import dr.evomodel.substmodel.CharacterTreePrior;
import dr.evomodel.substmodel.CharacterTreeStatistic;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The topology prior of a sampled character tree, and prior-only chains of the character-tree moves.
 *
 * @author Xiang Ji
 */
public class CharacterTreePriorTest extends TestCase {

    public void testTreeCounts() {
        assertEquals(4.0, Math.exp(logSum(CharacterTreePrior.logTreesByQ(3))), 1e-9);
        assertEquals(26.0, Math.exp(logSum(CharacterTreePrior.logTreesByQ(4))), 1e-9);
        assertEquals(236.0, Math.exp(logSum(CharacterTreePrior.logTreesByQ(5))), 1e-9);
        assertEquals(2752.0, Math.exp(logSum(CharacterTreePrior.logTreesByQ(6))), 1e-8);
        assertEquals(12818912.0, Math.exp(logSum(CharacterTreePrior.logTreesByQ(9))), 1e-3);
        double[] s6 = CharacterTreePrior.logTreesByQ(6);
        double[] expected = {1, 56, 490, 1260, 945};
        for (int q = 0; q < 5; q++) {
            assertEquals(expected[q], Math.exp(s6[q]), 1e-9 * expected[q]);
        }
        double[] s61 = CharacterTreePrior.logTreesByQ(61);
        assertEquals("N(61) = 2.964e105", 2.964, Math.exp(logSum(s61) - 105 * Math.log(10)), 5e-4);
    }

    private static double logSum(double[] values) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : values) max = Math.max(max, v);
        double sum = 0.0;
        for (double v : values) sum += Math.exp(v - max);
        return max + Math.log(sum);
    }

    public void testPriorValues() {
        Parameter keep = new Parameter.Default(1);
        CharacterTreeModel tree = new CharacterTreeModel("characterTree", "((A,G),(C,T));", Nucleotides.INSTANCE, true, keep);
        CharacterTreePrior uniform = new CharacterTreePrior(tree, CharacterTreePrior.Type.UNIFORM, 1.0);
        CharacterTreePrior classes = new CharacterTreePrior(tree, CharacterTreePrior.Type.RESOLUTION_CLASS, 1.0);
        CharacterTreePrior polytomy = new CharacterTreePrior(tree, CharacterTreePrior.Type.POLYTOMY, 2.0);
        assertEquals(-Math.log(26), uniform.getLogLikelihood(), 1e-14);
        assertEquals(-Math.log(3 * 15), classes.getLogLikelihood(), 1e-14);
        assertEquals(Math.log(0.3846153846153846 / 15), polytomy.getLogLikelihood(), 1e-14);
        for (int i = 0; i < 3; i++) {
            if (tree.getSkeleton().getRoot().getNumber() != 4 + i) keep.setParameterValue(i, 0.0);
        }
        assertEquals("(A,C,G,T);", tree.getContractedNewick());
        assertEquals("star: 1/26 shared by its 15 resolutions", -Math.log(26) - Math.log(15), uniform.getLogLikelihood(), 1e-14);
        assertEquals(-Math.log(3) - Math.log(15), classes.getLogLikelihood(), 1e-14);
        assertEquals(Math.log(0.10256410256410256 / 15), polytomy.getLogLikelihood(), 1e-14);
        double sum = 0.0;
        for (String t : CharacterTreeStatistic.treeList(Nucleotides.INSTANCE)) {
            CharacterTreeModel m = new CharacterTreeModel("t", t, Nucleotides.INSTANCE, true, new Parameter.Default(1));
            sum += Math.exp(new CharacterTreePrior(m, CharacterTreePrior.Type.POLYTOMY, 2.0).logTreeProbability(
                    m.getContractedInternalNodeCount() - 1));
        }
        assertEquals(1.0, sum, 1e-14);
    }

    private static final class Chain {
        final CharacterTreeModel tree;
        final Parameter rates;
        final ExponentialDistributionModel slab = new ExponentialDistributionModel(new Parameter.Default(1, 1.0));
        final CharacterTreePrior prior;
        final SimpleMCMCOperator[] operators;
        final int rootSlot;

        Chain(CharacterTreeModel tree, CharacterTreePrior.Type type, boolean spr) {
            this.tree = tree;
            final int internal = tree.getSkeleton().getInternalNodeCount();
            this.rates = new Parameter.Default(internal, 1.0);
            this.prior = new CharacterTreePrior(tree, type, 1.0);
            this.rootSlot = tree.getSkeleton().getRoot().getNumber() - tree.getExternalNodeCount();
            SimpleMCMCOperator topology = spr ? new CharacterTreeSPR(tree, 1.0) : new CharacterTreeNNI(tree, 1.0);
            this.operators = new SimpleMCMCOperator[]{topology, new CharacterTreeKeepFlip(tree, rates, slab, 0.5, 0.1, 1.0)};
        }

        double logPosterior() {
            double sum = prior.getLogLikelihood();
            for (int i = 0; i < rates.getDimension(); i++) {
                if (i != rootSlot) sum += slab.logPdf(rates.getParameterValue(i));
            }
            return sum;
        }

        Map<String, Integer> run(int iterations, int thin) {
            Map<String, Integer> counts = new HashMap<String, Integer>();
            double current = logPosterior();
            for (int it = 0; it < iterations; it++) {
                tree.storeModelState();
                rates.storeParameterValues();
                final int which = MathUtils.nextInt(3);
                double hastings;
                if (which < 2) {
                    hastings = operators[which].doOperation();
                } else {
                    int slot;
                    do {
                        slot = MathUtils.nextInt(rates.getDimension());
                    } while (slot == rootSlot);
                    final double s = 0.75 + MathUtils.nextDouble() * (1.0 / 0.75 - 0.75);
                    rates.setParameterValue(slot, rates.getParameterValue(slot) * s);
                    hastings = -Math.log(s);
                }
                final double proposed = logPosterior();
                if (Math.log(MathUtils.nextDouble()) < proposed - current + hastings) {
                    current = proposed;
                    tree.acceptModelState();
                    rates.acceptParameterValues();
                } else {
                    tree.restoreModelState();
                    rates.restoreParameterValues();
                }
                if (it % thin == 0) {
                    final String t = tree.getContractedNewick();
                    final Integer c = counts.get(t);
                    counts.put(t, c == null ? 1 : c + 1);
                }
            }
            return counts;
        }
    }

    private static double chiSquare(Map<String, Integer> counts, Map<String, Double> target) {
        double n = 0;
        for (int c : counts.values()) n += c;
        double stat = 0.0;
        for (Map.Entry<String, Double> e : target.entrySet()) {
            final Integer c = counts.get(e.getKey());
            final double expected = n * e.getValue();
            final double observed = c == null ? 0 : c;
            stat += (observed - expected) * (observed - expected) / expected;
        }
        return stat;
    }

    public void testPriorOnlyChainsAtFourStates() {
        MathUtils.setSeed(20261114);
        List<String> trees = CharacterTreeStatistic.treeList(Nucleotides.INSTANCE);
        for (CharacterTreePrior.Type type : new CharacterTreePrior.Type[]{CharacterTreePrior.Type.UNIFORM,
                CharacterTreePrior.Type.RESOLUTION_CLASS}) {
            CharacterTreeModel tree = new CharacterTreeModel("characterTree", "((A,G),(C,T));", Nucleotides.INSTANCE, true,
                    new Parameter.Default(1));
            Chain chain = new Chain(tree, type, false);
            Map<String, Integer> counts = chain.run(600000, 5);
            Map<String, Double> target = new HashMap<String, Double>();
            for (String t : trees) {
                CharacterTreeModel m = new CharacterTreeModel("t", t, Nucleotides.INSTANCE, true, new Parameter.Default(1));
                target.put(t, Math.exp(chain.prior.logTreeProbability(m.getContractedInternalNodeCount() - 1)));
            }
            assertEquals(type + " visits every tree", 26, counts.size());
            final double stat = chiSquare(counts, target);
            assertTrue(type + " chi-square " + stat + " on 25 df (0.001 critical value 52.6)", stat < 52.6);
        }
    }

    public void testPriorOnlyChainWithSubtreePruneRegraft() {
        MathUtils.setSeed(1114);
        String[] codes = {"a", "b", "c", "d", "e"};
        GeneralDataType type = new GeneralDataType(codes);
        CharacterTreeModel tree = new CharacterTreeModel("characterTree", "((a,b),((c,d),e));", type, true, new Parameter.Default(1));
        Chain chain = new Chain(tree, CharacterTreePrior.Type.UNIFORM, true);
        Map<String, Integer> counts = chain.run(1500000, 5);
        Map<String, Double> target = new HashMap<String, Double>();
        for (String t : CharacterTreeStatistic.treeList(type)) {
            target.put(t, 1.0 / 236);
        }
        assertEquals(236, counts.size());
        final double stat = chiSquare(counts, target);
        assertTrue("chi-square " + stat + " on 235 df (0.001 critical value 309.0)", stat < 309.0);
    }
}
