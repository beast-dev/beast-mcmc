/*
 * HierarchicalSubstitutionModelTest.java
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

import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Xiang Ji
 */
public class HierarchicalSubstitutionModelTest extends TestCase {

    private static final int S = 4;

    private static class Built {
        HierarchicalSubstitutionModel model;
        Parameter rates;
        Parameter frequencies;
        TreeModel tree;
    }

    private static String tipSet(Tree tree, NodeRef node) {
        if (tree.isExternal(node)) {
            return tree.getTaxonId(node.getNumber());
        }
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tree.getChildCount(node); i++) {
            sb.append(tipSet(tree, tree.getChild(node, i)));
        }
        char[] chars = sb.toString().toCharArray();
        java.util.Arrays.sort(chars);
        return new String(chars);
    }

    private static Built build(String newick, Map<String, Double> rateByClade, double[] pi) throws Exception {
        Tree imported = new NewickImporter(newick).importTree(null);
        TreeModel tree = new DefaultTreeModel("characterTree", imported);
        Parameter rates = new Parameter.Default(tree.getInternalNodeCount(), 0.0);
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(S + i);
            Double rate = rateByClade.get(tipSet(tree, node));
            rates.setParameterValue(i, rate == null ? 0.0 : rate);
        }
        Parameter frequencies = new Parameter.Default(pi);
        FrequencyModel frequencyModel = new FrequencyModel(Nucleotides.INSTANCE, frequencies);
        Built built = new Built();
        built.model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, tree), frequencyModel);
        built.rates = rates;
        built.frequencies = frequencies;
        built.tree = tree;
        return built;
    }

    private static Map<String, Double> agct(double root, double ag, double ct) {
        Map<String, Double> map = new HashMap<>();
        map.put("ACGT", root);
        map.put("AG", ag);
        map.put("CT", ct);
        return map;
    }

    private static double[] tn93(double[] pi, double lambdaRoot, double lambdaAG, double lambdaCT, double t) {
        final double piR = pi[0] + pi[2], piY = pi[1] + pi[3];
        final double beta = lambdaRoot;
        final double[] classPi = {piR, piY};
        final double[] classAlpha = {lambdaRoot + lambdaAG / piR, lambdaRoot + lambdaCT / piY};
        final int[] cls = {0, 1, 0, 1};
        double[] p = new double[S * S];
        for (int i = 0; i < S; i++) {
            for (int j = 0; j < S; j++) {
                if (cls[i] != cls[j]) {
                    p[i * S + j] = pi[j] * (1.0 - Math.exp(-beta * t));
                } else {
                    int c = cls[i];
                    double piC = classPi[c], piO = classPi[1 - c];
                    double rateC = piC * classAlpha[c] + piO * beta;
                    p[i * S + j] = pi[j] + pi[j] * (piO / piC) * Math.exp(-beta * t)
                            + ((i == j ? 1.0 : 0.0) - pi[j] / piC) * Math.exp(-rateC * t);
                }
            }
        }
        return p;
    }

    private static double[] fromEigen(EigenDecomposition ed, double t) {
        double[] evec = ed.getEigenVectors(), ievc = ed.getInverseEigenVectors(), eval = ed.getEigenValues();
        double[] p = new double[S * S];
        for (int i = 0; i < S; i++)
            for (int j = 0; j < S; j++) {
                double sum = 0.0;
                for (int k = 0; k < S; k++)
                    sum += evec[i * S + k] * Math.exp(eval[k] * t) * ievc[k * S + j];
                p[i * S + j] = sum;
            }
        return p;
    }

    private static double maxAbsDiff(double[] a, double[] b) {
        double d = 0.0;
        for (int i = 0; i < a.length; i++)
            d = Math.max(d, Math.abs(a[i] - b[i]));
        return d;
    }

    public void testWorkedExample() throws Exception {
        Built b = build("((A:1,G:1):1,(C:1,T:1):1);", agct(4, 2, 3), new double[]{0.25, 0.25, 0.25, 0.25});
        double[] p = new double[S * S];
        b.model.getTransitionProbabilities(0.1, p);
        final double d = 0.691986, dCT = 0.665873, ag = 0.143174, ct = 0.169287, tv = 0.082420;
        double[] expected = {d, tv, ag, tv, tv, dCT, tv, ct, ag, tv, d, tv, tv, ct, tv, dCT};
        assertEquals(0.0, maxAbsDiff(p, expected), 1e-6);
    }

    public void testTN93AndStateOrder() throws Exception {
        double[] pi = {0.1, 0.2, 0.3, 0.4};
        for (String newick : new String[]{"((A:1,G:1):1,(C:1,T:1):1);", "((T:1,C:1):1,(G:1,A:1):1);"}) {
            Built b = build(newick, agct(1.3, 0.7, 2.1), pi);
            double[] p = new double[S * S];
            b.model.getTransitionProbabilities(0.37, p);
            assertEquals(newick, 0.0, maxAbsDiff(p, tn93(pi, 1.3, 0.7, 2.1, 0.37)), 1e-14);
        }
    }

    public void testEigenDecompositionAndGenerator() throws Exception {
        double[] pi = {0.1, 0.2, 0.3, 0.4};
        String[] trees = {"((A:1,G:1):1,(C:1,T:1):1);", "((T:1,C:1):1,(G:1,A:1):1);", "(A:1,C:1,G:1,T:1);",
                "(((A:1,C:1):1,G:1):1,T:1);"};
        for (String newick : trees) {
            Map<String, Double> rates = new HashMap<>();
            rates.put("ACGT", 1.1);
            rates.put("AG", 0.6);
            rates.put("CT", 1.7);
            rates.put("AC", 0.9);
            rates.put("ACG", 0.4);
            Built b = build(newick, rates, pi);
            for (double t : new double[]{0.0, 0.05, 0.8, 3.0}) {
                double[] p = new double[S * S];
                b.model.getTransitionProbabilities(t, p);
                assertEquals(newick + " t=" + t, 0.0, maxAbsDiff(p, fromEigen(b.model.getEigenDecomposition(), t)), 1e-13);
            }
            double[] q = new double[S * S];
            b.model.getInfinitesimalMatrix(q);
            EigenDecomposition ed = b.model.getEigenDecomposition();
            double[] qFromEigen = new double[S * S];
            for (int i = 0; i < S; i++)
                for (int j = 0; j < S; j++)
                    for (int k = 0; k < S; k++)
                        qFromEigen[i * S + j] += ed.getEigenVectors()[i * S + k] * ed.getEigenValues()[k]
                                * ed.getInverseEigenVectors()[k * S + j];
            assertEquals(newick + " Q", 0.0, maxAbsDiff(q, qFromEigen), 1e-13);
        }
    }

    public void testModelChangesPropagate() throws Exception {
        Built b = build("((A:1,G:1):1,(C:1,T:1):1);", agct(1.0, 1.0, 1.0), new double[]{0.25, 0.25, 0.25, 0.25});
        double[] before = b.model.getEigenDecomposition().getEigenValues().clone();
        final boolean[] fired = {false};
        b.model.addModelListener(new dr.inference.model.ModelListener() {
            public void modelChangedEvent(dr.inference.model.Model model, Object object, int index) {
                fired[0] = true;
            }
            public void modelRestored(dr.inference.model.Model model) {
            }
        });
        b.rates.setParameterValue(0, 5.0);
        assertTrue("rate change reaches listeners", fired[0]);
        fired[0] = false;
        b.frequencies.setParameterValue(0, 0.4);
        b.frequencies.setParameterValue(1, 0.1);
        assertTrue("frequency change reaches listeners", fired[0]);
        assertTrue("eigen system recomputed", maxAbsDiff(before, b.model.getEigenDecomposition().getEigenValues()) > 0.1);
    }

    public void testBundle() throws Exception {
        Built b = build("((A:1,G:1):1,(C:1,T:1):1);", agct(4, 2, 3), new double[]{0.25, 0.25, 0.25, 0.25});
        final int n = S - 1;
        double[] bundle = new double[b.model.getBundleLength(1)];
        assertEquals(1 + 6 * n + n + S, bundle.length);
        b.model.getBundleHeader(bundle);
        b.model.fillBundleBlock(0.1, bundle, b.model.getBundleHeaderLength());
        assertEquals(65536.0 + n, bundle[0]);

        boolean[] seen = new boolean[2 * S - 1];
        for (int t = 0; t < n; t++) {
            int p = (int) bundle[1 + 3 * t];
            assertTrue(p >= S && !seen[p]);
            for (int c = 1; c <= 2; c++) {
                int child = (int) bundle[1 + 3 * t + c];
                assertTrue("children before parents", child < S || seen[child]);
            }
            seen[p] = true;
            for (int c = 0; c < 3; c++)
                assertEquals(bundle[1 + 3 * t + c], bundle[1 + 3 * n + 3 * (n - 1 - t) + c]);
        }
        assertEquals((double) (2 * S - 2), bundle[1 + 3 * (n - 1)]);

        final int block = 1 + 6 * n;
        double[] hByParent = new double[2 * S - 1];
        for (int m = 0; m < n; m++)
            hByParent[S + m] = bundle[block + m];
        assertEquals(0.329680, hByParent[2 * S - 2], 1e-6);
        double[] internalH = {bundle[block], bundle[block + 1]};
        java.util.Arrays.sort(internalH);
        assertEquals(0.243017, internalH[0], 1e-6);
        assertEquals(0.347469, internalH[1], 1e-6);
        double[] prZ0 = {0.548812, 0.496585, 0.548812, 0.496585};
        for (int x = 0; x < S; x++)
            assertEquals(prZ0[x], bundle[block + n + x], 1e-6);
    }
}
