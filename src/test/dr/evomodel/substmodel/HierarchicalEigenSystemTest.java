/*
 * HierarchicalEigenSystemTest.java
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301 USA
 *
 */

package test.dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.NewickImporter;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Both eigensystems represent the same generator and closed-form transition probabilities.
 */
public class HierarchicalEigenSystemTest extends TestCase {

    private static class Built {
        HierarchicalSubstitutionModel model;
        Parameter rates;
        Parameter frequencies;
    }

    private static Built build(String newick, DataType dataType, double[] pi, boolean normalized) throws Exception {
        TreeModel tree = new DefaultTreeModel("characterTree", new NewickImporter(newick).importTree(null));
        Built built = new Built();
        built.rates = new Parameter.Default(tree.getInternalNodeCount(), 1.0);
        built.frequencies = new Parameter.Default(pi.clone());
        FrequencyModel frequencyModel = new FrequencyModel(dataType, built.frequencies);
        built.model = new HierarchicalSubstitutionModel("hierarchical",
                new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(
                        built.rates, tree, frequencyModel, normalized));
        return built;
    }

    private static double entry(EigenDecomposition eigen, int n, int i, int j, double t, boolean generator) {
        double result = 0.0;
        for (int k = 0; k < n; k++) {
            double lambda = eigen.getEigenValues()[k];
            result += eigen.getEigenVectors()[i * n + k] * (generator ? lambda : Math.exp(t * lambda))
                    * eigen.getInverseEigenVectors()[k * n + j];
        }
        return result;
    }

    private static void assertReconstruction(HierarchicalSubstitutionModel model, double t, boolean full) {
        int n = model.getDataType().getStateCount();
        double[] q = new double[n * n], p = new double[n * n];
        model.getInfinitesimalMatrix(q);
        model.getTransitionProbabilities(t, p);
        EigenDecomposition eigen = model.getEigenDecomposition();
        for (int i = 0; i < n; i++) {
            int[] columns;
            if (full) {
                columns = new int[n];
                for (int j = 0; j < n; j++) columns[j] = j;
            } else {
                columns = new int[]{i, (i + 1) % n, (37 * i + 17) % n};
            }
            for (int j : columns) {
                String where = model.getEigenSystemType() + " S=" + n + " [" + i + "," + j + "]";
                assertEquals(where + " Q", q[i * n + j], entry(eigen, n, i, j, t, true), 2e-11);
                assertEquals(where + " P at " + t, p[i * n + j], entry(eigen, n, i, j, t, false), 2e-11);
            }
        }
        double[] after = new double[n * n];
        model.getInfinitesimalMatrix(after);
        assertTrue("decomposition must not mutate the generator", Arrays.equals(q, after));
    }

    public void testBothModesForStateOrderLadderAndRepeatedEigenvalues() throws Exception {
        String[] trees = {"((A:1,G:1):1,(C:1,T:1):1);", "((T:1,C:1):1,(G:1,A:1):1);",
                "(((A:1,C:1):1,G:1):1,T:1);", "(T:1,G:1,C:1,A:1);"};
        for (boolean normalized : new boolean[]{false, true}) {
            for (String tree : trees) {
                Built built = build(tree, Nucleotides.INSTANCE, new double[]{0.1, 0.2, 0.3, 0.4}, normalized);
                for (HierarchicalSubstitutionModel.EigenSystemType type : HierarchicalSubstitutionModel.EigenSystemType.values()) {
                    built.model.setEigenSystemType(type);
                    for (double t : new double[]{0.0, 0.01, 0.03, 0.1, 1.9}) {
                        assertReconstruction(built.model, t, true);
                    }
                }
            }
        }
    }

    private static String balanced(String[] codes, int from, int to) {
        if (to - from == 1) return codes[from] + ":1";
        int middle = from + (to - from) / 2;
        return "(" + balanced(codes, from, middle) + "," + balanced(codes, middle, to) + "):1";
    }

    public void testBenchmarkStateCountsPreserveExactUniformFrequencies() throws Exception {
        for (int n : new int[]{4, 8, 16, 20, 32, 61, 128, 256}) {
            String[] codes = new String[n];
            for (int i = 0; i < n; i++) codes[i] = "state" + i;
            double[] pi = new double[n];
            Arrays.fill(pi, 1.0 / n);
            for (boolean normalized : n >= 128 ? new boolean[]{false, true} : new boolean[]{true}) {
                Built built = build(balanced(codes, 0, n) + ";", new GeneralDataType(codes), pi, normalized);
                for (HierarchicalSubstitutionModel.EigenSystemType type : HierarchicalSubstitutionModel.EigenSystemType.values()) {
                    built.model.setEigenSystemType(type);
                    // Full reconstructions are cubic; sample entries for the larger benchmark state spaces.
                    assertReconstruction(built.model, 0.1, n <= 32);
                    for (int i = 0; i < n; i++) {
                        assertEquals(type + " S=" + n + " normalized=" + normalized + " frequency " + i,
                                1.0 / n, built.frequencies.getParameterValue(i), 0.0);
                    }
                }
            }
        }
    }

    public void testCacheInvalidatedByModeRatesFrequenciesAndRestore() throws Exception {
        Built built = build("((A:1,G:1):1,(C:1,T:1):1);", Nucleotides.INSTANCE,
                new double[]{0.1, 0.2, 0.3, 0.4}, true);
        HierarchicalSubstitutionModel model = built.model;
        EigenDecomposition analytical = model.getEigenDecomposition();
        assertSame(analytical, model.getEigenDecomposition());
        model.setEigenSystemType(HierarchicalSubstitutionModel.EigenSystemType.NUMERICAL);
        EigenDecomposition numerical = model.getEigenDecomposition();
        assertNotSame("changing the eigen method rebuilds the cache", analytical, numerical);
        assertSame(numerical, model.getEigenDecomposition());

        for (HierarchicalSubstitutionModel.EigenSystemType type : HierarchicalSubstitutionModel.EigenSystemType.values()) {
            model.setEigenSystemType(type);
            EigenDecomposition before = model.getEigenDecomposition();
            built.rates.setParameterValue(0, built.rates.getParameterValue(0) + 0.7);
            EigenDecomposition rateChanged = model.getEigenDecomposition();
            assertNotSame("rate change in " + type, before, rateChanged);
            assertSame(rateChanged, model.getEigenDecomposition());
            built.frequencies.setParameterValue(0, built.frequencies.getParameterValue(0) + 0.01);
            built.frequencies.setParameterValue(3, built.frequencies.getParameterValue(3) - 0.01);
            EigenDecomposition frequencyChanged = model.getEigenDecomposition();
            assertNotSame("frequency change in " + type, rateChanged, frequencyChanged);
            assertSame(frequencyChanged, model.getEigenDecomposition());
            assertReconstruction(model, 0.1, true);

            double[] storedQ = new double[16];
            model.getInfinitesimalMatrix(storedQ);
            model.storeModelState();
            built.rates.setParameterValue(0, built.rates.getParameterValue(0) + 1.3);
            built.frequencies.setParameterValue(0, built.frequencies.getParameterValue(0) + 0.01);
            built.frequencies.setParameterValue(3, built.frequencies.getParameterValue(3) - 0.01);
            EigenDecomposition proposed = model.getEigenDecomposition();
            model.restoreModelState();
            EigenDecomposition restored = model.getEigenDecomposition();
            assertNotSame("restore invalidates the proposal in " + type, proposed, restored);
            assertSame(restored, model.getEigenDecomposition());
            double[] restoredQ = new double[16];
            model.getInfinitesimalMatrix(restoredQ);
            assertTrue("restore recovers the generator in " + type, Arrays.equals(storedQ, restoredQ));
            assertReconstruction(model, 0.1, true);
        }
    }
}
