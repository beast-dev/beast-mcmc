/*
 * HierarchicalEigenSystemBeagleTest.java
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

package test.dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.InstanceDetails;
import dr.evolution.datatype.GeneralDataType;
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
 * Upload both eigensystems to the same actual SSE instance and compare its pruning with
 * Java pruning using the model's closed-form transitions. Optional native tests skip when
 * SSE is unavailable; use -Dbeast.test.requireSSE=true to make that an error.
 */
public class HierarchicalEigenSystemBeagleTest extends TestCase {

    private static final double[] CATEGORY_RATES = {0.4, 1.8};
    private static final double[] CATEGORY_WEIGHTS = {0.35, 0.65};
    private static final double[] PATTERN_WEIGHTS = {3, 1, 7, 2, 1, 4, 5};
    private static final int[] MATRIX_INDICES = {0, 1, 2, 3, 4, 5};
    private static final int[] OPERATIONS = {
            4, Beagle.NONE, Beagle.NONE, 0, 0, 1, 1,
            5, Beagle.NONE, Beagle.NONE, 2, 2, 3, 3,
            6, Beagle.NONE, Beagle.NONE, 4, 4, 5, 5
    };
    private static final long REQUIRED_FLAGS = BeagleFlag.PROCESSOR_CPU.getMask()
            | BeagleFlag.PRECISION_DOUBLE.getMask() | BeagleFlag.VECTOR_SSE.getMask()
            | BeagleFlag.EIGEN_REAL.getMask();

    private static String balanced(String[] codes, int from, int to) {
        if (to - from == 1) return codes[from] + ":1";
        int middle = from + (to - from) / 2;
        return "(" + balanced(codes, from, middle) + "," + balanced(codes, middle, to) + "):1";
    }

    private static Beagle sseInstance(int n) {
        if (Boolean.getBoolean("java.only")) {
            if (Boolean.getBoolean("beast.test.requireSSE")) {
                fail("beast.test.requireSSE=true requires native BEAGLE; java.only=true selects the Java backend");
            }
            System.out.println("Skipping hierarchical eigensystem SSE test: java.only=true");
            return null;
        }
        Beagle beagle;
        try {
            beagle = BeagleFactory.loadBeagleInstance(4, 7, 4, n, PATTERN_WEIGHTS.length,
                    1, MATRIX_INDICES.length, CATEGORY_RATES.length, 0, null,
                    REQUIRED_FLAGS, REQUIRED_FLAGS);
        } catch (RuntimeException e) {
            if (Boolean.getBoolean("beast.test.requireSSE")) throw e;
            System.out.println("Skipping hierarchical eigensystem SSE test at S=" + n + ": " + e.getMessage());
            return null;
        }
        InstanceDetails details = beagle.getDetails();
        assertNotNull("native BEAGLE instance details", details);
        assertEquals("required CPU/double/SSE/real-eigen flags", REQUIRED_FLAGS,
                details.getFlags() & REQUIRED_FLAGS);
        assertEquals("dense eigen representation", 0L,
                details.getFlags() & (BeagleFlag.HIERARCHICAL_REPRESENTATION.getMask()
                        | BeagleFlag.SPECTRAL_REPRESENTATION.getMask()));
        String implementation = details.getImplementationName();
        if (implementation == null) implementation = "implementation name unavailable from this JNI library";
        System.out.println("Hierarchical eigensystem S=" + n + ": " + implementation
                + "; resource=" + details.getResourceNumber()
                + " (" + BeagleFlag.toString(details.getFlags()) + ")");
        return beagle;
    }

    private static int[][] states(int n) {
        return new int[][]{
                {0, 0, 1, 2, n - 1, 1, 0},
                {0, 1, 1, 2, n - 2, 2, n - 1},
                {0, 2, 0, 2, 1, 3, 0},
                {0, 3, 0, 3, 0, 1, n - 1}
        };
    }

    private static double[] reference(HierarchicalSubstitutionModel model, int[][] states, double branchLength) {
        int n = model.getDataType().getStateCount();
        double[] pi = model.getFrequencyModel().getFrequencies();
        double[] likelihood = new double[PATTERN_WEIGHTS.length];
        double[] transition = new double[n * n];
        double[] left = new double[n], right = new double[n];
        for (int category = 0; category < CATEGORY_RATES.length; category++) {
            model.getTransitionProbabilities(branchLength * CATEGORY_RATES[category], transition);
            for (int pattern = 0; pattern < likelihood.length; pattern++) {
                for (int i = 0; i < n; i++) {
                    left[i] = transition[i * n + states[0][pattern]]
                            * transition[i * n + states[1][pattern]];
                    right[i] = transition[i * n + states[2][pattern]]
                            * transition[i * n + states[3][pattern]];
                }
                double root = 0.0;
                for (int i = 0; i < n; i++) {
                    double leftSum = 0.0, rightSum = 0.0;
                    for (int j = 0; j < n; j++) {
                        leftSum += transition[i * n + j] * left[j];
                        rightSum += transition[i * n + j] * right[j];
                    }
                    root += pi[i] * leftSum * rightSum;
                }
                likelihood[pattern] += CATEGORY_WEIGHTS[category] * root;
            }
        }
        for (int pattern = 0; pattern < likelihood.length; pattern++) {
            likelihood[pattern] = Math.log(likelihood[pattern]);
        }
        return likelihood;
    }

    private static double weightedTotal(double[] siteLogLikelihoods) {
        double result = 0.0;
        for (int pattern = 0; pattern < siteLogLikelihoods.length; pattern++) {
            result += PATTERN_WEIGHTS[pattern] * siteLogLikelihoods[pattern];
        }
        return result;
    }

    public void testAnalyticalAndNumericalOnSameSSEBackend() throws Throwable {
        for (int n : new int[]{4, 8, 16, 20, 32, 61, 128, 256}) {
            Beagle beagle = sseInstance(n);
            if (beagle == null) continue;
            try {
                String[] codes = new String[n];
                for (int i = 0; i < n; i++) codes[i] = "state" + i;
                TreeModel tree = new DefaultTreeModel("characterTree",
                        new NewickImporter(balanced(codes, 0, n) + ";").importTree(null));
                Parameter rates = new Parameter.Default(tree.getInternalNodeCount(), 1.0);
                double[] pi = new double[n];
                Arrays.fill(pi, 1.0 / n);
                FrequencyModel frequency = new FrequencyModel(new GeneralDataType(codes), new Parameter.Default(pi));
                HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hierarchical",
                        new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(rates, tree, frequency, true));
                int[][] states = states(n);
                for (int tip = 0; tip < states.length; tip++) beagle.setTipStates(tip, states[tip]);
                beagle.setPatternWeights(PATTERN_WEIGHTS);
                beagle.setStateFrequencies(0, pi);
                beagle.setCategoryRates(CATEGORY_RATES);
                beagle.setCategoryWeights(0, CATEGORY_WEIGHTS);
                double maxTotalError = 0.0, maxSiteError = 0.0;

                // Reuse exactly the same buffers and native implementation across both modes and a model update.
                for (int update = 0; update < 2; update++) {
                    if (update == 1) {
                        rates.setParameterValue(0, 2.3);
                        rates.setParameterValue(tree.getRoot().getNumber() - n, 0.6);
                    }
                    for (double branchLength : new double[]{0.01, 0.03, 0.1}) {
                        double[] expectedSites = reference(model, states, branchLength);
                        double expected = weightedTotal(expectedSites);
                        double analytical = Double.NaN;
                        for (HierarchicalSubstitutionModel.EigenSystemType type
                                : HierarchicalSubstitutionModel.EigenSystemType.values()) {
                            model.setEigenSystemType(type);
                            EigenDecomposition eigen = model.getEigenDecomposition();
                            beagle.setEigenDecomposition(0, eigen.getEigenVectors(),
                                    eigen.getInverseEigenVectors(), eigen.getEigenValues());
                            double[] branchLengths = new double[MATRIX_INDICES.length];
                            Arrays.fill(branchLengths, branchLength);
                            beagle.updateTransitionMatrices(0, MATRIX_INDICES, null, null,
                                    branchLengths, MATRIX_INDICES.length);
                            beagle.updatePartials(OPERATIONS, 3, Beagle.NONE);
                            double[] observed = new double[1];
                            beagle.calculateRootLogLikelihoods(new int[]{6}, new int[]{0}, new int[]{0},
                                    new int[]{Beagle.NONE}, 1, observed);
                            String where = type + " S=" + n + " update=" + update + " b=" + branchLength;
                            assertEquals(where + " weighted log likelihood", expected, observed[0], 1e-9);
                            maxTotalError = Math.max(maxTotalError, Math.abs(expected - observed[0]));
                            double[] observedSites = new double[PATTERN_WEIGHTS.length];
                            beagle.getSiteLogLikelihoods(observedSites);
                            for (int pattern = 0; pattern < observedSites.length; pattern++) {
                                assertEquals(where + " pattern=" + pattern,
                                        expectedSites[pattern], observedSites[pattern], 1e-10);
                                maxSiteError = Math.max(maxSiteError,
                                        Math.abs(expectedSites[pattern] - observedSites[pattern]));
                            }
                            if (type == HierarchicalSubstitutionModel.EigenSystemType.ANALYTICAL) {
                                analytical = observed[0];
                            } else {
                                assertEquals(where + " vs analytical", analytical, observed[0], 1e-9);
                            }
                        }
                    }
                }
                System.out.println("S=" + n + " maximum absolute log-likelihood error: total="
                        + maxTotalError + ", pattern=" + maxSiteError);
            } finally {
                beagle.finalize();
            }
        }
    }
}
