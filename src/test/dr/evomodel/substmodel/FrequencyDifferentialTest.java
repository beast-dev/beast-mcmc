/*
 * FrequencyDifferentialTest.java
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

import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.HiddenNucleotides;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentiableSubstitutionModelUtil;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.MarkovModulatedSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.WrtFrequencyParameter;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.evomodel.substmodel.aminoacid.WAG;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Checks the generic derivative of the normalized generator Q~ = Q / beta with respect to stationary frequencies.
 */
public class FrequencyDifferentialTest extends TestCase {

    private static final int K = 4;
    private static final double[] PI = {0.1, 0.2, 0.3, 0.4};
    private static final double KAPPA = 2.7;

    private Parameter frequencies;
    private HKY hky;

    public FrequencyDifferentialTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        frequencies = new Parameter.Default(PI.clone());
        hky = new HKY(KAPPA, new FrequencyModel(Nucleotides.INSTANCE, frequencies));
    }

    public void testMatchesFiniteDifferenceOfGenerator() {
        final double h = 1E-6;
        for (int k = 0; k < K; k++) {
            double[][] analytic = toArray(hky.getInfinitesimalDifferentialMatrix(hky.factory(frequencies, k)));

            double[] plus = generatorWithFrequency(k, PI[k] + h);
            double[] minus = generatorWithFrequency(k, PI[k] - h);
            frequencies.setParameterValue(k, PI[k]);

            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    double numeric = (plus[i * K + j] - minus[i * K + j]) / (2 * h);
                    assertEquals("k=" + k + " (" + i + "," + j + ")", numeric, analytic[i][j], 1E-7);
                }
            }
        }
    }

    public void testFirstDimensionMatchesLegacyFreqA() throws Exception {
        @SuppressWarnings("unchecked")
        Class<? extends Enum> wrtClass = (Class<? extends Enum>)
                Class.forName("dr.evomodel.substmodel.nucleotide.HKY$WrtHKYModelParameter");
        WrtParameter freqA = (WrtParameter) Enum.valueOf(wrtClass, "FREQ_A");

        double[][] legacy = toArray(DifferentiableSubstitutionModelUtil.getInfinitesimalDifferentialMatrix(freqA, hky));
        double[][] generic = toArray(hky.getInfinitesimalDifferentialMatrix(hky.factory(frequencies, 0)));

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                assertEquals(legacy[i][j], generic[i][j], 1E-14);
            }
        }
    }

    public void testFrequencyWeightedSumIsMinusGenerator() {
        // Q is homogeneous of degree 1 and beta of degree 2 in pi, so sum_k pi_k dQ~/dpi_k = -Q~
        double[] q = new double[K * K];
        hky.getInfinitesimalMatrix(q);

        double[][] sum = new double[K][K];
        for (int k = 0; k < K; k++) {
            double[][] d = toArray(hky.getInfinitesimalDifferentialMatrix(hky.factory(frequencies, k)));
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    sum[i][j] += PI[k] * d[i][j];
                }
            }
        }
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                assertEquals(-q[i * K + j], sum[i][j], 1E-13);
            }
        }
    }

    public void testExactTransitionDifferentialMatchesFiniteDifference() {
        final double h = 1E-6;
        for (double time : new double[]{0.05, 0.7, 3.0}) {
            for (int k = 0; k < K; k++) {
                WrappedMatrix differential = hky.getInfinitesimalDifferentialMatrix(hky.factory(frequencies, k));
                double[] d = DifferentialMassProvider.Mode.EXACT.dispatch(time, hky, differential);

                double[] p = expm(generatorWithFrequency(k, PI[k]), time);
                double[] plus = expm(generatorWithFrequency(k, PI[k] + h), time);
                double[] minus = expm(generatorWithFrequency(k, PI[k] - h), time);
                frequencies.setParameterValue(k, PI[k]);

                for (int i = 0; i < K; i++) {
                    for (int j = 0; j < K; j++) {
                        double pd = 0.0;
                        for (int l = 0; l < K; l++) {
                            pd += p[i * K + l] * d[l * K + j];
                        }
                        double numeric = (plus[i * K + j] - minus[i * K + j]) / (2 * h);
                        assertEquals("t=" + time + " k=" + k + " (" + i + "," + j + ")", numeric, pd, 1E-7);
                    }
                }
            }
        }
    }

    public void testFactoryRange() {
        for (int k = 0; k < K; k++) {
            WrtParameter wrt = hky.factory(frequencies, k);
            assertTrue(wrt instanceof WrtFrequencyParameter);
            assertEquals(k, ((WrtFrequencyParameter) wrt).getDim());
        }
        try {
            hky.factory(frequencies, K);
            fail("dimension " + K + " must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    public void testUnsupportedModelRefused() {
        // only the listed classes are accepted; a subclass may change the generator
        Parameter subclassFrequencies = new Parameter.Default(PI.clone());
        HKY subclass = new HKY(KAPPA, new FrequencyModel(Nucleotides.INSTANCE, subclassFrequencies)) {
        };
        try {
            DifferentiableSubstitutionModelUtil.frequencyFactory(subclass, subclassFrequencies, 0);
            fail("frequency gradients must be refused for a class outside the allowlist");
        } catch (IllegalArgumentException expected) {
            // expected
        }
        assertNull(DifferentiableSubstitutionModelUtil.frequencyFactory(subclass, frequencies, 0));
    }

    public void testOtherModelsMatchFiniteDifferenceOfGenerator() {
        Parameter tn93Frequencies = new Parameter.Default(PI.clone());
        checkAgainstFiniteDifference(new TN93(new Parameter.Default(2.0), new Parameter.Default(4.0),
                new FrequencyModel(Nucleotides.INSTANCE, tn93Frequencies)), tn93Frequencies);

        Parameter gtrFrequencies = new Parameter.Default(PI.clone());
        checkAgainstFiniteDifference(new GTR(new Parameter.Default(new double[]{1.0, 2.0, 0.5, 1.5, 3.0, 1.0}),
                new FrequencyModel(Nucleotides.INSTANCE, gtrFrequencies)), gtrFrequencies);

        Parameter variableFrequencies = new Parameter.Default(PI.clone());
        checkAgainstFiniteDifference(new GTR(new Parameter.Default(1.0), new Parameter.Default(2.0),
                new Parameter.Default(0.5), new Parameter.Default(1.5), new Parameter.Default(3.0),
                new Parameter.Default(1.0), new FrequencyModel(Nucleotides.INSTANCE, variableFrequencies)),
                variableFrequencies);

        Parameter aaFrequencies = new Parameter.Default(skewedAminoAcidFrequencies());
        checkAgainstFiniteDifference(new EmpiricalAminoAcidModel(WAG.INSTANCE,
                new FrequencyModel(AminoAcids.INSTANCE, aaFrequencies)), aaFrequencies);
    }

    public void testGtrWithHkyRatesEqualsHky() {
        Parameter gtrFrequencies = new Parameter.Default(PI.clone());
        GTR gtr = new GTR(new Parameter.Default(new double[]{1.0, KAPPA, 1.0, 1.0, KAPPA, 1.0}),
                new FrequencyModel(Nucleotides.INSTANCE, gtrFrequencies));
        for (int k = 0; k < K; k++) {
            double[][] viaGtr = toArray(gtr.getInfinitesimalDifferentialMatrix(gtr.factory(gtrFrequencies, k)));
            double[][] viaHky = toArray(hky.getInfinitesimalDifferentialMatrix(hky.factory(frequencies, k)));
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    assertEquals("k=" + k + " (" + i + "," + j + ")", viaHky[i][j], viaGtr[i][j], 1E-14);
                }
            }
        }
    }

    public void testForeignParametersRefused() {
        Parameter foreign = new Parameter.Default(2.0);
        Parameter tn93Frequencies = new Parameter.Default(PI.clone());
        TN93 tn93 = new TN93(new Parameter.Default(2.0), new Parameter.Default(4.0),
                new FrequencyModel(Nucleotides.INSTANCE, tn93Frequencies));
        assertTrue(tn93.factory(tn93Frequencies, 0) instanceof WrtFrequencyParameter);
        try {
            tn93.factory(foreign, 0);
            fail("a parameter TN93 does not use must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("has no differential"));
        }

        EmpiricalAminoAcidModel wag = new EmpiricalAminoAcidModel(WAG.INSTANCE,
                new FrequencyModel(AminoAcids.INSTANCE, new Parameter.Default(skewedAminoAcidFrequencies())));
        try {
            wag.factory(foreign, 0);
            fail("empirical amino-acid models are differentiable only wrt their frequencies");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("only with respect to its own frequency parameter"));
        }
    }

    public void testLogRateFactoryRefusesForeignParameter() {
        Parameter logRates = new Parameter.Default(K * (K - 1), 0.0);
        LogRateSubstitutionModel model = new LogRateSubstitutionModel("logRate", Nucleotides.INSTANCE,
                new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(PI.clone())),
                new LogAdditiveCtmcRateProvider.DataAugmented.Basic("rates", logRates));
        assertNotNull(model.factory(logRates, 0));
        try {
            model.factory(frequencies, 0);
            fail("a parameter other than the log rates must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    public void testMarkovModulatedFactoryRefusesBaseFrequencies() {
        List<SubstitutionModel> baseModels = new ArrayList<SubstitutionModel>();
        baseModels.add(hky);
        baseModels.add(new HKY(1.5, new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(PI.clone()))));
        MarkovModulatedSubstitutionModel model = new MarkovModulatedSubstitutionModel("mm", baseModels,
                new Parameter.Default(2, 1.0),
                DataType.getRegisteredDataTypeByName(HiddenNucleotides.DESCRIPTION + "2"), null);
        try {
            model.factory(frequencies, 0);
            fail("base-model frequency gradients must be refused");
        } catch (IllegalArgumentException expected) {
            // expected
        }
    }

    // dQ~/dpi_k against central differences of the model's own normalized generator, at any state count
    private static void checkAgainstFiniteDifference(BaseSubstitutionModel model, Parameter pi) {
        final int stateCount = pi.getDimension();
        final double h = 1E-6;
        for (int k = 0; k < stateCount; k++) {
            WrappedMatrix analytic = ((DifferentiableSubstitutionModel) model)
                    .getInfinitesimalDifferentialMatrix(((DifferentiableSubstitutionModel) model).factory(pi, k));

            final double value = pi.getParameterValue(k);
            double[] plus = new double[stateCount * stateCount];
            double[] minus = new double[stateCount * stateCount];
            pi.setParameterValue(k, value + h);
            model.getInfinitesimalMatrix(plus);
            pi.setParameterValue(k, value - h);
            model.getInfinitesimalMatrix(minus);
            pi.setParameterValue(k, value);

            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    double numeric = (plus[i * stateCount + j] - minus[i * stateCount + j]) / (2 * h);
                    assertEquals(model.getModelName() + " k=" + k + " (" + i + "," + j + ")",
                            numeric, analytic.get(i, j), 1E-7 * Math.max(1.0, Math.abs(numeric)));
                }
            }
        }
    }

    private static double[] skewedAminoAcidFrequencies() {
        double[] pi = new double[AminoAcids.INSTANCE.getStateCount()];
        double sum = 0.0;
        for (int k = 0; k < pi.length; k++) {
            pi[k] = 1.0 + k;
            sum += pi[k];
        }
        for (int k = 0; k < pi.length; k++) {
            pi[k] /= sum;
        }
        return pi;
    }

    private double[] generatorWithFrequency(int k, double value) {
        frequencies.setParameterValue(k, value);
        double[] q = new double[K * K];
        hky.getInfinitesimalMatrix(q);
        return q;
    }

    private static double[][] toArray(WrappedMatrix matrix) {
        double[][] out = new double[K][K];
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                out[i][j] = matrix.get(i, j);
            }
        }
        return out;
    }

    // Scaling and squaring with a Taylor series; independent of the model's eigendecomposition
    private static double[] expm(double[] q, double time) {
        int squarings = 10;
        double scale = time / (1 << squarings);
        double[] a = new double[K * K];
        for (int i = 0; i < K * K; i++) {
            a[i] = q[i] * scale;
        }
        double[] result = identity();
        double[] term = identity();
        for (int n = 1; n <= 20; n++) {
            term = multiply(term, a);
            for (int i = 0; i < K * K; i++) {
                term[i] /= n;
                result[i] += term[i];
            }
        }
        for (int s = 0; s < squarings; s++) {
            result = multiply(result, result);
        }
        return result;
    }

    private static double[] identity() {
        double[] out = new double[K * K];
        for (int i = 0; i < K; i++) {
            out[i * K + i] = 1.0;
        }
        return out;
    }

    private static double[] multiply(double[] x, double[] y) {
        double[] out = new double[K * K];
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double sum = 0.0;
                for (int l = 0; l < K; l++) {
                    sum += x[i * K + l] * y[l * K + j];
                }
                out[i * K + j] = sum;
            }
        }
        return out;
    }
}
