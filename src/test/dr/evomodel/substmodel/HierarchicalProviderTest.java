/*
 * HierarchicalProviderTest.java
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
import dr.evolution.datatype.TwoStates;
import dr.evomodel.substmodel.AbstractHierarchicalRateProvider;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.EqualRatesProvider;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HierarchicalSubstitutionModel;
import dr.evomodel.substmodel.RelativeRatesProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.ModelListener;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Rate providers of the hierarchical substitution model against the dense models they represent.
 *
 * @author Xiang Ji
 */
public class HierarchicalProviderTest extends TestCase {

    private static final int S = 4;
    private static final String[] TIP_ORDERS = {"((A:1,G:1):1,(C:1,T:1):1);", "((T:1,C:1):1,(G:1,A:1):1);"};
    private static final double[][] FREQUENCIES = {{0.1, 0.2, 0.3, 0.4}, {0.25, 0.25, 0.25, 0.25}, {0.4, 0.1, 0.15, 0.35}};
    private static final double[] TIMES = {0.01, 0.1, 1.0, 5.0};

    private static double maxDiff(SubstitutionModel a, SubstitutionModel b) {
        final int n = a.getDataType().getStateCount();
        double[] x = new double[n * n], y = new double[n * n];
        a.getInfinitesimalMatrix(x);
        b.getInfinitesimalMatrix(y);
        double d = 0.0;
        for (int i = 0; i < x.length; i++)
            d = Math.max(d, Math.abs(x[i] - y[i]));
        for (double t : TIMES) {
            a.getTransitionProbabilities(t, x);
            b.getTransitionProbabilities(t, y);
            for (int i = 0; i < x.length; i++)
                d = Math.max(d, Math.abs(x[i] - y[i]));
        }
        return d;
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

    private static FrequencyModel frequencies(double[] pi) {
        return new FrequencyModel(Nucleotides.INSTANCE, new Parameter.Default(pi));
    }

    public void testHKYEqualsDense() {
        for (double[] pi : FREQUENCIES) {
            for (String newick : TIP_ORDERS) {
                for (double k : new double[]{0.5, 1.0, 2.5, 8.0}) {
                    FrequencyModel f = frequencies(pi);
                    Parameter kappa = new Parameter.Default(k);
                    HKY hky = new HKY(kappa, f);
                    HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                            RelativeRatesProvider.hky(kappa, f, AbstractHierarchicalRateProvider.characterTree(newick)));
                    assertEquals(newick + " kappa " + k, 0.0, maxDiff(model, hky), 1e-13);
                }
            }
        }
    }

    public void testTN93EqualsDense() {
        for (double[] pi : FREQUENCIES) {
            for (String newick : TIP_ORDERS) {
                for (double[] k : new double[][]{{0.5, 3.0}, {2.5, 2.5}, {8.0, 1.5}}) {
                    FrequencyModel f = frequencies(pi);
                    Parameter kappa1 = new Parameter.Default(k[0]), kappa2 = new Parameter.Default(k[1]);
                    TN93 tn93 = new TN93(kappa1, kappa2, f);
                    HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                            RelativeRatesProvider.tn93(kappa1, kappa2, f, AbstractHierarchicalRateProvider.characterTree(newick)));
                    assertEquals(newick + " kappa " + Arrays.toString(k), 0.0, maxDiff(model, tn93), 1e-13);
                }
            }
        }
    }

    public void testTiedGTREqualsDense() {
        FrequencyModel f = frequencies(FREQUENCIES[0]);
        Parameter rates = new Parameter.Default(new double[]{2.2, 1.4, 0.6, 1.4, 0.6, 0.6});
        GTR gtr = new GTR(rates, f);
        RelativeRatesProvider provider = RelativeRatesProvider.gtr(rates, f,
                AbstractHierarchicalRateProvider.characterTree("(((A:1,C:1):1,G:1):1,T:1);"), 1e-10);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", provider);
        assertEquals(0.0, maxDiff(model, gtr), 1e-13);
        assertEquals(0.628, provider.getNormalization(), 1e-14);
        assertEquals(0.3821656050955414, provider.getLambda("AC"), 1e-12);
        assertEquals(0.7643312101910829, provider.getLambda("ACG"), 1e-12);
        assertEquals(0.9554140127388535, provider.getLambda("ACGT"), 1e-12);
    }

    public void testXmlLiterals() {
        FrequencyModel f = frequencies(FREQUENCIES[0]);
        RelativeRatesProvider provider = RelativeRatesProvider.hky(new Parameter.Default(2.5), f);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", provider);
        assertEquals(1.03, provider.getNormalization(), 1e-15);
        assertEquals(0.5825242718446603, provider.getLambda("AG"), 1e-15);
        assertEquals(0.8737864077669906, provider.getLambda("CT"), 1e-15);
        assertEquals(0.9708737864077672, provider.getLambda("ACGT"), 1e-15);
        assertTrue(Arrays.equals(new String[]{"ACGT", "CT", "AG"}, provider.getCladeNames())
                || Arrays.equals(new String[]{"ACGT", "AG", "CT"}, provider.getCladeNames()));
        double[] lambdas = model.getBeagleModelLambdas();
        assertEquals(65539.0, lambdas[0]);
        assertEquals(0.5825242718446603, lambdas[1], 1e-15);
        assertEquals(0.8737864077669906, lambdas[2], 1e-15);
        assertEquals(0.9708737864077672, lambdas[3], 1e-15);
    }

    public void testEigenSystemAndBundle() {
        for (double k : new double[]{1.0, 2.5}) {
            FrequencyModel f = frequencies(FREQUENCIES[0]);
            RelativeRatesProvider provider = RelativeRatesProvider.hky(new Parameter.Default(k), f);
            HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", provider);
            double[] p = new double[S * S];
            for (double t : new double[]{0.0, 0.05, 0.8, 3.0}) {
                model.getTransitionProbabilities(t, p);
                double[] q = fromEigen(model.getEigenDecomposition(), t);
                for (int i = 0; i < S * S; i++)
                    assertEquals("kappa " + k + " t " + t, p[i], q[i], 1e-13);
            }

            TreeModel tree = provider.getTreeModel();
            Parameter lambdas = new Parameter.Default(tree.getInternalNodeCount());
            for (int i = 0; i < tree.getInternalNodeCount(); i++)
                lambdas.setParameterValue(i, provider.getRate(tree, S + i));
            HierarchicalSubstitutionModel reference = new HierarchicalSubstitutionModel("default",
                    new HierarchicalSubstitutionModel.HierarchicalRateProvider.Default(lambdas, tree), f);
            double[] categoryRates = {0.3, 1.0, 2.7};
            double[] a = new double[model.getBundleLength(3)], b = new double[reference.getBundleLength(3)];
            for (double t : new double[]{0.0, 0.05, 0.37, 1.9}) {
                model.fillBundle(t, categoryRates, a);
                reference.fillBundle(t, categoryRates, b);
                assertTrue("bundle kappa " + k + " t " + t, Arrays.equals(a, b));
            }
        }
    }

    public void testChangesAndRestore() {
        Parameter pi = new Parameter.Default(FREQUENCIES[0]);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, pi);
        Parameter kappa = new Parameter.Default(2.5);
        HKY hky = new HKY(kappa, f);
        RelativeRatesProvider provider = RelativeRatesProvider.hky(kappa, f);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", provider);
        final int[] fired = {0};
        model.addModelListener(new ModelListener() {
            public void modelChangedEvent(Model m, Object object, int index) {
                fired[0]++;
            }

            public void modelRestored(Model m) {
            }
        });

        double[] p0 = new double[S * S];
        model.getTransitionProbabilities(0.3, p0);
        double[] lambdas0 = model.getBeagleModelLambdas().clone();
        final int computed = provider.getComputeCount();
        for (int i = 0; i < 5; i++)
            model.getTransitionProbabilities(0.3, new double[S * S]);
        assertEquals("no recompute on reads", computed, provider.getComputeCount());

        model.storeModelState();
        fired[0] = 0;
        kappa.setParameterValue(0, 6.0);
        assertTrue("kappa change reaches the model", fired[0] >= 1);
        assertEquals(0.0, maxDiff(model, hky), 1e-13);
        assertEquals("one recompute per change", computed + 1, provider.getComputeCount());
        model.restoreModelState();
        double[] p = new double[S * S];
        model.getTransitionProbabilities(0.3, p);
        assertTrue("P restored", Arrays.equals(p0, p));
        assertTrue("BEAGLE lambdas restored", Arrays.equals(lambdas0, model.getBeagleModelLambdas()));

        fired[0] = 0;
        pi.setParameterValue(0, 0.15);
        pi.setParameterValue(3, 0.35);
        assertTrue("frequency change reaches the model", fired[0] >= 1);
        assertEquals(0.0, maxDiff(model, hky), 1e-13);

        Parameter kappa1 = new Parameter.Default(2.0), kappa2 = new Parameter.Default(3.0);
        TN93 tn93 = new TN93(kappa1, kappa2, f);
        HierarchicalSubstitutionModel tnModel = new HierarchicalSubstitutionModel("hrate",
                RelativeRatesProvider.tn93(kappa1, kappa2, f));
        kappa1.setParameterValue(0, 0.7);
        assertEquals(0.0, maxDiff(tnModel, tn93), 1e-13);
        kappa2.setParameterValue(0, 9.0);
        assertEquals(0.0, maxDiff(tnModel, tn93), 1e-13);
    }

    public void testRejections() {
        FrequencyModel f = frequencies(FREQUENCIES[0]);
        for (double k : new double[]{1.0, 2.0}) {
            try {
                RelativeRatesProvider.hky(new Parameter.Default(k), f,
                        AbstractHierarchicalRateProvider.characterTree("((A:1,C:1):1,(G:1,T:1):1);"));
                fail("HKY on ((A,C),(G,T)) must be refused");
            } catch (IllegalArgumentException e) {
                assertTrue(e.getMessage(), e.getMessage().contains("must have the clades {A,G} and {C,T}")
                        && e.getMessage().contains("{A,C}") && e.getMessage().contains("{G,T}"));
            }
        }
        try {
            RelativeRatesProvider.tn93(new Parameter.Default(2.0), new Parameter.Default(3.0), f,
                    AbstractHierarchicalRateProvider.characterTree("(((A:1,G:1):1,C:1):1,T:1);"));
            fail("TN93 on a caterpillar must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("must have the clades {A,G} and {C,T}"));
        }

        RelativeRatesProvider equal = RelativeRatesProvider.gtr(new Parameter.Default(6, 1.0), f,
                AbstractHierarchicalRateProvider.characterTree("((A:1,C:1):1,(G:1,T:1):1);"), 1e-10);
        assertEquals(0.0, equal.getLambda("AC"), 0.0);
        assertEquals(0.0, equal.getLambda("GT"), 0.0);

        try {
            RelativeRatesProvider.gtr(new Parameter.Default(new double[]{1.0, 2.0, 1.1, 1.0, 2.0, 1.0}), f,
                    AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE), 1e-10);
            fail("untied rates must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("not tree-structured") && e.getMessage().contains("ACGT"));
        }

        try {
            RelativeRatesProvider.hky(new Parameter.Default(2.0),
                    new FrequencyModel(TwoStates.INSTANCE, new Parameter.Default(new double[]{0.3, 0.7})));
            fail("kappa on two states must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("need a nucleotide frequency model"));
        }

        RelativeRatesProvider provider = RelativeRatesProvider.hky(new Parameter.Default(2.0), f);
        try {
            new HierarchicalSubstitutionModel("hrate", provider, frequencies(FREQUENCIES[0]));
            fail("a second FrequencyModel must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("share one FrequencyModel"));
        }
        try {
            provider.getLambda("AC");
            fail("AC is not a clade");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("{A,C} is not a clade"));
        }
    }

    public void testKappaBoundsAreTheGeneratorConstraints() {
        FrequencyModel f = frequencies(FREQUENCIES[0]);
        Parameter kappa = new Parameter.Default(0.5);
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", RelativeRatesProvider.hky(kappa, f));
        assertEquals(0.0, kappa.getBounds().getLowerLimit(0));
        assertTrue("kappa 0.5 is a valid generator", model.isValidGenerator());
        assertFalse("with negative clade rates", model.hasNonNegativeRates());
        kappa.setParameterValue(0, 0.0);
        assertTrue("kappa 0 is on the boundary", model.isValidGenerator());

        Parameter rates = new Parameter.Default(new double[]{2.2, 1.4, 0.6, 1.4, 0.6, 0.6});
        RelativeRatesProvider.gtr(rates, f, AbstractHierarchicalRateProvider.characterTree("(((A:1,C:1):1,G:1):1,T:1);"), 1e-10);
        for (int i = 0; i < 6; i++)
            assertEquals(0.0, rates.getBounds().getLowerLimit(i));
    }

    public void testOfExistingModel() {
        Parameter pi = new Parameter.Default(FREQUENCIES[0]);
        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, pi);
        HKY hky = new HKY(new Parameter.Default(2.5), f);
        RelativeRatesProvider provider = RelativeRatesProvider.of(hky, AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE));
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", provider);
        assertSame(f, model.getFrequencyModel());
        assertEquals(0.5825242718446603, provider.getLambda("AG"), 1e-15);
        assertEquals(0.0, maxDiff(model, hky), 1e-13);

        double[] p0 = new double[S * S];
        model.getTransitionProbabilities(0.3, p0);
        model.storeModelState();
        hky.setKappa(6.0);
        assertEquals("setKappa reaches the provider", 0.0, maxDiff(model, hky), 1e-13);
        model.restoreModelState();
        assertEquals("kappa restored with the borrowed model", 2.5, hky.getKappa());
        double[] p = new double[S * S];
        model.getTransitionProbabilities(0.3, p);
        assertTrue(Arrays.equals(p0, p));
        pi.setParameterValue(0, 0.15);
        pi.setParameterValue(3, 0.35);
        assertEquals("frequency change", 0.0, maxDiff(model, hky), 1e-13);
        hky.setKappa(0.5);
        assertEquals("kappa below 1", 0.0, maxDiff(model, hky), 1e-13);

        for (String newick : TIP_ORDERS) {
            Parameter kappa1 = new Parameter.Default(1.7), kappa2 = new Parameter.Default(3.1);
            TN93 tn93 = new TN93(kappa1, kappa2, f);
            HierarchicalSubstitutionModel tnModel = new HierarchicalSubstitutionModel("hrate",
                    RelativeRatesProvider.of(tn93, AbstractHierarchicalRateProvider.characterTree(newick)));
            assertEquals(newick, 0.0, maxDiff(tnModel, tn93), 1e-13);
            kappa2.setParameterValue(0, 0.4);
            assertEquals(newick + " kappa2 change", 0.0, maxDiff(tnModel, tn93), 1e-13);
        }

        try {
            RelativeRatesProvider.of(new GTR(new Parameter.Default(6, 1.0), f),
                    AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE));
            fail("GTR is not supported yet");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("supports HKYModel, TN93Model"));
        }
        try {
            RelativeRatesProvider.of(hky, AbstractHierarchicalRateProvider.characterTree("((A:1,C:1):1,(G:1,T:1):1);"));
            fail("HKY on ((A,C),(G,T)) must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("must have the clades {A,G} and {C,T}"));
        }
    }

    private static double[] ramp(int n, double step) {
        double[] pi = new double[n];
        double total = 0.0;
        for (int i = 0; i < n; i++)
            total += pi[i] = 1.0 + step * i;
        for (int i = 0; i < n; i++)
            pi[i] /= total;
        return pi;
    }

    private static GeneralSubstitutionModel equalRatesModel(dr.evolution.datatype.DataType type, FrequencyModel f) {
        final int n = type.getStateCount();
        if (n == 2) {
            return new GeneralSubstitutionModel("binary", type, f, new Parameter.Default(0), 0);
        }
        return new GeneralSubstitutionModel("equal", type, f, new Parameter.Default(n * (n - 1) / 2, 1.0), -1);
    }

    private static void checkEqualRates(String what, HierarchicalSubstitutionModel model, double[] pi, boolean normalised) {
        final int n = pi.length;
        double sumSquares = 0.0;
        for (double p : pi)
            sumSquares += p * p;
        double[] lambdas = model.getBeagleModelLambdas();
        assertEquals(what, 65536.0 + n - 1, lambdas[0]);
        int zeros = 0;
        for (int r = 1; r < n; r++) {
            if (lambdas[r] == 0.0)
                zeros++;
            else
                assertEquals(what + " root", normalised ? 1.0 / (1.0 - sumSquares) : 1.0, lambdas[r], 1e-14);
        }
        assertEquals(what + " zero clade rates", n - 2, zeros);
        double[] bundle = new double[model.getBundleLength(1)];
        model.getBundleHeader(bundle);
        model.fillBundleBlock(0.3, bundle, model.getBundleHeaderLength());
        int zeroH = 0;
        for (int m = 0; m < n - 1; m++)
            if (bundle[model.getBundleHeaderLength() + m] == 0.0)
                zeroH++;
        assertEquals(what + " H = 0 below the root", n - 2, zeroH);
    }

    public void testEqualRatesNucleotides() {
        for (double[] pi : FREQUENCIES) {
            for (String newick : new String[]{null, "(((A:1,C:1):1,G:1):1,T:1);"}) {
                FrequencyModel f = frequencies(pi);
                HKY hky = new HKY(new Parameter.Default(1.0), f);
                TreeModel tree = newick == null ? AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE)
                        : AbstractHierarchicalRateProvider.characterTree(newick);
                HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                        new EqualRatesProvider(tree, f, true));
                assertEquals(newick, 0.0, maxDiff(model, hky), 1e-13);
                checkEqualRates("nucleotides " + newick, model, pi, true);
            }
        }
        HierarchicalSubstitutionModel jc = new HierarchicalSubstitutionModel("hrate",
                new EqualRatesProvider(AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE),
                        frequencies(FREQUENCIES[1]), true));
        assertEquals(4.0 / 3.0, jc.getBeagleModelLambdas()[3], 1e-15);
    }

    public void testEqualRatesOtherDataTypes() {
        String[] codes = {"a", "b", "c", "d", "e", "f", "g"};
        dr.evolution.datatype.DataType[] types = {TwoStates.INSTANCE, new dr.evolution.datatype.GeneralDataType(codes),
                dr.evolution.datatype.AminoAcids.INSTANCE};
        for (dr.evolution.datatype.DataType type : types) {
            final int n = type.getStateCount();
            for (double step : new double[]{0.0, 0.37}) {
                double[] pi = ramp(n, step);
                Parameter piParameter = new Parameter.Default(pi);
                FrequencyModel f = new FrequencyModel(type, piParameter);
                GeneralSubstitutionModel dense = equalRatesModel(type, f);
                TreeModel tree = AbstractHierarchicalRateProvider.defaultCharacterTree(type);
                assertEquals(n, tree.getExternalNodeCount());
                HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate",
                        new EqualRatesProvider(tree, f, true));
                assertTrue(model.isBinary());
                assertEquals(type.getDescription(), 0.0, maxDiff(model, dense), 1e-13);
                checkEqualRates(type.getDescription(), model, pi, true);

                HierarchicalSubstitutionModel borrowed = new HierarchicalSubstitutionModel("hrate",
                        EqualRatesProvider.of(dense, tree));
                assertEquals(type.getDescription() + " ratesFrom", 0.0, maxDiff(borrowed, dense), 1e-13);
                piParameter.setParameterValue(0, piParameter.getParameterValue(0) + 0.01);
                piParameter.setParameterValue(n - 1, piParameter.getParameterValue(n - 1) - 0.01);
                assertEquals(type.getDescription() + " frequency change", 0.0, maxDiff(model, dense), 1e-13);
                assertEquals(type.getDescription() + " frequency change, ratesFrom", 0.0, maxDiff(borrowed, dense), 1e-13);
            }
        }

        dr.evolution.datatype.DataType seven = new dr.evolution.datatype.GeneralDataType(codes);
        double[] pi = ramp(7, 0.2);
        FrequencyModel f = new FrequencyModel(seven, new Parameter.Default(pi));
        GeneralSubstitutionModel unnormalised = equalRatesModel(seven, f);
        unnormalised.setNormalization(false);
        TreeModel caterpillar = AbstractHierarchicalRateProvider.characterTree("((((((a:1,b:1):1,c:1):1,d:1):1,e:1):1,f:1):1,g:1);");
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", new EqualRatesProvider(caterpillar, f, false));
        assertEquals("unnormalised, caterpillar", 0.0, maxDiff(model, unnormalised), 1e-13);
        checkEqualRates("unnormalised", model, pi, false);
        assertEquals("ratesFrom an unnormalised model", 0.0, maxDiff(new HierarchicalSubstitutionModel("hrate",
                EqualRatesProvider.of(unnormalised, caterpillar)), unnormalised), 1e-13);
    }

    public void testEqualRatesFromRejectsUnequalRates() {
        dr.evolution.datatype.DataType seven = new dr.evolution.datatype.GeneralDataType(
                new String[]{"a", "b", "c", "d", "e", "f", "g"});
        FrequencyModel f = new FrequencyModel(seven, new Parameter.Default(ramp(7, 0.1)));
        Parameter rates = new Parameter.Default(21, 1.0);
        GeneralSubstitutionModel dense = new GeneralSubstitutionModel("equal", seven, f, rates, -1);
        EqualRatesProvider provider = EqualRatesProvider.of(dense, AbstractHierarchicalRateProvider.defaultCharacterTree(seven));
        HierarchicalSubstitutionModel model = new HierarchicalSubstitutionModel("hrate", provider);
        model.getBeagleModelLambdas();
        rates.setParameterValue(4, 1.5);
        try {
            model.getBeagleModelLambdas();
            fail("unequal rates must be refused");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("does not have equal rates"));
        }
        rates.setParameterValue(4, 1.0);
        try {
            EqualRatesProvider.of(new GeneralSubstitutionModel("unequal", seven, f, new Parameter.Default(21, 2.0) {{
                setParameterValue(0, 0.5);
            }}, -1), AbstractHierarchicalRateProvider.defaultCharacterTree(seven));
            fail("unequal rates must be refused at construction");
        } catch (IllegalArgumentException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("does not have equal rates"));
        }
    }

    public void testDefaultCharacterTrees() {
        TreeModel amino = AbstractHierarchicalRateProvider.defaultCharacterTree(dr.evolution.datatype.AminoAcids.INSTANCE);
        assertEquals(20, amino.getExternalNodeCount());
        assertEquals(19, amino.getInternalNodeCount());
        int depth = 0;
        for (int i = 0; i < amino.getExternalNodeCount(); i++) {
            int d = 0;
            for (dr.evolution.tree.NodeRef node = amino.getExternalNode(i); !amino.isRoot(node); node = amino.getParent(node))
                d++;
            depth = Math.max(depth, d);
        }
        assertEquals(5, depth);
        TreeModel binary = AbstractHierarchicalRateProvider.defaultCharacterTree(TwoStates.INSTANCE);
        assertEquals(2, binary.getExternalNodeCount());
        TreeModel nucleotides = AbstractHierarchicalRateProvider.defaultCharacterTree(Nucleotides.INSTANCE);
        RelativeRatesProvider hky = RelativeRatesProvider.hky(new Parameter.Default(2.0), frequencies(FREQUENCIES[0]), nucleotides);
        assertEquals(3, hky.getCladeNames().length);
    }
}
