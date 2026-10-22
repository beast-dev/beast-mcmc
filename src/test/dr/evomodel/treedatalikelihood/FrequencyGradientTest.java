/*
 * FrequencyGradientTest.java
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
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.DifferentialMassProvider.Mode;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionParameterGradient;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionModelGradient;
import dr.evomodel.treedatalikelihood.discrete.HomogeneousSubstitutionModelGradient.Backend;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Analytic gradient wrt stationary frequencies (plan 5, FQ6 and FQ8) on the b1 oracle fixture:
 * 5 taxa, 44 sites, HKY kappa = 2.5, pi = (0.1, 0.2, 0.3, 0.4), strict clock 1.0. TN93 and GTR reuse the fixture,
 * which also checks their relative-rate gradients (kappa1/kappa2, exchangeabilities) on the adjoint.
 */
public class FrequencyGradientTest extends TestCase {

    private static final String[][] SEQUENCES = {
            {"A", "AAAAAACCCCGGGGGTTTTTTTAACAAAGTTACGARACCCGGTA"},
            {"B", "AAAAAACCCCGGGGGTTTTTTTAACGGCGTTATATAACCCGGCC"},
            {"C", "AAAAAACCCCGGGGGTTTTTTTGGTAAGGCCACGAA-CCCAATC"},
            {"D", "AAAAAACCCCGGGGGTTTTTTTGGTAATATTATATGGCCCGGCC"},
            {"E", "AAAAAACCCCGGGGGTTTTTTTGGTAAAATTGCGAGGTTTGGCA"},
    };
    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";
    private static final double[] PI = {0.1, 0.2, 0.3, 0.4};
    private static final int SITE_COUNT = 44;

    public FrequencyGradientTest(String name) {
        super(name);
    }

    public void testGammaLiterals() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4), PartialsRescalingScheme.NONE, false);
        assertEquals(-195.67140160238029, f.likelihood.getLogLikelihood(), 1E-10);
        assertArrayEquals(new double[]{209.49044382424393, 59.60081659045633, 26.154964286906672, 16.20375178298842},
                f.gradient(Backend.ADJOINT).getGradientLogDensity(), 1E-8);
    }

    public void testGammaInvariantLiteralsAndEngineParity() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4, 0.3),
                PartialsRescalingScheme.ALWAYS, false);
        assertEquals(-196.41894163112337, f.likelihood.getLogLikelihood(), 1E-10);

        double[] expected = {211.85358175396522, 61.314533744023336, 26.036086712810548, 15.564703851122388};
        double[] adjoint = f.gradient(Backend.ADJOINT).getGradientLogDensity();
        double[] perDimension = f.gradient(Backend.PER_DIMENSION).getGradientLogDensity();
        assertArrayEquals(expected, adjoint, 1E-8);
        assertArrayEquals(adjoint, perDimension, 1E-9);
    }

    public void testTangentCheckPasses() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4), PartialsRescalingScheme.DYNAMIC, false);
        HomogeneousSubstitutionModelGradient gradient = f.gradient(Backend.AUTO);
        assertEquals(1.0, gradient.getSimplexSum());
        String report = gradient.getReport(); // throws if the tangent check fails
        assertTrue(report.contains("engine: adjoint"));
    }

    public void testFrequencyWeightedSumMatchesClockDerivative() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false);
        double piG = assertFrequencyWeightedSumIdentity(f, f.gradient(Backend.ADJOINT).getGradientLogDensity());
        assertEquals(64.2441088627, piG, 1E-8);
    }

    public void testGtrWithHkyRatesMatchesHkyLiterals() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4), PartialsRescalingScheme.NONE, false,
                Model.GTR_HKY_RATES);
        assertEquals(-195.67140160238029, f.likelihood.getLogLikelihood(), 1E-10);

        double[] expected = {209.49044382424393, 59.60081659045633, 26.154964286906672, 16.20375178298842};
        assertArrayEquals(expected, f.gradient(Backend.ADJOINT).getGradientLogDensity(), 1E-8);
        assertArrayEquals(expected, f.gradient(Backend.PER_DIMENSION).getGradientLogDensity(), 1E-8);
    }

    public void testTn93AndGtrEnginesAgreeAndPassTangentCheck() throws Exception {
        for (Model model : new Model[]{Model.TN93, Model.GTR}) {
            Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4, 0.3),
                    PartialsRescalingScheme.ALWAYS, false, model);

            double[] adjoint = f.gradient(Backend.ADJOINT).getGradientLogDensity();
            double[] perDimension = f.gradient(Backend.PER_DIMENSION).getGradientLogDensity();
            assertArrayEquals(adjoint, perDimension, 1E-9);
            assertFrequencyWeightedSumIdentity(f, adjoint);

            String report = f.gradient(Backend.AUTO).getReport(); // throws if the tangent check fails
            assertTrue(model + ": " + report, report.contains("engine: adjoint"));
        }
    }

    public void testSeparateRootFrequencies() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, true);
        HomogeneousSubstitutionModelGradient branchOnly = f.gradient(Backend.ADJOINT);
        HomogeneousSubstitutionModelGradient rootOnly = new HomogeneousSubstitutionModelGradient(f.likelihood,
                f.rootFrequencies, Backend.ADJOINT, 1E-5);

        double[] root = rootOnly.getGradientLogDensity();
        double piRoot = 0.0;
        for (int k = 0; k < root.length; k++) {
            piRoot += f.rootFrequencies.getParameterValue(k) * root[k];
        }
        assertEquals(SITE_COUNT, piRoot, 1E-10);

        double[] rootTermOfBranchOnly = branchOnly.getRootTerm();
        for (double value : rootTermOfBranchOnly) {
            assertEquals(0.0, value, 0.0);
        }
        // both reports run the tangent check against finite differences
        branchOnly.getReport();
        rootOnly.getReport();
    }

    public void testKappaThroughAdjointMatchesPerDimension() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false);
        HomogeneousSubstitutionModelGradient adjoint = new HomogeneousSubstitutionModelGradient(f.likelihood,
                f.kappa, Backend.ADJOINT, 1E-3);
        HomogeneousSubstitutionParameterGradient reference = new HomogeneousSubstitutionParameterGradient(
                "kappa", f.likelihood, f.kappa, f.delegate, 0, Mode.EXACT, 1E-3);

        double viaAdjoint = adjoint.getGradientLogDensity()[0];
        assertEquals(3.277118846996384, viaAdjoint, 1E-9);
        assertEquals(reference.getGradientLogDensity()[0], viaAdjoint, 1E-9);
        assertTrue(Double.isNaN(adjoint.getSimplexSum()));
    }

    public void testAdjointFollowsTreeAndClockChanges() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false);
        HomogeneousSubstitutionModelGradient gradient = f.gradient(Backend.ADJOINT);
        gradient.getGradientLogDensity();

        f.clockRate.setParameterValue(0, 1.7);
        f.tree.setNodeHeight(f.tree.getRoot(), f.tree.getNodeHeight(f.tree.getRoot()) + 0.2);

        double[] sameObject = gradient.getGradientLogDensity();
        double[] fresh = new HomogeneousSubstitutionModelGradient(f.likelihood, f.frequencies, Backend.ADJOINT, 1E-3)
                .getGradientLogDensity();
        assertArrayEquals(fresh, sameObject, 1E-10);
        gradient.getReport(); // tangent check at the new state
    }

    public void testAdjointBackendNeedsPreOrder() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false,
                PreOrderSettings.getDefault());
        try {
            f.gradient(Backend.ADJOINT);
            fail("backend=adjoint must be refused without usePreOrder");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("usePreOrder"));
        }
    }

    public void testFrequencyParameterDetection() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, true);
        assertTrue(HomogeneousSubstitutionModelGradient.isFrequencyParameter(f.likelihood, f.frequencies));
        assertTrue(HomogeneousSubstitutionModelGradient.isFrequencyParameter(f.likelihood, f.rootFrequencies));
        assertFalse(HomogeneousSubstitutionModelGradient.isFrequencyParameter(f.likelihood, f.kappa));
    }

    public void testRateGradientsMatchFiniteDifferences() throws Exception {
        for (Model model : new Model[]{Model.TN93, Model.GTR}) {
            Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4, 0.3),
                    PartialsRescalingScheme.ALWAYS, false, model);
            for (Parameter rates : f.rateParameters) {
                assertTrue(HomogeneousSubstitutionModelGradient.isRateParameter(f.likelihood, rates));
                double[] adjoint = rateGradient(f, rates, Backend.ADJOINT);
                assertArrayEquals(finiteDifference(f, rates), adjoint, 1E-6);
                assertArrayEquals(adjoint, rateGradient(f, rates, Backend.PER_DIMENSION), 1E-9);
            }
        }
    }

    public void testGtrRateGradientIsScaleInvariant() throws Exception {
        // Q~ does not change when all exchangeabilities are scaled, so s.g = 0
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel", 0.5, 4), PartialsRescalingScheme.NONE, false,
                Model.GTR);
        Parameter rates = f.rateParameters.get(0);
        double[] g = rateGradient(f, rates, Backend.ADJOINT);
        double sg = 0.0;
        double scale = 0.0;
        for (int m = 0; m < g.length; m++) {
            sg += rates.getParameterValue(m) * g[m];
            scale += Math.abs(rates.getParameterValue(m) * g[m]);
        }
        assertEquals(0.0, sg, 1E-10 * scale);
    }

    public void testGtrWithHkyRatesMatchesHkyKappa() throws Exception {
        // dlogL/dkappa = dlogL/ds_AG + dlogL/ds_CT at s = (1, kappa, 1, 1, kappa, 1)
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false,
                Model.GTR_HKY_RATES);
        double[] g = rateGradient(f, f.rateParameters.get(0), Backend.ADJOINT);
        assertEquals(3.277118846996384, g[1] + g[4], 1E-9);
    }

    public void testPerDimensionKeepsUnnamedParametersApart() throws Exception {
        // TN93's kappas have no ids here; each must get its own per-dimension traits
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false, Model.TN93);
        double first = rateGradient(f, f.rateParameters.get(0), Backend.PER_DIMENSION)[0];
        double second = rateGradient(f, f.rateParameters.get(1), Backend.PER_DIMENSION)[0];
        assertEquals(rateGradient(f, f.rateParameters.get(0), Backend.ADJOINT)[0], first, 1E-9);
        assertEquals(rateGradient(f, f.rateParameters.get(1), Backend.ADJOINT)[0], second, 1E-9);
    }

    public void testRateParameterDetection() throws Exception {
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false, Model.TN93);
        assertFalse(HomogeneousSubstitutionModelGradient.isRateParameter(f.likelihood, f.frequencies));
        Fixture hky = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false);
        assertFalse(HomogeneousSubstitutionModelGradient.isRateParameter(hky.likelihood, hky.kappa));
    }

    private static double[] rateGradient(Fixture f, Parameter rates, Backend backend) {
        HomogeneousSubstitutionModelGradient gradient = new HomogeneousSubstitutionModelGradient(f.likelihood,
                rates, backend, 1E-5);
        assertTrue(Double.isNaN(gradient.getSimplexSum()));
        for (double value : gradient.getRootTerm()) {
            assertEquals(0.0, value, 0.0);
        }
        return gradient.getGradientLogDensity();
    }

    private static double[] finiteDifference(Fixture f, Parameter parameter) {
        final double h = 1E-6;
        double[] numeric = new double[parameter.getDimension()];
        for (int d = 0; d < numeric.length; d++) {
            final double value = parameter.getParameterValue(d);
            parameter.setParameterValue(d, value + h);
            double plus = f.likelihood.getLogLikelihood();
            parameter.setParameterValue(d, value - h);
            double minus = f.likelihood.getLogLikelihood();
            parameter.setParameterValue(d, value);
            numeric[d] = (plus - minus) / (2 * h);
        }
        return numeric;
    }

    // sum_k pi_k R_k = N and sum_k pi_k B_k = -r d log L / d r, so pi.g = N - r d log L / d r for any normalized Q
    private static double assertFrequencyWeightedSumIdentity(Fixture f, double[] g) {
        double piG = 0.0;
        for (int k = 0; k < g.length; k++) {
            piG += f.frequencies.getParameterValue(k) * g[k];
        }

        final double h = 1E-6;
        f.clockRate.setParameterValue(0, 1.0 + h);
        double plus = f.likelihood.getLogLikelihood();
        f.clockRate.setParameterValue(0, 1.0 - h);
        double minus = f.likelihood.getLogLikelihood();
        f.clockRate.setParameterValue(0, 1.0);
        double clockDerivative = (plus - minus) / (2 * h);

        assertEquals(SITE_COUNT - clockDerivative, piG, 1E-6 * Math.abs(piG));
        return piG;
    }

    private static void assertArrayEquals(double[] expected, double[] actual, double relativeTolerance) {
        assertEquals(expected.length, actual.length);
        for (int k = 0; k < expected.length; k++) {
            assertEquals("k=" + k, expected[k], actual[k], relativeTolerance * Math.max(1.0, Math.abs(expected[k])));
        }
    }

    private enum Model {HKY, TN93, GTR_HKY_RATES, GTR}

    private static class Fixture {

        final Parameter frequencies = new Parameter.Default(PI.clone());
        final Parameter kappa = new Parameter.Default(2.5);
        final Parameter clockRate = new Parameter.Default(1.0);
        final List<Parameter> rateParameters = new ArrayList<Parameter>();
        final Parameter rootFrequencies;
        final DefaultTreeModel tree;
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;

        Fixture(GammaSiteRateModel siteRateModel, PartialsRescalingScheme scheme, boolean separateRoot)
                throws Exception {
            this(siteRateModel, scheme, separateRoot, new PreOrderSettings(true, false, true, true), Model.HKY);
        }

        Fixture(GammaSiteRateModel siteRateModel, PartialsRescalingScheme scheme, boolean separateRoot, Model model)
                throws Exception {
            this(siteRateModel, scheme, separateRoot, new PreOrderSettings(true, false, true, true), model);
        }

        Fixture(GammaSiteRateModel siteRateModel, PartialsRescalingScheme scheme, boolean separateRoot,
                PreOrderSettings settings) throws Exception {
            this(siteRateModel, scheme, separateRoot, settings, Model.HKY);
        }

        Fixture(GammaSiteRateModel siteRateModel, PartialsRescalingScheme scheme, boolean separateRoot,
                PreOrderSettings settings, Model model) throws Exception {

            SimpleAlignment alignment = new SimpleAlignment();
            alignment.setDataType(Nucleotides.INSTANCE);
            for (String[] row : SEQUENCES) {
                alignment.addSequence(new Sequence(new Taxon(row[0]), row[1]));
            }
            Tree newick = new NewickImporter(NEWICK).importTree(null);
            tree = new DefaultTreeModel("treeModel", newick);

            SubstitutionModel substitutionModel = create(model, new FrequencyModel(Nucleotides.INSTANCE, frequencies));
            HomogeneousBranchModel branchModel;
            if (separateRoot) {
                rootFrequencies = new Parameter.Default(new double[]{0.4, 0.3, 0.2, 0.1});
                branchModel = new HomogeneousBranchModel(substitutionModel, new FrequencyModel(Nucleotides.INSTANCE, rootFrequencies));
            } else {
                rootFrequencies = null;
                branchModel = new HomogeneousBranchModel(substitutionModel);
            }

            delegate = new BeagleDataLikelihoodDelegate(tree, new SitePatterns(alignment, null, 0, -1, 1, true),
                    branchModel, siteRateModel, true, false, scheme, false, settings);
            likelihood = new TreeDataLikelihood(delegate, tree, new StrictClockBranchRates(clockRate));
        }

        private SubstitutionModel create(Model model, FrequencyModel frequencyModel) {
            switch (model) {
                case HKY:
                    return new HKY(kappa, frequencyModel);
                case TN93:
                    rateParameters.add(new Parameter.Default(2.0));
                    rateParameters.add(new Parameter.Default(4.0));
                    return new TN93(rateParameters.get(0), rateParameters.get(1), frequencyModel);
                case GTR_HKY_RATES:
                    rateParameters.add(new Parameter.Default(new double[]{1.0, 2.5, 1.0, 1.0, 2.5, 1.0}));
                    return new GTR(rateParameters.get(0), frequencyModel);
                case GTR:
                    rateParameters.add(new Parameter.Default(new double[]{1.0, 2.0, 0.5, 1.5, 3.0, 1.0}));
                    return new GTR(rateParameters.get(0), frequencyModel);
                default:
                    throw new IllegalArgumentException(model.name());
            }
        }

        HomogeneousSubstitutionModelGradient gradient(Backend backend) {
            return new HomogeneousSubstitutionModelGradient(likelihood, frequencies, backend, 1E-5);
        }
    }
}
