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
import dr.evomodel.substmodel.nucleotide.HKY;
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

/**
 * Analytic gradient wrt HKY stationary frequencies (plan 5, FQ6) on the b1 oracle fixture:
 * 5 taxa, 44 sites, HKY kappa = 2.5, pi = (0.1, 0.2, 0.3, 0.4), strict clock 1.0.
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
        // sum_k pi_k R_k = N and sum_k pi_k B_k = -r d log L / d r, so pi.g = N - r d log L / d r
        Fixture f = new Fixture(new GammaSiteRateModel("siteModel"), PartialsRescalingScheme.NONE, false);
        double[] g = f.gradient(Backend.ADJOINT).getGradientLogDensity();
        double piG = 0.0;
        for (int k = 0; k < PI.length; k++) {
            piG += PI[k] * g[k];
        }

        final double h = 1E-6;
        f.clockRate.setParameterValue(0, 1.0 + h);
        double plus = f.likelihood.getLogLikelihood();
        f.clockRate.setParameterValue(0, 1.0 - h);
        double minus = f.likelihood.getLogLikelihood();
        f.clockRate.setParameterValue(0, 1.0);
        double clockDerivative = (plus - minus) / (2 * h);

        assertEquals(SITE_COUNT - clockDerivative, piG, 1E-6 * Math.abs(piG));
        assertEquals(64.2441088627, piG, 1E-8);
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

    private static void assertArrayEquals(double[] expected, double[] actual, double relativeTolerance) {
        assertEquals(expected.length, actual.length);
        for (int k = 0; k < expected.length; k++) {
            assertEquals("k=" + k, expected[k], actual[k], relativeTolerance * Math.max(1.0, Math.abs(expected[k])));
        }
    }

    private static class Fixture {

        final Parameter frequencies = new Parameter.Default(PI.clone());
        final Parameter kappa = new Parameter.Default(2.5);
        final Parameter clockRate = new Parameter.Default(1.0);
        final Parameter rootFrequencies;
        final DefaultTreeModel tree;
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;

        Fixture(GammaSiteRateModel siteRateModel, PartialsRescalingScheme scheme, boolean separateRoot)
                throws Exception {
            this(siteRateModel, scheme, separateRoot, new PreOrderSettings(true, false, true, true));
        }

        Fixture(GammaSiteRateModel siteRateModel, PartialsRescalingScheme scheme, boolean separateRoot,
                PreOrderSettings settings) throws Exception {

            SimpleAlignment alignment = new SimpleAlignment();
            alignment.setDataType(Nucleotides.INSTANCE);
            for (String[] row : SEQUENCES) {
                alignment.addSequence(new Sequence(new Taxon(row[0]), row[1]));
            }
            Tree newick = new NewickImporter(NEWICK).importTree(null);
            tree = new DefaultTreeModel("treeModel", newick);

            HKY hky = new HKY(kappa, new FrequencyModel(Nucleotides.INSTANCE, frequencies));
            HomogeneousBranchModel branchModel;
            if (separateRoot) {
                rootFrequencies = new Parameter.Default(new double[]{0.4, 0.3, 0.2, 0.1});
                branchModel = new HomogeneousBranchModel(hky, new FrequencyModel(Nucleotides.INSTANCE, rootFrequencies));
            } else {
                rootFrequencies = null;
                branchModel = new HomogeneousBranchModel(hky);
            }

            delegate = new BeagleDataLikelihoodDelegate(tree, new SitePatterns(alignment, null, 0, -1, 1, true),
                    branchModel, siteRateModel, true, false, scheme, false, settings);
            likelihood = new TreeDataLikelihood(delegate, tree, new StrictClockBranchRates(clockRate));
        }

        HomogeneousSubstitutionModelGradient gradient(Backend backend) {
            return new HomogeneousSubstitutionModelGradient(likelihood, frequencies, backend, 1E-5);
        }
    }
}
