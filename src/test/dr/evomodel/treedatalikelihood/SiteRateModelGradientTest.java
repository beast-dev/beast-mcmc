/*
 * SiteRateModelGradientTest.java
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
import dr.evomodel.treedatalikelihood.discrete.DiscreteTraitBranchRateGradient;
import dr.evomodel.treedatalikelihood.discrete.SiteRateModelGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Analytic gradients wrt the gamma shape, the proportion of invariant sites and the relative rate on the b1 fixture
 * (5 taxa, 44 sites, HKY kappa = 2.5, pi = (0.1, 0.2, 0.3, 0.4), strict clock 1.0).
 */
public class SiteRateModelGradientTest extends TestCase {

    private static final String[][] SEQUENCES = {
            {"A", "AAAAAACCCCGGGGGTTTTTTTAACAAAGTTACGARACCCGGTA"},
            {"B", "AAAAAACCCCGGGGGTTTTTTTAACGGCGTTATATAACCCGGCC"},
            {"C", "AAAAAACCCCGGGGGTTTTTTTGGTAAGGCCACGAA-CCCAATC"},
            {"D", "AAAAAACCCCGGGGGTTTTTTTGGTAATATTATATGGCCCGGCC"},
            {"E", "AAAAAACCCCGGGGGTTTTTTTGGTAAAATTGCGAGGTTTGGCA"},
    };
    private static final String NEWICK = "((A:0.12,B:0.12):0.31,(C:0.25,(D:0.07,E:0.07):0.18):0.18);";
    private static final int SITE_COUNT = 44;
    private static final PartialsRescalingScheme[] SCHEMES = {
            PartialsRescalingScheme.NONE, PartialsRescalingScheme.ALWAYS, PartialsRescalingScheme.DYNAMIC};

    public SiteRateModelGradientTest(String name) {
        super(name);
    }

    public void testRelativeRateMatchesFiniteDifference() throws Exception {
        // the relative rate moves every category rate exactly in proportion, so this checks the rate term tightly
        for (PartialsRescalingScheme scheme : SCHEMES) {
            Fixture f = new Fixture(0.5, 0.3, 4, scheme, false);
            assertGradient(f, f.nu, 1E-6, 1E-7);
        }
    }

    public void testRelativeRateMatchesClockDerivative() throws Exception {
        // nu and the strict clock rate both scale every branch, so nu d log L / d nu = r d log L / d r
        Fixture f = new Fixture(0.5, 0.3, 4, PartialsRescalingScheme.NONE, false);
        double viaNu = f.nu.getParameterValue(0) * gradient(f, f.nu).getGradientLogDensity()[0];
        assertEquals(numeric(f, f.clockRate, 1E-6), viaNu, 1E-7 * Math.abs(viaNu));
    }

    public void testShapeMatchesFiniteDifference() throws Exception {
        // log L inherits ~1e-9 jitter in alpha from GammaFunction.incompleteGammaP (AS32, accuracy 1e-8), so a
        // difference of log L checks alpha only to ~1e-4; GammaQuantileDerivativeTest checks the quantile part tightly
        for (double alpha : new double[]{0.2, 0.5, 2.0, 20.0}) {
            assertGradient(new Fixture(alpha, Double.NaN, 4, PartialsRescalingScheme.NONE, false),
                    null, 1E-4, 1E-4);
        }
        for (PartialsRescalingScheme scheme : SCHEMES) {
            assertGradient(new Fixture(0.5, 0.3, 4, scheme, false), null, 1E-4, 1E-4);
        }
    }

    public void testInvariantMatchesFiniteDifference() throws Exception {
        for (PartialsRescalingScheme scheme : SCHEMES) {
            Fixture f = new Fixture(0.5, 0.3, 4, scheme, false);
            assertGradient(f, f.invariant, 1E-6, 1E-7);
        }
        Fixture noGamma = new Fixture(Double.NaN, 0.3, 0, PartialsRescalingScheme.NONE, false);
        assertGradient(noGamma, noGamma.invariant, 1E-6, 1E-7);
    }

    public void testCategoryWeightDerivativesSumToSiteCount() throws Exception {
        // L_s = sum_c w_c L_{s,c}, so sum_c w_c d log L / d w_c = sum_s w_s
        Fixture f = new Fixture(0.5, 0.3, 4, PartialsRescalingScheme.ALWAYS, false);
        f.likelihood.getLogLikelihood();
        double[] derivatives = new double[f.siteRateModel.getCategoryCount()];
        f.delegate.getCategoryWeightDerivatives(f.tree.getRoot().getNumber(), derivatives);
        double sum = 0.0;
        for (int c = 0; c < derivatives.length; c++) {
            sum += f.siteRateModel.getCategoryProportions()[c] * derivatives[c];
        }
        assertEquals(SITE_COUNT, sum, 1E-10);
    }

    public void testBranchRateGradientUnaffected() throws Exception {
        // the site-rate delegate shares the infinitesimal-matrix buffers with the branch-rate delegate
        Fixture f = new Fixture(0.5, 0.3, 4, PartialsRescalingScheme.NONE, true);
        DiscreteTraitBranchRateGradient branchGradient = new DiscreteTraitBranchRateGradient("branchRates",
                f.likelihood, f.delegate, f.branchRates, false);
        branchGradient.getGradientLogDensity();

        gradient(f, f.shape).getGradientLogDensity();
        double[] afterShape = branchGradient.getGradientLogDensity();

        f.branchRates.setParameterValue(2, 1.3); // a change the branch-rate delegate may not see as a model change
        gradient(f, f.shape).getGradientLogDensity();
        double[] afterChange = branchGradient.getGradientLogDensity();
        double[] numeric = numericAll(f, f.branchRates);
        for (int i = 0; i < numeric.length; i++) {
            assertEquals("branch " + i, numeric[i], afterChange[i], 1E-6 * Math.max(1.0, Math.abs(numeric[i])));
        }
        assertEquals(afterShape.length, afterChange.length);
    }

    public void testRefusals() throws Exception {
        Fixture f = new Fixture(0.5, Double.NaN, 4, PartialsRescalingScheme.NONE, false);
        try {
            gradient(f, new Parameter.Default(1.0));
            fail("a parameter outside the site model must be refused");
        } catch (IllegalArgumentException expected) {
            assertTrue(expected.getMessage().contains("not a differentiable parameter"));
        }
    }

    // a null parameter means the fixture's gamma shape
    private static void assertGradient(Fixture f, Parameter parameter, double relativeStep, double tolerance) {
        Parameter p = parameter == null ? f.shape : parameter;
        final double analytic = gradient(f, p).getGradientLogDensity()[0];
        final double numeric = numeric(f, p, relativeStep);
        assertEquals(p.getId() + " = " + p.getParameterValue(0), numeric, analytic,
                tolerance * Math.max(1.0, Math.abs(numeric)));
    }

    private static SiteRateModelGradient gradient(Fixture f, Parameter parameter) {
        return new SiteRateModelGradient(f.likelihood, parameter, 1E-4);
    }

    private static double numeric(Fixture f, Parameter parameter, double relativeStep) {
        return numericAt(f, parameter, 0, relativeStep);
    }

    private static double[] numericAll(Fixture f, Parameter parameter) {
        double[] out = new double[parameter.getDimension()];
        for (int i = 0; i < out.length; i++) {
            out[i] = numericAt(f, parameter, i, 1E-6);
        }
        return out;
    }

    private static double numericAt(Fixture f, Parameter parameter, int i, double relativeStep) {
        final double value = parameter.getParameterValue(i);
        final double h = relativeStep * Math.max(1.0, value);
        parameter.setParameterValue(i, value + h);
        final double plus = f.likelihood.getLogLikelihood();
        parameter.setParameterValue(i, value - h);
        final double minus = f.likelihood.getLogLikelihood();
        parameter.setParameterValue(i, value);
        return (plus - minus) / (2 * h);
    }

    private static class Fixture {

        final Parameter shape;
        final Parameter invariant;
        final Parameter nu = new Parameter.Default("nu", 1.3);
        final Parameter clockRate = new Parameter.Default("clockRate", 1.0);
        final Parameter branchRates;
        final GammaSiteRateModel siteRateModel;
        final DefaultTreeModel tree;
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;

        Fixture(double alpha, double pInv, int gammaCategories, PartialsRescalingScheme scheme,
                boolean arbitraryBranchRates) throws Exception {
            SimpleAlignment alignment = new SimpleAlignment();
            alignment.setDataType(Nucleotides.INSTANCE);
            for (String[] row : SEQUENCES) {
                alignment.addSequence(new Sequence(new Taxon(row[0]), row[1]));
            }
            Tree newick = new NewickImporter(NEWICK).importTree(null);
            tree = new DefaultTreeModel("treeModel", newick);

            shape = Double.isNaN(alpha) ? null : new Parameter.Default("alpha", alpha);
            invariant = Double.isNaN(pInv) ? null : new Parameter.Default("pInv", pInv);
            siteRateModel = new GammaSiteRateModel("siteModel", nu, 1.0, shape, gammaCategories, invariant);

            HKY hky = new HKY(new Parameter.Default(2.5), new FrequencyModel(Nucleotides.INSTANCE,
                    new Parameter.Default(new double[]{0.1, 0.2, 0.3, 0.4})));
            delegate = new BeagleDataLikelihoodDelegate(tree,
                    new SitePatterns(alignment, null, 0, -1, 1, true), new HomogeneousBranchModel(hky),
                    siteRateModel, true, false, scheme, false, new PreOrderSettings(true, true, true, true));

            BranchRateModel branchRateModel;
            if (arbitraryBranchRates) {
                branchRates = new Parameter.Default("branchRates", tree.getNodeCount() - 1, 1.0);
                for (int i = 0; i < branchRates.getDimension(); i++) {
                    branchRates.setParameterValue(i, 0.8 + 0.1 * i);
                }
                branchRateModel = new ArbitraryBranchRates(tree, branchRates,
                        new ArbitraryBranchRates.BranchRateTransform.None(), false);
            } else {
                branchRates = null;
                branchRateModel = new StrictClockBranchRates(clockRate);
            }
            likelihood = new TreeDataLikelihood(delegate, tree, branchRateModel);
        }
    }
}
