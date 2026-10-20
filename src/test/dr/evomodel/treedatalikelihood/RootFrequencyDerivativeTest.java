/*
 * RootFrequencyDerivativeTest.java
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

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

/**
 * Root term of the frequency gradient, R_k = sum_s w_s (sum_c w_c p_root,c,s,k) / L_s.
 */
public class RootFrequencyDerivativeTest extends TraceCorrelationAssert {

    private static final double[] PI = {0.1, 0.2, 0.3, 0.4};
    private static final double[] ROOT_PI = {0.35, 0.15, 0.2, 0.3};

    public RootFrequencyDerivativeTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testWeightedSumEqualsSiteCount() {
        for (boolean gamma : new boolean[]{false, true}) {
            Fixture f = new Fixture(PartialsRescalingScheme.NONE, gamma, false);
            double[] r = f.rootDerivatives();
            double sum = 0.0;
            for (int k = 0; k < r.length; k++) {
                sum += PI[k] * r[k];
            }
            double totalWeight = 0.0;
            for (double w : f.delegate.getPatternWeights()) {
                totalWeight += w;
            }
            assertEquals(alignment.getSiteCount(), totalWeight, 0.0);
            assertEquals("gamma=" + gamma, totalWeight, sum, 1E-12 * totalWeight);
        }
    }

    public void testInvariantUnderRescaling() {
        for (boolean gamma : new boolean[]{false, true}) {
            double[] reference = new Fixture(PartialsRescalingScheme.NONE, gamma, false).rootDerivatives();
            for (PartialsRescalingScheme scheme : new PartialsRescalingScheme[]{
                    PartialsRescalingScheme.ALWAYS, PartialsRescalingScheme.DYNAMIC}) {
                double[] r = new Fixture(scheme, gamma, false).rootDerivatives();
                for (int k = 0; k < r.length; k++) {
                    assertEquals(scheme + " gamma=" + gamma + " k=" + k, reference[k], r[k],
                            1E-10 * Math.abs(reference[k]));
                }
            }
        }
    }

    public void testSeparateRootFrequenciesMatchFiniteDifference() {
        // Q does not depend on the root frequencies, so R is the full gradient wrt them (off the simplex)
        final double h = 1E-6;
        for (boolean gamma : new boolean[]{false, true}) {
            Fixture f = new Fixture(PartialsRescalingScheme.ALWAYS, gamma, true);
            double[] r = f.rootDerivatives();
            for (int k = 0; k < r.length; k++) {
                f.rootFrequencies.setParameterValue(k, ROOT_PI[k] + h);
                double plus = f.likelihood.getLogLikelihood();
                f.rootFrequencies.setParameterValue(k, ROOT_PI[k] - h);
                double minus = f.likelihood.getLogLikelihood();
                f.rootFrequencies.setParameterValue(k, ROOT_PI[k]);
                f.likelihood.getLogLikelihood();

                assertEquals("gamma=" + gamma + " k=" + k, (plus - minus) / (2 * h), r[k], 1E-6 * Math.abs(r[k]));
            }
        }
    }

    public void testDoesNotChangeSiteLogLikelihoods() {
        Fixture f = new Fixture(PartialsRescalingScheme.ALWAYS, true, false);
        f.likelihood.getLogLikelihood();
        double[] before = f.delegate.getSiteLogLikelihoods();
        f.rootDerivatives();
        double[] after = f.delegate.getSiteLogLikelihoods();
        assertEquals(before.length, after.length);
        for (int s = 0; s < before.length; s++) {
            assertEquals(Double.doubleToLongBits(before[s]), Double.doubleToLongBits(after[s]));
        }
    }

    private class Fixture {

        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;
        final DefaultTreeModel tree;
        final Parameter rootFrequencies;

        Fixture(PartialsRescalingScheme scheme, boolean gammaWithInvariant, boolean separateRoot) {
            tree = new DefaultTreeModel(ladderTree());

            HKY hky = new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE, PI.clone()));
            HomogeneousBranchModel branchModel;
            if (separateRoot) {
                rootFrequencies = new Parameter.Default(ROOT_PI.clone());
                branchModel = new HomogeneousBranchModel(hky, new FrequencyModel(Nucleotides.INSTANCE, rootFrequencies));
            } else {
                rootFrequencies = null;
                branchModel = new HomogeneousBranchModel(hky);
            }

            GammaSiteRateModel siteRateModel = gammaWithInvariant
                    ? new GammaSiteRateModel("siteModel", 0.5, 4, 0.2)
                    : new GammaSiteRateModel("siteModel");
            SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

            delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel, siteRateModel,
                    true, false, scheme, false, PreOrderSettings.getDefault());
            likelihood = new TreeDataLikelihood(delegate, tree,
                    new StrictClockBranchRates(new Parameter.Default(1.0)));
        }

        double[] rootDerivatives() {
            likelihood.getLogLikelihood();
            double[] out = new double[Nucleotides.INSTANCE.getStateCount()];
            delegate.getRootFrequencyDerivatives(tree.getRoot().getNumber(), out);
            return out;
        }

        private SimpleTree ladderTree() {
            SimpleNode current = new SimpleNode();
            current.setTaxon(taxa[0]);
            current.setHeight(0.0);
            for (int i = 1; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.0);
                SimpleNode parent = new SimpleNode();
                parent.setHeight(0.05 * i);
                parent.addChild(current);
                parent.addChild(tip);
                current = parent;
            }
            return new SimpleTree(current);
        }
    }
}
