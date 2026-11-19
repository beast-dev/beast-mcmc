/*
 * CtmcClockRateGradientTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package test.dr.evomodel.treedatalikelihood;

import beagle.BeagleInfo;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.BeagleFunctionality;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode;
import dr.evomodel.treedatalikelihood.discrete.CtmcClockRateGradient;
import dr.evomodel.treedatalikelihood.discrete.LogCtmcRateGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * The gradient with respect to a strict clock rate, formed from the cross products that also give the gradient with
 * respect to the substitution rates, against central differences of the log-likelihood. The rate matrix is
 * asymmetric and normalized, and the gamma categories rescale the branches, which r d logL / d r = sum_ij Q_ij
 * (d logL / d Q)_ij must allow for.
 */
public class CtmcClockRateGradientTest extends TestCase {

    public CtmcClockRateGradientTest(String name) {
        super(name);
    }

    @Override
    public void runBare() throws Throwable {
        if (!BeagleFunctionality.IS_MATRIX_VECTOR_SUPPORTED()) {
            System.err.println("SKIPPED " + getClass().getSimpleName() + "." + getName() +
                    ": needs native BEAGLE 4.2.0 or later, found " + BeagleInfo.getVersion());
            return;
        }
        super.runBare();
    }

    public void testClockRateGradient() {

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(Nucleotides.INSTANCE);
        SimpleNode root = parent(
                parent(tip(alignment, "a", "ACGTA"), tip(alignment, "b", "ACGCA"), 0.3),
                parent(
                        parent(tip(alignment, "c", "GCTTA"), tip(alignment, "d", "GCTTC"), 0.2),
                        tip(alignment, "e", "TAGTC"), 0.7),
                1.1);
        DefaultTreeModel tree = new DefaultTreeModel(new SimpleTree(root));

        double[] logRates = new double[12];
        for (int k = 0; k < logRates.length; ++k) {
            logRates[k] = 0.3 * (k % 5) - 0.6;
        }
        LogRateSubstitutionModel model = new LogRateSubstitutionModel("model", Nucleotides.INSTANCE,
                new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.3, 0.2, 0.15, 0.35}),
                new LogAdditiveCtmcRateProvider.DataAugmented.Basic("logRates",
                        new Parameter.Default("logRates", logRates)));
        model.setNormalization(true);
        model.setScaleRatesByFrequencies(false);

        Parameter clockRate = new Parameter.Default("clock.rate", 0.8);
        StrictClockBranchRates branchRates = new StrictClockBranchRates(clockRate);

        BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree,
                new SitePatterns(alignment), new HomogeneousBranchModel(model),
                new GammaSiteRateModel("siteModel", 0.5, 4), false, false, PartialsRescalingScheme.DEFAULT, true,
                new PreOrderSettings(true, false, false, false, true, false, false));
        TreeDataLikelihood likelihood = new TreeDataLikelihood(delegate, tree, branchRates);

        // the substitution-rate gradient first, so that the clock-rate gradient reads the cross products it caches
        LogCtmcRateGradient rateGradient = new LogCtmcRateGradient("test", likelihood, delegate, model,
                ApproximationMode.EXACT_SPECTRAL, false);
        CtmcClockRateGradient clockGradient = new CtmcClockRateGradient("test", likelihood, delegate, branchRates,
                ApproximationMode.EXACT_SPECTRAL);

        for (double rate : new double[]{0.8, 2.5}) {
            clockRate.setParameterValue(0, rate);
            rateGradient.getGradientLogDensity();
            final double gradient = clockGradient.getGradientLogDensity()[0];

            final double h = 1E-5 * rate;
            clockRate.setParameterValue(0, rate + h);
            final double up = likelihood.getLogLikelihood();
            clockRate.setParameterValue(0, rate - h);
            final double down = likelihood.getLogLikelihood();
            clockRate.setParameterValue(0, rate);
            final double numeric = (up - down) / (2 * h);

            assertEquals("clock rate " + rate, numeric, gradient, 1E-6 * Math.max(1.0, Math.abs(numeric)));
        }
    }

    private static SimpleNode tip(SimpleAlignment alignment, String name, String states) {
        Taxon taxon = new Taxon(name);
        alignment.addSequence(new Sequence(taxon, states));
        SimpleNode tip = new SimpleNode();
        tip.setTaxon(taxon);
        tip.setHeight(0.0);
        return tip;
    }

    private static SimpleNode parent(SimpleNode left, SimpleNode right, double height) {
        SimpleNode parent = new SimpleNode();
        parent.setHeight(height);
        parent.addChild(left);
        parent.addChild(right);
        return parent;
    }
}
