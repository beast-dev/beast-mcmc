/*
 * EpochSubstitutionGradientTest.java
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
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GlmSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.AbstractLogAdditiveSubstitutionModelGradient.ApproximationMode;
import dr.evomodel.treedatalikelihood.discrete.FixedEffectSubstitutionModelGradient;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.distribution.LogLinearModel;
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Checks the substitution-model cross products of an epoch model with degree-2 nodes at the epoch transition times.
 * <p>
 * Each epoch has its own unnormalized GLM substitution model with an intercept, a design column of ones. Changing the
 * intercept of an epoch scales its generator, which is the same as scaling the length of every segment in the epoch.
 * So the first-order gradient with respect to the intercept is exact, and it is compared with central differences of
 * the log likelihood. This requires the cross products of each epoch to cover exactly the segments in that epoch.
 *
 * @author Marc A Suchard
 */
public class EpochSubstitutionGradientTest extends TraceCorrelationAssert {

    private static final int STATES = 4;
    private static final int RATES = STATES * (STATES - 1);

    public EpochSubstitutionGradientTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testInterceptGradient() {
        PartialsRescalingScheme[] schemes = {PartialsRescalingScheme.NONE, PartialsRescalingScheme.DYNAMIC,
                PartialsRescalingScheme.ALWAYS};

        for (int seed = 0; seed < 3; ++seed) {
            for (PartialsRescalingScheme scheme : schemes) {
                for (int categories : new int[]{1, 4}) {
                    checkInterceptGradient(new Fixture(seed, scheme, categories),
                            "seed " + seed + ", " + scheme.getText() + ", " + categories + " categories");
                }
            }
        }
    }

    private void checkInterceptGradient(Fixture f, String label) {

        final double h = 1E-5;
        for (int e = 0; e < f.gradients.size(); ++e) {
            final double[] gradient = f.gradients.get(e).getGradientLogDensity();
            assertEquals(2, gradient.length); // the intercept, then the covariate

            final Parameter intercept = f.intercepts.get(e);
            final double value = intercept.getParameterValue(0);

            intercept.setParameterValue(0, value + h);
            final double up = f.likelihood.getLogLikelihood();
            intercept.setParameterValue(0, value - h);
            final double down = f.likelihood.getLogLikelihood();
            intercept.setParameterValue(0, value);

            final double numerical = (up - down) / (2 * h);
            assertTrue(label + ", epoch " + e + " is not used", numerical != 0.0);
            assertEquals(label + ", epoch " + e, numerical, gradient[0],
                    1E-5 * Math.max(1.0, Math.abs(numerical)));
        }
    }

    private class Fixture {

        final TreeDataLikelihood likelihood;
        final List<Parameter> intercepts = new ArrayList<Parameter>();
        final List<FixedEffectSubstitutionModelGradient> gradients =
                new ArrayList<FixedEffectSubstitutionModelGradient>();

        Fixture(long seed, PartialsRescalingScheme scheme, int categories) {
            Random random = new Random(seed);

            DefaultTreeModel tree = new DefaultTreeModel(randomTree(random));

            final int epochCount = 3 + random.nextInt(3);
            final double[] times = new double[epochCount - 1];
            double time = 0.0;
            for (int i = 0; i < times.length; ++i) {
                time += 0.05 + 0.2 * random.nextDouble();
                times[i] = time;
            }

            double[] ones = new double[RATES];
            Arrays.fill(ones, 1.0);

            List<GlmSubstitutionModel> glms = new ArrayList<GlmSubstitutionModel>();
            for (int e = 0; e < epochCount; ++e) {
                double[] frequencies = new double[STATES];
                double sum = 0.0;
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] = 0.5 + random.nextDouble();
                    sum += frequencies[j];
                }
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] /= sum;
                }

                double[] covariate = new double[RATES];
                for (int k = 0; k < RATES; ++k) {
                    covariate[k] = random.nextGaussian();
                }

                Parameter intercept = new Parameter.Default("intercept." + e, 1, 0.5 * random.nextGaussian());
                Parameter effect = new Parameter.Default("effect." + e, 1, 0.3 * random.nextGaussian());

                LogLinearModel glm = new LogLinearModel(null);
                glm.addIndependentParameter(intercept, new DesignMatrix("ones." + e,
                        new Parameter[]{new Parameter.Default(ones)}, false), null);
                glm.addIndependentParameter(effect, new DesignMatrix("covariate." + e,
                        new Parameter[]{new Parameter.Default(covariate)}, false), null);

                GlmSubstitutionModel model = new GlmSubstitutionModel("glm." + e, Nucleotides.INSTANCE,
                        new FrequencyModel(Nucleotides.INSTANCE, frequencies), glm);
                model.setNormalization(false);

                glms.add(model);
                intercepts.add(intercept);
            }

            EpochBranchModel branchModel = new EpochBranchModel(tree,
                    new ArrayList<SubstitutionModel>(glms), new Parameter.Default(times));
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);

            SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);

            Parameter rates = new Parameter.Default(tree.getNodeCount() - 1, 1.0);
            for (int i = 0; i < rates.getDimension(); ++i) {
                rates.setParameterValue(i, 0.5 + random.nextDouble());
            }
            ArbitraryBranchRates branchRates = new ArbitraryBranchRates(tree, rates,
                    ArbitraryBranchRates.make(false, false, false), false);

            BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel,
                    siteRateModel, false, false, scheme, false,
                    new PreOrderSettings(true, false, false, false, false, false, true));

            likelihood = new TreeDataLikelihood(delegate, tree, branchRates);

            for (GlmSubstitutionModel model : glms) {
                gradients.add(new FixedEffectSubstitutionModelGradient("test", likelihood, delegate, model,
                        ApproximationMode.FIRST_ORDER));
            }
        }

        private SimpleTree randomTree(Random random) {
            List<SimpleNode> active = new ArrayList<SimpleNode>();
            for (int i = 0; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.0);
                active.add(tip);
            }

            while (active.size() > 1) {
                SimpleNode left = active.remove(random.nextInt(active.size()));
                SimpleNode right = active.remove(random.nextInt(active.size()));

                SimpleNode parent = new SimpleNode();
                parent.setHeight(Math.max(left.getHeight(), right.getHeight()) + 0.02 + 0.2 * random.nextDouble());
                parent.addChild(left);
                parent.addChild(right);
                active.add(parent);
            }

            return new SimpleTree(active.get(0));
        }
    }
}
