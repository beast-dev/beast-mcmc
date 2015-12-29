/*
 * OrnsteinUhlenbeckPriorLikelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.coalescent;

import dr.evomodelxml.coalescent.OrnsteinUhlenbeckPriorLikelihoodParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.NormalDistribution;

/**
 * Ornstein-Uhlenbeck prior.
 * <p/>
 * A diffusion process - Basically a correlated sequence of normally distributed values, where correlation is time
 * dependent.
 * <p/>
 * <p/>
 * Very experimental and only slightly tested at this time
 *
 * @author Joseph Heled
 *         Date: 25/04/2008
 */

// It should be a model since it may be the only user of parameter sigma

public class OrnsteinUhlenbeckPriorLikelihood extends AbstractModelLikelihood {
    private final Parameter mean;
    private final Parameter sigma;
    private final Parameter lambda;
    private final boolean logSpace;
    private final boolean normalize;
    private Parameter data;
    private Parameter times;
    private ParametricDistributionModel popMeanPrior = null;
    private VariableDemographicModel m = null;

    private OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                             boolean logSpace, boolean normalize) {
        super(OrnsteinUhlenbeckPriorLikelihoodParser.OU);

        this.logSpace = logSpace;
        this.normalize = normalize;

        this.mean = mean;
        if (mean != null) {
            addVariable(mean);
        }

        this.sigma = sigma;
        addVariable(sigma);

        this.lambda = lambda;
        if (lambda != null) {
            addVariable(lambda);
        }
    }

    public OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                            VariableDemographicModel demographicModel, boolean logSpace, boolean normalize,
                                            ParametricDistributionModel popMeanPrior) {
        this(mean, sigma, lambda, logSpace, normalize);
        this.m = demographicModel;

        this.data = this.times = null;
        this.popMeanPrior = popMeanPrior;
    }

    public OrnsteinUhlenbeckPriorLikelihood(Parameter mean, Parameter sigma, Parameter lambda,
                                            Parameter dataParameter, Parameter timesParameter, boolean logSpace, boolean normalize) {
        this(mean, sigma, lambda, logSpace, normalize);
        dataParameter.addParameterListener(this);
        timesParameter.addParameterListener(this);
        this.data = dataParameter;
        this.times = timesParameter;
    }

    // log of normal distribution coefficient.
    final private double logNormalCoef = -0.5 * Math.log(2 * Math.PI);

    // A specialized version where everything is normalized. Time is normalized to 1. Data moved to mean zero and rescaled
    // according to time. Lambda is 0.5. The prior on mean is added.

    private double reNormalize(VariableDemographicModel m) {
        final double[] tps = m.getDemographicFunction().allTimePoints();
        // get a copy since we re-scale data
        final double[] vals = m.getPopulationValues().getParameterValues();

        assert !logSpace : "not implemented yet";

        final double len = tps[tps.length - 1];

        // compute mean
        double popMean = 0;
        // todo not correct when using midpoints
        if (m.getType() == VariableDemographicModel.Type.LINEAR) {
            for (int k = 0; k < tps.length; ++k) {
                final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0));
                popMean += dt * (vals[k + 1] + vals[k]);
            }
            popMean /= (2 * len);
        } else {
            for (int k = 0; k < tps.length; ++k) {
                final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0));
                popMean += dt * vals[k];
            }
            popMean /= len;
        }

        // Normalize to time interval = 1 and mean = 0
        final double sigma = this.sigma.getStatisticValue(0) / Math.sqrt(len);
        final double lam = 0.5 * len;
        for (int k = 0; k < vals.length; ++k) {
            vals[k] = (vals[k] - popMean) / len;
        }

        // optimized version of the code in getLogLikelihood.
        // get factors out when possible. logpdf of a normal is -x^2/2, when mean is 0
        double ll = 0.0;

        for (int k = 0; k < tps.length; ++k) {
            final double dt = (tps[k] - (k > 0 ? tps[k - 1] : 0)) / len;

            final double a = Math.exp(-lam * dt);
            final double d = (vals[k + 1] - vals[k] * a);

            final double c = 1 - a * a;
            ll += (d * d / c) - 0.5 * Math.log(c);
        }

        final double f2 = (2 * lam) / (sigma * sigma);
        ll = tps.length * (logNormalCoef - Math.log(sigma)) + ll * f2 / -2;

        if (popMeanPrior != null) {
            ll += popMeanPrior.logPdf(popMean);
        } else {
            // default Jeffreys
            ll -= Math.log(popMean);
        }
        return ll;
    }

    public double getLogLikelihood() {
        if (lastValue > 0) {
            return lastValue;
        }
        double logL;

        if (normalize) {
            assert m != null;
            logL = reNormalize(m);
        } else {

//        final double[] tps = m != null ? m.getDemographicFunction().allTimePoints() : times.getParameterValues();
//
//        final double[] vals = m != null ? m.getPopulationValues().getParameterValues() : data.getParameterValues();

            final double[] tps = times.getParameterValues();
            final double[] vals = data.getParameterValues();

            if (logSpace) {
                for (int k = 0; k < vals.length; ++k) {
                    vals[k] = Math.log(vals[k]);
                }
            }

            final double lambda = this.lambda.getStatisticValue(0);
            final double mean = this.mean.getStatisticValue(0);

            double sigma = this.sigma.getStatisticValue(0);
            if (normalize) {
                // make the process have a SD of sigma
                sigma *= Math.sqrt(2 * lambda);
            }

            logL = NormalDistribution.logPdf(vals[0], mean, sigma);

            final double xScale = -lambda * (normalize ? 1.0 / tps[tps.length - 1] : 1.0);

            for (int k = 0; k < tps.length; ++k) {
                double dt = tps[k] - (k > 0 ? tps[k - 1] : 0);
                double a = Math.exp(xScale * dt);
                double den = sigma * Math.sqrt((1 - a * a) / (2 * lambda));
                double z = (vals[k + 1] - (vals[k] * a + mean * (1 - a))) / den;
                logL += NormalDistribution.logPdf(z, 0, 1);
            }
        }

        lastValue = logL;

        return logL;
    }

    public void makeDirty() {
        lastValue = -1;
    }

    // simply saves last value
    double lastValue = -1;

    public Model getModel() {
        return this;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }
}