/*
 * BaselineIncrementField.java
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

package dr.math.distributions;

import dr.inference.distribution.RandomField;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.model.GradientProvider;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * @author Marc Suchard
 * @author Yucai Shao
 * @author Andy Magee
 */
public class BaselineIncrementField extends RandomFieldDistribution
        implements BayesianBridgeStatisticsProvider {

    public static final String TYPE = "BaselineIncrementField";

    private final Distribution baseline;
    private final Distribution increments;

    private final GradientProvider baselineGradient;
    private final GradientProvider incrementGradient;

    private final BayesianBridgeStatisticsProvider bayesianBridge;

    public BaselineIncrementField(String name,
                                  Distribution baseline,
                                  Distribution increments,
                                  RandomField.WeightProvider weights) {
        super(name);

        this.baseline = baseline;
        this.increments = increments;

        if (baseline instanceof Model) {
            addModel((Model) baseline);
        }
        if (increments instanceof Model) {
            addModel((Model) increments);
        }

        baselineGradient = (baseline instanceof GradientProvider) ?
                (GradientProvider) baseline : null;

        incrementGradient = (increments instanceof GradientProvider) ?
                (GradientProvider) increments : null;

        bayesianBridge = (increments instanceof BayesianBridgeStatisticsProvider) ?
                (BayesianBridgeStatisticsProvider) baseline : null;

        if (weights != null) {
            throw new IllegalArgumentException("Unsure how weights influence this field");
        }
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
//        if (model == baseline || model == increments) {
            // Do nothing
            // TODO do we need a fireModelChangedEvent()?
//        }
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new IllegalArgumentException("Unknown variable");
    }

    @Override
    public double getCoefficient(int i) {
        throw new RuntimeException("Should not be called");
    }

    @Override
    public Parameter getGlobalScale() {
        if (bayesianBridge != null) {
            return bayesianBridge.getGlobalScale();
        } else {
            throw new IllegalArgumentException("Not a Bayesian bridge");
        }
    }

    @Override
    public Parameter getLocalScale() {
        if (bayesianBridge != null) {
            return bayesianBridge.getLocalScale();
        } else {
            throw new IllegalArgumentException("Not a Bayesian bridge");
        }
    }

    @Override
    public Parameter getExponent() {
        if (bayesianBridge != null) {
            return bayesianBridge.getExponent();
        } else {
            throw new IllegalArgumentException("Not a Bayesian bridge");
        }
    }

    @Override
    public Parameter getSlabWidth() {
        if (bayesianBridge != null) {
            return bayesianBridge.getSlabWidth();
        } else {
            throw new IllegalArgumentException("Not a Bayesian bridge");
        }
    }

    @Override
    public int getDimension() {
        if (bayesianBridge != null) {
            return bayesianBridge.getDimension();
        } else {
            throw new IllegalArgumentException("Not a Bayesian bridge");
        }
    }

    @Override
    public double[] getGradientLogDensity(Object x) {
        double[] field = (double[]) x;

        double[] sub = new double[field.length - 1];
        System.arraycopy(field, 1, sub, 0, field.length - 1);

        double[] baselineGrad = baselineGradient.getGradientLogDensity(new double[]{field[0]});
        double[] incrementGrad = incrementGradient.getGradientLogDensity(sub);

        double[] gradient = new double[field.length];
        gradient[0] = baselineGrad[0];
        System.arraycopy(incrementGrad, 0, gradient, 1, field.length - 1);

        return gradient;
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(double[] x) {

        double logPdf = baseline.logPdf(x[0]);
        for (int i = 1; i < x.length; ++i) {
            logPdf += increments.logPdf(x[i]);
        }

        return logPdf;
    }

    @Override
    public double[][] getScaleMatrix() { throw new RuntimeException("Not yet implemented");}

    @Override
    public double[] getMean() { throw new RuntimeException("Not yet implemented"); }

    @Override
    public String getType() { return TYPE; }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getIncrement(int i, Parameter field) {
        return field.getParameterValue(i + 1);
    }
}
