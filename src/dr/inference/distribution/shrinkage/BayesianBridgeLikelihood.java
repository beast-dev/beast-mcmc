/*
 * MultivariateNormalDistributionModel.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution.shrinkage;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.*;

import static dr.inferencexml.distribution.shrinkage.BayesianBridgeLikelihoodParser.BAYESIAN_BRIDGE;

/**
 * A model for scaled-mixture-of-normals distributed data.
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Xiang Ji
 */

public class BayesianBridgeLikelihood extends AbstractModelLikelihood implements
        BayesianBridgeStatisticsProvider, PriorPreconditioningProvider,
        GradientWrtParameterProvider, HessianWrtParameterProvider {

    public BayesianBridgeLikelihood(Parameter coefficients,
                                    BayesianBridgeDistributionModel distribution) {

        super(BAYESIAN_BRIDGE);

        this.coefficients = coefficients;
        this.distribution = distribution;
        this.dim = coefficients.getDimension();

        addModel(distribution);
        addVariable(coefficients);
    }

    public Parameter getGlobalScale() { return distribution.getGlobalScale(); }

    public Parameter getExponent() { return distribution.getExponent(); }

    public Parameter getLocalScale() {return distribution.getLocalScale(); }

    public Parameter getSlabWidth() {return distribution.getSlabWidth(); }

    public double getCoefficient(int i) { return coefficients.getParameterValue(i); }

    @Override
    public double getLogLikelihood() {
        return distribution.logPdf(coefficients.getParameterValues());
    }

    @Override
    public double[] getGradientLogDensity() {
        return distribution.gradientLogPdf(coefficients.getParameterValues());
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return distribution.getDiagonalHessianLogDensity(coefficients.getParameterValues());
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return coefficients;
    }

    @Override
    public final void makeDirty() {
        // Do nothing
    }

    @Override
    public final void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    @Override
    public final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
       // no intermediates need to be recalculated...
    }

    @Override
    public void storeState() {
       // Do nothing
    }

    @Override
    public void restoreState() {
        // Do nothing
    }

    @Override
    public void acceptState() {
    } // no additional state needs accepting

    @Override
    public double getStandardDeviation(int index) {
        if (distribution instanceof PriorPreconditioningProvider) {
            return ((PriorPreconditioningProvider) distribution).getStandardDeviation(index);
        } else {
            throw new RuntimeException("Not a prior conditioner");
        }
    }

    private final Parameter coefficients;
    private final BayesianBridgeDistributionModel distribution;
    private final int dim;
}
