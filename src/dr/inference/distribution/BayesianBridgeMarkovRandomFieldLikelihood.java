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

package dr.inference.distribution;

import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.BayesianBridgeStatisticsProvider;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.*;
import dr.util.InverseFirstOrderFiniteDifferenceTransform;
import dr.xml.Reportable;

import static dr.inferencexml.distribution.shrinkage.BayesianBridgeLikelihoodParser.BAYESIAN_BRIDGE;

public class BayesianBridgeMarkovRandomFieldLikelihood extends AbstractModelLikelihood implements
        BayesianBridgeStatisticsProvider, PriorPreconditioningProvider,
        GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    public BayesianBridgeMarkovRandomFieldLikelihood(Parameter variables,
                                                                BayesianBridgeDistributionModel bridge,
                                                                ParametricDistributionModel firstElementDistribution,
                                                                InverseFirstOrderFiniteDifferenceTransform transform) {

        super(BAYESIAN_BRIDGE);

        this.variables = variables;
        this.bridge = bridge;
        this.firstElementDistribution = firstElementDistribution;
        this.transform = transform;
        this.dim = variables.getDimension();

        addModel(bridge);
        addModel(firstElementDistribution);
        addVariable(variables);
    }

    private double[][] getUnconstrainedValues() {
        double[] concatenated = transform.inverse(variables.getParameterValues(), 0, dim);

        double[][] values = new double[2][];
        values[0] = new double[1];
        values[0][0] = concatenated[0];
        values[1] = new double[dim - 1];
        for (int i = 0; i < dim - 1; i++) {
            values[1][i] = concatenated[i+1];
        }
        return values;
    }

//    private double[] getFirstElementVariableValues() {
//        double[] vals = new double[1];
//        vals[0] = variables.getParameterValue(0);
//        return vals;
//    }
//
//    private double[] getBridgeVariableValues() {
//        double[] vals = new double[dim - 1];
//        for (int i = 0; i < dim - 1; i++) {
//            vals[i] = variables.getParameterValue(i + 1);
//        }
//        return vals;
//    }

    @Override
    public double getLogLikelihood() {
        double[][] unconstrained = getUnconstrainedValues();
        double logPdf = 0.0;
        logPdf += firstElementDistribution.logPdf(unconstrained[0]);
        logPdf += bridge.logPdf(unconstrained[1]);
        logPdf -= transform.getLogJacobian(transform.inverse(variables.getParameterValues(),0,dim));
        return logPdf;
    }

    @Override
    public double[] getGradientLogDensity() {
        double[][] transformedVariables = getUnconstrainedValues();
        double[] grad = new double[dim];

        grad[0] = ((GradientProvider)firstElementDistribution).getGradientLogDensity(transformedVariables[0])[0];

        double[] bridgeGrad = bridge.getGradientLogDensity(transformedVariables[1]);
        for (int i = 0; i < dim - 1; i++) {
            grad[i + 1] = bridgeGrad[i];
        }

        return transform.updateGradientLogDensity(grad, transform.inverse(variables.getParameterValues(), 0, dim), 0, dim);

    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
//        return distribution.getDiagonalHessianLogDensity(coefficients.getParameterValues());
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    /*********************
    * BayesianBridgeStatisticsProvider interface, so we can use a BayesianBridgeShrinkageOperator
    *********************/
    @Override
    public double getCoefficient(int i) {
        // This treats the first, non-Bridge, variable as a Bridge Coefficient
        // There must therefore be a mask on this variable to prevent it contributing to the Gibbs sampler
        return transform.inverse(variables.getParameterValues(),0,dim)[i];
    }

    @Override
    public Parameter getGlobalScale() {
        return bridge.getGlobalScale();
    }

    @Override
    public Parameter getLocalScale() {
        return bridge.getLocalScale();
    }

    @Override
    public Parameter getExponent() {
        return bridge.getExponent();
    }

    @Override
    public Parameter getSlabWidth() {
        return bridge.getSlabWidth();
    }

    @Override
    // This may cause problems since the Bridge's dimension is smaller than the model's...
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
        return variables;
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


    /*********************
     * PriorPreconditioningProvider interface
     *********************/
    @Override
    // NB: Prior preconditioning only makes sense if the HMC operator is using a firstOrderFiniteDifferenceTransform
    public double getStandardDeviation(int index) {
        if (bridge instanceof PriorPreconditioningProvider && firstElementDistribution instanceof PriorPreconditioningProvider) {
            if (index == 0) {
                return ((PriorPreconditioningProvider) firstElementDistribution).getStandardDeviation(index);
            } else {
                return ((PriorPreconditioningProvider) bridge).getStandardDeviation(index - 1);
            }
        } else {
            throw new RuntimeException("Not a prior conditioner");
        }
    }

    private final Parameter variables;
    private final BayesianBridgeDistributionModel bridge;
    private final ParametricDistributionModel firstElementDistribution;
    private final int dim;
    private final InverseFirstOrderFiniteDifferenceTransform transform;

    /*********************
     * Reportable interface
     *********************/
    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1e-6, 1E-3);
    }
}
