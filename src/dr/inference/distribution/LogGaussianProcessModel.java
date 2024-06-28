/*
 * LogGaussianProcessModel.java
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

package dr.inference.distribution;

import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.inference.model.*;
import org.apache.commons.math.linear.*;

import java.util.List;


/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

// TODO this class is not really a LogLinearModel nor GeneralizedLinearModel; need to disassociate.
// TODO disassocation requires refactoring
// TODO substitute LogLinearModel with GeneralizedAdditiveModel
public class LogGaussianProcessModel extends LogLinearModel
        implements LogAdditiveCtmcRateProvider.DataAugmented {

    //private final Parameter kernelParameter;
    private final List<DesignMatrix> designMatrices;
    //private final Parameter designMatrix;
    private final Parameter gaussianNoise;

    private boolean precisionKnown;
    private boolean vInvKnown;

    private RealMatrix vInvReal; // changes ......
    private RealMatrix predictiveKernel;
    private RealMatrix predictiveInvVariance; //changes ....
    private double logDetinvV;
    private final GaussianProcessKernel kernel;

    private final Parameter realizedField;

    // Constructor
    public LogGaussianProcessModel(Parameter dependentParameter,
                                   Parameter gaussianNoise,
                                   Parameter kernelParameter,
                                   List<DesignMatrix> designMatrices) {
        super(dependentParameter);
        this.gaussianNoise = gaussianNoise;
        addVariable(gaussianNoise);
        //this.kernelParameter = kernelParameter;
        //addVariable(kernelParameter);

        this.realizedField = dependentParameter;
        addVariable(realizedField);

        this.designMatrices = designMatrices;
        for (DesignMatrix matrix : designMatrices) {
            addVariable(matrix);
        }

        precisionKnown = false;

        this.kernel = new GaussianProcessKernel("name", kernelParameter);
        addModel(kernel);
    }

    public Parameter getFieldParameter() { return realizedField; }

    @Override
    public Parameter getLogRateParameter() { return realizedField; }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == kernel) {
            precisionKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        if (variable == kernelParameter) {
//            precisionKnown = false;
//        } else
        if (variable == dependentParam) {
            precisionKnown = false; //  TODO CHECK
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    public static class GaussianProcessKernel extends AbstractModel {

        private final Parameter parameter;
        public GaussianProcessKernel(String name,
                                     Parameter parameter) {
            super(name);
            this.parameter = parameter;

            addVariable(parameter);
        }

        double getCorrelation(double x, double y) {
            double sigma = parameter.getParameterValue(1);
            double length = parameter.getParameterValue(2);
            return x * y;
//            sigma * Math.exp(-Math.pow(x-y,2)/(2*Math.pow(length, 2)));
//
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            if (variable == parameter) {
                fireModelChanged(parameter);
            }
        }

        @Override
        protected void storeState() {
        }

        @Override
        protected void restoreState() {
        }

        @Override
        protected void acceptState() {
        }
    }


    @Override
    protected double calculateLogLikelihood() {
        double[] rates = dependentParam.getParameterValues();
        RealMatrix ratesReal = MatrixUtils.createColumnRealMatrix(rates);

        if (!precisionKnown) {
            double noiseVar = gaussianNoise.getParameterValue(0);
            precision(noiseVar);
            logDetinvV = Math.log(predictiveInvVariance.getDeterminant());
            precisionKnown = true;
        }

        // TODO this seems incorrect (the GP mean is 0, no? ... and there's no predictive variance)
//        RealMatrix muReal = predictiveKernel.multiply(vInvReal).multiply(ratesReal);
//        RealMatrix ratesCenteredReal = ratesReal.subtract(muReal);
//        double exponent = ratesCenteredReal.transpose().multiply(predictiveInvVariance).multiply(ratesCenteredReal).getEntry(0, 0);
        double exponent = ratesReal.transpose().multiply(vInvReal).multiply(ratesReal).getEntry(0,0);
        
        // Compute the log likelihood
        double logLikelihood = -0.5 * (N * Math.log(2 * Math.PI) - logDetinvV) - 0.5 * exponent ;
        return logLikelihood;
    }

    public double[][] kernel() {
       // double[] kernelParameters = kernelParameter.getParameterValues();
        //double sigma = kernelParameters[1];
        //double length = kernelParameters[2];
        int P = designMatrices.size();
        double[][] K = new double[N][N];
        for (int m = 0; m < P; m++) {
            double[] X = designMatrices.get(m).getParameterValues();
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    K[i][j] += kernel.getCorrelation(X[i], X[j]); // X[i] * X[j];
                    //K[i][j] += sigma * Math.exp(-Math.pow(X[i]-X[j],2)/(2*Math.pow(length, 2)));
                }
            }
        }
        return K;
    }

    public void precision(double noiseVar) {
        // kernel from starting predictors
        double[][] kDouble = kernel();
        RealMatrix predictorKernel = MatrixUtils.createRealMatrix(kDouble);

        // error noise variance matrix
        double[][] noise = new double[N][N];
        for (int i = 0; i < N; i++) {
            noise[i][i] = noiseVar;
        }

        RealMatrix noiseReal = MatrixUtils.createRealMatrix(noise);
        vInvReal = predictorKernel.add(noiseReal).inverse();

        // kernel from additional predictors in the posterior
        predictiveKernel = MatrixUtils.createRealMatrix(kDouble); /// CHANGE FOR NEW VARS TODO

        double[][] k2Double = kernel(); //change when adding new x's /// CHANGE FOR NEW VARS TODO
        RealMatrix k2 = MatrixUtils.createRealMatrix(k2Double);

        RealMatrix predictiveVarianceReal = k2.subtract(predictorKernel.multiply(vInvReal).multiply(predictiveKernel.transpose()));
        predictiveInvVariance = predictiveVarianceReal.inverse();
    }

    @Override
    public double[] getXBeta() {
        // compute the mean and then exponentiate
        final int fieldDim = dependentParam.getDimension();
        double[] rates = new double[fieldDim];

        // here we just exponentiate the log-mean rate into the actual mean rates
        for (int i = 0; i < fieldDim; i++) {
            rates[i] = Math.exp(dependentParam.getParameterValue(i));
        }
        return rates;
    }

}