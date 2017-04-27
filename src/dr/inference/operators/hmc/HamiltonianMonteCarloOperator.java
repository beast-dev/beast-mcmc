/*
 * HMCOperator.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.operators.hmc;

import dr.inference.model.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.math.distributions.NormalDistribution;

/**
 * @author Max Tolkoff
 */
public class HamiltonianMonteCarloOperator extends AbstractCoercableOperator {
    GradientWrtParameterProvider derivative;
    Parameter parameter;
    double stepSize;
    int nSteps;
    NormalDistribution drawDistribution;


    public HamiltonianMonteCarloOperator(CoercionMode mode, double weight, GradientWrtParameterProvider derivative, Parameter parameter, double stepSize, int nSteps, double drawVariance) {
        super(mode);
        setWeight(weight);
        this.derivative = derivative;
        this.parameter = parameter;
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "HMC Operator";
    }

    @Override
    public double doOperation() { //throws OperatorFailedException {
        double functionalStepSize = stepSize;

        double[] HMCDerivative = derivative.getGradientLogDensity(); /* Sign change */
        double[] momentum = new double[HMCDerivative.length];
        for (int i = 0; i < momentum.length; i++) {
            momentum[i] = (Double) drawDistribution.nextRandom();
        }

        double prop=0;
        for (int i = 0; i < momentum.length ; i++) {
            prop += momentum[i] * momentum[i] / (2 * Math.pow(drawDistribution.getSD(), 2));
        }


        for (int i = 0; i < momentum.length; i++) {
                momentum[i] = momentum[i] + functionalStepSize / 2 * HMCDerivative[i];  /* Sign change */
        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j < momentum.length; j++) {
                    parameter.setParameterValue(j, parameter.getParameterValue(j) + functionalStepSize * momentum[j] / (Math.pow(drawDistribution.getSD() , 2)));
            }
//            parameter.fireParameterChangedEvent();

            HMCDerivative = derivative.getGradientLogDensity(); /* Sign change */

            if(i != nSteps){

                for (int j = 0; j < momentum.length; j++) {
                    momentum[j] = momentum[j] + functionalStepSize / 2 * HMCDerivative[j];  /* Sign change */
                }
            }
        }

        for (int i = 0; i < momentum.length; i++) {
                momentum[i] = momentum[i] + functionalStepSize / 2 * HMCDerivative[i];  /* Sign change */
        }

        double res=0;
        for (int i = 0; i <momentum.length ; i++) {
            res += momentum[i] * momentum[i] / (2 * Math.pow(drawDistribution.getSD(), 2));
        }
        return prop - res;
    }

    @Override
    public double getCoercableParameter() {
        return Math.log(stepSize);
    }

    @Override
    public void setCoercableParameter(double value) {
        stepSize = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return stepSize;
    }
}
