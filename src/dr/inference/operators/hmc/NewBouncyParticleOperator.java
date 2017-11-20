/*
 * NewHamiltonianMonteCarloOperator.java
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

import dr.evolution.tree.Tree;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;

/**
 * @author Zhenyu Zhang
 * @author Marc A. Suchard
 */

public class NewBouncyParticleOperator extends SimpleMCMCOperator {

    final GradientWrtParameterProvider gradientProvider;
    final MatrixParameterInterface precisionMatrix;
    // Something that returns the conditional distribution
    final double[][] sigma0;

//    protected double stepSize;
//    protected final int nSteps;
//    final NormalDistribution drawDistribution;
//    final LeapFrogEngine leapFrogEngine;


    public NewBouncyParticleOperator(CoercionMode mode, double weight,
                                     GradientWrtParameterProvider gradientProvider,
                                     MatrixParameterInterface precisionMatrix,
                                     double[][] sigma0) {
//        super(mode);

        this.gradientProvider = gradientProvider;
        this.precisionMatrix = precisionMatrix;
        this.sigma0 = sigma0;

        setWeight(weight);

    }

//    public NewBouncyParticleOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
//                                     Parameter parameter, Transform transform,
//                                     double stepSize, int nSteps, double drawVariance) {
////        super(mode);
////        setWeight(weight);
////        setTargetAcceptanceProbability(0.8); // Stan default
////
////        this.gradientProvider = gradientProvider;
////        this.stepSize = stepSize;
////        this.nSteps = nSteps;
////        this.drawDistribution = new NormalDistribution(0, Math.sqrt(drawVariance));
////        this.leapFrogEngine = (transform != null ?
////                new LeapFrogEngine.WithTransform(parameter, transform) :
////                new LeapFrogEngine.Default(parameter));
//
//    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Bouncy Particle operator";
    }

//    static double getScaledDotProduct(final double[] momentum,
//                                      final double sigmaSquared) {
//        double total = 0.0;
//        for (double m : momentum) {
//            total += m * m;
//        }
//
//        return total / (2 * sigmaSquared);
//    }

//    static double[] drawInitialMomentum(final NormalDistribution distribution, final int dim) {
//        double[] momentum = new double[dim];
//        for (int i = 0; i < dim; i++) {
//            momentum[i] = (Double) distribution.nextRandom();
//        }
//        return momentum;
//    }

    @Override
    public double doOperation() {
//        return leapFrog();
        return 0.0;
    }

    private long count = 0;

    private static final boolean DEBUG = false;

//    private double leapFrog() {

//        if (DEBUG) {
//            if (count % 5 == 0) {
//                System.err.println(stepSize);
//            }
//            ++count;
//        }
//
//        final int dim = gradientProvider.getDimension();
//
//        final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();
//
//        final double[] momentum = drawInitialMomentum(drawDistribution, dim);
//        final double[] position = leapFrogEngine.getInitialPosition();
//
//        final double prop = getScaledDotProduct(momentum, sigmaSquared) +
//                leapFrogEngine.getParameterLogJacobian();
//
//        leapFrogEngine.updateMomentum(position, momentum,
//                gradientProvider.getGradientLogDensity(), stepSize / 2);
//
//
//        if (DEBUG) {
//            System.err.println("nSteps = " + nSteps);
//        }
//
//        for (int i = 0; i < nSteps; i++) { // Leap-frog
//
//            leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared);
//
//            if (i < (nSteps - 1)) {
//                leapFrogEngine.updateMomentum(position, momentum,
//                        gradientProvider.getGradientLogDensity(), stepSize);
//            }
//        }
//
//        leapFrogEngine.updateMomentum(position, momentum,
//                gradientProvider.getGradientLogDensity(), stepSize / 2);
//
//        final double res = getScaledDotProduct(momentum, sigmaSquared) +
//                leapFrogEngine.getParameterLogJacobian();
//
//        return prop - res; //hasting ratio
//    }

//    @Override
//    public double getCoercableParameter() {
//        return Math.log(stepSize);
//    }
//
//    @Override
//    public void setCoercableParameter(double value) {
//        if (DEBUG) {
//            System.err.println("Setting coercable paramter: " + getCoercableParameter() + " -> " + value);
//        }
//        stepSize = Math.exp(value);
//    }
//
//    @Override
//    public double getRawParameter() {
//        return 0;
//    }
}
