
/*
 * NoUTurnOperator.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.GeneralOperator;
import dr.inference.operators.GibbsOperator;
import dr.math.MathUtils;
import dr.util.Transform;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Zhenyu Zhang
 */

public class NoUTurnOperator extends HamiltonianMonteCarloOperator implements GeneralOperator, GibbsOperator {

    private final int dim = gradientProvider.getDimension();
    private final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();
    final double kappa = 0.75;
    final double t0 = 10;
    final double gamma = 0.05;
    final int adaptLength = 10;
    final double delta = 0.5;
    double deltaMax = 1000;
    int maxDepth = 100; // TODO, the last todo, make an option for these.
    OneNut nut = new OneNut();
    OneNut previousNut = new OneNut();
    private double initialStepSize;

    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Transform transform, double stepSize, int nSteps, double drawVariance) {
        super(mode, weight, gradientProvider, parameter, transform, stepSize, nSteps, drawVariance);
    }

    @Override
    public String getOperatorName() {
        return "NUTS operator";
    }

    @Override
    public double doOperation(Likelihood likelihood) {

        long m = getCount() + 1;

        if (m == 1) {

            previousNut.position = leapFrogEngine.getInitialPosition();
            initialStepSize = 0.0625;
//            final double[] initialPosition = Arrays.copyOf(leapFrogEngine.getInitialPosition(), dim);
//            initialStepSize = findReasonableStepSize(initialPosition);
//todo : in findReasonableStepSize() the leapfrog.update will also update gradientProvider.. Maybe in the function we can replace it with a deep copy?
            previousNut.stepSize = initialStepSize;
            previousNut.logStepSize = Math.log(previousNut.stepSize);
            previousNut.logStepSizeAve = 0;
            previousNut.h = 0;
            previousNut.mu = Math.log(10 * initialStepSize);
        }

        nut = nutsOneStep(m, previousNut);

        previousNut = new OneNut(nut);

        leapFrogEngine.setParameter(nut.position); // zy: not sure if necessary

        return 233;
    }

    private OneNut nutsOneStep(long m, OneNut previousNut) {

        if (previousNut == null) {
            System.exit(3);
        }

        OneNut returnNut = new OneNut(previousNut);

        double[] momentum = drawInitialMomentum(drawDistribution, dim);
        double sliceU = MathUtils.nextDouble() * Math.exp(getJointProbability(gradientProvider, momentum));
        double[] positionMinus = Arrays.copyOf(previousNut.position, dim);
        double[] positionPlus = Arrays.copyOf(previousNut.position, dim);
        double[] momentumMinus = Arrays.copyOf(momentum, dim);
        double[] momentumPlus = Arrays.copyOf(momentum, dim);

        int j = 0;
        int n = 1;
        boolean s = true;

        BinaryTree temp = new BinaryTree();

        while (s) {

            double random = MathUtils.nextDouble();

            if (random < 0.5) {

                temp = buildTree(positionMinus, momentumMinus, sliceU, -1, j, previousNut.stepSize, momentum);

                positionMinus = Arrays.copyOf(temp.positionMinus, dim);
                momentumMinus = Arrays.copyOf(temp.momentumMinus, dim);

            } else {

                temp = buildTree(positionPlus, momentumPlus, sliceU, 1, j, previousNut.stepSize, momentum);

                positionPlus = Arrays.copyOf(temp.positionPlus, dim);
                momentumPlus = Arrays.copyOf(temp.momentumPlus, dim);
            }

            if (temp.flagContinue) {

                if (MathUtils.nextDouble() < temp.numNodes / n) { //todo : (BUG) here the position may or may not be updated, but the gradientProvider is already updated.
                    returnNut.position = Arrays.copyOf(temp.positionFinal, dim);
                }

            }

            n += temp.numNodes;
            s = computeStopCriterion(temp.flagContinue, positionPlus, positionMinus, momentumPlus, momentumMinus);
            j++;
            if (j > maxDepth) {
                System.err.println("reach maximum tree depth");
                System.exit(-1);
            }
        }
        if (m <= adaptLength) { 

            returnNut.h = (1 - 1 / (m + t0)) * previousNut.h + 1 / (m + t0) * (delta - (temp.alpha / temp.nAlpha));
            returnNut.logStepSize = previousNut.mu - Math.sqrt(m) / gamma * returnNut.h;
            returnNut.logStepSizeAve = Math.pow(m, -kappa) * returnNut.logStepSize + (1 - Math.pow(m, -kappa)) * previousNut.logStepSizeAve;

            returnNut.stepSize = Math.exp(returnNut.logStepSize);

        } else {
            returnNut.stepSize = Math.exp(returnNut.logStepSizeAve);
        }
        returnNut.stepSize = initialStepSize; //todo: fix dual averaging part (gave too small step size). Now just use the step size from findReasonableStepSize().
        return returnNut;
    }

    private BinaryTree buildTree(double[] position, double[] momentum, double sliceU, int direction, int j, double
            stepSize, double[] initialMomentum) {

        if (j == 0) {

            BinaryTree tree = new BinaryTree();

            double logJointProbBefore = getJointProbability(gradientProvider, initialMomentum);

            // "one frog jump!"
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);
            leapFrogEngine.updatePosition(position, momentum, direction * stepSize, sigmaSquared);
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);


            double logJointProbAfter = getJointProbability(gradientProvider, momentum);
            
            tree.momentumMinus = Arrays.copyOf(momentum, dim);
            tree.positionMinus = Arrays.copyOf(position, dim);
            tree.momentumPlus = Arrays.copyOf(momentum, dim);
            tree.positionPlus = Arrays.copyOf(position, dim);
            tree.positionFinal = Arrays.copyOf(position, dim);

            tree.alpha = Math.min(1, Math.exp(logJointProbAfter - logJointProbBefore));

            tree.nAlpha = 1;

            tree.numNodes = (sliceU <= Math.exp(logJointProbAfter) ? 1 : 0);

            tree.flagContinue = (sliceU <= Math.exp(deltaMax + logJointProbAfter) ? true : false);

            return tree;

        } else {

            BinaryTree newNodes;
            BinaryTree tree = buildTree(position, momentum, sliceU, direction, j - 1, stepSize, initialMomentum);

            double[] positionMinus = Arrays.copyOf(tree.positionMinus, dim);
            double[] positionPlus = Arrays.copyOf(tree.positionPlus, dim);
            double[] momentumMinus = Arrays.copyOf(tree.momentumMinus, dim);
            double[] momentumPlus = Arrays.copyOf(tree.momentumPlus, dim);

            if (tree.flagContinue) {
                if (direction == -1) {

                    newNodes = buildTree(tree.positionMinus, tree.momentumMinus, sliceU, direction, j - 1, stepSize, initialMomentum);

                    positionMinus = Arrays.copyOf(newNodes.positionMinus, dim);
                    momentumMinus = Arrays.copyOf(newNodes.momentumMinus, dim);
                } else {

                    newNodes = buildTree(tree.positionPlus, tree.momentumPlus, sliceU, direction, j - 1, stepSize, initialMomentum);

                    positionPlus = Arrays.copyOf(newNodes.positionPlus, dim);
                    momentumPlus = Arrays.copyOf(newNodes.momentumPlus, dim);
                }
                double randomNum = MathUtils.nextDouble();

                if (randomNum < (double) newNodes.numNodes / (tree.numNodes + newNodes.numNodes)) {
                    tree.positionFinal = Arrays.copyOf(newNodes.positionFinal, dim);
                }

                tree.flagContinue = computeStopCriterion(newNodes.flagContinue, positionPlus, positionMinus, momentumPlus, momentumMinus);
                tree.positionPlus = Arrays.copyOf(positionPlus, dim);
                tree.positionMinus = Arrays.copyOf(positionMinus, dim);
                tree.momentumMinus = Arrays.copyOf(momentumMinus, dim);
                tree.momentumPlus = Arrays.copyOf(momentumPlus, dim);

                //values for dual averaging use..
                tree.alpha += newNodes.alpha;
                tree.nAlpha += newNodes.nAlpha;
                tree.numNodes += newNodes.numNodes;
            }
            return tree;
        }
    }

    private double findReasonableStepSize(double[] position) {

        double stepSize = 1;
        final double[] momentum = drawInitialMomentum(drawDistribution, dim);
        int count = 1;

        double probBefore = getJointProbability(gradientProvider, momentum);

        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
        leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared);
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);

        double probAfter = getJointProbability(gradientProvider, momentum);

        double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);

        double probRatio = Math.exp(probAfter - probBefore);

        while (Math.pow(probRatio, a) > Math.pow(2, -a)) {

            probBefore = probAfter;

            //"one frog jump!"
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
            leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared); //zy: after "updateposition"
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);

            probAfter = getJointProbability(gradientProvider, momentum);

            probRatio = Math.exp(probAfter - probBefore);

            stepSize = Math.pow(2, a) * stepSize;

            count++;

            if (count > 100) {
                System.out.println("cannot find a reasonable stepsize in 100 iterations!");
                System.exit(-100);
                break;
            }
        }

        return stepSize;
    }

    private boolean computeStopCriterion(boolean flagContinue, double[] positionPlus, double[] positionMinus, double[] momentumPlus, double[] momentumMinus) {
        boolean flag = flagContinue && getDotProduct(sumArray(positionPlus, positionMinus, false),
                momentumMinus) >= 0 &&
                getDotProduct(sumArray(positionPlus, positionMinus, false), momentumPlus) >= 0;
        return flag;
    }

    private static double getDotProduct(double[] x, double[] y) {
        final int dim = x.length;
        double total = 0.0;
        for (int i = 0; i < dim; i++) {
            total += x[i] * y[i];
        }
        return total;
    }

    private double[] sumArray(double[] a, double[] b, boolean sum) {
        double result[] = new double[dim];
        if (sum) {
            for (int i = 0; i < dim; i++) {
                result[i] = a[i] + b[i];
            }
        } else {
            for (int i = 0; i < dim; i++) {
                result[i] = a[i] - b[i];
            }
        }
        return result;
    }

    private double getJointProbability(GradientWrtParameterProvider gradientProvider, double[] momentum) {

        double logjointprob;

        if (gradientProvider != null && momentum != null) {
            logjointprob = gradientProvider.getLikelihood().getLogLikelihood() - getScaledDotProduct(momentum, sigmaSquared);
            return logjointprob;

        } else if (gradientProvider == null) {
            System.err.println("null in get joint prob (likelihood)!!");
            System.exit(-99);
            return -99;
        } else {
            System.err.println("null in get joint prob (momentum)!!");
            System.exit(-99);
            return -99;
        }

    }

    @Override
    public int getStepCount() {
        return (int) getCount();
    }

    public class BinaryTree {
        double[] positionMinus;
        double[] momentumMinus;
        double[] positionPlus;
        double[] momentumPlus;
        double[] positionFinal;
        int numNodes;
        boolean flagContinue;
        double alpha;
        int nAlpha;

        public BinaryTree() {
        }
    }

    public class OneNut {
        double[] position;
        double stepSize;
        double logStepSize;
        double logStepSizeAve;
        double h;
        double mu;

        public OneNut() {
        }

        public OneNut(OneNut copy) {

            this.position = java.util.Arrays.copyOf(copy.position, dim);
            this.stepSize = copy.stepSize;
            this.logStepSize = copy.logStepSize;
            this.logStepSizeAve = copy.logStepSizeAve;
            this.h = copy.h;
            this.mu = copy.mu;
            
            
        }

    }

}


