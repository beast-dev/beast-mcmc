
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

//    private final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();

    private class Options {
        private double kappa = 0.75;
        private double t0 = 10.0;
        private double gamma = 0.05;
        private double delta = 0.5;
        private double deltaMax = 1000.0;

        private int maxDepth = 100;
        private int adaptLength = 10;
    }

    // TODO Magic numbers; pass as options
    private final Options options = new Options();

    Step nut = new Step();
    Step previousNut = new Step();

//    private double initialStepSize;

    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Transform transform, double stepSize, int nSteps, double drawVariance) {
        super(mode, weight, gradientProvider, parameter, transform, stepSize, nSteps, drawVariance);

        // Initialization
        System.err.println("RUNNING");

    }

    class StepSize {
        final double initialStepSize;
        double stepSize;
        double logStepSize;
        double averageLogStepSize;
        double h;
        double mu;

        StepSize(double initialStepSize) {
            this.initialStepSize = initialStepSize;
            this.stepSize = initialStepSize;
            this.logStepSize = Math.log(stepSize);
            this.averageLogStepSize = 0;
            this.h = 0;
            this.mu = Math.log(10 * initialStepSize); // TODO Magic number?
        }

        public void update(long m, double alpha, double nAlpha, Options options) {
            
            if (m <= options.adaptLength) {

                h = (1 - 1 / (m + options.t0)) * h + 1 / (m + options.t0) * (options.delta - (alpha / nAlpha));
                logStepSize = mu - Math.sqrt(m) / options.gamma * h;
                averageLogStepSize = Math.pow(m, -options.kappa) * logStepSize +
                        (1 - Math.pow(m, -options.kappa)) * averageLogStepSize;

                stepSize = Math.exp(logStepSize);

            } else {

                stepSize = Math.exp(averageLogStepSize); // TODO No necessary
            }

            stepSize = initialStepSize; //todo: fix dual averaging part (gave too small step size). Now just use the step size from findReasonableStepSize().

        }
    }

    private StepSize stepSizeInformation;

    @Override
    public String getOperatorName() {
        return "No-UTurn-Sampler operator";
    }

    @Override
    public double doOperation(Likelihood likelihood) {

        final double[] initialPosition = leapFrogEngine.getInitialPosition();
        final double initialLogLikelihood = gradientProvider.getLikelihood().getLogLikelihood();

        if (stepSizeInformation == null) { // First call
            final double initialStepSize = 0.0625;
//            final double initialStepSize = findReasonableStepSize(initialPosition);
            stepSizeInformation = new StepSize(initialStepSize);

            // TODO Debug
            leapFrogEngine.setParameter(initialPosition);
            final double testLogLikelihood = gradientProvider.getLikelihood().getLogLikelihood();
            assert (testLogLikelihood == initialLogLikelihood);
        }

        previousNut.position = initialPosition; // TODO Remove
        Step step = takeOneStep(getCount() + 1, previousNut); // TODO Remove previousNut

        leapFrogEngine.setParameter(step.position);  // zy: not sure if necessary

        return 0.0;
    }

    private Step takeOneStep(long m, Step previousNut) {

        assert (previousNut != null);

        Step returnNut = new Step(previousNut);

        final double[] momentum = drawInitialMomentum(drawDistribution, dim);

        final double logSliceU = Math.log(MathUtils.nextDouble()) + getJointProbability(gradientProvider, momentum);

        double[] positionMinus = Arrays.copyOf(previousNut.position, dim);
        double[] positionPlus = Arrays.copyOf(previousNut.position, dim);
        double[] momentumMinus = Arrays.copyOf(momentum, dim);
        double[] momentumPlus = Arrays.copyOf(momentum, dim);

        int j = 0;
        int n = 1;
        boolean growTree = true;  // Variable `s` in paper

        BinaryTree temp = new BinaryTree();
        while (growTree) {

            double random = MathUtils.nextDouble();

            if (random < 0.5) {

                temp = buildTree(positionMinus, momentumMinus, logSliceU, -1, j, stepSizeInformation.stepSize, momentum);

                positionMinus = Arrays.copyOf(temp.positionMinus, dim);
                momentumMinus = Arrays.copyOf(temp.momentumMinus, dim);

            } else {

                temp = buildTree(positionPlus, momentumPlus, logSliceU, 1, j, stepSizeInformation.stepSize, momentum);

                positionPlus = Arrays.copyOf(temp.positionPlus, dim);
                momentumPlus = Arrays.copyOf(temp.momentumPlus, dim);
            }

            if (temp.flagContinue) {

                if (MathUtils.nextDouble() < temp.numNodes / n) { //todo : (BUG) here the position may or may not be updated, but the gradientProvider is already updated.
                    returnNut.position = Arrays.copyOf(temp.positionFinal, dim);
                }

            }

            n += temp.numNodes;

            growTree = computeStopCriterion(temp.flagContinue, positionPlus, positionMinus, momentumPlus, momentumMinus);
            j++;

            if (j > options.maxDepth) {
                throw new RuntimeException("Reach maximum tree depth"); // TODO Handle more gracefully
            }
        }

        stepSizeInformation.update(m, temp.alpha, temp.nAlpha, options);

        return returnNut;
    }

    private BinaryTree buildTree(double[] position, double[] momentum, double logSliceU, int direction, int j, double
            stepSize, double[] initialMomentum) {

        if (j == 0) {

            BinaryTree tree = new BinaryTree();

            double logJointProbBefore = getJointProbability(gradientProvider, initialMomentum);

            // "one frog jump!"
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);
            leapFrogEngine.updatePosition(position, momentum, direction * stepSize, 1.0);
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);


            double logJointProbAfter = getJointProbability(gradientProvider, momentum);
            
            tree.momentumMinus = Arrays.copyOf(momentum, dim);
            tree.positionMinus = Arrays.copyOf(position, dim);
            tree.momentumPlus = Arrays.copyOf(momentum, dim);
            tree.positionPlus = Arrays.copyOf(position, dim);
            tree.positionFinal = Arrays.copyOf(position, dim);

            tree.alpha = Math.min(1, Math.exp(logJointProbAfter - logJointProbBefore));

            tree.nAlpha = 1;

//            tree.numNodes = (sliceU <= Math.exp(logJointProbAfter) ? 1 : 0);
            tree.numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

//            tree.flagContinue = (sliceU <= Math.exp(options.deltaMax + logJointProbAfter) ? true : false);
            tree.flagContinue = (logSliceU <= options.deltaMax + logJointProbAfter ? true : false);

            return tree;

        } else {

            BinaryTree newNodes;
            BinaryTree tree = buildTree(position, momentum, logSliceU, direction, j - 1, stepSize, initialMomentum);

            double[] positionMinus = Arrays.copyOf(tree.positionMinus, dim);
            double[] positionPlus = Arrays.copyOf(tree.positionPlus, dim);
            double[] momentumMinus = Arrays.copyOf(tree.momentumMinus, dim);
            double[] momentumPlus = Arrays.copyOf(tree.momentumPlus, dim);

            if (tree.flagContinue) {
                if (direction == -1) {

                    newNodes = buildTree(tree.positionMinus, tree.momentumMinus, logSliceU, direction, j - 1, stepSize, initialMomentum);

                    positionMinus = Arrays.copyOf(newNodes.positionMinus, dim);
                    momentumMinus = Arrays.copyOf(newNodes.momentumMinus, dim);
                } else {

                    newNodes = buildTree(tree.positionPlus, tree.momentumPlus, logSliceU, direction, j - 1, stepSize, initialMomentum);

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

    private double findReasonableStepSize(double[] position) { // TODO Still needs to be reviewed

        double stepSize = 1;
        final double[] momentum = drawInitialMomentum(drawDistribution, dim);
        int count = 1;

        double probBefore = getJointProbability(gradientProvider, momentum);

        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
        leapFrogEngine.updatePosition(position, momentum, stepSize, 1.0);
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);

        double probAfter = getJointProbability(gradientProvider, momentum);

        double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);

        double probRatio = Math.exp(probAfter - probBefore);

        while (Math.pow(probRatio, a) > Math.pow(2, -a)) {

            probBefore = probAfter;

            //"one frog jump!"
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
            leapFrogEngine.updatePosition(position, momentum, stepSize, 1.0); //zy: after "updateposition"
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);

            probAfter = getJointProbability(gradientProvider, momentum);
            probRatio = Math.exp(probAfter - probBefore);  // TODO Work on log-scale

            stepSize = Math.pow(2, a) * stepSize;
            count++;

            int maxTries = 100; // TODO magic number
            if (count > maxTries) {
                throw new RuntimeException("Cannot find a reasonable stepsize in " + maxTries + " iterations");
            }
        }

        return stepSize;
    }

    private static boolean computeStopCriterion(boolean flagContinue,
                                                double[] positionPlus, double[] positionMinus,
                                                double[] momentumPlus, double[] momentumMinus) {

        double[] positionDifference =  subtractArray(positionPlus, positionMinus);

        boolean flag = flagContinue &&
                        getDotProduct(positionDifference, momentumMinus) >= 0 &&
                        getDotProduct(positionDifference, momentumPlus) >= 0; // TODO check signs
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

    private static double[] subtractArray(double[] a, double[] b) {

        assert (a.length == b.length);
        final int dim = a.length;

        double result[] = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    private static double[] sumArray(double[] a, double[] b) {

        assert (a.length == b.length);
        final int dim = a.length;

        double result[] = new double[dim];

        for (int i = 0; i < dim; i++) {
            result[i] = a[i] + b[i];
        }

        return result;
    }

    private static double getJointProbability(GradientWrtParameterProvider gradientProvider, double[] momentum) {

        assert (gradientProvider != null);
        assert (momentum != null);

        final double logJointProb = gradientProvider.getLikelihood().getLogLikelihood()
                - getScaledDotProduct(momentum, 1.0);
        return logJointProb;
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

    public class Step {
        double[] position;
//        double stepSize;
//        double logStepSize;
//        double logStepSizeAve;
//        double h;
//        double mu;

        public Step() {
        }

        public Step(double[] inPosition) { // TODO Remove
            this.position = Arrays.copyOf(inPosition, inPosition.length);
        }

        public Step(Step copy) {

            this.position = java.util.Arrays.copyOf(copy.position, dim);
//            this.stepSize = copy.stepSize;
//            this.logStepSize = copy.logStepSize;
//            this.logStepSizeAve = copy.logStepSizeAve;
//            this.h = copy.h;
//            this.mu = copy.mu;
        }

    }

}


