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

    private class Options { //TODO: these values might be adjusted for dual averaging.
        private double kappa = 0.75;
        private double t0 = 10.0;
        private double gamma = 0.05;
        private double targetAcceptRate = 0.9;
        private double logProbErrorTol = 100.0;
        private double muFactor = 10.0;

        private int findMax = 100;
        private int maxHeight = 10;
        private int adaptLength = 1000;
    }

    // TODO Magic numbers; pass as options
    private final Options options = new Options();

    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Transform transform, double stepSize, int nSteps, double drawVariance) {
        super(mode, weight, gradientProvider, parameter, transform, stepSize, nSteps, drawVariance, 0.0);
    }

    @Override
    protected InstabilityHandler getDefaultInstabilityHandler() {
        return InstabilityHandler.IGNORE;
    }

    private class StepSize {
        final double initialStepSize;
        double stepSize;
        double logStepSize;
        double averageLogStepSize;
        double h;
        double mu;

        private StepSize(double initialStepSize) {
            this.initialStepSize = initialStepSize;
            this.stepSize = initialStepSize;
            this.logStepSize = Math.log(stepSize);
            this.averageLogStepSize = 0;
            this.h = 0;
            this.mu = Math.log(options.muFactor * initialStepSize);
        }

        private void update(long m, double cumAcceptProb, double numAcceptProbStates, Options options) {

            if (m <= options.adaptLength) {

                h = (1 - 1 / (m + options.t0)) * h + 1 / (m + options.t0) * (options.targetAcceptRate - (cumAcceptProb / numAcceptProbStates));
                logStepSize = mu - Math.sqrt(m) / options.gamma * h;
                averageLogStepSize = Math.pow(m, -options.kappa) * logStepSize +
                        (1 - Math.pow(m, -options.kappa)) * averageLogStepSize;
                stepSize = Math.exp(logStepSize);
            }
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

        if (stepSizeInformation == null) {
            stepSizeInformation = findReasonableStepSize(initialPosition);

            final double testLogLikelihood = gradientProvider.getLikelihood().getLogLikelihood();
            assert (testLogLikelihood == initialLogLikelihood);
            assert (Arrays.equals(leapFrogEngine.getInitialPosition(), initialPosition));
        }

        double[] position = takeOneStep(getCount() + 1, initialPosition);
        leapFrogEngine.setParameter(position);

        return 0.0;
    }

    private double[] takeOneStep(long m, double[] initialPosition) {

        double[] endPosition = Arrays.copyOf(initialPosition, initialPosition.length);
        final double[] initialMomentum = drawInitialMomentum(drawDistribution, dim);

        final double initialJointDensity = getJointProbability(gradientProvider, initialMomentum);

        double logSliceU = Math.log(MathUtils.nextDouble()) + initialJointDensity;

        TreeState trajectoryTree = new TreeState(initialPosition, initialMomentum, 1, true);
            // Trajectory of Hamiltonian dynamics endowed with a binary tree structure.

        int height = 0;

        while (trajectoryTree.flagContinue) {

            double[] tmp = updateTrajectoryTree(trajectoryTree, height, logSliceU, initialJointDensity);
            if (tmp != null) {
                endPosition = tmp;
            }

            height++;

            if (height > options.maxHeight) {
                trajectoryTree.flagContinue = false;
                //throw new RuntimeException("Reach maximum tree height"); // TODO Handle more gracefully
            }
        }

        stepSizeInformation.update(m, trajectoryTree.cumAcceptProb, trajectoryTree.numAcceptProbStates, options);

        return endPosition;
    }

    private double[] updateTrajectoryTree(TreeState trajectoryTree, int depth, double logSliceU, double initialJointDensity) {

        double[] endPosition = null;

        final double uniform1 = MathUtils.nextDouble();
        int direction = (uniform1 < 0.5) ? -1 : 1;

        TreeState nextTrajectoryTree = buildTree(
                trajectoryTree.getPosition(direction), trajectoryTree.getMomentum(direction),
                direction, logSliceU, depth, stepSizeInformation.stepSize, initialJointDensity);

        if (nextTrajectoryTree.flagContinue) {

            final double uniform = MathUtils.nextDouble();
            final double acceptProb = (double) nextTrajectoryTree.numNodes / (double) trajectoryTree.numNodes;
            if (uniform < acceptProb) {
                endPosition = nextTrajectoryTree.getSample();
            }
        }

        trajectoryTree.mergeNextTree(nextTrajectoryTree, direction, true);

        return endPosition;
    }

    private TreeState buildTree(double[] position, double[] momentum, int direction,
                                double logSliceU, int height, double stepSize, double initialJointDensity) {

        if (height == 0) {
            return buildBaseCase(position, momentum, direction, logSliceU, stepSize, initialJointDensity);
        } else {
            return buildRecursiveCase(position, momentum, direction, logSliceU, height, stepSize, initialJointDensity);
        }
    }

    private void handleInstability() {
        throw new RuntimeException("Numerical instability; need to handle"); // TODO
    }

    private TreeState buildBaseCase(double[] inPosition, double[] inMomentum, int direction,
                                    double logSliceU, double stepSize, double initialJointDensity) {

        // Make deep copy of position and momentum
        double[] position = Arrays.copyOf(inPosition, inPosition.length);
        double[] momentum = Arrays.copyOf(inMomentum, inMomentum.length);

        leapFrogEngine.setParameter(position);

        // "one frog jump!"
        try {
            doLeap(position, momentum, direction * stepSize);
        } catch (NumericInstabilityException e) {
            handleInstability();
        }

        double logJointProbAfter = getJointProbability(gradientProvider, momentum);

        final int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

        final boolean flagContinue = (logSliceU < options.logProbErrorTol + logJointProbAfter);

        // Values for dual-averaging
        final double acceptProb = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
        final int numAcceptProbStates = 1;

        leapFrogEngine.setParameter(inPosition);

        return new TreeState(position, momentum, numNodes, flagContinue, acceptProb, numAcceptProbStates);
    }

    private TreeState buildRecursiveCase(double[] inPosition, double[] inMomentum, int direction,
                                         double logSliceU, int height, double stepSize, double initialJointDensity) {

        TreeState subtree = buildTree(inPosition, inMomentum, direction, logSliceU,
                height - 1, // Recursion
                stepSize, initialJointDensity);

        if (subtree.flagContinue) {

            TreeState nextSubtree = buildTree(subtree.getPosition(direction), subtree.getMomentum(direction), direction,
                    logSliceU, height - 1, stepSizeInformation.stepSize, initialJointDensity);

            subtree.mergeNextTree(nextSubtree, direction);

        }
        return subtree;
    }

    private void doLeap(final double[] position,
                        final double[] momentum,
                        final double stepSize) throws NumericInstabilityException {
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
        leapFrogEngine.updatePosition(position, momentum, stepSize, 1.0);
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
    }

    private StepSize findReasonableStepSize(double[] initialPosition) {

        double stepSize = 1;
        double[] momentum = drawInitialMomentum(drawDistribution, dim);
        int count = 1;

        double[] position = Arrays.copyOf(initialPosition, dim);

        double probBefore = getJointProbability(gradientProvider, momentum);

        try {
            doLeap(position, momentum, stepSize);
        } catch (NumericInstabilityException e) {
            handleInstability();
        }

        double probAfter = getJointProbability(gradientProvider, momentum);

        double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);

        double probRatio = Math.exp(probAfter - probBefore);

        while (Math.pow(probRatio, a) > Math.pow(2, -a)) {

            probBefore = probAfter;

            //"one frog jump!"
            try {
                doLeap(position, momentum, stepSize);
            } catch (NumericInstabilityException e) {
                handleInstability();
            }

            probAfter = getJointProbability(gradientProvider, momentum);
            probRatio = Math.exp(probAfter - probBefore);

            stepSize = Math.pow(2, a) * stepSize;
            count++;


            if (count > options.findMax) {
                throw new RuntimeException("Cannot find a reasonable step-size in " + options.findMax + " iterations");
            }
        }

        leapFrogEngine.setParameter(initialPosition);

        return new StepSize(stepSize);
    }

    private static boolean computeStopCriterion(boolean flagContinue, TreeState state) {
        return computeStopCriterion(flagContinue,
                state.getPosition(1), state.getPosition(-1),
                state.getMomentum(1), state.getMomentum(-1));
    }

    private static boolean computeStopCriterion(boolean flagContinue,
                                                double[] positionPlus, double[] positionMinus,
                                                double[] momentumPlus, double[] momentumMinus) {

        double[] positionDifference = subtractArray(positionPlus, positionMinus);

        return flagContinue &&
                getDotProduct(positionDifference, momentumMinus) >= 0 &&
                getDotProduct(positionDifference, momentumPlus) >= 0;
    }

    private static double getDotProduct(double[] x, double[] y) {

        assert (x.length == y.length);
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

    private  double getJointProbability(GradientWrtParameterProvider gradientProvider, double[] momentum) {

        assert (gradientProvider != null);
        assert (momentum != null);

        return gradientProvider.getLikelihood().getLogLikelihood() - getScaledDotProduct(momentum, 1.0)
                - leapFrogEngine.getParameterLogJacobian();
    }

    private class TreeState {

        private TreeState(double[] position, double[] moment,
                         int numNodes, boolean flagContinue) {
            this(position, moment, numNodes, flagContinue, 0.0, 0);
        }

        private TreeState(double[] position, double[] moment,
                         int numNodes, boolean flagContinue,
                         double cumAcceptProb, int numAcceptProbStates) {
            this.position = new double[3][];
            this.momentum = new double[3][];

            for (int i = 0; i < 3; ++i) {
                this.position[i] = position;
                this.momentum[i] = moment;
            }

            // Recursion variables
            this.numNodes = numNodes;
            this.flagContinue = flagContinue;

            // Dual-averaging variables
            this.cumAcceptProb = cumAcceptProb;
            this.numAcceptProbStates = numAcceptProbStates;
        }

        private double[] getPosition(int direction) {
            return position[getIndex(direction)];
        }

        private double[] getMomentum(int direction) {
            return momentum[getIndex(direction)];
        }

        private double[] getSample() {
            /*
            Returns a state chosen uniformly from the acceptable states along a hamiltonian dynamics trajectory tree.
            The sample is updated recursively while building trees.
            */
            return position[getIndex(0)];
        }

        private void setPosition(int direction, double[] position) {
            this.position[getIndex(direction)] = position;
        }

        private void setMomentum(int direction, double[] momentum) {
            this.momentum[getIndex(direction)] = momentum;
        }

        private void setSample(double[] position) { setPosition(0, position); }

        private int getIndex(int direction) { // valid directions: -1, 0, +1
            assert (direction >= -1 && direction <= 1);
            return direction + 1;
        }

        private void mergeNextTree(TreeState nextTree, int direction) {

            setPosition(direction, nextTree.getPosition(direction));
            setMomentum(direction, nextTree.getMomentum(direction));

            updateSample(nextTree);

            numNodes += nextTree.numNodes;
            flagContinue = computeStopCriterion(nextTree.flagContinue, this);

            cumAcceptProb += nextTree.cumAcceptProb;
            numAcceptProbStates += nextTree.numAcceptProbStates;
        }

        private void mergeNextTree(TreeState nextTree, int direction, boolean skipSample) {

            setPosition(direction, nextTree.getPosition(direction));
            setMomentum(direction, nextTree.getMomentum(direction));

            if (!skipSample) {
                updateSample(nextTree);
            }

            numNodes += nextTree.numNodes;
            flagContinue = computeStopCriterion(nextTree.flagContinue, this);

            cumAcceptProb += nextTree.cumAcceptProb;
            numAcceptProbStates += nextTree.numAcceptProbStates;
        }

        private void updateSample(TreeState nextTree) {
            double uniform = MathUtils.nextDouble();
            if (nextTree.numNodes > 0
                    && uniform < ((double) nextTree.numNodes / (double) (numNodes + nextTree.numNodes))) {
                setSample(nextTree.getSample());
            }
        }

        final private double[][] position;
        final private double[][] momentum;

        private int numNodes;
        private boolean flagContinue;

        private double cumAcceptProb;
        private int numAcceptProbStates;
    }
}


