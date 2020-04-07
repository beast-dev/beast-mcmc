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
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.GeneralOperator;
import dr.inference.operators.GibbsOperator;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.WrappedVector;
import dr.util.Transform;

import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Zhenyu Zhang
 */

public class NoUTurnOperator extends HamiltonianMonteCarloOperator implements GeneralOperator, GibbsOperator { //todo: NUTS should not be subclass of HMC

    private final int dim = gradientProvider.getDimension();

    class Options {
        private double logProbErrorTol = 100.0;
        private int findMax = 100;
        private int maxHeight = 10;
    }

    private final Options options = new Options();
    
    public NoUTurnOperator(AdaptationMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Transform transform, Parameter mask,
                           HamiltonianMonteCarloOperator.Options runtimeOptions,
                           MassPreconditioner.Type preconditioningType,
                           ReversibleHMCProvider reversibleHMCprovider) {
        super(mode, weight, gradientProvider, parameter, transform, mask, runtimeOptions, preconditioningType);
        this.reversibleHMCProvider = reversibleHMCprovider;
    }

    @Override
    protected InstabilityHandler getDefaultInstabilityHandler() {
        return InstabilityHandler.IGNORE;
    }

    private StepSize stepSizeInformation;

    @Override
    public String getOperatorName() {
        return "No-UTurn-Sampler operator";
    }

    @Override
    public double doOperation(Likelihood likelihood) {

        if (shouldCheckGradient()) {
            checkGradient(likelihood);
        }

        final double[] initialPosition = leapFrogEngine.getInitialPosition();

        if (stepSizeInformation == null) {
            stepSizeInformation = findReasonableStepSize(initialPosition, super.stepSize);
        }

        double[] position = takeOneStep(getCount() + 1, initialPosition);
        leapFrogEngine.setParameter(position);

        return 0.0;
    }

    private double[] takeOneStep(long m, double[] initialPosition) {

        double[] endPosition = Arrays.copyOf(initialPosition, initialPosition.length);
//        final double[][] mass = massProvider.getMass();
        //final WrappedVector initialMomentum = mask(preconditioning.drawInitialMomentum(), mask);
        final WrappedVector initialMomentum = mask(reversibleHMCProvider.drawMomentum(), mask);

        final double initialJointDensity = getJointProbability(gradientProvider, initialMomentum);

        double logSliceU = Math.log(MathUtils.nextDouble()) + initialJointDensity;

        TreeState trajectoryTree = new TreeState(initialPosition, initialMomentum.getBuffer(), 1, true);
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
            }
        }

        stepSizeInformation.update(m, trajectoryTree.cumAcceptProb, trajectoryTree.numAcceptProbStates);

        return endPosition;
    }

    private double[] updateTrajectoryTree(TreeState trajectoryTree, int depth, double logSliceU, double initialJointDensity) {

        double[] endPosition = null;

        final double uniform1 = MathUtils.nextDouble();
        int direction = (uniform1 < 0.5) ? -1 : 1;

        TreeState nextTrajectoryTree = buildTree(
                trajectoryTree.getPosition(direction), trajectoryTree.getMomentum(direction),
                direction, logSliceU, depth, stepSizeInformation.getStepSize(), initialJointDensity);

        if (nextTrajectoryTree.flagContinue) {

            final double uniform = MathUtils.nextDouble();
            final double acceptProb = (double) nextTrajectoryTree.numNodes / (double) trajectoryTree.numNodes;
            if (uniform < acceptProb) {
                endPosition = nextTrajectoryTree.getSample();
            }
        }

        trajectoryTree.mergeNextTree(nextTrajectoryTree, direction);

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


    private TreeState buildBaseCase(double[] inPosition, double[] inMomentum, int direction,
                                    double logSliceU, double stepSize, double initialJointDensity) {

        // Make deep copy of position and momentum
        WrappedVector position = new WrappedVector.Raw(Arrays.copyOf(inPosition, inPosition.length));
        //double[] position = Arrays.copyOf(inPosition, inPosition.length);
        WrappedVector momentum = new WrappedVector.Raw(Arrays.copyOf(inMomentum, inMomentum.length));

        leapFrogEngine.setParameter(position.getBuffer());

        // "one reversibleHMC integral
        reversibleHMCProvider.reversiblePositionUpdate(position, momentum, direction, stepSize);
//        try {
//            doLeap(position, momentum, direction * stepSize);
//        } catch (NumericInstabilityException e) {
//            handleInstability();
//        }

        double logJointProbAfter = getJointProbability(gradientProvider, momentum);

        final int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

        final boolean flagContinue = (logSliceU < options.logProbErrorTol + logJointProbAfter);

        // Values for dual-averaging
        final double acceptProb = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
        final int numAcceptProbStates = 1;

        leapFrogEngine.setParameter(inPosition);

        return new TreeState(position.getBuffer(), momentum.getBuffer(), numNodes, flagContinue, acceptProb, numAcceptProbStates);
    }

    private TreeState buildRecursiveCase(double[] inPosition, double[] inMomentum, int direction,
                                         double logSliceU, int height, double stepSize, double initialJointDensity) {

        TreeState subtree = buildTree(inPosition, inMomentum, direction, logSliceU,
                height - 1, // Recursion
                stepSize, initialJointDensity);

        if (subtree.flagContinue) {

            TreeState nextSubtree = buildTree(subtree.getPosition(direction), subtree.getMomentum(direction), direction,
                    logSliceU, height - 1, stepSizeInformation.getStepSize(), initialJointDensity);

            subtree.mergeNextTree(nextSubtree, direction);

        }
        return subtree;
    }


    private StepSize findReasonableStepSize(double[] initialPosition, double forcedInitialStepSize) {

        if (forcedInitialStepSize != 0) {
            return new StepSize(forcedInitialStepSize);
        } else {
            double stepSize = 0.1;
//        final double[] mass = massProvider.getMass();
            WrappedVector momentum = preconditioning.drawInitialMomentum();
            int count = 1;

            WrappedVector position = new WrappedVector.Raw(Arrays.copyOf(initialPosition, dim));

            double probBefore = getJointProbability(gradientProvider, momentum);
            reversibleHMCProvider.reversiblePositionUpdate(position, momentum, 1, stepSize);
//            try {
//                doLeap(position, momentum, stepSize);
//            } catch (NumericInstabilityException e) {
//                handleInstability();
//            }

            double probAfter = getJointProbability(gradientProvider, momentum);

            double a = ((probAfter - probBefore) > Math.log(0.5) ? 1 : -1);

            double probRatio = Math.exp(probAfter - probBefore);

            while (Math.pow(probRatio, a) > Math.pow(2, -a)) {

                probBefore = probAfter;
                reversibleHMCProvider.reversiblePositionUpdate(position, momentum, 1, stepSize);
                //"one frog jump!"
//                try {
//                    doLeap(position, momentum, stepSize);
//                } catch (NumericInstabilityException e) {
//                    handleInstability();
//                }

                probAfter = getJointProbability(gradientProvider, momentum);
                probRatio = Math.exp(probAfter - probBefore);

                stepSize = Math.pow(2, a) * stepSize;
                count++;

                if (count > options.findMax) {
                    throw new RuntimeException("Cannot find a reasonable step-size in " + options.findMax + " " +
                            "iterations");
                }
            }
            leapFrogEngine.setParameter(initialPosition);
            return new StepSize(stepSize);
        }
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

        double[] result = new double[dim];
        for (int i = 0; i < dim; i++) {
            result[i] = a[i] - b[i];
        }

        return result;
    }

    private  double getJointProbability(GradientWrtParameterProvider gradientProvider, WrappedVector momentum) {

        assert (gradientProvider != null);
        assert (momentum != null);

        return gradientProvider.getLikelihood().getLogLikelihood() - getKineticEnergy(momentum)
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

    private ReversibleHMCProvider reversibleHMCProvider;
}


