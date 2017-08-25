
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
    private  double tempStepsize = 0;


    private class Options { //TODO: these values might be adjusted for dual averaging.
        private double kappa = 0.75;
        private double t0 = 10.0;
        private double gamma = 0.05;
        private double delta = 0.3;
        private double deltaMax = 1000.0;

        private int findMax = 100;
        private int maxDepth = 100;
        private int adaptLength = 1000;
    }


    private final Options options = new Options();

    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Transform transform, double stepSize, int nSteps, double drawVariance) {
        super(mode, weight, gradientProvider, parameter, transform, stepSize, nSteps, drawVariance);


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
            this.mu = Math.log(10 * initialStepSize); // TODO Magic number? (for dual averaging use. leave for later)
        }

        public void update(long m, double alpha, double nAlpha, Options options) {
            
            if (m <= options.adaptLength) {
                h = (1 - 1 / (m + options.t0)) * h + 1 / (m + options.t0) * (options.delta - (alpha / nAlpha));
                logStepSize = mu - Math.sqrt(m) / options.gamma * h;
                averageLogStepSize = Math.pow(m, -options.kappa) * logStepSize +
                        (1 - Math.pow(m, -options.kappa)) * averageLogStepSize;
                stepSize = Math.exp(logStepSize);
                tempStepsize = stepSize;

            }

           //stepSize = 0.05; //todo: use to hand tune the stepsize.
            //stepSize = tempStepsize * (MathUtils.nextDouble()*0.2 + 0.9);//todo: randomize the stepsize each iteration.
            //System.err.println("m is " + m + " stepsize is " + stepSize);
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
        final double[] initialPOSITION = Arrays.copyOf(initialPosition,dim);

        if (stepSizeInformation == null) { // First call
            final double initialStepSize = findReasonableStepSize(initialPosition,options.findMax);
            stepSizeInformation = new StepSize(initialStepSize);


            final double testLogLikelihood = gradientProvider.getLikelihood().getLogLikelihood();
            assert (testLogLikelihood == initialLogLikelihood);
        }

        double[] position = takeOneStepNew(getCount() + 1, initialPOSITION);
        leapFrogEngine.setParameter(position);
        //System.err.println("position is" + Arrays.toString(position) + "jacobian " + leapFrogEngine.getParameterLogJacobian());

        return 0.0;
    }


    private double[] takeOneStepNew(long m, double[] initialPosition) {

        double[] endPosition = Arrays.copyOf(initialPosition, initialPosition.length);
        final double[] initialMomentum = drawInitialMomentum(drawDistribution, dim);

        final double initialJointDensity = getJointProbability(gradientProvider, initialMomentum);
        double logSliceU = Math.log(MathUtils.nextDouble()) + initialJointDensity;

        TreeState root = new TreeState(initialPosition, initialMomentum, 1, true);
        int j = 0;

        while (root.flagContinue) {


            double[] tmp = updateRoot(root, j, m,logSliceU, initialJointDensity);
            if (tmp != null) {
                endPosition = tmp;
            }

            j++;

            if (j > options.maxDepth) {
                throw new RuntimeException("Reach maximum tree depth"); // TODO Handle more gracefully
            }
        }
        stepSizeInformation.update(m, root.alpha, root.nAlpha, options);

        return endPosition;
    }

    private double[] updateRoot(TreeState root, int j, long m,double logSliceU, double initialJointDensity) { //todo no m here

        double[] endPosition = null;

        final double uniform1 = MathUtils.nextDouble();
        int direction = (uniform1 < 0.5) ? -1 : 1;

        TreeState node = buildTreeNew(m,root.getPosition(direction), root.getMomentum(direction), direction,
                            logSliceU, j, stepSizeInformation.stepSize, initialJointDensity);

        root.setPosition(direction, node.getPosition(direction));
        root.setMomentum(direction, node.getMomentum(direction));


        if (node.flagContinue) {
            double uniform2 = MathUtils.nextDouble();

            double p = (double) node.numNodes / (double)root.numNodes;

            if (uniform2 < p) {

                endPosition = node.getPosition(0);

            }
        }

        // Recursion
        root.numNodes += node.numNodes;
        root.flagContinue = computeStopCriterion(node.flagContinue, root);

        // Dual-averaging
        root.alpha += node.alpha;
        root.nAlpha += node.nAlpha;

        return endPosition;
    }

    private TreeState buildTreeNew(long m, double[] position, double[] momentum, int direction,
                                   double logSliceU, int j, double stepSize, double initialJointDensity) {


        if (j == 0) {

            return buildBaseCase( m, position, momentum, direction, logSliceU, j, stepSize, initialJointDensity);
        } else {

            return buildRecursiveCase(m, position, momentum, direction, logSliceU, j, stepSize, initialJointDensity);
        }
    }

    private TreeState buildBaseCase(long m,double[] inPosition, double[] inMomentum, int direction,
                                    double logSliceU, int j, double stepSize, double initialJointDensity) {


        // Make deep copy of position and momentum
        double[] position = Arrays.copyOf(inPosition, inPosition.length);
        double[] momentum = Arrays.copyOf(inMomentum, inMomentum.length);

        // "one frog jump!"
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);
        leapFrogEngine.updatePosition(position, momentum, direction * stepSize, 1.0);
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);

        double logJointProbAfter = getJointProbability(gradientProvider, momentum);

        final int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);

        final boolean flagContinue = (logSliceU < options.deltaMax + logJointProbAfter);

        // Values for dual-averaging
        final double alpha = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity));
        final int nAlpha = 1;

        leapFrogEngine.setParameter(inPosition);
        return new TreeState(position, momentum, numNodes, flagContinue, alpha, nAlpha);
    }

    private TreeState buildRecursiveCase(long m,double[] inPosition, double[] inMomentum, int direction,
                                         double logSliceU, int j, double stepSize, double initialJointDensity) {
        
        TreeState node = buildTreeNew(m,inPosition, inMomentum, direction, logSliceU,
                j - 1, // Recursion
                stepSize, initialJointDensity);

        if (node.flagContinue) {

            TreeState child = buildTreeNew(m,node.getPosition(direction), node.getMomentum(direction), direction,
                      logSliceU, j - 1, stepSizeInformation.stepSize, initialJointDensity);

            node.setPosition(direction, child.getPosition(direction));
            node.setMomentum(direction, child.getMomentum(direction));

            double uniform = MathUtils.nextDouble();
            if (child.numNodes > 0
                    && uniform <  ((double) child.numNodes / (double) (node.numNodes + child.numNodes))) {

                    node.setPosition(0, child.getPosition(0));
            }

            node.numNodes += child.numNodes;
            node.flagContinue = computeStopCriterion(child.flagContinue, node);

            node.alpha += child.alpha;
            node.nAlpha += child.nAlpha;

        }
        return node;
    }


    private double findReasonableStepSize(double[] position, int findmax) {

        double stepSize = 1;
        double[] momentum = drawInitialMomentum(drawDistribution, dim);
        int count = 1;

        double[] initialPosition = Arrays.copyOf(position,dim);

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
            probRatio = Math.exp(probAfter - probBefore);

            stepSize = Math.pow(2, a) * stepSize;
            count++;


            if (count > findmax) {
                throw new RuntimeException("Cannot find a reasonable stepsize in " + findmax + " iterations");
            }
        }

        leapFrogEngine.setParameter(initialPosition);


        return stepSize;
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

        boolean flag = flagContinue &&
                        getDotProduct(positionDifference, momentumMinus) >= 0 &&
                        getDotProduct(positionDifference, momentumPlus) >= 0;
        return flag;
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

    private double getJointProbability(GradientWrtParameterProvider gradientProvider, double[] momentum) {

        assert (gradientProvider != null);
        assert (momentum != null);
        return gradientProvider.getLikelihood().getLogLikelihood() - getScaledDotProduct(momentum, 1.0) - leapFrogEngine.getParameterLogJacobian();

    }

    private class TreeState {

        public TreeState(double[] position, double[] moment) {
            this(position, moment, 1, true);
        }

        public TreeState(double[] position, double[] moment,
                         int numNodes, boolean flagContinue) {
            this(position, moment, numNodes, flagContinue, 0.0, 0);
        }

        public TreeState(double[] position, double[] moment,
                         int numNodes, boolean flagContinue,
                         double alpha, int nAlpha) {
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
            this.alpha = alpha;
            this.nAlpha = nAlpha;
        }

        public double[] getPosition(int direction) {
            return position[getIndex(direction)];
        }

        public double[] getMomentum(int direction) {
            return momentum[getIndex(direction)];
        }

        public void setPosition(int direction, double[] position) {
            this.position[getIndex(direction)] = position;
        }

        public void setMomentum(int direction, double[] momentum) {
            this.momentum[getIndex(direction)] = momentum;
        }

        private int getIndex(int direction) {
            return direction + 1;
        }

        final private double[][] position;
        final private double[][] momentum;

        private int numNodes;
        private boolean flagContinue;

        private double alpha;
        private int nAlpha;
    }

}


