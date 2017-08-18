
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

import java.util.Scanner; //todo delete these two imports.
import java.io.FileNotFoundException;
import java.io.*;
/**
 * @author Marc A. Suchard
 * @author Zhenyu Zhang
 */

public class NoUTurnOperator extends HamiltonianMonteCarloOperator implements GeneralOperator, GibbsOperator { //todo gibbs?

    private final int dim = gradientProvider.getDimension();
    double[] logus = readFile1("logu.txt");
    int DEBUG = 0;

    double[] dir1 = readFile1("dir1.txt");

    double[] dir2 = readFile1("dir2.txt");

//    private final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();

    private class Options {
        private double kappa = 0.75;
        private double t0 = 10.0;
        private double gamma = 0.05;
        private double delta = 0.8;
        private double deltaMax = 1000.0;

        private int maxDepth = 10;
        private int adaptLength = 100;
    }

    // TODO Magic numbers; pass as options
    private final Options options = new Options();

//    Step nut = new Step();
//    Step previousNut = new Step();

//    private double initialStepSize;

    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, Transform transform, double stepSize, int nSteps, double drawVariance) {
        super(mode, weight, gradientProvider, parameter, transform, stepSize, nSteps, drawVariance);

        // Initialization
        System.err.println("RUNNING");

    }



//    private double[] readFile(String fileName) {
//        int i = 0;
//        double[] logu = new double[101];
//        try {
//            File file = new File(fileName);
//            Scanner scanner = new Scanner(file);
//            while (scanner.hasNextLine()) {
//                logu[i] = Double.parseDouble(scanner.nextLine());
////                System.out.println(scanner.nextLine());
//                i++;
//            }
//            scanner.close();
//        } catch (FileNotFoundException e) {
//            e.printStackTrace();
//        }
//        return logu;
//    }

    private double[] readFile1(String fileName) {
        double[] logu = new double[1020];
        try {
            Scanner scanner = new Scanner(new File(fileName));
            int i = 0;
            while (scanner.hasNextDouble()) {
                logu[i++] = scanner.nextDouble();
            }scanner.close();
        }catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return logu;
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
            //System.err.println("m is " + m + "stepsize is " + stepSize);

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
//            System.err.println("logu length is " + length);
            
            //System.err.println("initial prob is " + initialLogLikelihood);
//            final double initialStepSize = findReasonableStepSize(initialPosition);
            final double initialStepSize = 1;// TODO Magic number

            //System.err.println("initial stepsize is " + initialStepSize);

            stepSizeInformation = new StepSize(initialStepSize);


            //leapFrogEngine.setParameter(initialPosition); // Restart at initial possible todo: not working ...


            // TODO For debug purposes
            final double testLogLikelihood = gradientProvider.getLikelihood().getLogLikelihood();
          //  System.err.println("test prob is " + testLogLikelihood);
            assert (testLogLikelihood == initialLogLikelihood);



        }


//        double[] position = takeOneStep(getCount() + 1, initialPosition);

        double[] position = takeOneStepNew(getCount() + 1, initialPosition);

       // leapFrogEngine.setParameter(position);
        final double[] tmp_momentum = drawInitialMomentum(drawDistribution, dim);
        leapFrogEngine.updatePosition(position, tmp_momentum, 0.0, 1.0);

        return 0.0;
    }

//    private double[] takeOneStep(long m, double[] initialPosition) {
//
//        double[] endPosition = Arrays.copyOf(initialPosition, initialPosition.length);
//        final double[] initialMomentum = drawInitialMomentum(drawDistribution, dim);
//
//        final double logSliceU = Math.log(MathUtils.nextDouble()) + getJointProbability(gradientProvider, initialMomentum);
//
//        double[] positionMinus = Arrays.copyOf(initialPosition, dim); // TODO Why are copies necessary?
//        double[] positionPlus = Arrays.copyOf(initialPosition, dim);
//        double[] momentumMinus = Arrays.copyOf(initialMomentum, dim);
//        double[] momentumPlus = Arrays.copyOf(initialMomentum, dim);
//
//        int j = 0;
//
//        int n = 1;
//
//        boolean growTree = true;  // Variable `s` in paper
//
//        TreeNode child = new TreeNode(); // TODO Do no initialize here
//        while (growTree) {
//
//            double uniform = MathUtils.nextDouble();
//            int direction = (uniform < 0.5) ? -1 : 1;
//
////            child = buildTree(state.getPosition(direction), state.getMomentum(direction),
////                    logSliceU, direction, j, stepSizeInformation.stepSize, initialMomentum);
//
//            if (uniform < 0.5) {
//
//                child = buildTree(positionMinus, momentumMinus, logSliceU, -1, j, stepSizeInformation.stepSize, initialMomentum);
//
//                positionMinus = Arrays.copyOf(child.positionMinus, dim);
//                momentumMinus = Arrays.copyOf(child.momentumMinus, dim);
//
//            } else {
//
//                child = buildTree(positionPlus, momentumPlus, logSliceU, 1, j, stepSizeInformation.stepSize, initialMomentum);
//
//                positionPlus = Arrays.copyOf(child.positionPlus, dim);
//                momentumPlus = Arrays.copyOf(child.momentumPlus, dim);
//            }
//
//            if (child.flagContinue) {
//
//                if (MathUtils.nextDouble() < (double) (child.numNodes / n)) { //todo : (BUG) here the position may or may not be updated, but the gradientProvider is already updated.
//                    endPosition = Arrays.copyOf(child.positionFinal, dim);
//                }
//
//            }
//
//            n += child.numNodes;
//
//            growTree = computeStopCriterion(child.flagContinue, positionPlus, positionMinus, momentumPlus, momentumMinus);
//            j++;
////            state = child;
//
//            if (j > options.maxDepth) {
//                throw new RuntimeException("Reach maximum tree depth"); // TODO Handle more gracefully
//            }
//        }
//
//        stepSizeInformation.update(m, child.alpha, child.nAlpha, options);
//
//        return endPosition;
//    }


    private double[] takeOneStepNew(long m, double[] initialPosition) {

//        double[] endPosition = null;
        //if (m == 53) {DEBUG = 1;}

        double[] endPosition = Arrays.copyOf(initialPosition, initialPosition.length); // todo: check if it should be null
        final double[] initialMomentum = drawInitialMomentum(drawDistribution, dim);
        initialMomentum[0] = 0.3;

        // TODO Could save a likelihood evaluation by computing likelihood before guessing step size
        final double initialJointDensity = getJointProbability(gradientProvider, initialMomentum);

        double logSliceU = Math.log(MathUtils.nextDouble()) + initialJointDensity;
        logSliceU = initialJointDensity - 1;
        int index = ((int) m )- 1;
        logSliceU = logus[index];
//        System.err.println("m is " + m + "logu is" + logus[index]);

        TreeState root = new TreeState(initialPosition, initialMomentum, 1, true);

        int j = 0;

        while (root.flagContinue) {


            double[] tmp = updateRoot(root, j, m,logSliceU, initialJointDensity);
            if (tmp != null) {
                endPosition = tmp;
            }

            j++;
//            assert (j<=2);
            if (j > 2) {
                System.err.println("j is " + j);
            }
            if (j > options.maxDepth) {
                throw new RuntimeException("Reach maximum tree depth"); // TODO Handle more gracefully
            }
        }

        stepSizeInformation.update(m, root.alpha, root.nAlpha, options);
        stepSizeInformation.stepSize = 1; // todo

        return endPosition;
    }

    private double[] updateRoot(TreeState root, int j, long m,double logSliceU, double initialJointDensity) { //todo no m here

        double[] endPosition = null;

        final double uniform1 = MathUtils.nextDouble();
        int direction = (uniform1 < 0.5) ? -1 : 1;

        if(j == 0) {
            int index = ((int) m )- 1;
            direction = dir1[index]<0? -1:1;
        } else {
            int index = ((int) m )- 1;
            direction = dir2[index]<0? -1:1;
        }


        TreeState node = buildTreeNew(m,root.getPosition(direction), root.getMomentum(direction), direction,
                            logSliceU, j, stepSizeInformation.stepSize, initialJointDensity);

        root.setPosition(direction, node.getPosition(direction));
        root.setMomentum(direction, node.getMomentum(direction));


        if (node.flagContinue) {
          //  System.err.println("flagcintinue");
            double uniform2 = MathUtils.nextDouble();
            uniform2 = 0.3;

//            assert((node.numNodes / root.numNodes ) >= 1 || (node.numNodes / root.numNodes ) == 0);
//
//          //if (node.numNodes / root.numNodes > 1 ) {
//                System.err.println("j is "+ j + "n1 = " + node.numNodes + " n2 = " + root.numNodes);
//           //}
            double p = (double) node.numNodes / (double)root.numNodes;
            //assert( p == 0 || p >= 1);

            if(DEBUG ==1 ) {
                   System.err.println("n1 is " + node.numNodes + " n2 is " + root.numNodes +" p is "+p);

            }

            if (uniform2 < p) {

                 if(m == 53){
                     System.err.println("node.numNodes is " + node.numNodes + "the root.numnode is " + root.numNodes);
                     System.err.println("node.getPOSITION IS " + java.util.Arrays.toString(node.getPosition(0))
                             + "now j is " + j);
                 }
                endPosition = node.getPosition(0);
                if (DEBUG==1) {
                    System.err.println("node.getPOSITION IS " + java.util.Arrays.toString(node.getPosition(0))
                            + "now j is " + j);
                }
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

//        double logJointProbBefore = getJointProbability(gradientProvider, initialMomentum);

        // Make deep copy of position and momentum
        double[] position = Arrays.copyOf(inPosition, inPosition.length);
        double[] momentum = Arrays.copyOf(inMomentum, inMomentum.length);



        // "one frog jump!"
        leapFrogEngine.updatePosition(position, momentum, 0.0, 1.0);//todo : correct?

        double logJointProbBefore = getJointProbability(gradientProvider, momentum);

        if(DEBUG == 1) {

            System.err.println("Parameters into the jump is \n position" + Arrays.toString(position) +
                    "\n momentum " + Arrays.toString(momentum) + "\n gradient" + Arrays.toString(gradientProvider.getGradientLogDensity())

                    + "\n stepsize is " + direction * stepSize / 2);
        }

        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);
        leapFrogEngine.updatePosition(position, momentum, direction * stepSize, 1.0);
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);

        double logJointProbAfter = getJointProbability(gradientProvider, momentum);

        if(DEBUG == 1) {
            System.err.println("the position after one jump is " + Arrays.toString(position));
            System.err.println("the momentum after one jump is " + Arrays.toString(momentum));
            System.err.println("the logJointProbBefore is " + initialJointDensity);
            System.err.println("the logJointProbAfter is " + logJointProbAfter +
                    "\n *********************************");
        }
        // Values for recursion
       // System.err.println("the logsliceU is " +  logSliceU );
        final int numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);
       // System.err.println("the logJointProbBefore is " +  logJointProbBefore );



        final boolean flagContinue = (logSliceU < options.deltaMax + logJointProbAfter);

        // Values for dual-averaging
        final double alpha = Math.min(1.0, Math.exp(logJointProbAfter - initialJointDensity)); // TODO Is this current?
        final int nAlpha = 1;

        leapFrogEngine.updatePosition(inPosition, inMomentum, 0.0, 1.0);//todo : correct?
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
            uniform = 0.3;
            if (child.numNodes > 0
                    && uniform <  ((double) child.numNodes / (double) (node.numNodes + child.numNodes))) {
//                if (m == 3){
//                    System.err.println("\n the before position is " + java.util.Arrays.toString(node.getPosition(0)));
//                }
                node.setPosition(0, child.getPosition(0));

//                if (m == 3){
//                    System.err.println("\n the after position is " + java.util.Arrays.toString(node.getPosition(0)));
//                }
            }

            node.numNodes += child.numNodes;
            node.flagContinue = computeStopCriterion(child.flagContinue, node);

            node.alpha += child.alpha;
            node.nAlpha += child.nAlpha;

        }
        return node;
    }

//    private TreeNode buildTree(double[] position, double[] momentum, double logSliceU, int direction, int j, double
//            stepSize, double[] initialMomentum) {
//
//        if (j == 0) {
//
//            TreeNode tree = new TreeNode();
//
//            double logJointProbBefore = getJointProbability(gradientProvider, initialMomentum);
//
//            // "one frog jump!"
//            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);
//            leapFrogEngine.updatePosition(position, momentum, direction * stepSize, 1.0);
//            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), direction * stepSize / 2);
//
//            double logJointProbAfter = getJointProbability(gradientProvider, momentum);
//
//            tree.momentumMinus = Arrays.copyOf(momentum, dim);
//            tree.positionMinus = Arrays.copyOf(position, dim);
//            tree.momentumPlus = Arrays.copyOf(momentum, dim);
//            tree.positionPlus = Arrays.copyOf(position, dim);
//            tree.positionFinal = Arrays.copyOf(position, dim);
//
//            tree.numNodes = (logSliceU <= logJointProbAfter ? 1 : 0);
//            tree.flagContinue = (logSliceU <= options.deltaMax + logJointProbAfter ? true : false);
//
//            tree.alpha = Math.min(1, Math.exp(logJointProbAfter - logJointProbBefore));
//            tree.nAlpha = 1;
//
//            return tree;
//
//        } else {
//
//            TreeNode newNodes;
//            TreeNode tree = buildTree(position, momentum, logSliceU, direction, j - 1, stepSize, initialMomentum);
//
//            double[] positionMinus = Arrays.copyOf(tree.positionMinus, dim);
//            double[] positionPlus = Arrays.copyOf(tree.positionPlus, dim);
//            double[] momentumMinus = Arrays.copyOf(tree.momentumMinus, dim);
//            double[] momentumPlus = Arrays.copyOf(tree.momentumPlus, dim);
//
//            if (tree.flagContinue) {
//
//                if (direction == -1) {
//
//                    newNodes = buildTree(tree.positionMinus, tree.momentumMinus, logSliceU, direction, j - 1, stepSize, initialMomentum);
//
//                    positionMinus = Arrays.copyOf(newNodes.positionMinus, dim);
//                    momentumMinus = Arrays.copyOf(newNodes.momentumMinus, dim);
//                } else {
//
//                    newNodes = buildTree(tree.positionPlus, tree.momentumPlus, logSliceU, direction, j - 1, stepSize, initialMomentum);
//
//                    positionPlus = Arrays.copyOf(newNodes.positionPlus, dim);
//                    momentumPlus = Arrays.copyOf(newNodes.momentumPlus, dim);
//                }
//                double randomNum = MathUtils.nextDouble();
//
//                if (randomNum < (double) newNodes.numNodes / (tree.numNodes + newNodes.numNodes)) {
//                    tree.positionFinal = Arrays.copyOf(newNodes.positionFinal, dim);
//                }
//
//                tree.positionPlus = Arrays.copyOf(positionPlus, dim);
//                tree.positionMinus = Arrays.copyOf(positionMinus, dim);
//                tree.momentumMinus = Arrays.copyOf(momentumMinus, dim);
//                tree.momentumPlus = Arrays.copyOf(momentumPlus, dim);
//
//                tree.numNodes += newNodes.numNodes;
//                tree.flagContinue = computeStopCriterion(newNodes.flagContinue, positionPlus, positionMinus, momentumPlus, momentumMinus);
//
//                //values for dual averaging use..
//                tree.alpha += newNodes.alpha;
//                tree.nAlpha += newNodes.nAlpha;
//
//            }
//            return tree;
//        }
//    }

    private double findReasonableStepSize(double[] position) { // TODO Still needs to be reviewed

        double stepSize = 1;
        final double[] momentum = drawInitialMomentum(drawDistribution, dim);
        int count = 1;

        double[] initialposi = Arrays.copyOf(position,dim);

        double probBefore = getJointProbability(gradientProvider, momentum);
        System.err.println(probBefore);

        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);
        leapFrogEngine.updatePosition(position, momentum, stepSize, 1.0);
        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize / 2);

        double probAfter = getJointProbability(gradientProvider, momentum);
        System.err.println(probAfter);

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

        leapFrogEngine.updatePosition(initialposi, momentum, 0.0, 1.0); //todo: is it a way to return back to the initial position?

        System.err.println("final prob from findStepSize " + gradientProvider.getLikelihood().getLogLikelihood());
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

//    private static double[] sumArray(double[] a, double[] b) {
//
//        assert (a.length == b.length);
//        final int dim = a.length;
//
//        double result[] = new double[dim];
//
//        for (int i = 0; i < dim; i++) {
//            result[i] = a[i] + b[i];
//        }
//
//        return result;
//    }

    private static double getJointProbability(GradientWrtParameterProvider gradientProvider, double[] momentum) {

        assert (gradientProvider != null);
        assert (momentum != null);

        return gradientProvider.getLikelihood().getLogLikelihood() - getScaledDotProduct(momentum, 1.0);
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

        private int getIndex(int direction) { // valid directions: -1, 0, +1
            return direction + 1; // TODO Check bounds;
        }

        final private double[][] position;
        final private double[][] momentum;

        private int numNodes;
        private boolean flagContinue;

        private double alpha;
        private int nAlpha;
    }

//    public class TreeNode {
//        double[] positionMinus;
//        double[] momentumMinus;
//
//        double[] positionPlus;
//        double[] momentumPlus;
//
//        double[] positionFinal;
//
//        int numNodes;
//        boolean flagContinue;
//        double alpha;
//        int nAlpha;
//
//        public TreeNode() {
//        }
//
////        public TreeNode(double[] positionMinus, double[] momentumMinus,
////                        double[] positionPlus, double[] momentumPlus) {
////            this.positionMinus = positionMinus;
////            this.momentumMinus = momentumMinus;
////            this.positionPlus = positionPlus;
////            this.momentumPlus = momentumPlus;
////        }
////
////        public TreeNode(double[] initialPosition, double[] initialMomentum) {
////            this.positionMinus = initialPosition;
////            this.momentumMinus = initialMomentum;
////            this.positionPlus = initialPosition;
////            this.momentumMinus = initialMomentum;
////        }
//
//    }
}


