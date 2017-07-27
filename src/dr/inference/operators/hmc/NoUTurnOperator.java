
/**
 * @author Marc A. Suchard
 */

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

import dr.math.distributions.UniformDistribution;
import jebl.math.Random;

/**
 * @author Marc A. Suchard
 */

public class NoUTurnOperator extends HamiltonianMonteCarloOperator implements GeneralOperator {

    private final int dim = gradientProvider.getDimension();
    private final int chainLength = 1000000;
    private final double sigmaSquared = drawDistribution.getSD() * drawDistribution.getSD();
    private final Likelihood likelihood;
    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
                           Parameter parameter, double stepSize, int nSteps, double drawVariance, Likelihood likelihood)
    {
        super(mode, weight, gradientProvider, parameter, null, stepSize, nSteps, drawVariance);
        this.likelihood = likelihood;
    }

    @Override
    public String getOperatorName() {
        return "no U turn HMC operator";
    }

    @Override
    public double doOperation(Likelihood likelihood) {
        return leapFrog();
        // TODO Implement NUTS algorithm
    }

    @Override
    public double doOperation() {
        throw new IllegalStateException("Should not be called.");
    } //zy: when throw?
/**************************************************** BEGIN ***************************************************************************************************************/
    private double[][] nuts(double[] position,Likelihood likelihood,GradientWrtParameterProvider gradientProvider,double targetAcceptanceRatio,
                             int maxDepth){
        double[][] result = new double[chainLength][dim];
        OneNut nut = new OneNut();
        for (int i = 0;i < chainLength;i++){
            nut = nutsOneStep(position,i,likelihood,gradientProvider,nut,targetAcceptanceRatio,maxDepth);
            position = nut.position;
//            result[i][] = position;
        }
        return result;
    }

    private OneNut nutsOneStep(double[] position,int iter,Likelihood likelihood, GradientWrtParameterProvider gradientProvider,OneNut nut,
                               double targetAcceptance, int maxdepth){
        final double kappa = 0.75;
        final double t0 = 10;
        final double gamma = 0.05;
        final int adaptLength = 50;
        double h = 0;
        double log_eps;
        OneNut returnnut = new OneNut();

        if(iter == 1){
            returnnut.stepSize = findReasonableEpsilon(position);
            returnnut.stepSizeAve = 1;
            returnnut.mu = Math.log(10*returnnut.stepSize);
            returnnut.h = 0;
        } else {
            returnnut = nut;
        }

        double[] momentum = drawInitialMomentum(drawDistribution, dim);
        double sliceU = 0.0; //Math.random()*Math.exp(likelihood - getScaledDotProduct(momentum, sigmaSquared));
        double[] positionMinus = position;
        double[] positionPlus = position;
        double[] momentumMinus = momentum;
        double[] momentumPlus = momentum;
        int j=0;
        int n=1;
        boolean s = true;
        BinaryTree temp = new BinaryTree();
        while(s){
            if (Math.random()<0.5){ 
                temp = buildTree(positionMinus,momentumMinus,sliceU,-1,j,stepSize,position,momentum);
                positionMinus = temp.positionMinus;
                momentumMinus = temp.momentumMinus;
            } else {
                temp = buildTree(positionPlus,momentumPlus,sliceU,-1,j,stepSize,position,momentum);
                positionPlus = temp.positionPlus;
                momentumPlus = temp.momentumPlus;
            }
            if (temp.flagError){
                if (Math.random() < temp.numNodes/n){
                    position = temp.positionFinal;
                }
            }
            n += temp.numNodes;
            s = computeStopCriterion(temp.flagError,positionPlus,positionMinus,momentumPlus,momentumMinus);
            j++;
            if(j>maxdepth){
                System.err.println("reach maximum tree depth");
                break;
            }
        }
        if (iter < adaptLength){
            h = (1 - 1/(iter + t0))*h + 1/(iter + t0) * (targetAcceptance - temp.alpha / temp.nAlpha);
            log_eps = returnnut.mu - Math.sqrt(iter)/gamma * h;
            returnnut.stepSizeAve = Math.exp(Math.pow(iter,-kappa)*log_eps + (1 - Math.pow(iter,-kappa)) * Math.log(returnnut.stepSizeAve));
            returnnut.stepSize = Math.exp(log_eps);
        } else {
            returnnut.stepSize = returnnut.stepSizeAve;
        }
        return returnnut;
    }


    private BinaryTree buildTree(double[] position, double[] momentum, double sliceU, int direction, int j, double
            stepSize,
                                 double[] initialPosition, double[] initialMomentum) {
        BinaryTree tree = new BinaryTree();
        double[] positionMinus;
        double[] momentumMinus;
        double[] positionPlus;
        double[] momentumPlus;

        double logprobBefore = 0.0; //likelihood(initialPosition) + getScaledDotProduct(initialMomentum, sigmaSquared);
        if (j == 0) {
            leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize);
            leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared);

            double logprobAfter = 0.0; //likelihood(position) + getScaledDotProduct(momentum, sigmaSquared);
            tree.momentumMinus = momentum;
            tree.positionMinus = position;
            tree.momentumPlus = momentum;
            tree.positionPlus = position;
            tree.positionFinal = position;
//            tree.numNodes = Math.log(sliceU) <= likelihood - getScaledDotProduct(momentum, sigmaSquared) ? 1 : 0;
//            tree.flagError = Math.log(sliceU) < 1000 + likelihood - getScaledDotProduct(momentum, sigmaSquared) ? true : false;
            tree.alpha = Math.min(1, Math.exp(logprobAfter / logprobBefore));
            tree.nAlpha = 1;

            return tree;
        } else {
            BinaryTree newNodes;
            tree = buildTree(position, momentum, sliceU, direction, j - 1, stepSize, initialPosition, initialMomentum);
            positionMinus = tree.positionMinus;
            momentumMinus = tree.momentumMinus;
            positionPlus = tree.positionPlus;
            momentumPlus = tree.momentumPlus;
            while (tree.flagError) {
                if (direction == -1) {
                    newNodes = buildTree(tree.positionMinus, tree.momentumMinus, sliceU, direction, j - 1, stepSize,
                            initialPosition, initialMomentum);
                    positionMinus = newNodes.positionMinus;
                    momentumMinus = newNodes.momentumMinus;
                } else {
                    newNodes = buildTree(tree.positionPlus, tree.momentumPlus, sliceU, direction, j - 1, stepSize,
                            initialPosition, initialMomentum);
                    positionPlus = newNodes.positionPlus;
                    momentumPlus = newNodes.momentumPlus;
                }
                double randomNum = Math.random();
                if (randomNum < (double) newNodes.numNodes / (tree.numNodes + newNodes.numNodes)) {
                    tree.positionFinal = newNodes.positionFinal;
                }
                tree.alpha += newNodes.alpha;
                tree.nAlpha += newNodes.nAlpha;
                tree.numNodes += newNodes.numNodes;

                tree.flagError = computeStopCriterion(newNodes.flagError,newNodes.positionPlus,newNodes.positionMinus,newNodes.momentumPlus,newNodes.momentumMinus);
                tree.positionPlus = positionPlus;
                tree.positionMinus = positionMinus;
                tree.momentumMinus = momentumMinus;
                tree.momentumPlus = momentumPlus;
            }
            return tree;
        }

    }

    private double findReasonableEpsilon(double[] position) {
        double stepSize = 1;
        int count = 1;
        final double[] momentum = drawInitialMomentum(drawDistribution, dim);

        double probBefore = 0.0; //Math.exp(likelihood + getScaledDotProduct(momentum, sigmaSquared));

        leapFrogEngine.updateMomentum(position, momentum, gradientProvider.getGradientLogDensity(), stepSize);
        leapFrogEngine.updatePosition(position, momentum, stepSize, sigmaSquared); //zy: after "updateposition"
        // the interested parameters changed.

        double probAfter = 0.0; //Math.exp(likelihood + getScaledDotProduct(momentum, sigmaSquared));

        double probRatio = probAfter / probBefore;

        double a = (probRatio > 0.5 ? 1 : -1);

        while (Math.pow(probRatio, a) > Math.pow(2, -a)) {
            stepSize = Math.pow(2, a) * stepSize;
            count++;
            if (count > 100) {
                System.out.println("cannot find a reasonable epsilon in 100 iterations!");
                break;
            }
        }

        return stepSize;
    }

    private boolean computeStopCriterion(boolean flagError, double[] positionPlus,double[] positionMinus,double[] momentumPlus,double[] momentumMinus) {
        boolean flag = flagError && getDotProduct(sumArray(positionPlus, positionMinus, false),
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

    public class BinaryTree {
        double[] positionMinus;
        double[] momentumMinus;
        double[] positionPlus;
        double[] momentumPlus;
        double[] positionFinal;
        int numNodes;
        boolean flagError;
        double alpha;
        int nAlpha;

        public BinaryTree() {
        }
    }

    public class OneNut{
        double[] position;
        double stepSize;
        double stepSizeAve;
        double h;
        double mu;

        public OneNut(){};
    }
}


//package dr.inference.operators.hmc;
//
//        import dr.inference.hmc.GradientWrtParameterProvider;
//        import dr.inference.model.Likelihood;
//        import dr.inference.model.Parameter;
//        import dr.inference.operators.CoercionMode;
//        import dr.inference.operators.GeneralOperator;
//
///**
// * @author Marc A. Suchard
// */
//


//public class NoUTurnOperator extends HamiltonianMonteCarloOperator implements GeneralOperator {
//
//    public NoUTurnOperator(CoercionMode mode, double weight, GradientWrtParameterProvider gradientProvider,
//                           Parameter parameter, double stepSize, int nSteps, double drawVariance) {
//        super(mode, weight, gradientProvider, parameter, null, stepSize, nSteps, drawVariance);
//    }
//
//    @Override
//    public double doOperation(Likelihood likelihood) {
//        return leapFrog(); // TODO Implement NUTS algorithm
//    }
//
//    @Override
//    public double doOperation() {
//        throw new IllegalStateException("Should not be called.");
//    }
//
//    private void buildTree() {
//
//    }
//
//    private double findReasonableEpsilon() {
//        return 0.0;
//    }
//
//    private double computeStopCriterion() {
//        return 0.0;
//    }
//}

