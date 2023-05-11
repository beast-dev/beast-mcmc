/*
 * SmoothSkygridLikelihood.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.Reportable;
import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math.analysis.integration.RombergIntegrator;
import org.apache.commons.math.analysis.integration.UnivariateRealIntegrator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Xiang Ji
 * @author Yuwei Bao
 * @author Marc A. Suchard
 */
public class SmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Citable, Reportable {

    private final List<TreeModel> trees;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;
    private final Parameter smoothRate;
    private final SmoothSkygridPopulationSizeInverse populationSizeInverse;
    private final OldSmoothLineageCount lineageCount;

    private final GlobalSigmoidSmoothFunction smoothFunction;

    private final List<BigFastTreeIntervals> intervalsList;

    public SmoothSkygridLikelihood(String name,
                                   List<TreeModel> trees,
                                   Parameter logPopSizeParameter,
                                   Parameter gridPointParameter,
                                   Parameter smoothRate) {
        super(name);
        this.trees = trees;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;
        this.smoothRate = smoothRate;
        this.smoothFunction = new GlobalSigmoidSmoothFunction();
        this.populationSizeInverse = new SmoothSkygridPopulationSizeInverse(logPopSizeParameter, gridPointParameter, smoothFunction, smoothRate);
        this.lineageCount = new OldSmoothLineageCount(trees.get(0), smoothFunction, smoothRate);
        intervalsList = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            intervalsList.add(new BigFastTreeIntervals(trees.get(i)));
            addModel(intervalsList.get(i));
        }

        for (TreeModel tree:trees) {
            addModel(tree);
        }
        addVariable(logPopSizeParameter);
        addVariable(gridPointParameter);
        addVariable(smoothRate);
    }

    @Override
    public String getReport() {
        return "smoothSkygrid(" + getLogLikelihood() + ")";
    }

    public Tree getTree(int nt) {
        return trees.get(nt);
    }

    class OldSmoothLineageCount {

        private final Tree tree;
        private final GlobalSigmoidSmoothFunction smoothFunction;

        private final Parameter smoothRate;

        OldSmoothLineageCount(Tree tree, GlobalSigmoidSmoothFunction smoothFunction, Parameter smoothRate) {
            this.tree = tree;
            this.smoothFunction = smoothFunction;
            this.smoothRate = smoothRate;
        }

        double getLineageCount(double time) {
            double sum = 0;
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                sum += smoothFunction.getSmoothValue(time, tree.getNodeHeight(tree.getNode(i)), 0.0, 1.0, smoothRate.getParameterValue(0));
            }
            for (int i = tree.getExternalNodeCount(); i < tree.getNodeCount(); i++) {
                sum += smoothFunction.getSmoothValue(time, tree.getNodeHeight(tree.getNode(i)), 0.0, -1.0, smoothRate.getParameterValue(0));
            }
            return sum;
        }
    }

    class SmoothSkygridPopulationSizeInverse {

        private final Parameter logPopSizeParameter;
        private final Parameter gridPointParameter;
        private final GlobalSigmoidSmoothFunction smoothFunction;
        private final Parameter smoothRate;

        SmoothSkygridPopulationSizeInverse(Parameter logPopSizeParameter,
                                           Parameter gridPointParameter,
                                           GlobalSigmoidSmoothFunction smoothFunction,
                                           Parameter smoothRate) {
            this.logPopSizeParameter = logPopSizeParameter;
            this.gridPointParameter = gridPointParameter;
            this.smoothRate = smoothRate;
            this.smoothFunction = smoothFunction;
        }

        double getPopulationSizeInverse(double time) {
            double sum = 0;
            for(int i = 0; i < gridPointParameter.getDimension(); i++) {
                double increment = smoothFunction.getSmoothValue(time, gridPointParameter.getParameterValue(i),
                        i == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(0)) : 0.0,
                        i == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(1)) : Math.exp(-logPopSizeParameter.getParameterValue(i + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(i)),
                        smoothRate.getParameterValue(0));
                sum += increment;
            }
            return sum;
        }
    }

    @Override
    public Type getUnits() {
        return null;
    }

    @Override
    public void setUnits(Type units) {

    }

    public double[] getGradientWrtNodeHeight() {
        assert(trees.size() == 1);
        Tree tree = trees.get(0);
        BigFastTreeIntervals intervals = intervalsList.get(0);

        double[] gradient = new double[tree.getInternalNodeCount()];

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            gradient[i] += getLogSmoothPopulationSizeInverseDerivative(tree.getNodeHeight(node), tree.getNodeHeight(tree.getRoot()))
                    / getSmoothPopulationSizeInverse(tree.getNodeHeight(node), tree.getNodeHeight(tree.getRoot()));
        }

        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());
        final double firstInversePopulationSize = Math.exp(-logPopSizeParameter.getParameterValue(0));

        double rootDerivative = 0;

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);

            if (!tree.isRoot(node)) {
                final int intervalNum = intervals.getIntervalIndexForNode(node.getNumber());
                final double lineageCountDifference = getLineageCountDifference(intervalNum, intervals);
                gradient[i] += - firstInversePopulationSize * lineageCountDifference * smoothFunction.getSingleIntegrationDerivative(startTime, endTime, tree.getNodeHeight(node), smoothRate.getParameterValue(0));
            }

        }


            for (int i = 0; i < intervals.getIntervalCount(); i++) {
                NodeRef node = tree.getNode(i);
                final int intervalNum = intervals.getIntervalIndexForNode(node.getNumber());
                final double lineageCountDifference = getLineageCountDifference(intervalNum, intervals);
                rootDerivative += -lineageCountDifference * smoothFunction.getSingleIntegrationDerivativeWrtEndTime(endTime, tree.getNodeHeight(node), smoothRate.getParameterValue(0));
            }

        gradient[tree.getRoot().getNumber() - tree.getExternalNodeCount()] += firstInversePopulationSize * rootDerivative;


        int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {

            final NodeRef node = tree.getNode(i + tree.getExternalNodeCount());

            final int intervalIndex = intervals.getIntervalIndexForNode(node.getNumber());
            final double lineageCountDifference = getLineageCountDifference(intervalIndex, intervals);
            double thisDerivative = 0;

            for (int j = 0; j < maxGridIndex + 1; j++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
                final double gridTime = gridPointParameter.getParameterValue(j);
                thisDerivative += (nextPopSizeInverse - currentPopSizeInverse) *
                        smoothFunction.getPairProductIntegrationDerivative(startTime, endTime, tree.getNodeHeight(node), gridTime, smoothRate.getParameterValue(0));
            }

            gradient[i] -= lineageCountDifference * thisDerivative;

        }

        rootDerivative = 0.0;
        for (int i = 0; i <  intervals.getIntervalCount() + 1; i++) { //

            final double lineageCountDifference = getLineageCountDifference(i, intervals);
            final double intervalTime = intervals.getIntervalTime(i);

            for (int j = 0; j < maxGridIndex + 1; j++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
                final double gridTime = gridPointParameter.getParameterValue(j);
                rootDerivative -= lineageCountDifference * (nextPopSizeInverse - currentPopSizeInverse) *
                        smoothFunction.getPairProductIntegrationDerivativeWrtEndTime(startTime, endTime, intervalTime, gridTime, smoothRate.getParameterValue(0));
            }
        }
        gradient[tree.getRoot().getNumber() - tree.getExternalNodeCount()] += rootDerivative;

        return gradient;
    }

    private double getLineageCountDifference(int intervalIndex, BigFastTreeIntervals intervals) {
        if (intervalIndex == 0) {
            return ((double) intervals.getLineageCount(0) * (intervals.getLineageCount(0) - 1)) / 2.0;
        } else if (intervalIndex == intervals.getIntervalCount()) {
            return -((double) intervals.getLineageCount(intervalIndex - 1) * (intervals.getLineageCount(intervalIndex - 1) - 1)) / 2.0;
        } else {
            return ((double) intervals.getLineageCount(intervalIndex) * (intervals.getLineageCount(intervalIndex) - 1)
            - intervals.getLineageCount(intervalIndex - 1) * (intervals.getLineageCount(intervalIndex - 1) - 1)) / 2.0;
        }
    }

    private boolean DEBUG = true;
    private final static UnivariateRealIntegrator integrator = new RombergIntegrator();

    protected double calculateLogLikelihood() {
        assert(trees.size() == 1);
        if (!likelihoodKnown) {
            TreeModel tree = trees.get(0);
            final double startTime = 0;
            final double endTime = tree.getNodeHeight(tree.getRoot());
            final int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

            double[] tmpA = new double[tree.getNodeCount()];
            double[] tmpB = new double[tree.getNodeCount()];
            double[] tmpC = new double[tree.getNodeCount()];
            double[] tmpD = new double[maxGridIndex];
            double[] tmpE = new double[maxGridIndex];
            double[] tmpF = new double[maxGridIndex];
            double[] tmpLineageEffect = new double[tree.getNodeCount()];
            double[] tmpTimes = new double[tree.getNodeCount()];
            int[] tmpCounts = new int[tree.getNodeCount()];

            NodeRef[] nodes = new NodeRef[tree.getNodeCount()];
            System.arraycopy(tree.getNodes(), 0, nodes, 0, tree.getNodeCount());
            Arrays.parallelSort(nodes, (a, b) -> Double.compare(tree.getNodeHeight(a), tree.getNodeHeight(b)));

            double lastTime = tree.getNodeHeight(nodes[0]);
            double currentLineageEffect = getLineageCountEffect(tree, 0);
            int currentCount = 1;
            int index = 0;
            tmpTimes[index] = lastTime;
            for (int i = 1; i < nodes.length; i++) {
                NodeRef node = nodes[i];
                final double time = tree.getNodeHeight(node);
                if (time == lastTime) {
                    currentCount++;
                    currentLineageEffect += getLineageCountEffect(tree, node.getNumber());
                } else {
                    tmpLineageEffect[index] = currentLineageEffect;
                    tmpCounts[index] = currentCount;
                    index++;
                    tmpTimes[index] = time;
                    currentCount = 1;
                    currentLineageEffect = getLineageCountEffect(tree, node.getNumber());
                    lastTime = time;
                }
            }
            tmpLineageEffect[index] = currentLineageEffect;
            tmpCounts[index] = currentCount;
            final int uniqueTimes = index + 1;

            for (int i = 0; i < uniqueTimes; i++) {
                final double timeI = tmpTimes[i];
                double sum = 0;
                for (int j = 0; j < uniqueTimes; j++) {
                    if (j != i) {
                        final double timeJ = tmpTimes[j];
                        final double lineageCountEffect = tmpLineageEffect[j];
                        final double thisInverse = smoothFunction.getInverseOneMinusExponential(timeJ - timeI, smoothRate.getParameterValue(0));
                        sum += lineageCountEffect * thisInverse;
                    }
                }
                tmpA[i] = sum;
            }

            for (int i = 0; i < uniqueTimes; i++) {
                final double timeI = tmpTimes[i];
                double sum = 0;
                for (int k = 0; k < maxGridIndex; k++) {
                    final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                    final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                    final double gridTime = gridPointParameter.getParameterValue(k);
                    final double thisInverse = smoothFunction.getInverseOneMinusExponential(gridTime - timeI, smoothRate.getParameterValue(0));
                    sum += (nextPopSizeInverse - currentPopSizeInverse) * thisInverse;
                }
                tmpB[i] = sum;
            }

            for (int i = 0; i < uniqueTimes; i++) {
                final double timeI = tmpTimes[i];
                final double logDiff = smoothFunction.getLogOnePlusExponential(timeI - endTime, smoothRate.getParameterValue(0)) -
                        smoothFunction.getLogOnePlusExponential(timeI - startTime, smoothRate.getParameterValue(0));
                tmpC[i] = logDiff;
            }

            for (int k = 0; k < maxGridIndex; k++) {
                final double gridTime = gridPointParameter.getParameterValue(k);
                tmpD[k] = smoothFunction.getLogOnePlusExponential(gridTime - endTime, smoothRate.getParameterValue(0)) -
                        smoothFunction.getLogOnePlusExponential(gridTime - startTime, smoothRate.getParameterValue(0));
                double sum = 0;
                double quadraticSum = 0;
                for (int i = 0; i < uniqueTimes; i++) {
                    final double timeI = tmpTimes[i];
                    final double lineageCountEffect = tmpLineageEffect[i];
                    final double tmp = smoothFunction.getInverseOneMinusExponential(timeI - gridTime, smoothRate.getParameterValue(0)) *
                            lineageCountEffect;
                    sum += tmp;
                    quadraticSum += tmp * tmp;
                }
                tmpE[k] = sum;
                tmpF[k] = sum * sum - quadraticSum;
            }

            double tripleIntegrationSum = 0;
            double lineageEffectSqaredSum = 0;
            for (int i = 0; i < uniqueTimes; i++) {
                final double lineageCountEffect = tmpLineageEffect[i];
                lineageEffectSqaredSum += lineageCountEffect * lineageCountEffect;
                tripleIntegrationSum += lineageCountEffect * tmpA[i] * tmpB[i] * tmpC[i];
            }
            tripleIntegrationSum *= 2;

            for (int k = 0; k < maxGridIndex; k++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                tripleIntegrationSum += (nextPopSizeInverse - currentPopSizeInverse) * tmpF[k] * tmpD[k];
            }

            tripleIntegrationSum /= -smoothRate.getParameterValue(0) * 2;
            tripleIntegrationSum += -0.5 * (1 - lineageEffectSqaredSum)
                    * (Math.exp(-logPopSizeParameter.getParameterValue(maxGridIndex)) - Math.exp(-logPopSizeParameter.getParameterValue(0)))
                    * (endTime - startTime);

            double tripleWithQuadraticIntegrationSum = 0;
            final double commonFirstTermMultiplier = (Math.exp(-logPopSizeParameter.getParameterValue(maxGridIndex)) - Math.exp(-logPopSizeParameter.getParameterValue(0))) * (endTime - startTime);
            for (int i = 0; i < uniqueTimes; i++) {
                final double lineageCountEffect = tmpLineageEffect[i] * tmpLineageEffect[i];
                final double timeI = tmpTimes[i];
                double thisResult =  commonFirstTermMultiplier;
                final double commonSecondTermMultiplier = smoothFunction.getInverseOnePlusExponential(timeI - startTime, smoothRate.getParameterValue(0))
                        - smoothFunction.getInverseOnePlusExponential(timeI - endTime, smoothRate.getParameterValue(0));
                for (int k = 0; k < maxGridIndex; k++) {
                    final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                    final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                    final double gridTime = gridPointParameter.getParameterValue(k);
                    final double inverse = smoothFunction.getInverseOneMinusExponential(gridTime - timeI, smoothRate.getParameterValue(0));
                    thisResult += (nextPopSizeInverse - currentPopSizeInverse) / smoothRate.getParameterValue(0)
                            * (inverse * commonSecondTermMultiplier + (2.0 - inverse) * inverse * tmpC[i] +
                            (1 - inverse) * (1 - inverse) * tmpD[k]);
                }
                thisResult *= lineageCountEffect;
                tripleWithQuadraticIntegrationSum += thisResult;
            }
            tripleWithQuadraticIntegrationSum *= -0.5;

//            if (DEBUG) {
//                double checkTripleIntegrationSum = 0;
//                double tripleWithQuadraticIntegrationCheck = 0;
//                for (int i = 0; i < uniqueTimes; i++) {
//                    final double lineageCountEffectI = tmpLineageEffect[i];
//                    final double timeI = tmpTimes[i];
//                    for (int j = 0; j < uniqueTimes; j++) {
//                        final double lineageCountEffectJ = tmpLineageEffect[j];
//                        final double timeJ = tmpTimes[j];
//                        if (j != i) {
//                            for (int k = 0; k < maxGridIndex; k++) {
//                                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
//                                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
//                                final double gridTime = gridPointParameter.getParameterValue(k);
//                                checkTripleIntegrationSum += (nextPopSizeInverse - currentPopSizeInverse) * lineageCountEffectI * lineageCountEffectJ
//                                        * smoothFunction.getTripleProductIntegration(startTime, endTime, timeI, timeJ, gridTime, smoothRate.getParameterValue(0));
//                            }
//                        } else {
//                            for (int k = 0; k < maxGridIndex; k++) {
//                                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
//                                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
//                                final double gridTime = gridPointParameter.getParameterValue(k);
//                                tripleWithQuadraticIntegrationCheck += (nextPopSizeInverse - currentPopSizeInverse) * lineageCountEffectI* lineageCountEffectI
//                                        * smoothFunction.getTripleProductWithQuadraticIntegration(startTime, endTime, timeI, gridTime, smoothRate.getParameterValue(0));
//                            }
//                        }
//                    }
//                }
//                checkTripleIntegrationSum *= -0.5;
//                tripleWithQuadraticIntegrationCheck *= -0.5;
//            }

//            double tripleIntegrationCheck = 0;
//            for (int i = 0; i < tree.getNodeCount(); i++) {
//                final double timeI = tree.getNodeHeight(tree.getNode(i));
//                final double lineageCountEffectI = getLineageCountEffect(tree, i);
//                for (int j = 0; j < tree.getNodeCount(); j++) {
//                    final double timeJ = tree.getNodeHeight(tree.getNode(j));
//                    final double lineageCountEffectJ = getLineageCountEffect(tree, j);
//                    for (int k = 0; k < maxGridIndex; k++) {
//                        final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
//                        final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
//                        final double gridTime = gridPointParameter.getParameterValue(k);
//                        if (timeI == timeJ) {
//                            final double thisResult = (nextPopSizeInverse - currentPopSizeInverse) * lineageCountEffectI* lineageCountEffectI
//                                    * smoothFunction.getTripleProductWithQuadraticIntegration(startTime, endTime, timeI, gridTime, smoothRate.getParameterValue(0));
//                            tripleIntegrationCheck += thisResult;
//                            final double testResult = (nextPopSizeInverse - currentPopSizeInverse) * lineageCountEffectI* lineageCountEffectI
//                                    * checkNumericIntegration(startTime, endTime, smoothFunction,
//                                    timeI, timeJ, gridTime, smoothRate.getParameterValue(0));
//
//
//
//                        } else {
//                            final double thisResult = (nextPopSizeInverse - currentPopSizeInverse) * lineageCountEffectI * lineageCountEffectJ
//                                    * smoothFunction.getTripleProductIntegration(startTime, endTime, timeI, timeJ, gridTime, smoothRate.getParameterValue(0));
//                            tripleIntegrationCheck += thisResult;
//                            final double testResult = (nextPopSizeInverse - currentPopSizeInverse) * lineageCountEffectI* lineageCountEffectJ
//                                    * checkNumericIntegration(startTime, endTime, smoothFunction,
//                                    timeI, timeJ, gridTime, smoothRate.getParameterValue(0));
//                        }
//                    }
//                }
//            }
//            tripleIntegrationCheck *= -0.5;



            double firstDoubleIntegrationOffDiagonalSum = 0;
            double firstDoubleIntegrationDiagonalSum = 0;
            for (int i = 0; i < uniqueTimes; i++) {
                final double lineageCountEffect = tmpLineageEffect[i];
                final double timeI = tmpTimes[i];
                firstDoubleIntegrationOffDiagonalSum += lineageCountEffect * tmpA[i] * tmpC[i];
                firstDoubleIntegrationDiagonalSum += lineageCountEffect * lineageCountEffect
                        * smoothFunction.getQuadraticIntegration(startTime, endTime, timeI, smoothRate.getParameterValue(0));
            }
            firstDoubleIntegrationOffDiagonalSum /= smoothRate.getParameterValue(0);
            firstDoubleIntegrationOffDiagonalSum += 0.5 * (1 - lineageEffectSqaredSum) * (endTime - startTime);

            final double firstDoubleIntegrationSum = -(firstDoubleIntegrationDiagonalSum * 0.5 + firstDoubleIntegrationOffDiagonalSum) * Math.exp(-logPopSizeParameter.getParameterValue(0));



//            if (DEBUG) {
//                double checkFirstDoubleIntegrationSum = 0;
//                for (int i = 0; i < uniqueTimes; i++) {
//                    final double lineageCountEffectI = tmpLineageEffect[i];
//                    final double timeI = tmpTimes[i];
//                    for (int j = 0; j < uniqueTimes; j++) {
//                        final double lineageCountEffectJ = tmpLineageEffect[j];
//                        final double timeJ = tmpTimes[j];
//                        if (j != i) {
//                            checkFirstDoubleIntegrationSum += lineageCountEffectI * lineageCountEffectJ * smoothFunction.getPairProductIntegration(startTime, endTime, timeI, timeJ, smoothRate.getParameterValue(0));
//                        } else {
//                        checkFirstDoubleIntegrationSum += lineageCountEffectI * lineageCountEffectJ * smoothFunction.getQuadraticIntegration(startTime, endTime, timeI, smoothRate.getParameterValue(0));
//                        }
//                    }
//                }
//                checkFirstDoubleIntegrationSum *= -0.5 * Math.exp(-logPopSizeParameter.getParameterValue(0));
//            }




            double secondDoubleIntegrationSum = 0;
            for (int i = 0; i < uniqueTimes; i++) {
                secondDoubleIntegrationSum += 0.5 * tmpB[i] * tmpC[i] * tmpLineageEffect[i];
            }

            for (int k = 0; k < maxGridIndex; k++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                secondDoubleIntegrationSum += 0.5 * tmpE[k] * tmpD[k] * (nextPopSizeInverse - currentPopSizeInverse);
            }

            secondDoubleIntegrationSum /= smoothRate.getParameterValue(0);
            secondDoubleIntegrationSum += 0.5 * (endTime - startTime) * (Math.exp(-logPopSizeParameter.getParameterValue(maxGridIndex)) - Math.exp(-logPopSizeParameter.getParameterValue(0)));

//            if (DEBUG) {
//                double checkSecondDoubleIntegration = 0;
//                for (int i = 0; i < uniqueTimes; i++) {
//                    final double lineageCountEffectI = tmpLineageEffect[i];
//                    final double timeI = tmpTimes[i];
//                    for (int k = 0; k < maxGridIndex; k++) {
//                        final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
//                        final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
//                        final double gridTime = gridPointParameter.getParameterValue(k);
//                        checkSecondDoubleIntegration += lineageCountEffectI * (nextPopSizeInverse - currentPopSizeInverse) * smoothFunction.getPairProductIntegration(startTime, endTime, timeI, gridTime, smoothRate.getParameterValue(0));
//                    }
//                }
//                checkSecondDoubleIntegration *= 0.5;
//            }



            double singleIntegration = 0;
            for (int i = 0; i < uniqueTimes; i++) {
                final double timeI = tmpTimes[i];
                final double lineageCountEffectI = tmpLineageEffect[i];
                singleIntegration += lineageCountEffectI * smoothFunction.getSingleIntegration(startTime, endTime, timeI, smoothRate.getParameterValue(0));
            }
            singleIntegration *= 0.5 * Math.exp(-logPopSizeParameter.getParameterValue(0));

//            try {
//                double numeric = getNumericIntegration(tree, smoothFunction, logPopSizeParameter, gridPointParameter, smoothRate.getParameterValue(0), startTime, endTime);
//                numeric *= -0.5 * Math.exp(-logPopSizeParameter.getParameterValue(0));
//                System.err.println("here");
//            } catch (Exception e) {
//                throw new RuntimeException(e);
//            }

            double logPopulationSizeInverse = 0;
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
                logPopulationSizeInverse += Math.log(getSmoothPopulationSizeInverse(tree.getNodeHeight(node), tree.getNodeHeight(tree.getRoot())));
            }

            logLikelihood = logPopulationSizeInverse + singleIntegration + firstDoubleIntegrationSum + secondDoubleIntegrationSum + tripleIntegrationSum + tripleWithQuadraticIntegrationSum;

            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public static double getNumericIntegration(TreeModel tree, GlobalSigmoidSmoothFunction smoothFunction,
                                               Parameter logPopSizeParameter, Parameter gridPointParameter,
                                               double smoothRate,
                                               double startTime, double endTime) throws Exception {
        UnivariateRealFunction f = v -> getTestFunction(logPopSizeParameter, gridPointParameter, tree, smoothFunction, endTime, smoothRate, v);
        return integrator.integrate(f, startTime, endTime);
    }

    public static double checkNumericIntegration(double startTime, double endTime, GlobalSigmoidSmoothFunction smoothFunction,
                                                 double timeI, double timeJ, double gridTime,
                                                 double smoothRate) {
        UnivariateRealFunction f = v -> smoothFunction.getSmoothValue(v, timeI, 0, 1, smoothRate)
                * smoothFunction.getSmoothValue(v, timeJ, 0, 1, smoothRate)
                * smoothFunction.getSmoothValue(v, gridTime, 0, 1, smoothRate);
        try {
            return integrator.integrate(f, startTime, endTime);
        } catch (ConvergenceException e) {
            throw new RuntimeException(e);
        } catch (FunctionEvaluationException e) {
            throw new RuntimeException(e);
        }
    }

    public static double getTestFunction(Parameter logPopSizeParameter, Parameter gridPointParameter,
                                         TreeModel tree,
                                         GlobalSigmoidSmoothFunction smoothFunction, double endTime,
                                         double smoothRate, double t) {
        final int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);
        double smoothedPopulationInverse = Math.exp(-logPopSizeParameter.getParameterValue(0));
//        double smoothedPopulationInverse = 0;
        for (int k = 0; k < maxGridIndex; k++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
            final double gridTime = gridPointParameter.getParameterValue(k);
            smoothedPopulationInverse += (nextPopSizeInverse - currentPopSizeInverse) * smoothFunction.getSmoothValue(t, gridTime, 0, 1, smoothRate);
        }

        double quadraticSum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            quadraticSum += smoothFunction.getSmoothValue(t, tree.getNodeHeight(tree.getNode(i)), 0, 1, smoothRate)
                    * smoothFunction.getSmoothValue(t, tree.getNodeHeight(tree.getNode(i)), 0, 1, smoothRate);
        }

        double smoothedLineageCount = getSmoothLineageCount(tree, smoothFunction, smoothRate, t);


        return smoothedLineageCount * (smoothedLineageCount - 1);
    }

    private static double getSmoothLineageCount(TreeModel tree, GlobalSigmoidSmoothFunction smoothFunction, double smoothRate, double t) {
        double smoothedLineageCount = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double lineageCount = tree.isExternal(tree.getNode(i)) ? 1.0 : -1.0;
            smoothedLineageCount += lineageCount  * smoothFunction.getSmoothValue(t, tree.getNodeHeight(tree.getNode(i)), 0, 1, smoothRate);
        }
        return smoothedLineageCount;
    }

    private double getLineageCountEffect(Tree tree, int node) {
        if (tree.isExternal(tree.getNode(node))) {
            return 1;
        } else {
            return -1;
        }
    }

    protected double calculateLogLikelihood2() {
        assert(trees.size() == 1);
        if (!intervalsKnown) {
            for(IntervalList intervalList : intervalsList){
                intervalList.calculateIntervals();
            }
            intervalsKnown = true;
        }
        Tree tree = trees.get(0);
        BigFastTreeIntervals intervals = intervalsList.get(0);
        double logPopulationSizeInverse = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
//            logPopulationSizeInverse += Math.log(populationSizeInverse.getPopulationSizeInverse(tree.getNodeHeight(node)));
            logPopulationSizeInverse += Math.log(getSmoothPopulationSizeInverse(tree.getNodeHeight(node), tree.getNodeHeight(tree.getRoot())));
        }

        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());

        final double firstInversePopulationSize = Math.exp(-logPopSizeParameter.getParameterValue(0));
        double singleSigmoidIntegralSums =  getLineageCountDifference(0, intervals)
                * smoothFunction.getSingleIntegration(startTime, endTime, intervals.getIntervalTime(0), smoothRate.getParameterValue(0));
        for (int i = 1; i < intervals.getIntervalCount() + 1; i++) {
            final double lineageCountDifference = getLineageCountDifference(i, intervals);
            final double smoothIntegral = smoothFunction.getSingleIntegration(startTime, endTime, intervals.getIntervalTime(i), smoothRate.getParameterValue(0));
            singleSigmoidIntegralSums += lineageCountDifference * smoothIntegral;
        }
        singleSigmoidIntegralSums *= firstInversePopulationSize;

        int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

        double pairSigmoidIntegralSums = 0;

        for (int j = 0; j < maxGridIndex + 1; j++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
            final double gridTime = gridPointParameter.getParameterValue(j);
            pairSigmoidIntegralSums += getLineageCountDifference(0, intervals) * (nextPopSizeInverse - currentPopSizeInverse) *
                    smoothFunction.getPairProductIntegration(startTime, endTime, intervals.getIntervalTime(0), gridTime, smoothRate.getParameterValue(0));
        }

        for (int i = 1; i < intervals.getIntervalCount() + 1; i++) {

            final double lineageCountDifference = getLineageCountDifference(i, intervals);
            final double intervalTime = intervals.getIntervalTime(i);

            for (int j = 0; j < maxGridIndex + 1; j++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
                final double gridTime = gridPointParameter.getParameterValue(j);
                pairSigmoidIntegralSums += lineageCountDifference * (nextPopSizeInverse - currentPopSizeInverse) *
                        smoothFunction.getPairProductIntegration(startTime, endTime, intervalTime, gridTime, smoothRate.getParameterValue(0));
            }
        }

        return logPopulationSizeInverse - singleSigmoidIntegralSums - pairSigmoidIntegralSums;
    }

    public static int getMaxGridIndex(Parameter gridPointParameter, double endTime) {
        int maxGridIndex = gridPointParameter.getDimension() - 1;
        while (gridPointParameter.getParameterValue(maxGridIndex) > endTime && maxGridIndex > 0) {
            maxGridIndex--;
        }
        return maxGridIndex + 1;
    }

    private double getSmoothPopulationSizeInverse(double t, double endTime) {
        int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

        double populationSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(0));

        for (int j = 0; j < maxGridIndex + 1; j++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
            final double gridTime = gridPointParameter.getParameterValue(j);
            populationSizeInverse += (nextPopSizeInverse - currentPopSizeInverse) *
                    smoothFunction.getSmoothValue(t, gridTime, 0, 1, smoothRate.getParameterValue(0));
        }

        return populationSizeInverse;
    }

    private double getLogSmoothPopulationSizeInverseDerivative(double t, double endTime) {
        int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

        double derivative = 0;

        for (int j = 0; j < maxGridIndex + 1; j++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(j + 1));
            final double gridTime = gridPointParameter.getParameterValue(j);
            derivative += (nextPopSizeInverse - currentPopSizeInverse) *
                    smoothFunction.getDerivative(t, gridTime, 0, 1, smoothRate.getParameterValue(0));
        }

        return derivative;
    }

    private double oldCalculateLogLikelihood() {
        assert(trees.size() == 1);
        Tree tree = trees.get(0);
        double logPopulationSizeInverse = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
            logPopulationSizeInverse += Math.log(populationSizeInverse.getPopulationSizeInverse(tree.getNodeHeight(node)));
        }
        double integralBit = 0;
        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double stepLocation1 = tree.getNodeHeight(tree.getNode(i));
            final double preStepValue1 = 0;
            final double postStepValue1 = i < tree.getExternalNodeCount() ? 1 : -1;
            for (int j = 0; j < tree.getNodeCount(); j++) {
                final double stepLocation2 = tree.getNodeHeight(tree.getNode(j));
                final double preStepValue2 = j == 0 ? -1 : 0;
                final double postStepValue2 = (j < tree.getExternalNodeCount() ? 1 : -1) + preStepValue2;
                for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                    final double stepLocation3 = gridPointParameter.getParameterValue(k);
                    final double preStepValue3 = k == 0 ? Math.exp(-logPopSizeParameter.getParameterValue(0)) : 0;
                    final double postStepValue3 = k == 0? Math.exp(-logPopSizeParameter.getParameterValue(1)) :
                            Math.exp(-logPopSizeParameter.getParameterValue(k + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(k));
                    final double analytic = -0.5 * smoothFunction.getTripleProductIntegration(startTime, endTime,
                            stepLocation1, preStepValue1, postStepValue1,
                            stepLocation2, preStepValue2, postStepValue2,
                            stepLocation3, preStepValue3, postStepValue3,
                            smoothRate.getParameterValue(0));
                    integralBit += analytic;
                }
            }
        }
        return logPopulationSizeInverse + integralBit;
    }


    public static double getReciprocalPopSizeInInterval(double time, OldSmoothLineageCount lineageCount,
                                                        SmoothSkygridPopulationSizeInverse populationSizeInverse) {
        return 0.5 * lineageCount.getLineageCount(time) * (lineageCount.getLineageCount(time) - 1) * populationSizeInverse.getPopulationSizeInverse(time);
    }

    public OldSmoothLineageCount getLineageCount() {
        return lineageCount;
    }

    public SmoothSkygridPopulationSizeInverse getPopulationSizeInverse() {
        return populationSizeInverse;
    }

    @Override
    public int getNumberOfCoalescentEvents() {
        int sum = 0;
        for (Tree tree : trees) {
            sum += tree.getInternalNodeCount();
        }
        return sum;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        throw new RuntimeException("Not yet implemented.");
    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Differentiable skygrid coalescent";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(CommonCitations.GILL_2013_IMPROVING,
                new Citation(
                        new Author[] {
                                new Author( "Y", "Bao"),
                                new Author("MA", "Suchard"),
                                new Author( "X", "Ji"),
                        },
                        Citation.Status.IN_PREPARATION
                )
        );
    }

}
