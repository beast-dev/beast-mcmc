/*
 * SmoothSkygridLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.coalescent.smooth;

import dr.evolution.coalescent.IntervalList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.Reportable;

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

        this.tmpA = new double[trees.get(0).getNodeCount()];
        this.tmpB = new double[trees.get(0).getNodeCount()];
        this.tmpC = new double[trees.get(0).getNodeCount()];
        this.tmpADerivOverS = new double[trees.get(0).getNodeCount()];
        this.tmpBDerivOverS = new double[trees.get(0).getNodeCount()];
        this.tmpCDerivOverS = new double[trees.get(0).getNodeCount()];
        this.tmpD = new double[gridPointParameter.getDimension()];
        this.tmpE = new double[gridPointParameter.getDimension()];
        this.tmpF = new double[gridPointParameter.getDimension()];
        this.tmpLineageEffect = new double[trees.get(0).getNodeCount()];
        this.tmpTimes = new double[trees.get(0).getNodeCount()];
        this.tmpCounts = new int[trees.get(0).getNodeCount()];
        this.tmpSumsKnown = false;

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

    private double[] tmpA;
    private double[] tmpADerivOverS;
    private double[] tmpB;
    private double[] tmpBDerivOverS;
    private double[] tmpC;
    private double[] tmpCDerivOverS;
    private double[] tmpD;
    private double[] tmpE;
    private double[] tmpF;
    private double[] tmpLineageEffect;
    private double[] tmpTimes;
    private int[] tmpCounts;
    private int uniqueTimes;
    private boolean tmpSumsKnown;

    private void calculateTmpSums() {
        if (!tmpSumsKnown) {
            TreeModel tree = trees.get(0);
            final double startTime = 0;
            final double endTime = tree.getNodeHeight(tree.getRoot());
            final int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

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
            uniqueTimes = index + 1;

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
            tmpSumsKnown = true;
        }
    }

    private void calculateTmpSumDerivatives() {
        if (!tmpSumsKnown) {
            calculateTmpSums();
        }

        TreeModel tree = trees.get(0);
        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());
        final int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

        for (int i = 0; i < uniqueTimes; i++) {
            final double timeI = tmpTimes[i];
            double sum = 0;
            for (int j = 0; j < uniqueTimes; j++) {
                if (j != i) {
                    final double timeJ = tmpTimes[j];
                    final double lineageCountEffect = tmpLineageEffect[j];
                    final double thisInverse = smoothFunction.getInverseOneMinusExponential(timeJ - timeI, smoothRate.getParameterValue(0));
                    sum += lineageCountEffect * thisInverse * (1 - thisInverse);
                }
            }
            tmpADerivOverS[i] = - sum;
        }

        for (int i = 0; i < uniqueTimes; i++) {
            final double timeI = tmpTimes[i];
            double sum = 0;
            for (int k = 0; k < maxGridIndex; k++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                final double gridTime = gridPointParameter.getParameterValue(k);
                final double thisInverse = smoothFunction.getInverseOneMinusExponential(gridTime - timeI, smoothRate.getParameterValue(0));
                sum += (nextPopSizeInverse - currentPopSizeInverse) * thisInverse * (1 - thisInverse);
            }
            tmpBDerivOverS[i] = -sum;
        }

        for (int i = 0; i < uniqueTimes; i++) {
            final double timeI = tmpTimes[i];
            tmpCDerivOverS[i] = smoothFunction.getSingleIntegrationDerivative(startTime, endTime, timeI, smoothRate.getParameterValue(0));
        }
    }

    protected double calculateLogLikelihood() {
        assert(trees.size() == 1);
        if (!likelihoodKnown) {
            TreeModel tree = trees.get(0);
            final double startTime = 0;
            final double endTime = tree.getNodeHeight(tree.getRoot());
            final int maxGridIndex = getMaxGridIndex(gridPointParameter, endTime);

            calculateTmpSums();

            double lineageEffectSquaredSum = 0;
            for (int i = 0; i < uniqueTimes; i++) {
                lineageEffectSquaredSum += tmpLineageEffect[i] * tmpLineageEffect[i];
            }

            double tripleIntegrationSum = getTripleIntegration(startTime, endTime, maxGridIndex, lineageEffectSquaredSum);

            double doubleIntegrationSum = getDoubleIntegration(startTime, endTime, maxGridIndex, lineageEffectSquaredSum);

            final double singleIntegration = getSingleIntegration(startTime, endTime);

            double logPopulationSizeInverse = 0;
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getNode(tree.getExternalNodeCount() + i);
                logPopulationSizeInverse += Math.log(getSmoothPopulationSizeInverse(tree.getNodeHeight(node), tree.getNodeHeight(tree.getRoot())));
            }

            logLikelihood = logPopulationSizeInverse + singleIntegration + doubleIntegrationSum + tripleIntegrationSum;

            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    private double[] getGradientWrtNodeHeightNew() {
        if (!likelihoodKnown) {
            calculateLogLikelihood();
        }
        TreeModel tree = trees.get(0);
        final double startTime = 0;
        final double endTime = tree.getNodeHeight(tree.getRoot());
        double[] gradient = new double[tree.getInternalNodeCount()];
        getGradientWrtNodeHeightFromSingleIntegration(startTime, endTime, gradient);

        double lineageEffectSquaredSum = 0;
        for (int i = 0; i < uniqueTimes; i++) {
            lineageEffectSquaredSum += tmpLineageEffect[i] * tmpLineageEffect[i];
        }
        getGradientWrtNodeHeightFromDoubleIntegration(startTime, endTime, getMaxGridIndex(gridPointParameter, endTime), gradient);

        getGradientWrtNodeHeightFromTripleIntegration(startTime, endTime, getMaxGridIndex(gridPointParameter, endTime), gradient);
        return gradient;
    }

    double getTripleIntegration(double startTime, double endTime, int maxGridIndex, double lineageEffectSquaredSum) {
        double tripleIntegrationSum = 0;
        for (int i = 0; i < uniqueTimes; i++) {
            final double lineageCountEffect = tmpLineageEffect[i];
            tripleIntegrationSum += lineageCountEffect * tmpA[i] * tmpB[i] * tmpC[i];
        }
        tripleIntegrationSum *= 2;

        for (int k = 0; k < maxGridIndex; k++) {
            final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
            final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
            tripleIntegrationSum += (nextPopSizeInverse - currentPopSizeInverse) * tmpF[k] * tmpD[k];
        }

        tripleIntegrationSum /= -smoothRate.getParameterValue(0) * 2;
        tripleIntegrationSum += -0.5 * (1 - lineageEffectSquaredSum)
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
        return tripleIntegrationSum + tripleWithQuadraticIntegrationSum;
    }

    private void getGradientWrtNodeHeightFromTripleIntegration(double startTime, double endTime, int maxGridIndex,
                                                               double[] gradient) {
        for (int i = 0; i < uniqueTimes; i++) {
            final double lineageCountEffect = tmpLineageEffect[i];
            final double timeI = tmpTimes[i];
            gradient[i] += lineageCountEffect * (tmpADerivOverS[i] * tmpB[i] * tmpC[i] + tmpA[i] * tmpBDerivOverS[i] * tmpC[i] + tmpA[i] * tmpB[i] * tmpCDerivOverS[i]);
            for (int k = 0; k < maxGridIndex; k++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                final double gridTime = gridPointParameter.getParameterValue(k);
                final double tmpEInverse = smoothFunction.getInverseOneMinusExponential(timeI - gridTime, smoothRate.getParameterValue(0));

                gradient[i] += (nextPopSizeInverse - currentPopSizeInverse) * tmpD[k] * (tmpE[k] - lineageCountEffect * tmpEInverse ) * tmpEInverse * (1 - tmpEInverse) * lineageCountEffect;
            }


            final double startTimeInverse = smoothFunction.getInverseOnePlusExponential(timeI - startTime, smoothRate.getParameterValue(0));
            final double endTimeInverse = smoothFunction.getInverseOnePlusExponential(timeI - startTime, smoothRate.getParameterValue(0));
            final double commonSecondTermMultiplier = startTimeInverse - endTimeInverse;
            final double commonSecondTermMultiplierDerivativeOverS = - startTimeInverse * (1 - startTimeInverse) + endTimeInverse * (1 - endTimeInverse);

            for (int k = 0; k < maxGridIndex; k++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                final double gridTime = gridPointParameter.getParameterValue(k);
                final double inverse = smoothFunction.getInverseOneMinusExponential(gridTime - timeI, smoothRate.getParameterValue(0));
                final double inverseDerivativeOverS = -inverse * (1 - inverse);
                gradient[i] += (nextPopSizeInverse - currentPopSizeInverse)
                        * (inverseDerivativeOverS * commonSecondTermMultiplier + inverse * commonSecondTermMultiplierDerivativeOverS +
                        2 * (1 - inverse) * inverseDerivativeOverS * tmpC[i] + (2.0 - inverse) * inverse * tmpCDerivOverS[i] +
                        2 * (1 - inverse) * (-inverseDerivativeOverS) * tmpD[k]);
            }
        }

    }

    double getDoubleIntegration(double startTime, double endTime, int maxGridIndex, double lineageEffectSquaredSum) {
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
        firstDoubleIntegrationOffDiagonalSum += 0.5 * (1 - lineageEffectSquaredSum) * (endTime - startTime);

        final double firstDoubleIntegrationSum = -(firstDoubleIntegrationDiagonalSum * 0.5 + firstDoubleIntegrationOffDiagonalSum) * Math.exp(-logPopSizeParameter.getParameterValue(0));

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

        return firstDoubleIntegrationSum + secondDoubleIntegrationSum;

    }

    private void getGradientWrtNodeHeightFromDoubleIntegration(double startTime, double endTime, int maxGridIndex,
                                                               double[] gradient) {
        final double firstPopSize = Math.exp(-logPopSizeParameter.getParameterValue(0));
        for (int i = 0; i < uniqueTimes; i++) {
            final double lineageCountEffect = tmpLineageEffect[i];
            final double timeI = tmpTimes[i];
            //firstDoubleIntegrationOffDiagonalSum += lineageCountEffect * tmpA[i] * tmpC[i];
            gradient[i] += -lineageCountEffect * (tmpA[i] * tmpCDerivOverS[i] + tmpADerivOverS[i] * tmpC[i]) * firstPopSize;

            //firstDoubleIntegrationDiagonalSum
            gradient[i] += lineageCountEffect * lineageCountEffect
                    * (smoothFunction.getSingleIntegrationDerivative(startTime, endTime, timeI, smoothRate.getParameterValue(0))
                    + (smoothFunction.getDerivative(timeI, endTime, 0, 1, smoothRate.getParameterValue(0))
                    - smoothFunction.getDerivative(timeI, startTime, 0, 1, smoothRate.getParameterValue(0)) / smoothRate.getParameterValue(0))
                    ) * -0.5 * firstPopSize;

            gradient[i] += 0.5 * tmpLineageEffect[i] * (tmpB[i] * tmpCDerivOverS[i] + tmpBDerivOverS[i] * tmpC[i]);

            for (int k = 0; k < maxGridIndex; k++) {
                final double currentPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k));
                final double nextPopSizeInverse = Math.exp(-logPopSizeParameter.getParameterValue(k + 1));
                final double gridTime = gridPointParameter.getParameterValue(k);
                final double tmpEInverse = smoothFunction.getInverseOneMinusExponential(timeI - gridTime, smoothRate.getParameterValue(0));
                gradient[i] += 0.5 * tmpD[k] * (nextPopSizeInverse - currentPopSizeInverse) * tmpEInverse * (1 - tmpEInverse) * lineageCountEffect;
            }
        }
    }

    private double getSingleIntegration(double startTime, double endTime) {
        double singleIntegration = 0;
        for (int i = 0; i < uniqueTimes; i++) {
            final double timeI = tmpTimes[i];
            final double lineageCountEffectI = tmpLineageEffect[i];
            singleIntegration += lineageCountEffectI * smoothFunction.getSingleIntegration(startTime, endTime, timeI, smoothRate.getParameterValue(0));
        }
        singleIntegration *= 0.5 * Math.exp(-logPopSizeParameter.getParameterValue(0));
        return singleIntegration;
    }

    private void getGradientWrtNodeHeightFromSingleIntegration(double startTime, double endTime, double[] gradient) {
        for (int i = 0; i < uniqueTimes; i++) {
            final double timeI = tmpTimes[i];
            final double lineageCountEffectI = tmpLineageEffect[i];
            gradient[i] += lineageCountEffectI * smoothFunction.getSingleIntegrationDerivative(startTime, endTime, timeI, smoothRate.getParameterValue(0))
            * 0.5 * Math.exp(-logPopSizeParameter.getParameterValue(0));
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        super.handleModelChangedEvent(model, object, index);
        tmpSumsKnown = false;
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
