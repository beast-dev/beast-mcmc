/*
 * NewSmoothSkygridLikelihood.java
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

import dr.evolution.tree.NodeRef;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A likelihood function for a smooth skygrid coalescent process that nicely works with the newer tree intervals
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class NewSmoothSkygridLikelihood extends AbstractCoalescentLikelihood implements Reportable {

    private final List<TreeModel> trees;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;
    private final Parameter smoothRate;
    private final List<BigFastTreeIntervals> intervals;
    private final double magicSmallThreshold = 1E-4; // XJ: about 1-hour difference if time in years
    private boolean cacheKnown = false;

    public NewSmoothSkygridLikelihood(String name,
                                      List<TreeModel> trees,
                                      Parameter logPopSizeParameter,
                                      Parameter gridPointParameter,
                                      Parameter smoothRate) {
        super(name);
        this.trees = trees;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;
        this.smoothRate = smoothRate;
        this.intervals = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            intervals.add(new BigFastTreeIntervals(trees.get(i)));
            addModel(intervals.get(i));
            addModel(trees.get(i));
        }

        addVariable(logPopSizeParameter);
        addVariable(gridPointParameter);
        addVariable(smoothRate);

        if (trees.size() > 1) {
            throw new IllegalArgumentException("Not so many trees yet");
        }

        this.tmpA = new double[trees.get(0).getNodeCount()];
        this.tmpB = new double[trees.get(0).getNodeCount()];
        this.tmpC = new double[trees.get(0).getNodeCount()];
        this.tmpD = new double[gridPointParameter.getDimension()];
        this.tmpE = new double[gridPointParameter.getDimension()];
        this.tmpF = new double[gridPointParameter.getDimension()];
        this.uniqueNodeTimes = new double[trees.get(0).getNodeCount()];
        this.sumLineageEffects = new double[trees.get(0).getNodeCount()];
        this.uniqueTimeIndexForGrid = new int[gridPointParameter.getDimension()];
        this.gridIndexUniqueTime = new int[trees.get(0).getNodeCount()];
    }

    private double getSmoothRate() {
        return smoothRate.getParameterValue(0);
    }

    private double getDoubleSigmoidIntegral(double ti, double tj, double t) {
        final double s = getSmoothRate();
        if (Math.abs(ti - tj) < magicSmallThreshold) {
            return t + (GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti - t, s) + GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(t - ti, s)) / s;
        } else {
            return t + (GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti - t, s) * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(tj - ti, s)
                    + GlobalSigmoidSmoothFunction.getLogOnePlusExponential(tj - t, s) * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - tj, s)) / s;
        }
    }

    private double getCompleteDoubleSigmoidIntegral(double ti, double tj, double t) {
        return getDoubleSigmoidIntegral(ti, tj, t) - getDoubleSigmoidIntegral(ti, tj, 0);
    }


    private double[] tmpA;
    private double[] tmpB;
    private double[] tmpC;
    private double[] tmpD;
    private double[] tmpE;
    private double[] tmpF;
    private double firstDoubleIntegralExtra;
    private double secondDoubleIntegralExtra;

    private void cacheTmps() {
        if (!cacheKnown) {
            sortNodeTimes();
            TreeModel tree = trees.get(0);
            final double rootTime = tree.getNodeHeight(tree.getRoot());
            final double s = getSmoothRate();

            firstDoubleIntegralExtra = 0;
            secondDoubleIntegralExtra = 0;
            for (int i = 0; i < numberUniqueNodeTimes; i++) {
                final double ti = uniqueNodeTimes[i];
                final double gi = sumLineageEffects[i];
                double sum = 0;
                for (int j = 0; j < numberUniqueNodeTimes; j++) {
                    if (j != i) {
                        final double gj = sumLineageEffects[j];
                        final double tj = uniqueNodeTimes[j];
                        sum += 2 * gj * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(tj - ti, s);
                    }
                }
                tmpA[i] = sum;
                tmpC[i] = gi * (GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti - rootTime, s) - GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti, s));

                sum = 0;
                for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                    final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                    final double xk = gridPointParameter.getParameterValue(k);
                    if (k != gridIndexUniqueTime[i]) {
                        sum += popSizeInverseDifference * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(xk - ti, s);
                    }
                }
                tmpB[i] = sum;
            }

            for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                final double xk = gridPointParameter.getParameterValue(k);
                final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                tmpD[k] = popSizeInverseDifference * (GlobalSigmoidSmoothFunction.getLogOnePlusExponential(xk - rootTime, s) - GlobalSigmoidSmoothFunction.getLogOnePlusExponential(xk, s));
                double sum = 0;
                double squaredSum = 0;
                for (int i = 0; i < numberUniqueNodeTimes; i++) {
                    final double ti = uniqueNodeTimes[i];
                    final double gi = sumLineageEffects[i];
                    if (i != uniqueTimeIndexForGrid[k]) {
                        final double thisInverse = gi * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - xk, s);
                        sum += thisInverse;
                        squaredSum += thisInverse * thisInverse;
                    }
                }
                tmpE[k] = sum;
                tmpF[k] = squaredSum;
            }

            cacheKnown = true;
        }
    }

    private double getFirstDoubleIntegral() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        final double s = getSmoothRate();
        double sum = 0;
        for (int i = 0; i < numberUniqueNodeTimes; i++) {
            sum += tmpA[i] * tmpC[i];
        }
        sum /= s;
        final double result = sum + rootTime * (2 - 2 * tree.getExternalNodeCount());
        return result;
    }

    private double getFirstDoubleIntegralBruteForce() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double gi = getLineageEffect(i);
            final double ti = tree.getNodeHeight(tree.getNode(i));
            for (int j = 0; j < tree.getNodeCount(); j++) {
                final double gj = getLineageEffect(j);
                final double tj = tree.getNodeHeight(tree.getNode(j));
                if (i != j) {
                    sum += gi * gj * getCompleteDoubleSigmoidIntegral(ti, tj, rootTime);
                }
            }
        }
        return sum;
    }

    private double getFirstDoubleIntegralApproximate() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double gi = getLineageEffect(i);
            final double ti = tree.getNodeHeight(tree.getNode(i));
            for (int j = 0; j < tree.getNodeCount(); j++) {
                final double gj = getLineageEffect(j);
                final double tj = tree.getNodeHeight(tree.getNode(j));
                if (i != j) {
                    sum += gi * gj * (rootTime - ((ti > tj) ? ti : tj));
                }
            }
        }
        return sum;
    }

    private double getSecondDoubleIntegral() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        final double s = getSmoothRate();
        double sum = 0;
        for (int i = 0; i < numberUniqueNodeTimes; i++) {
            final double gi = sumLineageEffects[i];
            final double ti = uniqueNodeTimes[i];
            sum += gi * tmpC[i] + gi * gi * (GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(rootTime - ti, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(-ti, s));
        }
        return sum / s + rootTime * tree.getNodeCount();
    }

    private double getSecondDoubleIntegralBruteForce() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double ti = tree.getNodeHeight(tree.getNode(i));
            sum += getCompleteDoubleSigmoidIntegral(ti, ti, rootTime);
        }
        return sum;
    }

    private double getSecondDoubleIntegralApproximate() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            sum += rootTime - tree.getNodeHeight(tree.getNode(i));
        }
        return sum;
    }

    private double getThirdDoubleIntegral() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        final double s = getSmoothRate();
        double sum = 0;
        for (int i = 0; i < numberUniqueNodeTimes; i++) {
            sum += tmpC[i] * tmpB[i];
        }
        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            sum += tmpD[k] * tmpE[k];
            if (uniqueTimeIndexForGrid[k] > -1) {
                final double xk = gridPointParameter.getParameterValue(k);
                final double gi = sumLineageEffects[uniqueTimeIndexForGrid[k]];
                final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                sum += gi * popSizeInverseDifference * (GlobalSigmoidSmoothFunction.getLogOnePlusExponential(xk - rootTime, s) - GlobalSigmoidSmoothFunction.getLogOnePlusExponential(xk, s)
                        + GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(rootTime - xk, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(- xk, s));
            }
        }
        sum /= s;
        sum += secondDoubleIntegralExtra + (Math.exp(-logPopSizeParameter.getParameterValue(gridPointParameter.getDimension())) - Math.exp(-logPopSizeParameter.getParameterValue(0))) * rootTime;
        return sum;
    }

    private double getThirdDoubleIntegralBruteForce() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double ti = tree.getNodeHeight(tree.getNode(i));
            final double gi = getLineageEffect(i);
            for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                final double xk = gridPointParameter.getParameterValue(k);
                final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                sum += gi * popSizeInverseDifference * getCompleteDoubleSigmoidIntegral(ti, xk, rootTime);
            }
        }
        return sum;
    }

    private double getThirdDoubleIntegralApproximate() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double ti = tree.getNodeHeight(tree.getNode(i));
            final double gi = getLineageEffect(i);
            for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                final double xk = gridPointParameter.getParameterValue(k);
                final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                sum += gi * popSizeInverseDifference * (rootTime - ((ti > xk) ? ti : xk));
            }
        }
        return sum;
    }

    private double getPopSizeInverseDifference(int k) {
        return Math.exp(-logPopSizeParameter.getParameterValue(k + 1)) - Math.exp(-logPopSizeParameter.getParameterValue(k));
    }

    private double getSingleSigmoidIntegral(double ti, double t) {
        final double s = getSmoothRate();
        final double exponential = Math.exp(s * (t - ti));
        if (Double.isInfinite(exponential)) {
            return t - ti;
        } else {
            return Math.log(1 + exponential) / s;
        }
    }

    private double getAllSingleIntegralApproximate() {
        TreeModel tree = trees.get(0);
        final double rootHeight = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            sum -= getLineageEffect(i) * (rootHeight - tree.getNodeHeight(tree.getNode(i)));
        }
        return sum;
    }

    private double getCompleteSingleSigmoidIntegral(double ti, double t) {
        double result = getSingleSigmoidIntegral(ti, t) - getSingleSigmoidIntegral(ti, 0);
        return result;
    }

    private double getLineageEffect(int node) {
        TreeModel tree = trees.get(0);
        NodeRef nodeRef = tree.getNode(node);
        return tree.isExternal(nodeRef) ? 1 : -1;
    }

    private double getAllSingleIntegrals() {
        final double s = getSmoothRate();
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int i = 0; i < numberUniqueNodeTimes; i++) {
            sum -= sumLineageEffects[i] * getCompleteSingleSigmoidIntegral(uniqueNodeTimes[i], rootTime);
        }
        return sum;
    }

    private double getTripleIntegralFragment(double t0, double t1, double t2, double t) {
        final double s = getSmoothRate();
        final double numerator = GlobalSigmoidSmoothFunction.getLogOnePlusExponential(t0 - t, s);
        final double denominatorInverse = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(t1 - t0, s) *
                GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(t2 - t0, s);
        if ((isLimitingCase(t0, t) || isLimitingCase(t1, t) || isLimitingCase(t2, t)) && t0 == t) {
                return 0;
        } else {
            return numerator * denominatorInverse;
        }
    }

    private double getCubicIntegral(double t0, double t) {
        final double s = getSmoothRate();
        final double exponent = Math.exp(s * (t - t0));
        if (isLimitingCase(t0, t) || isLimitingCase(t, t0)) {
            return t > t0 ? t - t0 : 0;
        } else {
            return ((3 + 4 * exponent)/2/(1+exponent)/(1+exponent) + GlobalSigmoidSmoothFunction.getLogOnePlusExponential(t - t0, s)) / s;
        }
    }

    private boolean isLimitingCase(double t0, double t) {
        return Math.exp(getSmoothRate() * (t0 - t)) == 0;
    }

    private double getTripleIntegralWithQuadratic(double ti, double xk, double t) {
        final double s = getSmoothRate();
        final double oneOverOnePlusTiMinusT = GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti - t, s);
        final double oneOverOneMinusXkMinusTi = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(xk - ti, s);
        final double oneOverOneMinusTiMinusXk = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - xk, s);
        final boolean limitS = isLimitingCase(ti, t) || isLimitingCase(xk, t);
        if (limitS) {
            return t - oneOverOnePlusTiMinusT * oneOverOneMinusXkMinusTi / s
                    + (oneOverOneMinusXkMinusTi + oneOverOneMinusXkMinusTi * oneOverOneMinusTiMinusXk) * (ti < t ? 0 : ti - t)
                    + oneOverOneMinusTiMinusXk * oneOverOneMinusTiMinusXk * (xk < t ? 0 : xk - t);
        }
        return t + (-oneOverOnePlusTiMinusT * oneOverOneMinusXkMinusTi
                + (oneOverOneMinusXkMinusTi + oneOverOneMinusXkMinusTi * oneOverOneMinusTiMinusXk) * GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti - t, s)
                + oneOverOneMinusTiMinusXk * oneOverOneMinusTiMinusXk * GlobalSigmoidSmoothFunction.getLogOnePlusExponential(xk - t, s)) / s;
    }

    private double getTripleIntegral(double ti, double tj, double xk, double t) {
        final double s = getSmoothRate();

        if (Math.abs(ti - tj) < magicSmallThreshold && Math.abs(ti - xk) < magicSmallThreshold && Math.abs(tj - xk) < magicSmallThreshold) {
            return getCubicIntegral(ti, t);
        } else if (Math.abs(ti - tj) < magicSmallThreshold) {
            return getTripleIntegralWithQuadratic(ti, xk, t);
        } else if (Math.abs(ti - xk) < magicSmallThreshold) {
            return getTripleIntegralWithQuadratic(ti, tj, t);
        } else if (Math.abs(tj - xk) < magicSmallThreshold) {
            return getTripleIntegralWithQuadratic(tj, ti, t);
        } else {
            return (t + (getTripleIntegralFragment(ti, tj, xk, t) + getTripleIntegralFragment(tj, xk, ti, t) + getTripleIntegralFragment(xk, ti, tj, t)) / s);
        }
    }

    private double getCompleteTripleIntegral(double ti, double tj, double xk, double t) {
        return getTripleIntegral(ti, tj, xk, t) - getTripleIntegral(ti, tj, xk, 0);
    }

    private int numberUniqueNodeTimes;
    private double[] uniqueNodeTimes;
    private double[] sumLineageEffects;
    private int[] uniqueTimeIndexForGrid;
    private int[] gridIndexUniqueTime;

    private void sortNodeTimes() {
        TreeModel tree = trees.get(0);

        NodeRef[] nodes = new NodeRef[tree.getNodeCount()];
        System.arraycopy(tree.getNodes(), 0, nodes, 0, tree.getNodeCount());
        Arrays.parallelSort(nodes, (a, b) -> Double.compare(tree.getNodeHeight(a), tree.getNodeHeight(b)));

        double lastTime = tree.getNodeHeight(nodes[0]);
        int index = 0;
        uniqueNodeTimes[index] = lastTime;
        sumLineageEffects[index] = getLineageEffect(nodes[0].getNumber());
        for (int i = 1; i < nodes.length; i++) {
            NodeRef node = nodes[i];
            final double time = tree.getNodeHeight(node);
            if (Math.abs(time - lastTime) < magicSmallThreshold) {
                sumLineageEffects[index] += getLineageEffect(node.getNumber());
            } else {
                index++;
                uniqueNodeTimes[index] = time;
                sumLineageEffects[index] = getLineageEffect(node.getNumber());
                lastTime = time;
            }
        }
        numberUniqueNodeTimes = index + 1;

        Arrays.fill(uniqueTimeIndexForGrid, -1);
        Arrays.fill(gridIndexUniqueTime, -1);
        index = 0;
        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            while(uniqueNodeTimes[index] < gridPointParameter.getParameterValue(k) && index < numberUniqueNodeTimes) {
                index++;
            }
            if (Math.abs(uniqueNodeTimes[index] - gridPointParameter.getParameterValue(k)) < magicSmallThreshold) {
                uniqueTimeIndexForGrid[k] = index;
                gridIndexUniqueTime[index] = k;
            }
        }
    }


    private double getTripleIntegral() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());

        final double firstTripleIntegral = getFirstTripleIntegral();
        final double tripleWithSquareIntegral = getTripleIntegralWithSquares();
        final double allTs = rootTime * (Math.exp(-logPopSizeParameter.getParameterValue(gridPointParameter.getDimension())) - Math.exp(-logPopSizeParameter.getParameterValue(0)));

        return firstTripleIntegral + tripleWithSquareIntegral + allTs;

    }


    private double getFirstTripleIntegral() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        final double s = getSmoothRate();
        double sum = 0;
        for (int i = 0; i < numberUniqueNodeTimes; i++) {
            sum += tmpA[i] * tmpB[i] * tmpC[i];
        }

        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            sum += (tmpE[k] * tmpE[k] - tmpF[k]) * tmpD[k];
        }

        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            final double popSizeInverseDifference = getPopSizeInverseDifference(k);
            double thisExtra = 0;
            if (uniqueTimeIndexForGrid[k] > -1) {
                final double ti = uniqueNodeTimes[uniqueTimeIndexForGrid[k]];
                final double gi = sumLineageEffects[uniqueTimeIndexForGrid[k]];
                final double inverseOnePlusExpTiMinusT = GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti - rootTime, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti, s);
                thisExtra -= gi * inverseOnePlusExpTiMinusT * 0.5 * tmpA[uniqueTimeIndexForGrid[k]];

                double firstSumOverJ = 0;
                double secondSumOverJ = 0;
                for (int j = 0; j < numberUniqueNodeTimes; j++) {
                    if (j != uniqueTimeIndexForGrid[k]) {
                        final double gj = sumLineageEffects[j];
                        final double tj = uniqueNodeTimes[j];
                        final double inverseOneMinusExpTiMinusTj = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - tj, s);
                        final double inverseOneMinusExpTjMinusTi = 1 - inverseOneMinusExpTiMinusTj;
                        firstSumOverJ += inverseOneMinusExpTiMinusTj * inverseOneMinusExpTiMinusTj * tmpC[j];
                        secondSumOverJ += gj * (inverseOneMinusExpTjMinusTi + inverseOneMinusExpTjMinusTi * inverseOneMinusExpTiMinusTj);
                    }
                }
                thisExtra += gi * firstSumOverJ + secondSumOverJ * tmpC[uniqueTimeIndexForGrid[k]];
                sum += thisExtra * popSizeInverseDifference;
            }
        }
        sum /= s;


        return sum;
    }

    private double getTripleIntegralWithSquares() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        final double s = getSmoothRate();
        double sum = 0;
        for (int i = 0; i < numberUniqueNodeTimes; i++) {
            final double ti = uniqueNodeTimes[i];
            final double gi = sumLineageEffects[i];
            final double inverseOnePlusExpTiMinusT = GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti - rootTime, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti, s);
            sum -= gi * gi * inverseOnePlusExpTiMinusT * tmpB[i];
            double sumOverK = 0;

            for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                if (k != gridIndexUniqueTime[i]) {
                    final double xk = gridPointParameter.getParameterValue(k);
                    final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                    final double inverseOneMinusExpTiMinusXk = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - xk, s);
                    final double inverseOneMinusExpXkMinusTi = 1 - inverseOneMinusExpTiMinusXk;

                    sumOverK += popSizeInverseDifference * (inverseOneMinusExpXkMinusTi + inverseOneMinusExpXkMinusTi * inverseOneMinusExpXkMinusTi);
                }
            }
            sum += sumOverK * gi * tmpC[i];
        }

        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            final double xk = gridPointParameter.getParameterValue(k);
            double sumOverI = 0;
            for (int i = 0; i < numberUniqueNodeTimes; i++) {
                final double ti = uniqueNodeTimes[i];
                final double gi = sumLineageEffects[i];
                if (i != uniqueTimeIndexForGrid[k]) {
                    final double inverseOneMinusExpTiMinusXk = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - xk, s);
                    sumOverI += gi * gi * inverseOneMinusExpTiMinusXk * inverseOneMinusExpTiMinusXk;
                }
            }
            sum += sumOverI * tmpD[k];
        }

        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            if (uniqueTimeIndexForGrid[k] > -1) {
                final double ti = gridPointParameter.getParameterValue(k);
                final double gi = sumLineageEffects[uniqueTimeIndexForGrid[k]];
                final double popSizeInverseDifference = getPopSizeInverseDifference(k);

                sum += gi * gi * popSizeInverseDifference * (getCompleteTripleIntegral(ti, ti, ti, rootTime) - rootTime * s);

            }
        }

        sum /= s;

        return sum;
    }

    private double getTripleIntegralBruteForce() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            final double xk = gridPointParameter.getParameterValue(k);
            final double popSizeInverseDifference = getPopSizeInverseDifference(k);
            for (int i = 0; i < tree.getNodeCount(); i++) {
                final double ti = tree.getNodeHeight(tree.getNode(i));
                final double gi = getLineageEffect(i);
                for (int j = 0; j < tree.getNodeCount(); j++) {
                    final double gj = getLineageEffect(j);
                    final double tj = tree.getNodeHeight(tree.getNode(j));
                        final double newAmount = popSizeInverseDifference * gi * gj * getCompleteTripleIntegral(ti, tj, xk, rootTime);
                        sum += newAmount;
                }
            }
        }
        return sum;
    }

    private double getTripleIntegralApproximate() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        double sum = 0;
        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            final double xk = gridPointParameter.getParameterValue(k);
            final double popSizeInverseDifference = getPopSizeInverseDifference(k);
            for (int i = 0; i < tree.getNodeCount(); i++) {
                final double ti = tree.getNodeHeight(tree.getNode(i));
                final double gi = getLineageEffect(i);
                double cutoff = ti > xk ? ti : xk;
                for (int j = 0; j < tree.getNodeCount(); j++) {
                    final double gj = getLineageEffect(j);
                    final double tj = tree.getNodeHeight(tree.getNode(j));
                    final double integralApproximation = rootTime - (tj > cutoff ? tj : cutoff);
                    if (integralApproximation > 0) {
                        sum += popSizeInverseDifference * gi * gj * integralApproximation;
                    }
                }
            }
        }
        return sum;
    }

    @Override
    protected double calculateLogLikelihood() {
        if (!cacheKnown) {
            cacheTmps();
        }
        return 0;
    }

    @Override
    public Type getUnits() {
        return null;
    }

    @Override
    public void setUnits(Type units) {

    }

    @Override
    public int getNumberOfCoalescentEvents() {
        return 0;
    }

    @Override
    public double getCoalescentEventsStatisticValue(int i) {
        return 0;
    }

    final private static boolean DEBUG = true;

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append("smoothSkygrid(" + getLogLikelihood() + ")").append("\n");
        if (DEBUG) {
            sb.append("=========================\n");
            sb.append("smooth rate = " + getSmoothRate() + "\n");
            sb.append("Single integral = " + getAllSingleIntegrals() + "\n");
            sb.append("Single integral approximate (s -> Inf) = " + getAllSingleIntegralApproximate() + "\n");
            sb.append("First double integral = " + getFirstDoubleIntegral() + "\n");
            sb.append("First double integral brute force = " + getFirstDoubleIntegralBruteForce() + "\n");
            sb.append("First double integral approximate = " + getFirstDoubleIntegralApproximate() + "\n");

            sb.append("Second double integral = " + getSecondDoubleIntegral() + "\n");
            sb.append("Second double integral brute force = " + getSecondDoubleIntegralBruteForce() + "\n");
            sb.append("Second double integral approximate = " + getSecondDoubleIntegralApproximate() + "\n");

            sb.append("First + Second double integral = " + (getFirstDoubleIntegral() + getSecondDoubleIntegral()) + "\n");
            sb.append("First + Second double integral brute force = " + (getFirstDoubleIntegralBruteForce() + getSecondDoubleIntegralBruteForce()) + "\n");
            sb.append("First + Second double integral approximate = " + (getFirstDoubleIntegralApproximate() + getSecondDoubleIntegralApproximate()) + "\n");


            sb.append("Third double integral = " + getThirdDoubleIntegral() + "\n");
            sb.append("Third double integral brute force = " + getThirdDoubleIntegralBruteForce() + "\n");
            sb.append("Third double integral approximate = " + getThirdDoubleIntegralApproximate() + "\n");

            sb.append("Triple integral = " + getTripleIntegral() + "\n");
            sb.append("Triple integral brute force = " + getTripleIntegralBruteForce() + "\n");
            sb.append("Triple integral approximate = " + getTripleIntegralApproximate() + "\n");


        }

        return sb.toString();
    }
}
