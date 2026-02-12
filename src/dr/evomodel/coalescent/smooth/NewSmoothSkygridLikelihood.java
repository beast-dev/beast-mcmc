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
            TreeModel tree = trees.get(0);
            final double rootTime = tree.getNodeHeight(tree.getRoot());
            final double s = getSmoothRate();

            firstDoubleIntegralExtra = 0;
            secondDoubleIntegralExtra = 0;
            for (int i = 0; i < tree.getNodeCount(); i++) {
                final double ti = tree.getNodeHeight(tree.getNode(i));
                final double gi = getLineageEffect(i);
                final double singleExtra = (GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(rootTime - ti, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(-ti, s)) / s;
                double sum = 0;
                for (int j = 0; j < tree.getNodeCount(); j++) {
                    if (j != i) {
                        final double gj = getLineageEffect(j);
                        final double tj = tree.getNodeHeight(tree.getNode(j));
                        if (Math.abs(ti - tj) < magicSmallThreshold) {
                            sum += gj;
                            firstDoubleIntegralExtra += gj * singleExtra;
                        } else {
                            sum += 2 * gj * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(tj - ti, s);
                        }
                    }
                }
                tmpA[i] = sum;
                tmpC[i] = gi * (GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti - rootTime, s) - GlobalSigmoidSmoothFunction.getLogOnePlusExponential(ti, s));

                sum = 0;
                final double secondDoubleExtra = (GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(rootTime - ti, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(- ti, s)) / s;
                for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                    final double popSizeInverseDifference = getPopSizeInverseDifference(k);
                    final double xk = gridPointParameter.getParameterValue(k);
                    if (Math.abs(xk - ti) < magicSmallThreshold) {
                        sum += 0.5 * popSizeInverseDifference;
                        secondDoubleIntegralExtra += gi * popSizeInverseDifference * secondDoubleExtra;
                    } else {
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
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    final double ti = tree.getNodeHeight(tree.getNode(i));
                    final double gi = getLineageEffect(i);
                    if (Math.abs(ti - xk) < magicSmallThreshold) {
                        sum += 0.5 * gi;
                    } else {
                        sum += gi * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(ti - xk, s);
                    }
                }
                tmpE[k] = sum;
            }

            cacheKnown = true;
        }
    }

    private double getFirstDoubleIntegral() {
        TreeModel tree = trees.get(0);
        final double rootTime = tree.getNodeHeight(tree.getRoot());
        final double s = getSmoothRate();
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double gi = getLineageEffect(i);
            sum += tmpA[i] * tmpC[i];
        }
        sum /= s;
        final double result = sum + firstDoubleIntegralExtra + rootTime * (2 - 2 * tree.getExternalNodeCount());
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
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double gi = getLineageEffect(i);
            final double ti = tree.getNodeHeight(tree.getNode(i));
            sum += gi * tmpC[i] + (GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(rootTime - ti, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(-ti, s));
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
        for (int i = 0; i < tree.getNodeCount(); i++) {
            sum += tmpC[i] * tmpB[i];
        }
        for (int k = 0; k < gridPointParameter.getDimension(); k++) {
            sum += tmpD[k] * tmpE[k];
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
        for (int i = 0; i < tree.getNodeCount(); i++) {
            final double lineageEffect = getLineageEffect(i);
            sum -= lineageEffect * getCompleteSingleSigmoidIntegral(tree.getNodeHeight(tree.getNode(i)), rootTime);
        }
        return sum;
    }

    private double getTripleIntegralFragment(double t0, double t1, double t2, double t) {
        final double s = getSmoothRate();
        final double numerator = GlobalSigmoidSmoothFunction.getLogOnePlusExponential(t0 - t, s);
        final double denominatorInverse = GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(t1 - t0, s) *
                GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(t2 - t0, s);
        if ((isLimitingCase(t0, t) || isLimitingCase(t1, t) || isLimitingCase(t2, t)) && t0 == t) {
//            if (t0 >= t)
//                return (t0 - t) * s * denominatorInverse;
//            else
                return 0;
        } else {
            return numerator * denominatorInverse;
        }
    }

    private double getCubicIntegral(double t0, double t) {
        final double s = getSmoothRate();
        final double exponent = Math.exp(s * (t - t0));
        if (isLimitingCase(t0, t)) {
            return t0 - t;
        } else {
            return ((3 + 4 * exponent)/2/(1+exponent)/(1+exponent) + GlobalSigmoidSmoothFunction.getLogOnePlusExponential(t0 - t, s)) / s;
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
            return t + (-oneOverOnePlusTiMinusT * oneOverOneMinusXkMinusTi
                    + (oneOverOneMinusXkMinusTi + oneOverOneMinusXkMinusTi * oneOverOneMinusTiMinusXk) * (ti < t ? 0 : ti - t) * s
                    + oneOverOneMinusTiMinusXk * oneOverOneMinusTiMinusXk * (xk < t ? 0 : xk - t) * s) / s;
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

    private double getFirstTripleIntegralBruteForce() {
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
                    if (i != j) {
                        final double newAmount = popSizeInverseDifference * gi * gj * getCompleteTripleIntegral(ti, tj, xk, rootTime);
                        sum += newAmount;
                    }
                }
            }
        }
        return sum;
    }

    private double getFirstTripleIntegralApproximate() {
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
                    if (i != j) {
                        final double integralApproximation = (rootTime - (tj > cutoff ? tj : cutoff)) < 0 ? 0 : (rootTime - (tj > cutoff ? tj : cutoff));
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

            sb.append("Third double integral = " + getThirdDoubleIntegral() + "\n");
            sb.append("Third double integral brute force = " + getThirdDoubleIntegralBruteForce() + "\n");
            sb.append("Third double integral approximate = " + getThirdDoubleIntegralApproximate() + "\n");

            sb.append("First triple integral brute force = " + getFirstTripleIntegralBruteForce() + "\n");
            sb.append("First triple integral approximate = " + getFirstTripleIntegralApproximate() + "\n");


        }

        return sb.toString();
    }
}
