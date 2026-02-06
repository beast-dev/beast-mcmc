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
                final double singleExtra = (GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti - rootTime, s) - GlobalSigmoidSmoothFunction.getInverseOnePlusExponential(ti, s)) / s;
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
                    final double popSizeInverseDifference = Math.exp(-logPopSizeParameter.getParameterValue(k)) - Math.exp(-logPopSizeParameter.getParameterValue(k+ 1));
                    final double xk = gridPointParameter.getParameterValue(k);
                    if (Math.abs(xk - ti) < magicSmallThreshold) {
                        sum += 0.5 * popSizeInverseDifference;
                        secondDoubleIntegralExtra += secondDoubleExtra;
                    } else {
                        sum += popSizeInverseDifference * GlobalSigmoidSmoothFunction.getInverseOneMinusExponential(xk - ti, s);
                    }
                }
                tmpB[i] = sum;
            }

            for (int k = 0; k < gridPointParameter.getDimension(); k++) {
                final double xk = gridPointParameter.getParameterValue(k);
                final double popSizeInverseDifference = Math.exp(-logPopSizeParameter.getParameterValue(k)) - Math.exp(-logPopSizeParameter.getParameterValue(k+ 1));
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

        }

        return sb.toString();
    }
}
