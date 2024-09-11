/*
 * SkyGlideLikelihood.java
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

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MachineAccuracy;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * A likelihood function for a piece-wise linear log population size coalescent process that nicely works with the newer tree intervals
 *
 * @author Mathieu Fourment
 * @author Erick Matsen
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SkyGlideLikelihood extends AbstractModelLikelihood implements Reportable {

    private final List<TreeModel> trees;

    private final List<BigFastTreeIntervals> intervals;
    private final Parameter logPopSizeParameter;
    private final Parameter gridPointParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;

    public SkyGlideLikelihood(String name,
                              List<TreeModel> trees,
                              Parameter logPopSizeParameter,
                              Parameter gridPointParameter) {
        super(name);
        this.trees = trees;
        this.logPopSizeParameter = logPopSizeParameter;
        this.gridPointParameter = gridPointParameter;
        this.intervals = new ArrayList<>();
        for (int i = 0; i < trees.size(); i++) {
            BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals(trees.get(i));
            this.intervals.add(treeIntervals);
            addModel(treeIntervals);
        }
        addVariable(logPopSizeParameter);
    }

    public List<TreeModel> getTrees() {
        return trees;
    }

    public BigFastTreeIntervals getIntervals(int treeIndex) {
        return intervals.get(treeIndex);
    }

    public TreeModel getTree(int treeIndex) {
        return trees.get(treeIndex);
    }

    @Override
    public String getReport() {
        return "skyGlideLikelihood(" + getLogLikelihood() + ")";
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {
        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            double lnL = 0;
            for (int i = 0; i < trees.size(); i++) {
                lnL += getSingleTreeLogLikelihood(i);
            }
            logLikelihood = lnL;
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public Parameter getLogPopSizeParameter() {
        return logPopSizeParameter;
    }

    public double[] getGradientWrtLogPopulationSize() {
        double[] gradient = new double[logPopSizeParameter.getDimension()];

        for (int index = 0; index < trees.size(); index++) {
            BigFastTreeIntervals interval = intervals.get(index);
            Tree thisTree = trees.get(index);
            int currentGridIndex = 0;
            for (int i = 0; i < interval.getIntervalCount(); i++) {
                final int lineageCount = interval.getLineageCount(i);
                int[] nodeIndices = interval.getNodeNumbersForInterval(i);
                final double intervalStart = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[0]));
                final double intervalEnd = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[1]));

                if (intervalStart != intervalEnd) {
                    int[] gridIndices = getGridPoints(currentGridIndex, intervalStart, intervalEnd);
                    final int firstGridIndex = gridIndices[0];
                    final int lastGridIndex = gridIndices[1];


                    if (firstGridIndex == lastGridIndex) {
                        updateIntervalGradientWrtLogPopSize(intervalStart, intervalEnd, firstGridIndex, lineageCount, gradient);
                    } else {
                        updateIntervalGradientWrtLogPopSize(intervalStart, gridPointParameter.getParameterValue(firstGridIndex), firstGridIndex, lineageCount, gradient);
                        currentGridIndex = firstGridIndex;
                        while(currentGridIndex + 1 < lastGridIndex) {
                            updateIntervalGradientWrtLogPopSize(gridPointParameter.getParameterValue(currentGridIndex), gridPointParameter.getParameterValue(currentGridIndex + 1), currentGridIndex + 1, lineageCount, gradient);
                            currentGridIndex++;
                        }
                        updateIntervalGradientWrtLogPopSize(gridPointParameter.getParameterValue(currentGridIndex), intervalEnd, currentGridIndex + 1, lineageCount, gradient);
                    }
                    currentGridIndex = lastGridIndex;
                }
            }
            updateSingleTreePopulationInverseGradientWrtLogPopSize(index, gradient);
        }

        return gradient;
    }

    public double[] getDiagonalHessianLogDensityWrtLogPopSize() {
        double[] diagonalHessian = new double[logPopSizeParameter.getDimension()];

        for (int index = 0; index < trees.size(); index++) {
            BigFastTreeIntervals interval = intervals.get(index);
            Tree thisTree = trees.get(index);
            int currentGridIndex = 0;
            for (int i = 0; i < interval.getIntervalCount(); i++) {
                final int lineageCount = interval.getLineageCount(i);
                int[] nodeIndices = interval.getNodeNumbersForInterval(i);
                final double intervalStart = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[0]));
                final double intervalEnd = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[1]));

                if (intervalStart != intervalEnd) {
                    int[] gridIndices = getGridPoints(currentGridIndex, intervalStart, intervalEnd);
                    final int firstGridIndex = gridIndices[0];
                    final int lastGridIndex = gridIndices[1];


                    if (firstGridIndex == lastGridIndex) {
                        updateIntervalDiagonalHessianWrtLogPopSize(intervalStart, intervalEnd, firstGridIndex, lineageCount, diagonalHessian);
                    } else {
                        updateIntervalDiagonalHessianWrtLogPopSize(intervalStart, gridPointParameter.getParameterValue(firstGridIndex), firstGridIndex, lineageCount, diagonalHessian);
                        currentGridIndex = firstGridIndex;
                        while(currentGridIndex + 1 < lastGridIndex) {
                            updateIntervalDiagonalHessianWrtLogPopSize(gridPointParameter.getParameterValue(currentGridIndex), gridPointParameter.getParameterValue(currentGridIndex + 1), currentGridIndex + 1, lineageCount, diagonalHessian);
                            currentGridIndex++;
                        }
                        updateIntervalDiagonalHessianWrtLogPopSize(gridPointParameter.getParameterValue(currentGridIndex), intervalEnd, currentGridIndex + 1, lineageCount, diagonalHessian);
                    }
                    currentGridIndex = lastGridIndex;
                }
            }
        }

        return diagonalHessian;
    }

    private void updateIntervalDiagonalHessianWrtLogPopSize(double intervalStart, double intervalEnd, int gridIndex,
                                                            int lineageCount, double[] diagonalHessian) {
        final double slope = getGridSlope(gridIndex);
        final double intercept = getGridIntercept(gridIndex);
        final double lineageMultiplier = -0.5 * lineageCount * (lineageCount - 1);
        assert(slope != 0 || intercept != 0);
        final double realSmall = getMagicUnderFlowBound(slope);
        if (intervalStart != intervalEnd) {
            final double expIntervalStart = Math.exp(-slope * intervalStart);
            final double expIntervalEnd = Math.exp(-slope * intervalEnd);
            final double expIntercept = Math.exp(-intercept);

            final double thisGridTime = gridIndex < gridPointParameter.getDimension() ? gridPointParameter.getParameterValue(gridIndex) : 0;
            final double lastGridTime = gridIndex == 0 ? 0 : gridPointParameter.getParameterValue(gridIndex - 1);

            final double secondDerivativeWrtIntercept = getLinearInverseIntegral(intervalStart, intervalEnd, gridIndex);
            final double secondDerivativeWrtSlope = Math.abs(slope) < realSmall ? expIntercept
                    * (intervalEnd * intervalEnd * intervalEnd - intervalStart * intervalStart * intervalStart) / 3 :
                    expIntercept * (-2 / slope / slope * (intervalEnd * expIntervalEnd - intervalStart * expIntervalStart)
                    +(intervalStart * intervalStart * expIntervalStart - intervalEnd * intervalEnd * expIntervalEnd) / slope
                    +2 / slope / slope / slope * (expIntervalStart - expIntervalEnd));
            final double derivativeWrtSlope = Math.abs(slope) < realSmall ? expIntercept * (intervalStart * intervalStart - intervalEnd * intervalEnd) / 2 :
                    expIntercept * ((intervalEnd * expIntervalEnd - intervalStart * expIntervalStart) / slope - (expIntervalStart - expIntervalEnd) / slope / slope);
            final double secondDerivativeWrtInterceptSlope = -derivativeWrtSlope;

            final double partialInterceptPartialFirstLogPopSize = thisGridTime > 0 ? thisGridTime / (thisGridTime - lastGridTime) : 1;
            final double partialInterceptPartialSecondLogPopSize = thisGridTime > 0 ? - lastGridTime / (thisGridTime - lastGridTime) : 0;
            final double partialSlopePartialFirstLogPopSize = thisGridTime > 0 ? - 1 / (thisGridTime - lastGridTime) : 0;
            final double partialSLopePartialSecondLogPopSize = thisGridTime > 0 ? 1 / (thisGridTime - lastGridTime) : 0;

            diagonalHessian[gridIndex] += lineageMultiplier * (secondDerivativeWrtIntercept * partialInterceptPartialFirstLogPopSize * partialInterceptPartialFirstLogPopSize
                    + 2 * secondDerivativeWrtInterceptSlope * partialSlopePartialFirstLogPopSize * partialInterceptPartialFirstLogPopSize
                    + secondDerivativeWrtSlope * partialSlopePartialFirstLogPopSize * partialSlopePartialFirstLogPopSize);
            if (gridIndex < gridPointParameter.getDimension()) {
                diagonalHessian[gridIndex + 1] += lineageMultiplier * (secondDerivativeWrtIntercept * partialInterceptPartialSecondLogPopSize * partialInterceptPartialSecondLogPopSize
                        + 2 * secondDerivativeWrtInterceptSlope * partialSLopePartialSecondLogPopSize * partialInterceptPartialSecondLogPopSize
                        + secondDerivativeWrtSlope * partialSLopePartialSecondLogPopSize * partialSLopePartialSecondLogPopSize);
            }
        }
    }

    public enum NodeHeightDerivativeType {
        GRADIENT {
            @Override
            double getNodeHeightDerivative(double intercept, double slope, double time, double lineageMultiplier) {
                return lineageMultiplier * Math.exp(-intercept - slope * time);
            }

            @Override
            void updateSingleTreePopulationInverseGradientWrtNodeHeight(SkyGlideLikelihood likelihood, int treeIndex, double[] derivatives) {

                int currentGridIndex = 0;
                BigFastTreeIntervals interval = likelihood.getIntervals(treeIndex);
                TreeModel tree = likelihood.getTree(treeIndex);

                for (int i = 0; i < interval.getIntervalCount(); i++) {
                    if (interval.getIntervalType(i) == IntervalType.COALESCENT) {
                        final double time = interval.getIntervalTime(i + 1);
                        final int nodeIndex = interval.getNodeNumbersForInterval(i)[1];
                        currentGridIndex = likelihood.getGridIndex(time, currentGridIndex);
                        final double slope = likelihood.getGridSlope(currentGridIndex);
                        derivatives[nodeIndex - tree.getExternalNodeCount()] -= slope;
                    }
                }
            }
        },
        DIAGONAL_HESSIAN {
            @Override
            double getNodeHeightDerivative(double intercept, double slope, double time, double lineageMultiplier) {
                return - lineageMultiplier * Math.exp(-intercept - slope * time) * slope;
            }

            @Override
            void updateSingleTreePopulationInverseGradientWrtNodeHeight(SkyGlideLikelihood likelihood, int treeIndex, double[] derivatives) {

            }
        };
        abstract double getNodeHeightDerivative(double intercept, double slope, double time, double lineageMultiplier);
        abstract void updateSingleTreePopulationInverseGradientWrtNodeHeight(SkyGlideLikelihood likelihood, int treeIndex, double[] derivatives);
    }

    public double[] getGradientWrtNodeHeight(int treeIndex) {
        return getDerivativeWrtNodeHeight(treeIndex, NodeHeightDerivativeType.GRADIENT);
    }

    public double[] getDiagonalHessianWrtNodeHeight(int treeIndex) {
        return getDerivativeWrtNodeHeight(treeIndex, NodeHeightDerivativeType.DIAGONAL_HESSIAN);
    }

    public double[] getDerivativeWrtNodeHeight(int treeIndex, NodeHeightDerivativeType derivativeType) {

        BigFastTreeIntervals interval = intervals.get(treeIndex);
        Tree thisTree = trees.get(treeIndex);
        double[] gradient = new double[thisTree.getInternalNodeCount()];

        int currentGridIndex = 0;
        double tmp = 0;
        double numSameHeightNodes = 1;
        for (int i = 0; i < interval.getIntervalCount(); i++) {
            final int lineageCount = interval.getLineageCount(i);
            int[] nodeIndices = interval.getNodeNumbersForInterval(i);
            final double intervalStart = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[0]));
            final double intervalEnd = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[1]));

            if (!(thisTree.isExternal(thisTree.getNode(nodeIndices[0])) && thisTree.isExternal(thisTree.getNode(nodeIndices[1])))) {
                int[] gridIndices = getGridPoints(currentGridIndex, intervalStart, intervalEnd);
                final int firstGridIndex = gridIndices[0];
                final int lastGridIndex = gridIndices[1];

                if (intervalStart == intervalEnd) {
                    if (interval.getIntervalType(i) == IntervalType.COALESCENT)
                        numSameHeightNodes++;
                } else {
                    final double firstGridSlope = getGridSlope(firstGridIndex);
                    final double firstGridIntercept = getGridIntercept(firstGridIndex);

                    final double lastGridSlope = getGridSlope(lastGridIndex);
                    final double lastGridIntercept = getGridIntercept(lastGridIndex);

                    final double lineageMultiplier = 0.5 * lineageCount * (lineageCount - 1);
                    if (!thisTree.isExternal(thisTree.getNode(nodeIndices[0]))) {
                        tmp += derivativeType.getNodeHeightDerivative(firstGridIntercept, firstGridSlope, intervalStart, lineageMultiplier);
                    }

                    int count = 0;
                    int j = 0;
                    while(numSameHeightNodes - count > 0 && tmp != 0) {
                        final int nodeIndex = interval.getNodeNumbersForInterval(i - j)[0];
                        if (!thisTree.isExternal(thisTree.getNode(nodeIndex))) {
                            count++;
                            gradient[nodeIndex - thisTree.getExternalNodeCount()] += tmp / numSameHeightNodes;
                        }
                        j++;
                    }
                    numSameHeightNodes = 1;
                    tmp = 0;

                    if (interval.getIntervalType(i) == IntervalType.COALESCENT) {
                        tmp = -derivativeType.getNodeHeightDerivative(lastGridIntercept, lastGridSlope, intervalEnd, lineageMultiplier);
                    }
                }
                currentGridIndex = lastGridIndex;
            }
        }
        int count = 0;
        int j = 0;
        while(numSameHeightNodes - count > 0 && tmp != 0) {
            final int nodeIndex = interval.getNodeNumbersForInterval(interval.getIntervalCount() - 1 - j)[1];
            if (!thisTree.isExternal(thisTree.getNode(nodeIndex))) {
                count++;
                gradient[nodeIndex - thisTree.getExternalNodeCount()] += tmp / numSameHeightNodes;
            }
            j++;
        }

        derivativeType.updateSingleTreePopulationInverseGradientWrtNodeHeight(this, treeIndex, gradient);

        return gradient;
    }


    public double getSingleTreeLogLikelihood(int index) {
        BigFastTreeIntervals interval = intervals.get(index);
        Tree thisTree = trees.get(index);
        int currentGridIndex = 0;
        double lnL = 0;
        for (int i = 0; i < interval.getIntervalCount(); i++) {
            final int lineageCount = interval.getLineageCount(i);
            int[] nodeIndices = interval.getNodeNumbersForInterval(i);
            final double intervalStart = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[0]));
            final double intervalEnd = thisTree.getNodeHeight(thisTree.getNode(nodeIndices[1]));

            if (intervalStart != intervalEnd) {
                int[] gridIndices = getGridPoints(currentGridIndex, intervalStart, intervalEnd);
                final int firstGridIndex = gridIndices[0];
                final int lastGridIndex = gridIndices[1];
                double sum = 0;

                if (firstGridIndex == lastGridIndex) {
                    sum += getLinearInverseIntegral(intervalStart, intervalEnd, firstGridIndex);
                } else {
                    sum += getLinearInverseIntegral(intervalStart, gridPointParameter.getParameterValue(firstGridIndex), firstGridIndex);
                    currentGridIndex = firstGridIndex;
                    while(currentGridIndex + 1 < lastGridIndex) {
                        sum += getLinearInverseIntegral(gridPointParameter.getParameterValue(currentGridIndex), gridPointParameter.getParameterValue(currentGridIndex + 1), currentGridIndex + 1);
                        currentGridIndex++;
                    }
                    sum += getLinearInverseIntegral(gridPointParameter.getParameterValue(currentGridIndex), intervalEnd, currentGridIndex + 1);
                }
                currentGridIndex = lastGridIndex;
                lnL -= 0.5 * lineageCount * (lineageCount - 1) * sum;
            }
        }
        lnL += getSingleTreePopulationInverseLogLikelihood(index);
        return lnL;
    }

    private double getSingleTreePopulationInverseLogLikelihood(int index) {
        int currentGridIndex = 0;
        double lnL = 0;
        BigFastTreeIntervals interval = intervals.get(index);

        for (int i = 0; i < interval.getIntervalCount(); i++) {
            if (interval.getIntervalType(i) == IntervalType.COALESCENT) {
                final double time = interval.getIntervalTime(i + 1);
                currentGridIndex = getGridIndex(time, currentGridIndex);
                lnL -= getLogPopulationSize(time, currentGridIndex);
            }
        }

        return lnL;
    }

    private void updateSingleTreePopulationInverseGradientWrtLogPopSize(int index, double[] gradient) {
        int currentGridIndex = 0;
        BigFastTreeIntervals interval = intervals.get(index);

        for (int i = 0; i < interval.getIntervalCount(); i++) {
            if (interval.getIntervalType(i) == IntervalType.COALESCENT) {
                final double time = interval.getIntervalTime(i + 1);
                currentGridIndex = getGridIndex(time, currentGridIndex);
                updateLogPopSizeDerivative(time, currentGridIndex, gradient);
            }
        }
    }

    private double getLogPopulationSize(double time, int gridIndex) {
        final double slope = getGridSlope(gridIndex);
        final double intercept = getGridIntercept(gridIndex);
        return intercept + slope * time;
    }

    private void updateLogPopSizeDerivative(double time, int gridIndex, double[] gradient) {
        updateGridSlopeDerivativeWrtLogPopSize(gridIndex, gradient, -time);
        updateGridInterceptDerivativeWrtLogPopSize(gridIndex, gradient, -1.0);
    }

    private int getGridIndex(double time, int startGridIndex) {
        int index = startGridIndex;
        while (index < gridPointParameter.getDimension() && gridPointParameter.getParameterValue(index) < time) {
            index++;
        }
        return index;
    }

    private double getLinearInverseIntegral(double start, double end, int gridIndex) {
        final double slope = getGridSlope(gridIndex);
        final double intercept = getGridIntercept(gridIndex);
        assert(slope != 0 || intercept != 0);
        if (start == end) {
            return 0;
        }

        final double realSmall = getMagicUnderFlowBound(slope);

        if (Math.abs(slope) < realSmall) {
            return Math.exp(-intercept) * (end - start);
        } else {
            return Math.exp(-intercept) * (Math.exp(-slope * start) - Math.exp(-slope * end)) / slope;
        }
    }

    private double getMagicUnderFlowBound(double slope) { // TODO: arbitrary magic bound
        return MachineAccuracy.SQRT_EPSILON*(Math.abs(slope) + 1.0);
    }

    private void updateIntervalGradientWrtLogPopSize(double intervalStart, double intervalEnd, int gridIndex, int lineageCount,
                                                     double[] gradient) {
        final double slope = getGridSlope(gridIndex);
        final double intercept = getGridIntercept(gridIndex);
        final double lineageMultiplier = -0.5 * lineageCount * (lineageCount - 1);
        assert(slope != 0 || intercept != 0);
        final double realSmall = getMagicUnderFlowBound(slope);
        if (intervalStart != intervalEnd) {
            final double slopeMultiplier = Math.abs(slope) < realSmall ? Math.exp(-intercept) * (intervalStart * intervalStart - intervalEnd * intervalEnd) / 2
                    : Math.exp(-intercept) * ( (-intervalStart * Math.exp(-slope * intervalStart) + intervalEnd * Math.exp(-slope * intervalEnd))
                    - (Math.exp(-slope * intervalStart) - Math.exp(-slope * intervalEnd)) / slope) / slope;
            final double interceptMultiplier = Math.abs(slope) < realSmall ? (intervalEnd - intervalStart) * (-Math.exp(-intercept))
                    : Math.exp(-intercept) * (-Math.exp(-slope * intervalStart ) + Math.exp(-slope * intervalEnd )) / slope;

                updateGridInterceptDerivativeWrtLogPopSize(gridIndex, gradient, lineageMultiplier * interceptMultiplier);
                updateGridSlopeDerivativeWrtLogPopSize(gridIndex, gradient, lineageMultiplier * slopeMultiplier);
        }
    }


    private double getGridSlope(int gridIndex) {
        if (gridIndex == gridPointParameter.getDimension()) {
            return 0;
        }
        final double thisGridTime = gridPointParameter.getParameterValue(gridIndex);
        final double lastGridTime = gridIndex == 0 ? 0 : gridPointParameter.getParameterValue(gridIndex - 1);
        return (logPopSizeParameter.getParameterValue(gridIndex + 1) - logPopSizeParameter.getParameterValue(gridIndex))
                / (thisGridTime - lastGridTime);
    }

    private void updateGridSlopeDerivativeWrtLogPopSize(int gridIndex, double[] gradient, double multiplier) {
        if (gridIndex != gridPointParameter.getDimension()) {
            final double thisGridTime = gridPointParameter.getParameterValue(gridIndex);
            final double lastGridTime = gridIndex == 0 ? 0 : gridPointParameter.getParameterValue(gridIndex - 1);

            gradient[gridIndex + 1] += multiplier / (thisGridTime - lastGridTime);
            gradient[gridIndex] -= multiplier / (thisGridTime - lastGridTime);
        }
    }

    private double getGridIntercept(int gridIndex) {
        if (gridIndex == gridPointParameter.getDimension() || gridIndex == 0) {
            return logPopSizeParameter.getParameterValue(gridIndex);
        }

        final double thisGridTime = gridPointParameter.getParameterValue(gridIndex);
        final double lastGridTime = gridPointParameter.getParameterValue(gridIndex - 1);

        return (thisGridTime * logPopSizeParameter.getParameterValue(gridIndex)
                - lastGridTime * logPopSizeParameter.getParameterValue(gridIndex + 1)) /
                (thisGridTime - lastGridTime);
    }

    private void updateGridInterceptDerivativeWrtLogPopSize(int gridIndex, double[] gradient, double multiplier) {
        if (gridIndex == gridPointParameter.getDimension() || gridIndex == 0) {
            gradient[gridIndex] += multiplier;
        } else {
            final double thisGridTime = gridPointParameter.getParameterValue(gridIndex);
            final double lastGridTime = gridPointParameter.getParameterValue(gridIndex - 1);

            final double firstDerivative = thisGridTime / (thisGridTime - lastGridTime) * multiplier;
            final double secondDerivative = -lastGridTime / (thisGridTime - lastGridTime) * multiplier;

            gradient[gridIndex] += firstDerivative;
            gradient[gridIndex + 1] += secondDerivative;
        }
    }

    private int[] getGridPoints(int startGridIndex, double startTime, double endTime) { // return smallest grid index that is higher than each time
        int i = startGridIndex;
        while (i < gridPointParameter.getDimension() && gridPointParameter.getParameterValue(i) < startTime) {
            i++;
        }
        int firstGridIndex = i;

        while (i < gridPointParameter.getDimension() && gridPointParameter.getParameterValue(i) < endTime) {
            i++;
        }
        int lastGridIndex = i;
        return new int[]{firstGridIndex, lastGridIndex};
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }
}
