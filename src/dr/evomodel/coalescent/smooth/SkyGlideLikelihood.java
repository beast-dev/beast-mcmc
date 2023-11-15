/*
 * SkyGlideLikelihood.java
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

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
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
    }

    @Override
    public String getReport() {
        return "skyGlideLikelihood(" + getLogLikelihood() + ")";
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

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
        double lnL = 0;
        for (int i = 0; i < trees.size(); i++) {
            lnL += getSingleTreeLogLikelihood(i);
        }
        return lnL;
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
                    if (firstGridIndex < gridPointParameter.getDimension() - 1) {
                        sum += getLinearInverseIntegral(intervalStart, intervalEnd, firstGridIndex);
                    } else {
                        sum += getLinearInverseIntegral(intervalStart, intervalEnd, gridPointParameter.getDimension() - 1);
                    }
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

    private double getLogPopulationSize(double time, int gridIndex) {
        final double slope = getGridSlope(gridIndex);
        final double intercept = getGridIntercept(gridIndex);
        return intercept + slope * time;
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

        if (slope == 0) {
            return Math.exp(-intercept) * (end - start);
        } else {
            return (Math.exp(-(slope * start + intercept)) - Math.exp(-(slope * end + intercept))) / slope;
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

    }
}
