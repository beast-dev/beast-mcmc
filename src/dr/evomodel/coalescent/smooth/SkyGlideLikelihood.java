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

import dr.evolution.tree.NodeRef;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
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
            this.intervals.add(new BigFastTreeIntervals(trees.get(0)));
        }
    }

    @Override
    public String getReport() {
        return null;
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
        return null;
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
        int currentGridIndex = 0;
        double lnL = 0;
        for (int i = 0; i < interval.getIntervalCount() - 1; i++) {
            final double intervalStart = interval.getIntervalTime(i);
            final double intervalEnd = interval.getIntervalTime(i + 1);
            final int lineageCount = interval.getLineageCount(i);
            int[] gridIndices = getGridPoints(currentGridIndex, intervalStart, intervalEnd);
            final int firstGridIndex = gridIndices[0];
            final int lastGridIndex = gridIndices[1];
            if (firstGridIndex == Integer.MAX_VALUE) { // no grid points within interval
                lnL += 0.5 * lineageCount * (lineageCount - 1) * getLinearInverseIntegral(intervalStart, intervalEnd, currentGridIndex);
            } else {
                double sum = 0;
                sum += getLinearInverseIntegral(intervalStart, gridPointParameter.getParameterValue(firstGridIndex), currentGridIndex);
                currentGridIndex = firstGridIndex;
                while(currentGridIndex < lastGridIndex) {
                    sum += getLinearInverseIntegral(gridPointParameter.getParameterValue(currentGridIndex), gridPointParameter.getParameterValue(currentGridIndex + 1), currentGridIndex);
                    currentGridIndex++;
                }
                sum += getLinearInverseIntegral(gridPointParameter.getParameterValue(currentGridIndex), intervalEnd, currentGridIndex);
                lnL += 0.5 * lineageCount * (lineageCount - 1) * sum;
            }
        }
        lnL += getSingleTreePopulationInverseLogLikelihood(index);
        return lnL;
    }

    private double getSingleTreePopulationInverseLogLikelihood(int index) {
        TreeModel tree = trees.get(index);
        int currentGridIndex = 0;
        double lnL = 0;

        NodeRef[] nodes = new NodeRef[tree.getNodeCount()];
        System.arraycopy(tree.getNodes(), 0, nodes, 0, tree.getNodeCount());
        Arrays.parallelSort(nodes, (a, b) -> Double.compare(tree.getNodeHeight(a), tree.getNodeHeight(b)));

        for (NodeRef node : nodes) {
            final double time = tree.getNodeHeight(node);
            currentGridIndex = getGridIndex(time, currentGridIndex);
            lnL -= getLogPopulationSize(time, currentGridIndex);
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
        double gridStart = gridPointParameter.getParameterValue(startGridIndex);
        while (index < gridPointParameter.getDimension() - 1 && gridPointParameter.getParameterValue(index + 1) < time) {
            index++;
        }
        return index;
    }

    private double getLinearInverseIntegral(double start, double end, int gridIndex) {
        final double slope = getGridSlope(gridIndex);
        final double intercept = getGridIntercept(gridIndex);
        assert(slope != 0 || intercept != 0);
        if (slope == 0) {
//            return (end - start) / intercept;
            return Math.exp(-intercept) * (end - start);
        } else {
//            return (Math.log(slope * end + intercept) - Math.log(slope * start + intercept)) / slope;
            return (Math.exp(-(slope * start + intercept)) - Math.exp(-(slope * end + intercept))) / slope;
        }
    }

    private double getGridSlope(int gridIndex) {
        if (gridIndex == gridPointParameter.getDimension() - 1) {
            return 0;
        }
        return (logPopSizeParameter.getParameterValue(gridIndex + 1) - logPopSizeParameter.getParameterValue(gridIndex))
                / (gridPointParameter.getParameterValue(gridIndex + 1) - gridPointParameter.getParameterValue(gridIndex));
    }

    private double getGridIntercept(int gridIndex) {
        if (gridIndex == gridPointParameter.getDimension() - 1) {
            return logPopSizeParameter.getParameterValue(gridIndex);
        }
        return (gridPointParameter.getParameterValue(gridIndex + 1) * logPopSizeParameter.getParameterValue(gridIndex)
                - gridPointParameter.getParameterValue(gridIndex) * logPopSizeParameter.getParameterValue(gridIndex + 1)) /
                (gridPointParameter.getParameterValue(gridIndex + 1) - gridPointParameter.getParameterValue(gridIndex));
    }

    private int[] getGridPoints(int startGridIndex, double startTime, double endTime) {
        int firstGridIndex = Integer.MAX_VALUE;
        int lastGridIndex = -1;
        int i = startGridIndex;
        double time = gridPointParameter.getParameterValue(i);
        while (time < endTime) {
            if (time >= startTime) {
                if (firstGridIndex > i) firstGridIndex = i;
                if (lastGridIndex < i) lastGridIndex = i;
            }
            i++;
            time = gridPointParameter.getParameterValue(i);
        }
        return new int[]{firstGridIndex, lastGridIndex};
    }

    @Override
    public void makeDirty() {

    }
}
