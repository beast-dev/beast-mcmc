/*
 * CoalescentGradient.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.coalescent;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.Binomial;
import dr.util.ComparableDouble;
import dr.util.HeapSort;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CoalescentGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    private final CoalescentLikelihood likelihood;
    private final Parameter parameter;
    private final Tree tree;

    public CoalescentGradient(CoalescentLikelihood likelihood,
                              TreeModel tree) {
        this.likelihood = likelihood;
        this.tree = tree;
        this.parameter = new NodeHeightProxyParameter("NodeHeights", tree, true);
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        final double logLikelihood = likelihood.getLogLikelihood();
        double[] gradient = new double[tree.getInternalNodeCount()];

        if (logLikelihood == Double.NEGATIVE_INFINITY) {
            Arrays.fill(gradient, Double.NaN);
            return gradient;
        }

        int[] intervalIndices = new int[tree.getInternalNodeCount()];
        int[] nodeIndices = new int[tree.getInternalNodeCount()];
        double[] sortedValues = new double[tree.getInternalNodeCount()];
        getIntervalIndexForInternalNodes(intervalIndices, nodeIndices, sortedValues);

        IntervalList intervals = likelihood.getIntervalList();

        DemographicFunction demographicFunction = likelihood.getDemoModel().getDemographicFunction();

        int numSameHeightNodes = 1;
        double thisGradient = 0.0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getNode(tree.getExternalNodeCount() + nodeIndices[i]);
            final int lineageCount = intervals.getLineageCount(intervalIndices[nodeIndices[i]]);

            final double intensityGradient = demographicFunction.getIntensityGradient(tree.getNodeHeight(node));

            final double intervalLength = intervals.getInterval(intervalIndices[nodeIndices[i]]);

            thisGradient -= demographicFunction.getLogDemographicGradient(tree.getNodeHeight(node));
            if ( intervalLength != 0.0) {
                thisGradient -= Binomial.choose2(lineageCount) * intensityGradient;
            } else {
                numSameHeightNodes++;
            }

            if (!tree.isRoot(node) && intervals.getInterval(intervalIndices[nodeIndices[i]] + 1) != 0.0) {
                final int nextLineageCount = intervals.getLineageCount(intervalIndices[nodeIndices[i]] + 1);
                thisGradient += Binomial.choose2(nextLineageCount) * intensityGradient;

                for (int j = 0; j < numSameHeightNodes; j++) {
                    gradient[nodeIndices[i - j]] = thisGradient / ((double) numSameHeightNodes);
                }
                thisGradient = 0.0;
                numSameHeightNodes = 1;
            }
        }
        for (int j = 0; j < numSameHeightNodes; j++) {
            gradient[nodeIndices[tree.getInternalNodeCount() - j - 1]] = thisGradient / ((double) numSameHeightNodes);
        }

        return gradient;
    }

    private void getIntervalIndexForInternalNodes(int[] intervalIndices, int[] nodeIndices, double[] sortedValues) {
        double[] nodeHeights = new double[tree.getInternalNodeCount()];
        ArrayList<ComparableDouble> sortedInternalNodes = new ArrayList<ComparableDouble>();
        for (int i = 0; i < nodeIndices.length; i++) {
            final double nodeHeight = tree.getNodeHeight(tree.getNode(tree.getExternalNodeCount() + i));
            sortedInternalNodes.add(new ComparableDouble(nodeHeight));
            nodeHeights[i] = nodeHeight;
        }
        HeapSort.sort(sortedInternalNodes, nodeIndices);
        for (int i = 0; i < nodeIndices.length; i++) {
            sortedValues[i] = nodeHeights[nodeIndices[i]];
        }

        IntervalList intervals = likelihood.getIntervalList();
        int intervalIndex = 0;
        double finishTime = intervals.getInterval(intervalIndex);
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            while(intervalIndex < intervals.getIntervalCount() - 1 && sortedValues[i] - finishTime > realSmallNumber) {
                intervalIndex++;
                finishTime += intervals.getInterval(intervalIndex);
            }
            intervalIndices[nodeIndices[i]] = intervalIndex;
        }
    }

    private final double realSmallNumber = 1E-10;

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1E-2);
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "CoalescentGradient check");
    }
}
