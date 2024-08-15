/*
 * CoalescentGradient.java
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

package dr.evolution.coalescent;


import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.GradientProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.Binomial;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CoalescentGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    private final CoalescentLikelihood likelihood;
    private final Parameter parameter;
    private final Tree tree;
    private final GradientProvider provider;

    public enum Wrt {
        NODE_HEIGHTS,
        PARAMETER
    }

    public CoalescentGradient(CoalescentLikelihood likelihood,
                              TreeModel tree,
                              Parameter wrtParameter,
                              Wrt wrt,
                              double tolerance) {
        this.likelihood = likelihood;
        this.tree = tree;
        if (wrt == Wrt.NODE_HEIGHTS) {
            this.parameter = new NodeHeightProxyParameter("NodeHeights", tree, true);
            this.provider = new GradientProvider() {
                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return getGradientLogDensityWrtNodeHeights();
                }
            };
        } else {
            this.parameter = wrtParameter;
            this.provider = null; // TODO return gradient wrt parameter
        }
        this.tolerance = tolerance;
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

        if (likelihood.getPopulationSizeModel() != null) {
            throw new RuntimeException("Not yet implemented!");
        }

        return provider.getGradientLogDensity(null);
    }

    private double[] getGradientLogDensityWrtNodeHeights() {

        final double logLikelihood = likelihood.getLogLikelihood();
        double[] gradient = new double[tree.getInternalNodeCount()];

        if (logLikelihood == Double.NEGATIVE_INFINITY) {
            Arrays.fill(gradient, Double.NaN);
            return gradient;
        }

        IntervalList intervals = likelihood.getIntervalList();
        BigFastTreeIntervals bigFastTreeIntervals = (BigFastTreeIntervals) intervals;

        DemographicFunction demographicFunction = likelihood.getDemoModel().getDemographicFunction();

        int numSameHeightNodes = 1;
        double thisGradient = 0;
        for (int i = 0; i < bigFastTreeIntervals.getIntervalCount(); i++) {
            if (bigFastTreeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                final double time = bigFastTreeIntervals.getIntervalTime(i + 1);
                final int lineageCount = bigFastTreeIntervals.getLineageCount(i);
                final double kChoose2 = Binomial.choose2(lineageCount);
                final double intensityGradient = demographicFunction.getIntensityGradient(time);
                thisGradient += demographicFunction.getLogDemographicGradient(time);

                if (bigFastTreeIntervals.getInterval(i) != 0) {
                    thisGradient -= kChoose2 * intensityGradient;
                } else {
                    numSameHeightNodes++;
                }

                if ( i < bigFastTreeIntervals.getIntervalCount() - 1
                        && bigFastTreeIntervals.getInterval(i + 1) != 0) {

                    final int nextLineageCount = bigFastTreeIntervals.getLineageCount(i + 1);
                    thisGradient += Binomial.choose2(nextLineageCount) * intensityGradient;

                    for (int j = 0; j < numSameHeightNodes; j++) {
                        final int nodeIndex = bigFastTreeIntervals.getNodeNumbersForInterval(i - j)[1];
                        gradient[nodeIndex - tree.getExternalNodeCount()] = thisGradient / (double) numSameHeightNodes;
                    }

                    thisGradient = 0;
                    numSameHeightNodes = 1;
                }
            }
        }

        int j = numSameHeightNodes;
        int v = bigFastTreeIntervals.getIntervalCount() - 1;
        while(j > 0) {
            if (bigFastTreeIntervals.getIntervalType(v) == IntervalType.COALESCENT) {
                gradient[bigFastTreeIntervals.getNodeNumbersForInterval(v)[1] - tree.getExternalNodeCount()] = thisGradient / (double) numSameHeightNodes;
                j--;
            }
            v--;
        }

        return gradient;
    }

    private final double tolerance;

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, tolerance);
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "CoalescentGradient check");
    }
}
