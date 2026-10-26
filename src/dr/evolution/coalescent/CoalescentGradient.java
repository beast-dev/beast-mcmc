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
import dr.evolution.tree.NodeRef;
import dr.evomodel.coalescent.CoalescentLikelihood;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.demographicmodel.DifferentiableDemographicModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.GradientProvider;
import dr.inference.model.HessianProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.Binomial;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CoalescentGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    private final CoalescentLikelihood likelihood;
    private final Parameter parameter;
    private final Tree tree;
    private final HessianProvider provider;

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
        if (likelihood.getPopulationSizeModel() != null) {
            throw new IllegalArgumentException("Coalescent gradients are not yet implemented for a population size model");
        }
        if (wrt == Wrt.NODE_HEIGHTS) {
            this.parameter = new NodeHeightProxyParameter("NodeHeights", tree, true);
            this.provider = new HessianProvider() {
                @Override
                public double[] getDiagonalHessianLogDensity(Object x) {
                    return getNodeHeightDerivatives(true);
                }

                @Override
                public double[][] getHessianLogDensity(Object x) {
                    throw new RuntimeException("Not implemented yet");
                }

                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return getNodeHeightDerivatives(false);
                }
            };
        } else {
            checkDifferentiable(likelihood.getDemoModel(), wrtParameter);
            this.parameter = wrtParameter;
            this.provider = new HessianProvider() {
                @Override
                public double[] getDiagonalHessianLogDensity(Object x) {
                    throw new UnsupportedOperationException("Not yet implemented.");
                }

                @Override
                public double[][] getHessianLogDensity(Object x) {
                    throw new UnsupportedOperationException("Not yet implemented.");
                }

                @Override
                public int getDimension() {
                    return 1;
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return getGradientLogDensityWrtParameter();
                }
            };
        }
        this.wrt = wrt;
        this.tolerance = tolerance;
    }

    private static void checkDifferentiable(DemographicModel model, Parameter wrtParameter) {
        if (wrtParameter == null || wrtParameter.getDimension() != 1) {
            throw new IllegalArgumentException("The coalescent gradient needs a scalar demographic parameter");
        }
        if (!(model instanceof DifferentiableDemographicModel)
                || !((DifferentiableDemographicModel) model).isDifferentiableWrtParameter(wrtParameter)) {
            throw new IllegalArgumentException("No analytic gradient of " + model.getId() + " with respect to "
                    + wrtParameter.getId());
        }
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
    public double[] getDiagonalHessianLogDensity() {
        return provider.getDiagonalHessianLogDensity(null);
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new UnsupportedOperationException("Not yet implemented.");
    }

    @Override
    public double[] getGradientLogDensity() {
        return provider.getGradientLogDensity(null);
    }

    /*
     * Moving the m internal nodes at time t (m > 1 only for tied heights) lengthens the interval below them, with k
     * lineages, and shortens the one above, with k - m, and moves m of the -log N(t) terms, so
     * d log L / dt = -[C(k, 2) - C(k - m, 2)] I'(t) - m (log N)'(t), split equally among the m nodes. The sweep reads
     * the tree directly, so it needs no node-to-interval mapping and works with any interval list over the whole tree.
     */
    private double[] getNodeHeightDerivatives(boolean secondDerivative) {

        final double logLikelihood = likelihood.getLogLikelihood();
        double[] derivatives = new double[tree.getInternalNodeCount()];

        if (logLikelihood == Double.NEGATIVE_INFINITY) {
            Arrays.fill(derivatives, Double.NaN);
            return derivatives;
        }

        IntervalList intervals = likelihood.getIntervalList();
        checkCoversWholeTree(intervals);
        final double startTime = intervals.getStartTime();

        DemographicFunction demographicFunction = likelihood.getDemoModel().getDemographicFunction();

        // by height, with tips before internal nodes of the same height
        final int nodeCount = tree.getNodeCount();
        Integer[] order = new Integer[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (a, b) -> {
            int byHeight = Double.compare(tree.getNodeHeight(tree.getNode(a)), tree.getNodeHeight(tree.getNode(b)));
            return byHeight != 0 ? byHeight
                    : Boolean.compare(!tree.isExternal(tree.getNode(a)), !tree.isExternal(tree.getNode(b)));
        });

        final int externalCount = tree.getExternalNodeCount();
        int lineageCount = 0;
        int i = 0;
        while (i < nodeCount) {
            final NodeRef node = tree.getNode(order[i]);
            if (tree.isExternal(node)) {
                lineageCount++;
                i++;
                continue;
            }
            final double height = tree.getNodeHeight(node);
            int j = i + 1;
            while (j < nodeCount && !tree.isExternal(tree.getNode(order[j]))
                    && tree.getNodeHeight(tree.getNode(order[j])) == height) {
                j++;
            }
            final int tied = j - i;
            final double time = height - startTime;
            final double kChoose2Change = Binomial.choose2(lineageCount) - Binomial.choose2(lineageCount - tied);

            final double value;
            if (secondDerivative) {
                value = -kChoose2Change * demographicFunction.getIntensitySecondDerivative(time)
                        - tied * demographicFunction.getLogDemographicSecondDerivative(time);
            } else {
                // log L has -log N(t) at each coalescent time, and getLogDemographicGradient is d log N / dt
                value = -kChoose2Change * demographicFunction.getIntensityGradient(time)
                        - tied * demographicFunction.getLogDemographicGradient(time);
            }
            for (int k = i; k < j; k++) {
                derivatives[order[k] - externalCount] = value / tied;
            }
            lineageCount -= tied;
            i = j;
        }

        return derivatives;
    }

    private void checkCoversWholeTree(IntervalList intervals) {
        int coalescentCount = 0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                coalescentCount++;
            }
        }
        if (coalescentCount != tree.getInternalNodeCount()) {
            throw new UnsupportedOperationException("Node-height gradients need intervals over the whole tree "
                    + "(no included or excluded subtrees)");
        }
    }

    /*
     * log L = sum_i { -C(k_i, 2) [I(f_i) - I(s_i)] } - sum_{coalescent i} log N(f_i), over intervals [s_i, f_i], so
     * d log L / d theta = sum_n (d a_n / d theta) d log L / d a_n over the demographic function's arguments a_n.
     */
    private double[] getGradientLogDensityWrtParameter() {

        final double logLikelihood = likelihood.getLogLikelihood();

        if (logLikelihood == Double.NEGATIVE_INFINITY) {
            return new double[]{Double.NaN};
        }

        IntervalList intervals = likelihood.getIntervalList();
        DemographicModel demographicModel = likelihood.getDemoModel();
        DemographicFunction demographicFunction = demographicModel.getDemographicFunction();
        final double[] argumentDerivatives =
                ((DifferentiableDemographicModel) demographicModel).getArgumentDerivatives(parameter);

        double gradient = 0.0;
        double startTime = 0.0;
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            final double finishTime = startTime + intervals.getInterval(i);
            final double kChoose2 = Binomial.choose2(intervals.getLineageCount(i));
            final boolean coalescent = intervals.getIntervalType(i) == IntervalType.COALESCENT;

            for (int n = 0; n < argumentDerivatives.length; n++) {
                if (argumentDerivatives[n] == 0.0) {
                    continue;
                }
                double term = 0.0;
                if (finishTime != startTime) {
                    term -= kChoose2 * (demographicFunction.getIntensityArgumentGradient(n, finishTime)
                            - demographicFunction.getIntensityArgumentGradient(n, startTime));
                }
                if (coalescent) {
                    term -= demographicFunction.getLogDemographicArgumentGradient(n, finishTime);
                }
                gradient += argumentDerivatives[n] * term;
            }
            startTime = finishTime;
        }

        return new double[]{gradient};
    }

    private final Wrt wrt;
    private final double tolerance;

    @Override
    public String getReport() {
        // growth rates can be negative; node heights and sizes are positive
        final double lowerBound = wrt == Wrt.PARAMETER
                ? (parameter.getBounds() == null ? Double.NEGATIVE_INFINITY : parameter.getBounds().getLowerLimit(0)) : 0.0;
        return GradientWrtParameterProvider.getReportAndCheckForError(this, lowerBound, Double.POSITIVE_INFINITY, tolerance);
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "CoalescentGradient check");
    }
}
