/*
 * GMRFSkyrideGradient.java
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

package dr.evomodel.coalescent;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightTransform;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * Gradients of the GMRF Skyride (coalescent + field) wrt its log population sizes, precision and node heights.
 *
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class GMRFSkyrideGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    private final OldGMRFSkyrideLikelihood skyrideLikelihood;
    private final WrtParameter wrtParameter;
    private final Parameter parameter;
    private final Tree tree;
    private final NodeHeightTransform nodeHeightTransform;
    private final double tolerance;

    public GMRFSkyrideGradient(OldGMRFSkyrideLikelihood gmrfSkyrideLikelihood,
                               WrtParameter wrtParameter,
                               TreeModel tree,
                               NodeHeightTransform nodeHeightTransform) {
        this(gmrfSkyrideLikelihood, wrtParameter, tree, nodeHeightTransform, 1E-4);
    }

    public GMRFSkyrideGradient(OldGMRFSkyrideLikelihood gmrfSkyrideLikelihood,
                               WrtParameter wrtParameter,
                               TreeModel tree,
                               NodeHeightTransform nodeHeightTransform,
                               double tolerance) {
        if (gmrfSkyrideLikelihood instanceof GMRFMultilocusSkyrideLikelihood) {
            throw new IllegalArgumentException("These gradients are for the Skyride; for a Skygrid use " +
                    "wrtParameter=\"logPopulationSizes\", \"precision\" or \"nodeHeight\" (GMRFGradient)");
        }
        this.skyrideLikelihood = gmrfSkyrideLikelihood;
        this.wrtParameter = wrtParameter;
        this.tree = tree != null ? tree : gmrfSkyrideLikelihood.getTree();
        this.nodeHeightTransform = nodeHeightTransform;
        this.tolerance = tolerance;
        switch (wrtParameter) {
            case LOG_POPULATION_SIZES:
                this.parameter = skyrideLikelihood.getPopSizeParameter();
                break;
            case PRECISION:
                this.parameter = skyrideLikelihood.getPrecisionParameter();
                break;
            default:
                if (nodeHeightTransform == null) {
                    this.parameter = new NodeHeightProxyParameter("internalNodeHeights", (TreeModel) this.tree, true);
                } else {
                    this.parameter = nodeHeightTransform.getParameter();
                }
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return skyrideLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return wrtParameter.getGradientLogDensity(skyrideLikelihood, tree);
    }

    private MultivariateFunction numeric1 = new MultivariateFunction() {
        @Override
        public double evaluate(double[] argument) {

            if (nodeHeightTransform != null) {
                wrtParameter.update(nodeHeightTransform, argument);
            } else {
                for (int i = 0; i < parameter.getDimension(); i++) {
                    parameter.setParameterValueQuietly(i, argument[i]);
                }
                parameter.fireParameterChangedEvent();
            }

            skyrideLikelihood.makeDirty();
            return skyrideLikelihood.getLogLikelihood();
        }

        @Override
        public int getNumArguments() {
            return getParameter().getDimension();
        }

        @Override
        public double getLowerBound(int n) {
            return 0;
        }

        @Override
        public double getUpperBound(int n) {
            return Double.POSITIVE_INFINITY;
        }
    };

    @Override
    public String getReport() {
        if (nodeHeightTransform == null) {
            return GradientWrtParameterProvider.getReportAndCheckForError(this,
                    wrtParameter.getParameterLowerBound(), Double.POSITIVE_INFINITY, tolerance);
        }

        double[] savedValues = getParameter().getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(numeric1, getParameter().getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            getParameter().setParameterValue(i, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("analytic: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");
        return sb.toString();
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return wrtParameter.getDiagonalHessianLogDensity(skyrideLikelihood, getDimension());
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!");
    }

    public enum WrtParameter {

        LOG_POPULATION_SIZES {
            @Override
            double[] getGradientLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, Tree tree) {
                return skyrideLikelihood.getGradientWrtLogPopulationSize();
            }

            @Override
            double[] getDiagonalHessianLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, int dimension) {
                return skyrideLikelihood.getDiagonalHessianWrtLogPopulationSize();
            }

            @Override
            double getParameterLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }
        },

        PRECISION {
            @Override
            double[] getGradientLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, Tree tree) {
                return skyrideLikelihood.getGradientWrtPrecision();
            }

            @Override
            double[] getDiagonalHessianLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, int dimension) {
                return skyrideLikelihood.getDiagonalHessianWrtPrecision();
            }
        },

        COALESCENT_INTERVAL {
            @Override
            double[] getGradientLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, Tree tree) {
                // h_(k) = c_0 + ... + c_k, so d / d c_j = sum_{k >= j} d / d h_(k)
                double[] heightGradient = skyrideLikelihood.getGradientWrtSortedCoalescentHeights();
                double[] intervalGradient = new double[heightGradient.length];
                double accumulatedGradient = 0.0;
                for (int i = heightGradient.length - 1; i > -1; i--) {
                    accumulatedGradient += heightGradient[i];
                    intervalGradient[i] = accumulatedGradient;
                }
                return intervalGradient;
            }

            @Override
            void update(NodeHeightTransform nodeHeightTransform, double[] values) {
                nodeHeightTransform.inverse(values, 0, values.length);
            }
        },

        NODE_HEIGHTS {
            @Override
            double[] getGradientLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, Tree tree) {
                return mapToNodes(skyrideLikelihood.getGradientWrtSortedCoalescentHeights(), tree);
            }

            @Override
            void update(NodeHeightTransform nodeHeightTransform, double[] values) {
                nodeHeightTransform.transform(values, 0, values.length);
            }
        };

        abstract double[] getGradientLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, Tree tree);

        double[] getDiagonalHessianLogDensity(OldGMRFSkyrideLikelihood skyrideLikelihood, int dimension) {
            return new double[dimension];
        }

        double getParameterLowerBound() {
            return 0.0;
        }

        void update(NodeHeightTransform nodeHeightTransform, double[] values) {
            throw new IllegalArgumentException("No node-height transform for " + name());
        }

        /*
         * The j-th internal node in height order is the j-th coalescent time. Nodes at the same height share the
         * derivative of moving them together equally.
         */
        static double[] mapToNodes(double[] sortedGradient, Tree tree) {
            final int internalCount = tree.getInternalNodeCount();
            final int externalCount = tree.getExternalNodeCount();
            Integer[] order = new Integer[internalCount];
            for (int i = 0; i < internalCount; i++) {
                order[i] = externalCount + i;
            }
            Arrays.sort(order, (a, b) -> Double.compare(tree.getNodeHeight(tree.getNode(a)),
                    tree.getNodeHeight(tree.getNode(b))));

            double[] gradient = new double[internalCount];
            int i = 0;
            while (i < internalCount) {
                final double height = tree.getNodeHeight(tree.getNode(order[i]));
                int j = i;
                double sum = 0.0;
                while (j < internalCount && tree.getNodeHeight(tree.getNode(order[j])) == height) {
                    sum += sortedGradient[j];
                    j++;
                }
                for (int k = i; k < j; k++) {
                    gradient[order[k] - externalCount] = sum / (j - i);
                }
                i = j;
            }
            return gradient;
        }
    }
}
