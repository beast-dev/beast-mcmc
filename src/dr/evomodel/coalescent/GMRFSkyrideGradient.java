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

import dr.evomodel.treedatalikelihood.discrete.NodeHeightTransform;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class GMRFSkyrideGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    private final GMRFSkyrideLikelihood skyrideLikelihood;
    private final WrtParameter wrtParameter;
    private final Parameter parameter;
    private final OldAbstractCoalescentLikelihood.IntervalNodeMapping intervalNodeMapping;
    private final NodeHeightTransform nodeHeightTransform;

    public GMRFSkyrideGradient(GMRFSkyrideLikelihood gmrfSkyrideLikelihood,
                               WrtParameter wrtParameter,
                               Parameter parameter,
                               NodeHeightTransform nodeHeightTransform) {

        this.skyrideLikelihood = gmrfSkyrideLikelihood;
        this.intervalNodeMapping = skyrideLikelihood.getIntervalNodeMapping();
        this.wrtParameter = wrtParameter;
        this.nodeHeightTransform = nodeHeightTransform;
        this.parameter = parameter;
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
        return wrtParameter.getGradientLogDensity(skyrideLikelihood, intervalNodeMapping);
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
        double[] savedValues = getParameter().getParameterValues();
        double[] testGradient = NumericalDerivative.gradient(numeric1, getParameter().getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            getParameter().setParameterValue(i, savedValues[i]);
        }
        double[] testDiagonalHessian = NumericalDerivative.diagonalHessian(numeric1, getParameter().getParameterValues());
        for (int i = 0; i < savedValues.length; ++i) {
            getParameter().setParameterValue(i, savedValues[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("analytic: ").append(new dr.math.matrixAlgebra.Vector(getGradientLogDensity()));
        sb.append("\n");
        sb.append("numeric: ").append(new dr.math.matrixAlgebra.Vector(testGradient));
        sb.append("\n");

        sb.append("analytic diagonal Hessian: ").append(new dr.math.matrixAlgebra.Vector(getDiagonalHessianLogDensity()));
        sb.append("\n");
        sb.append("numeric diagonal Hessian: ").append(new dr.math.matrixAlgebra.Vector(testDiagonalHessian));
        sb.append("\n");

        return sb.toString();
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return new double[getDimension()];
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!");
    }

    public enum WrtParameter {

        COALESCENT_INTERVAL {
            @Override
            double[] getGradientLogDensity(GMRFSkyrideLikelihood skyrideLikelihood,
                                           OldAbstractCoalescentLikelihood.IntervalNodeMapping intervalNodeMapping) {
                double[] unSortedNodeHeightGradient = super.getGradientLogDensityWrtUnsortedNodeHeight(skyrideLikelihood);
                double[] intervalGradient = new double[unSortedNodeHeightGradient.length];
                double accumulatedGradient = 0.0;
                for (int i = unSortedNodeHeightGradient.length - 1; i > -1; i--) {
                    accumulatedGradient += unSortedNodeHeightGradient[i];
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
            double[] getGradientLogDensity(GMRFSkyrideLikelihood skyrideLikelihood,
                                           OldAbstractCoalescentLikelihood.IntervalNodeMapping intervalNodeMapping) {
                double[] unSortedNodeHeightGradient = getGradientLogDensityWrtUnsortedNodeHeight(skyrideLikelihood);
                return intervalNodeMapping.sortByNodeNumbers(unSortedNodeHeightGradient);
            }

            @Override
            void update(NodeHeightTransform nodeHeightTransform, double[] values) {
                nodeHeightTransform.transform(values, 0, values.length);
            }
        };

        abstract double[] getGradientLogDensity(GMRFSkyrideLikelihood skyrideLikelihood,
                                                OldAbstractCoalescentLikelihood.IntervalNodeMapping intervalNodeMapping);

        abstract void update(NodeHeightTransform nodeHeightTransform, double[] values);

        double[] getGradientLogDensityWrtUnsortedNodeHeight(GMRFSkyrideLikelihood skyrideLikelihood) {
            double[] unSortedNodeHeightGradient = new double[skyrideLikelihood.getCoalescentIntervalDimension()];
            double[] gamma = skyrideLikelihood.getPopSizeParameter().getParameterValues();

            int index = 0;
            for (int i = 0; i < skyrideLikelihood.getIntervalCount(); i++) {
                if (skyrideLikelihood.getIntervalType(i) == OldAbstractCoalescentLikelihood.CoalescentEventType.COALESCENT) {
                    double weight = -Math.exp(-gamma[index]) * skyrideLikelihood.getLineageCount(i) * (skyrideLikelihood.getLineageCount(i) - 1);
                    if (index < skyrideLikelihood.getCoalescentIntervalDimension() - 1 && i < skyrideLikelihood.getIntervalCount() - 1) {
                        weight -= -Math.exp(-gamma[index + 1]) * skyrideLikelihood.getLineageCount(i + 1) * (skyrideLikelihood.getLineageCount(i + 1) - 1);
                    }
                    unSortedNodeHeightGradient[index] = weight / 2.0;
                    index++;
                }
            }
            return unSortedNodeHeightGradient;
        }
    }
}
