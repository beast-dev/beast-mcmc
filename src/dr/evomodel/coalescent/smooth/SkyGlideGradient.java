/*
 * SkyGlideGradient.java
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

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.List;

/**
 * A likelihood function for a piece-wise linear log population size coalescent process that nicely works with the newer tree intervals
 *
 * @author Mathieu Fourment
 * @author Erick Matsen
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SkyGlideGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    private final SkyGlideLikelihood likelihood;

    private final Parameter parameter;

    private final WrtParameter wrtParameter;

    private final double tolerance;

    private int treeIndex = -1;

    public SkyGlideGradient(SkyGlideLikelihood likelihood,
                            Parameter parameter,
                            double tolerance) {
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.wrtParameter = factory(parameter);
        this.tolerance = tolerance;
    }

    private WrtParameter factory(Parameter parameter) {
        if (parameter == likelihood.getLogPopSizeParameter()) {
            return WrtParameter.LOG_POP_SIZE;
        } else if (parameter instanceof NodeHeightProxyParameter){
            List<TreeModel> trees = likelihood.getTrees();
            final TreeModel tree = ((NodeHeightProxyParameter) parameter).getTree();
            for (int i = 0; i < trees.size(); i++) {
                if (trees.get(i) == tree) {
                    treeIndex = i;
                    return WrtParameter.NODE_HEIGHT;
                }
            }
            throw new RuntimeException("Parameter not recognized.");
        } else {
            throw new RuntimeException("Parameter not recognized.");
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
    public double[] getGradientLogDensity() {
        return wrtParameter.getGradientLogDensity(likelihood, treeIndex);
    }

    @Override
    public String getReport() {
        String output = GradientWrtParameterProvider.getReportAndCheckForError(this, wrtParameter.getParameterLowerBound(), wrtParameter.getParameterUpperBound(), tolerance)
                + "\n" + HessianWrtParameterProvider.getReportAndCheckForError(this, wrtParameter == WrtParameter.NODE_HEIGHT ? null : tolerance);
        return output;
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return wrtParameter.getDiagonalHessianLogDensity(likelihood, treeIndex);
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented.");
    }

    public enum WrtParameter {
        LOG_POP_SIZE {
            @Override
            double[] getGradientLogDensity(SkyGlideLikelihood likelihood, int treeIndex) {
                return likelihood.getGradientWrtLogPopulationSize();
            }

            @Override
            double[] getDiagonalHessianLogDensity(SkyGlideLikelihood likelihood, int treeIndex) {
                return likelihood.getDiagonalHessianLogDensityWrtLogPopSize();
            }

            @Override
            double getParameterLowerBound() {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            double getParameterUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
        },
        NODE_HEIGHT {
            @Override
            double[] getGradientLogDensity(SkyGlideLikelihood likelihood, int treeIndex) {
                return likelihood.getGradientWrtNodeHeight(treeIndex);
            }

            @Override
            double[] getDiagonalHessianLogDensity(SkyGlideLikelihood likelihood, int treeIndex) {
                return likelihood.getDiagonalHessianWrtNodeHeight(treeIndex);
            }

            @Override
            double getParameterLowerBound() {
                return 0;
            }

            @Override
            double getParameterUpperBound() {
                return Double.POSITIVE_INFINITY;
            }
        };
        abstract double[] getGradientLogDensity(SkyGlideLikelihood likelihood, int treeIndex);
        abstract double[] getDiagonalHessianLogDensity(SkyGlideLikelihood likelihood, int treeIndex);
        abstract double getParameterLowerBound();
        abstract double getParameterUpperBound();
    }
}
