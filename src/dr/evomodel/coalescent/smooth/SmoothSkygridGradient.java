/*
 * SmoothSkygridGradient.java
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

import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author Yuwei Bao
 * @author Marc A. Suchard
 */
public class SmoothSkygridGradient implements GradientWrtParameterProvider, Reportable {

    private final SmoothSkygridLikelihood skygridLikelihood;

    private final WrtParameter wrtParameter;

    private final Parameter parameter;

    private final Double tolerance;

    public SmoothSkygridGradient(SmoothSkygridLikelihood skygridLikelihood,
                                 WrtParameter wrtParameter,
                                 Double tolerance) {
        this.skygridLikelihood = skygridLikelihood;
        this.wrtParameter = wrtParameter;
        this.tolerance = tolerance;
        this.parameter = wrtParameter.getParameter(skygridLikelihood);
    }

    @Override
    public Likelihood getLikelihood() {
        return skygridLikelihood;
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
        return wrtParameter.getGradientLogDensity(skygridLikelihood);
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, tolerance, 1E-3);
    }

    public enum WrtParameter {

        NODE_HEIGHT("nodeHeight") {
            @Override
            Parameter getParameter(SmoothSkygridLikelihood skygridLikelihood) {
                TreeModel treeModel = (TreeModel) skygridLikelihood.getTree(0);
                return new NodeHeightProxyParameter("allInternalNode", treeModel, true);
            }

            @Override
            double[] getGradientLogDensity(SmoothSkygridLikelihood likelihood) {
                return likelihood.getGradientWrtNodeHeight();
            }
        };

        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(SmoothSkygridLikelihood skygridLikelihood);

        abstract double[] getGradientLogDensity(SmoothSkygridLikelihood likelihood);

        private final String name;

        public static WrtParameter factory(String match) {
            for (WrtParameter type : WrtParameter.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }

    }

}
