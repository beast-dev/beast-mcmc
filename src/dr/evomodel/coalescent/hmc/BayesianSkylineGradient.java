/*
 * BayesianSkylineGradient.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.hmc;

import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BayesianSkylineGradient implements
        GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    private final BayesianSkylineLikelihood likelihood;
    private final WrtParameter wrtParameter;

    public BayesianSkylineGradient(BayesianSkylineLikelihood likelihood,
                                   WrtParameter wrtParameter) {
        this.likelihood = likelihood;
        this.wrtParameter = wrtParameter;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return wrtParameter.getParameter(likelihood);
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return wrtParameter.getGradientLogDensity(likelihood);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return wrtParameter.getDiagonalHessianLogDensity(likelihood);
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, wrtParameter.getParameterLowerBound(), Double.POSITIVE_INFINITY, null);
    }

    public enum WrtParameter {
        NODE_HEIGHT("nodeHeight") {
            @Override
            Parameter getParameter(BayesianSkylineLikelihood likelihood) {
                return null;
            }

            @Override
            double[] getGradientLogDensity(BayesianSkylineLikelihood likelihood) {
                return new double[0];
            }

            @Override
            double[] getDiagonalHessianLogDensity(BayesianSkylineLikelihood likelihood) {
                return new double[0];
            }

            @Override
            double getParameterLowerBound() {
                return 0;
            }

            @Override
            public void getWarning(BayesianSkylineLikelihood likelihood) {

            }
        };

        private final String name;

        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(BayesianSkylineLikelihood likelihood);

        abstract double[] getGradientLogDensity(BayesianSkylineLikelihood likelihood);

        abstract double[] getDiagonalHessianLogDensity(BayesianSkylineLikelihood likelihood);

        abstract double getParameterLowerBound();

        public abstract void getWarning(BayesianSkylineLikelihood likelihood);

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
