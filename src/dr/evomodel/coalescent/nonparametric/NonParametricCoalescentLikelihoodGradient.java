/*
 * NonParametricCoalescentLikelihoodGradient.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.evomodel.coalescent.nonparametric;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.Reportable;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 * @author Xiang Ji
 */

public class NonParametricCoalescentLikelihoodGradient extends AbstractModel implements
        GradientWrtParameterProvider, Reportable {

    private final NonParametricCoalescentLikelihood likelihood;
    private final Parameter parameter;
    private final GradientProvider provider;
    private boolean gradientKnown;
    private double[] gradient;

    public NonParametricCoalescentLikelihoodGradient(NonParametricCoalescentLikelihood likelihood,
                                                     Parameter parameter) {
        super("NonparametricCoalescentLikelihoodGradient");

        this.likelihood = likelihood;
        this.parameter = parameter;

        addVariable(parameter);
        addModel(likelihood);

        if (parameter == likelihood.getLogPopSizes()) {
            provider = new GradientProvider() {
                @Override
                public int getDimension() {
                    return parameter.getDimension();
                }

                @Override
                public double[] getGradientLogDensity(Object x) {
                    return likelihood.getGradientLogDensity(x);
                }
            };
        }
        else {
            throw new IllegalArgumentException("Not yet implemented");
//            provider = likelihood.getGradientWrt(parameter);
        }
        this.gradientKnown = false;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY, null);
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
//        if(!gradientKnown) {
            gradient = provider.getGradientLogDensity(parameter.getParameterValues());
//            gradientKnown = true;
//        }
        return gradient;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
//        if (model instanceof MultilocusNonparametricCoalescentLikelihood) {
//            gradientKnown = false;
//        } else {
//            throw new RuntimeException("Unknown object");
//        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        gradientKnown = false;
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
}
