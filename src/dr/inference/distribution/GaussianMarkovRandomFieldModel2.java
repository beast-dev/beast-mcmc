/*
 * MultivariateNormalDistributionModel.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.*;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.distributions.GaussianMarkovRandomField2;
import dr.math.distributions.GaussianProcessRandomGenerator;


/**
 * A class that acts as a model for gaussian random walk
 *
 * @author Marc Suchard
 * Pratyusa Datta
 */

public class GaussianMarkovRandomFieldModel2 extends AbstractModelLikelihood implements
        GradientWrtParameterProvider, HessianWrtParameterProvider {

    public GaussianMarkovRandomFieldModel2(Parameter coefficients,
                                          GaussianMarkovRandomField2 distribution) {
        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.coefficients = coefficients;
        this.distribution = distribution;
        this.dim = coefficients.getDimension();

        addModel(distribution);
        addVariable(coefficients);
    }


    public Parameter getincrementPrecision() { return distribution.getincrementPrecision(); }

    public Parameter getstart() { return distribution.getstart(); }
    public double[][] getScaleMatrix() {
        return distribution.getScaleMatrix();
    }


    public double[] getMean() {
        return distribution.getMean();
    }

    public String getType() {
        return distribution.getType();
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************
    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return coefficients;
    }

    @Override
    public final void makeDirty() {
        // Do nothing
    }

    @Override
    public final void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    @Override
    public final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // no intermediates need to be recalculated...
    }

    @Override
    public void storeState() {
        // Do nothing
    }

    @Override
    public void restoreState() {
        // Do nothing
    }

    @Override
    public void acceptState() {
    } // no additional state needs accepting


    @Override
    public int getDimension() {
        return dim;
    }


    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    public Parameter getIncrementPrecision() { return distribution.getincrementPrecision(); }

    public Parameter getStart() { return distribution.getstart(); }

    public double getLogLikelihood() {
        return distribution.logPdf(coefficients.getParameterValues());
    }

    public double[] getGradientLogDensity() {
        return distribution.gradLogPdf(coefficients.getParameterValues());
    }


    // **************************************************************
    // Private instance variables and functions
    // **************************************************************




    private final Parameter coefficients;
    private final GaussianMarkovRandomField2 distribution;
    private final int dim;


    @Override
    public double[] getDiagonalHessianLogDensity() {
        return new double[0];
    }

    @Override
    public double[][] getHessianLogDensity() {
        return new double[0][];
    }
}
