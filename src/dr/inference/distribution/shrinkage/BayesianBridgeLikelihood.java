/*
 * MultivariateNormalDistributionModel.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution.shrinkage;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;

import static dr.inferencexml.distribution.shrinkage.BayesianBridgeLikelihoodParser.BAYESIAN_BRIDGE;

/**
 * A model for scaled-mixture-of-normals distributed data.
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Xiang Ji
 */

public abstract class BayesianBridgeLikelihood extends AbstractModelLikelihood implements GradientWrtParameterProvider {

    BayesianBridgeLikelihood(Parameter coefficients,
                             Parameter globalScale,
                             Parameter exponent) {

        super(BAYESIAN_BRIDGE);

        this.coefficients = coefficients;
        this.globalScale = globalScale;
        this.exponent = exponent;

        this.dim = coefficients.getDimension();

        addVariable(coefficients);
        addVariable(globalScale);
        addVariable(exponent);
    }

    abstract double calculateLogLikelihood();

    abstract double[] calculateGradientLogDensity();

    public Parameter getGlobalScale() { return globalScale; }

    public Parameter getExponent() { return exponent; }

    public abstract Parameter getLocalScale();

    @Override
    public double getLogLikelihood() { return calculateLogLikelihood(); }

    @Override
    public double[] getGradientLogDensity() {
        return calculateGradientLogDensity();
    }

    @Override
    public int getDimension() {
        return dim;
    }

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
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
       // no intermediates need to be recalculated...
    }

    @Override
    protected void storeState() {
       // Do nothing
    }

    @Override
    protected void restoreState() {
        // Do nothing
    }

    @Override
    protected void acceptState() {
    } // no additional state needs accepting

    final Parameter coefficients;
    final Parameter globalScale;
    final Parameter exponent;

    final int dim;
}
