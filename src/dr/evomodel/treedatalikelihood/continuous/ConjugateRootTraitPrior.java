/*
 * ConjugateRootTraitPrior.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import static dr.evomodelxml.treedatalikelihood.ContinuousDataLikelihoodParser.CONJUGATE_ROOT_PRIOR;

/**
 * @author Marc A. Suchard
 */
public class ConjugateRootTraitPrior extends AbstractModel {

    private static final String PRIOR_SAMPLE_SIZE = AbstractMultivariateTraitLikelihood.PRIOR_SAMPLE_SIZE;

    private final Parameter mean;
    private final Parameter pseudoObservations;

    public ConjugateRootTraitPrior(Parameter mean,
                                   Parameter pseudoObservations) {
        super("ConjugateRootTraitPrior");
        this.mean = mean;
        this.pseudoObservations = pseudoObservations;

        addVariable(mean);
        addVariable(pseudoObservations);
    }

    public double[] getMean() {
        return mean.getParameterValues();
    }

    public Parameter getMeanParameter() {
        return mean;
    }

    public double getPseudoObservations() {
        return pseudoObservations.getParameterValue(0);
    }

    public static ConjugateRootTraitPrior parseConjugateRootTraitPrior(XMLObject xo,
                                                                       final int dim) throws XMLParseException {

        XMLObject cxo = xo.getChild(CONJUGATE_ROOT_PRIOR);

        Parameter meanParameter = (Parameter) cxo.getChild(MultivariateDistributionLikelihood.MVN_MEAN)
                .getChild(Parameter.class);

        if (meanParameter.getDimension() != dim) {
            throw new XMLParseException("Root prior mean dimension (" + meanParameter.getDimension() +
                    ") does not match trait diffusion dimension (" + dim + ")");
        }

        Parameter sampleSizeParameter = (Parameter) cxo.getChild(PRIOR_SAMPLE_SIZE).getChild(Parameter.class);

        return new ConjugateRootTraitPrior(meanParameter, sampleSizeParameter);
    }

    public static final XMLSyntaxRule[] rules = {
            new ElementRule(MultivariateDistributionLikelihood.MVN_MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(PRIOR_SAMPLE_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("No submodels");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged(variable);
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { }

    @Override
    protected void acceptState() { }
}
