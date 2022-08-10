/*
 * CompoundGradientParser.java
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

package dr.inferencexml.hmc;

import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.evomodel.treedatalikelihood.discrete.SequenceDistanceStatistic;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.CompoundDerivative;
import dr.inference.hmc.GradientWrtIncrement;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andy Magee
 */
public class GradientWrtIncrementParser extends AbstractXMLObjectParser {

    public final static String GRAD_INCREMENTS = "gradientWrtIncrements";
    public final static String INCREMENT_TRANSFORM = "incrementTransformType";

    @Override
    public String getParserName() {
        return GRAD_INCREMENTS;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GradientWrtParameterProvider gradient = (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM);
        GradientWrtIncrement.IncrementTransformType type = GradientWrtIncrement.IncrementTransformType.factory(ttype);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new GradientWrtIncrement(gradient, parameter, type);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class,"The gradient with respect to the (non-increment) parameters."),
            new ElementRule(Parameter.class,"The increments."),
            new StringAttributeRule(INCREMENT_TRANSFORM,"The transformation to get from the parameter to the increments, f in increments[i] = f(parameter[i]) - f(parameter[i - 1])")
    };

    @Override
    public String getParserDescription() {
        return "Takes a gradient with respect to a vector-parameter and returns the gradient with respect to the increments increment[i] = f(parameter[i]) - f(parameter[i - 1]).";
    }

    @Override
    public Class getReturnType() {
        return GradientWrtIncrement.class;
    }
}
