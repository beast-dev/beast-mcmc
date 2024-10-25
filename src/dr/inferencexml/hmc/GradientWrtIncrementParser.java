/*
 * GradientWrtIncrementParser.java
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

package dr.inferencexml.hmc;

import dr.evomodel.speciation.SpeciationLikelihoodGradient;
import dr.evomodel.treedatalikelihood.discrete.SequenceDistanceStatistic;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.CompoundDerivative;
import dr.inference.hmc.GradientWrtIncrement;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andy Magee
 * @author Yucai Shao
 */
public class GradientWrtIncrementParser extends AbstractXMLObjectParser {

    public final static String GRAD_INCREMENTS = "gradientWrtIncrements1D";
    public final static String INCREMENT_TRANSFORM = "incrementTransformType";

    @Override
    public String getParserName() {
        return GRAD_INCREMENTS;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GradientWrtParameterProvider gradient = (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        double upper = Double.POSITIVE_INFINITY;
        double lower = Double.NEGATIVE_INFINITY;
        if( xo.hasAttribute("upper") && xo.hasAttribute("lower")) {
            upper = xo.getDoubleAttribute("upper");
            lower = xo.getDoubleAttribute("lower");
        }

        Transform incrementTransform = null;
        String ttype = (String) xo.getAttribute(INCREMENT_TRANSFORM);
        if (ttype.equalsIgnoreCase("log")) {
            incrementTransform = Transform.LOG;
        } else if (ttype.equalsIgnoreCase("logit")) {
            incrementTransform = new Transform.ScaledLogitTransform(upper, lower);
        } else if (ttype.equalsIgnoreCase("none")) {
            incrementTransform = new Transform.NoTransform();
        } else {
            throw new RuntimeException("Invalid transform type");
        }

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new GradientWrtIncrement(gradient, parameter, incrementTransform);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GradientWrtParameterProvider.class,"The gradient with respect to the (non-increment) parameters."),
            new ElementRule(Parameter.class,"The increments."),
            new StringAttributeRule(INCREMENT_TRANSFORM,"The transformation to get from the parameter to the increments, f in increments[i] = f(parameter[i]) - f(parameter[i - 1]).",false)
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