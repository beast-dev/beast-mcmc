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

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.hmc.CompoundDerivative;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.NumericalGradient;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andy Magee
 */
public class NumericalGradientParser extends AbstractXMLObjectParser {

    public static final String NUM_GRAD = "numericalGradient";

    @Override
    public String getParserName() {
        return NUM_GRAD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Likelihood likelihood = (Likelihood)xo.getChild(Likelihood.class);
        Parameter parameter = (Parameter)xo.getChild(Parameter.class);
        return new NumericalGradient(likelihood, parameter);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Likelihood.class, false),
            new ElementRule(Parameter.class, false)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return NumericalGradient.class;
    }
}
