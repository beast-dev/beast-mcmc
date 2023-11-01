/*
 * MultivariateNormalDistributionModelParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.distribution;

import dr.inference.distribution.GaussianMarkovRandomFieldModel2;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.GaussianMarkovRandomField2;
import dr.xml.*;

public class GaussianMarkovRandomFieldParser2 extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "gaussianMarkovRandomField2";
    private static final String DIMENSION = "dim";
    private static final String PRECISION = "precision";
    private static final String START = "start";

    public String getParserName() {

        return NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter coefficients = (Parameter) xo.getChild(Parameter.class);

        int dim = coefficients.getDimension();

        XMLObject cxo = xo.getChild(PRECISION);
        Parameter incrementPrecision = (Parameter) cxo.getChild(Parameter.class);

        if (incrementPrecision.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Scale must be > 0.0");
        }

        cxo = xo.getChild(START);
        Parameter start = (Parameter) cxo.getChild(Parameter.class);



        return new GaussianMarkovRandomFieldModel2(coefficients, new GaussianMarkovRandomField2(dim, incrementPrecision, start));
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
//            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(PRECISION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(START,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return GaussianMarkovRandomFieldModel2.class;
    }

}
