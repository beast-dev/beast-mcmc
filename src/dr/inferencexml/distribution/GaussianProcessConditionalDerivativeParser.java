/*
 * GaussianProcessPredictionParser.java
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

package dr.inferencexml.distribution;

import dr.inference.model.Parameter;
import dr.math.distributions.gp.AdditiveGaussianProcessDistribution;
import dr.math.distributions.gp.GaussianProcessConditionalDerivative;
import dr.xml.*;

import java.util.Arrays;

public class GaussianProcessConditionalDerivativeParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianProcessConditionalDerivative";
    private static final String FIELD = "field";
    private static final String PRIORMEAN = "priorMean";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AdditiveGaussianProcessDistribution gp = (AdditiveGaussianProcessDistribution)
                xo.getChild(AdditiveGaussianProcessDistribution.class);

        Parameter field = (Parameter) xo.getChild(FIELD).getChild(Parameter.class);

        Parameter priorMean;
        if (xo.hasChildNamed(PRIORMEAN)) {
            priorMean = (Parameter) xo.getChild(PRIORMEAN).getChild(Parameter.class);
        } else {
                double[] zeros = new double[field.getDimension()];
                Arrays.fill(zeros, 0.0);
                priorMean = new Parameter.Default(zeros);
        }

        return new GaussianProcessConditionalDerivative(gp, field, priorMean);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(AdditiveGaussianProcessDistribution.class),
            new ElementRule(FIELD, Parameter.class),
            new ElementRule(PRIORMEAN, Parameter.class, "Specify the prior mean of the gp derivative", true)
    };

    public String getParserDescription() { // TODO update
        return "Describes the derivative of a gaussian process";
    }

    public Class getReturnType() { return GaussianProcessConditionalDerivative.class; }
}
