/*
 * NegativeBinomialDistributionModelParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.NegativeBinomialDistributionModel;
import dr.inference.distribution.PoissonDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Reads a Poisson distribution model from a DOM Document element.
 */
public class NegativeBinomialDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NEGATIVE_BINOMIAL_DISTRIBUTION_MODEL = "negativeBinomialDistributionModel";
    public static final String MEAN = "mean";
    public static final String ALPHA = "alpha";

    public String getParserName() {
        return NEGATIVE_BINOMIAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter meanParam, alphaParameter;

        XMLObject cxo = xo.getChild(MEAN);
        if (cxo.getChild(0) instanceof Parameter) {
            meanParam = (Parameter) cxo.getChild(Parameter.class);
        } else {
            meanParam = new Parameter.Default(cxo.getDoubleChild(0));
        }

        cxo = xo.getChild(ALPHA);
        if (cxo.getChild(0) instanceof Parameter) {
            alphaParameter = (Parameter) cxo.getChild(Parameter.class);
        } else {
            alphaParameter = new Parameter.Default(cxo.getDoubleChild(0));
        }

        return new NegativeBinomialDistributionModel(meanParam, alphaParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            ),
            new ElementRule(ALPHA,
                    new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(Parameter.class),
                                    new ElementRule(Double.class)
                            )}
            )
    };

    public String getParserDescription() {
        return "Describes a Negative Binomial distribution with a given mean and scale " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return NegativeBinomialDistributionModel.class;
    }

}
