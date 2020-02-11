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

import dr.inference.distribution.AutoRegressiveNormalDistributionModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AutoRegressiveNormalDistributionModelParser extends AbstractXMLObjectParser {

    public static final String NORMAL_DISTRIBUTION_MODEL = "autoRegressiveNormalDistributionModel";
    private static final String DIMENSION = "dim";
    private static final String MARGINAL_PRECISION = "scale";
    private static final String DECAY_PRECISION = "rho";

    public String getParserName() {
        return NORMAL_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);

        XMLObject cxo = xo.getChild(MARGINAL_PRECISION);
        Parameter marginal = (Parameter) cxo.getChild(Parameter.class);

        if (marginal.getParameterValue(0) <= 0.0) {
            throw new XMLParseException("Scale must be > 0.0");
        }

        cxo = xo.getChild(DECAY_PRECISION);
        Parameter decay = (Parameter) cxo.getChild(Parameter.class);

        if (Math.abs(decay.getParameterValue(0)) >= 1.0) {
            throw new XMLParseException("|Rho| must be < 1.0");
        }

        return new AutoRegressiveNormalDistributionModel(dim, marginal, decay);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(MARGINAL_PRECISION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(DECAY_PRECISION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };

    public String getParserDescription() {
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() {
        return AutoRegressiveNormalDistributionModel.class;
    }

}
