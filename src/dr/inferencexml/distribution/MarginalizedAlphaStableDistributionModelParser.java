/*
 * NormalDistributionModelParser.java
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

import dr.inference.distribution.MarginalizedAlphaStableDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class MarginalizedAlphaStableDistributionModelParser extends AbstractXMLObjectParser {

    private static final String MARGINALIZED_ALPHA_STABLE_DISTRIBUTION_MODEL = "marginalizedAlphaStableDistributionModel";
    private static final String SCALE = "scale";
    private static final String ALPHA = "alpha";

    public String getParserName() {
        return MARGINALIZED_ALPHA_STABLE_DISTRIBUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter scaleParam = parseParameter(xo, SCALE);
        Parameter alphaParam = parseParameter(xo, ALPHA);

        String name = xo.hasId() ? xo.getId() : MARGINALIZED_ALPHA_STABLE_DISTRIBUTION_MODEL;

        return new MarginalizedAlphaStableDistributionModel(name, scaleParam, alphaParam);
    }

    private static Parameter parseParameter(XMLObject xo, String name) throws XMLParseException {
        XMLObject cxo = xo.getChild(name);
        return (cxo.getChild(0) instanceof Parameter) ?
                (Parameter) cxo.getChild(Parameter.class) :
                new Parameter.Default(cxo.getDoubleChild(0));
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static ElementRule parameterRules(String name) {
        return new ElementRule(name,
                new XMLSyntaxRule[]{
                        new XORRule(
                                new ElementRule(Parameter.class),
                                new ElementRule(Double.class)
                        )}
        );
    }

    private final XMLSyntaxRule[] rules = {
            parameterRules(SCALE),
            parameterRules(ALPHA),
    };

    public String getParserDescription() {
        return null;
    }

    public Class getReturnType() {
        return MarginalizedAlphaStableDistributionModel.class;
    }
}
