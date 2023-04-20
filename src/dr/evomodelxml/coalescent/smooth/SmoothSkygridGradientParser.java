/*
 * SmoothSkygridGradientParser.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.coalescent.smooth;


import dr.evomodel.coalescent.smooth.SmoothSkygridGradient;
import dr.evomodel.coalescent.smooth.SmoothSkygridLikelihood;
import dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser;
import dr.xml.*;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SmoothSkygridGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "smoothSkygridGradient";

    private static final String WRT_PARAMETER = "wrtParameter";

    private static final String TOLERANCE = HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        SmoothSkygridLikelihood skygridLikelihood = (SmoothSkygridLikelihood) xo.getChild(SmoothSkygridLikelihood.class);

        String wrtParameterCase = (String) xo.getAttribute(WRT_PARAMETER);

        SmoothSkygridGradient.WrtParameter type = SmoothSkygridGradient.WrtParameter.factory(wrtParameterCase);

        double tolerance = xo.getAttribute(TOLERANCE, 1E-4);

        return new SmoothSkygridGradient(skygridLikelihood, type, tolerance);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SmoothSkygridLikelihood.class),
            AttributeRule.newStringRule(WRT_PARAMETER),
            AttributeRule.newDoubleRule(TOLERANCE, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SmoothSkygridGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
