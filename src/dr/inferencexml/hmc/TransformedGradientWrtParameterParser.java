/*
 * TransformedGradientWrtParameterParser.java
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

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.TransformedGradientWrtParameter;
import dr.inference.model.GradientProvider;
import dr.inference.model.Parameter;
import dr.inference.model.ReciprocalLikelihood;
import dr.inference.model.TransformedParameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */

public class TransformedGradientWrtParameterParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "transformedGradient";
    private static final String WRT = "wrt";
    private static final String JACOBIAN = "includeJacobian";
    private static final String INVERSE = "inverse";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GradientWrtParameterProvider gradient = (GradientWrtParameterProvider)
                xo.getChild(GradientWrtParameterProvider.class);

        TransformedParameter parameter = (TransformedParameter) xo.getChild(TransformedParameter.class);

        boolean includeJacobian = xo.getAttribute(JACOBIAN, true);
        boolean inverse = xo.getAttribute(INVERSE, false);

        ReciprocalLikelihood reciprocalLikelihood = (ReciprocalLikelihood) xo.getChild(ReciprocalLikelihood.class);

        if (xo.hasChildNamed(WRT)) {

            Parameter wrt = (Parameter) xo.getElementFirstChild(WRT);

            if (wrt != parameter.getUntransformedParameter()) {
                throw new XMLParseException("Mismatch between transformed and untransformed parameters");
            }
        }

        return new TransformedGradientWrtParameter(gradient, parameter, reciprocalLikelihood,
                includeJacobian, inverse);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new XORRule(
                        new ElementRule(GradientWrtParameterProvider.class),
                        new ElementRule(GradientProvider.class)
                ),
                new ElementRule(TransformedParameter.class),
                new ElementRule(ReciprocalLikelihood.class, true),
                new ElementRule(WRT, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                }, true),
                AttributeRule.newBooleanRule(JACOBIAN, true),
                AttributeRule.newBooleanRule(INVERSE, true),
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return GradientWrtParameterProvider.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}