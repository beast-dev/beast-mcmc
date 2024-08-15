/*
 * SkyGlideGradientParser.java
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

package dr.evomodelxml.coalescent.smooth;

import dr.evomodel.coalescent.smooth.SkyGlideGradient;
import dr.evomodel.coalescent.smooth.SkyGlideLikelihood;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser;
import dr.xml.*;


/**
 * @author Mathieu Fourment
 * @author Erick Matsen
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SkyGlideGradientParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "skyGlideGradient";
    private static final String GRADIENT_CHECK_TOLERANCE = HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        SkyGlideLikelihood likelihood = (SkyGlideLikelihood) xo.getChild(SkyGlideLikelihood.class);
        double gradientCheckTolerance = xo.getAttribute(GRADIENT_CHECK_TOLERANCE, 1E-3);

        return new SkyGlideGradient(likelihood, parameter, gradientCheckTolerance);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(SkyGlideLikelihood.class),
            new ElementRule(Parameter.class)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return SkyGlideGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
