/*
 * BayesianSkylineGradientParser.java
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

package dr.evomodelxml.coalescent;

import dr.evomodel.coalescent.BayesianSkylineLikelihood;
import dr.evomodel.coalescent.hmc.BayesianSkylineGradient;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BayesianSkylineGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "skylineGradient";
    private static final String WRT_PARAMETER = "wrtParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        BayesianSkylineLikelihood skylineLikelihood = (BayesianSkylineLikelihood) xo.getChild(BayesianSkylineLikelihood.class);
        String wrtParameterCase = (String) xo.getAttribute(WRT_PARAMETER);
        BayesianSkylineGradient.WrtParameter type = BayesianSkylineGradient.WrtParameter.factory(wrtParameterCase);

        return new BayesianSkylineGradient(skylineLikelihood, type);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(WRT_PARAMETER),
            new ElementRule(BayesianSkylineLikelihood.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BayesianSkylineGradient.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
