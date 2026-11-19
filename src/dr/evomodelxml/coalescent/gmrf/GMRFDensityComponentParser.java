/*
 * GMRFDensityComponentParser.java
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

package dr.evomodelxml.coalescent.gmrf;

import dr.evomodel.coalescent.gmrf.GMRFDensityComponent;
import dr.evomodel.coalescent.gmrf.GMRFSkyrideLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class GMRFDensityComponentParser extends AbstractXMLObjectParser {

    public String getParserName() {
        return GMRFDensityComponent.PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        GMRFSkyrideLikelihood skyride = (GMRFSkyrideLikelihood) xo.getChild(GMRFSkyrideLikelihood.class);
        String componentString = (String) xo.getAttribute(GMRFDensityComponent.COMPONENT);
        boolean returnCoalescent;
        if (componentString.compareToIgnoreCase(GMRFDensityComponent.COALESCENT_TAG) == 0) {
            returnCoalescent = true;
        } else if (componentString.compareToIgnoreCase(GMRFDensityComponent.GMRF_TAG) == 0) {
            returnCoalescent = false;
        } else {
            throw new XMLParseException("Unknown component of GMRF Skyride");
        }

        return new GMRFDensityComponent(skyride, returnCoalescent);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a single component (coalescent or field) of the GMRF Skyride.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(GMRFDensityComponent.COMPONENT),
            new ElementRule(GMRFSkyrideLikelihood.class),
    };
}
