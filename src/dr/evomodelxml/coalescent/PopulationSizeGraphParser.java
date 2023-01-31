/*
 * PopulationSizeGraphParser.java
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

import dr.evomodel.coalescent.PopulationSizeGraph;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.xml.*;

/**
 */
@Deprecated
public class PopulationSizeGraphParser extends AbstractXMLObjectParser {

    public static String POPGRAPH_STATISTIC = "popGraph";

    public String getParserName() {
        return POPGRAPH_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Object child = xo.getChild(0);

        if (child instanceof VariableDemographicModel) {
            final double dim = xo.getDoubleAttribute("time");
            return new PopulationSizeGraph((VariableDemographicModel) child, dim);
        }

        throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a statistic that is the population size at evenly spaced intervals over tree.";
    }

    public Class getReturnType() {
        return PopulationSizeGraph.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(VariableDemographicModel.class, 1, 1),
            AttributeRule.newDoubleRule("time")
    };
}
