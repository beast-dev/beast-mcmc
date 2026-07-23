/*
 * BooleanLikelihoodParser.java
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

package dr.inferencexml.model;

import dr.inference.model.BooleanLikelihood;
import dr.inference.model.BooleanStatistic;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

/**
 * Reads a distribution likelihood from a DOM Document element.
 */
public class BooleanLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BOOLEAN_LIKELIHOOD = "booleanLikelihood";

    public static final String DATA = "data";

    public String getParserName() { return BOOLEAN_LIKELIHOOD; }

    public Object parseXMLObject(XMLObject xo) {

        BooleanLikelihood likelihood = new BooleanLikelihood();

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof BooleanStatistic) {
                likelihood.addData( (BooleanStatistic)xo.getChild(i));
            }
        }

        return likelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A function that log likelihood of a set of boolean statistics. "+
                "If all the statistics are true then it returns 0.0 otherwise -infinity.";
    }

    public Class getReturnType() { return BooleanLikelihood.class; }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
        new ElementRule(BooleanStatistic.class, 1, Integer.MAX_VALUE )
    };

}
