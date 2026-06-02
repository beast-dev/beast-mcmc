/*
 * StrictClockBranchRatesParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class StrictClockBranchRatesParser extends AbstractXMLObjectParser {

    public static final String STRICT_CLOCK_BRANCH_RATES = "strictClockBranchRates";
    public static final String RATE = "rate";

    public String getParserName() {
        return STRICT_CLOCK_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);

        Logger.getLogger("dr.evomodel").info("\nUsing strict molecular clock model.");

        return new StrictClockBranchRates(rateParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a strict clock model. " +
                        "All branches have the same rate of molecular evolution.";
    }

    public Class getReturnType() {
        return StrictClockBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(RATE, Parameter.class, "The molecular evolutionary rate parameter", false),
    };
}
