/*
 * BranchRateTransformParser.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchRateTransformParser extends AbstractXMLObjectParser {
    public static final String BRANCH_RATE_TRANFORM="branchRateTransform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return ArbitraryBranchRatesParser.parseTransform(xo);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ArbitraryBranchRatesParser.SCALE, Parameter.class, "optional scale parameter", true),
                new ElementRule(ArbitraryBranchRatesParser.LOCATION, Parameter.class, "optional location parameter", true),
                AttributeRule.newBooleanRule(ArbitraryBranchRatesParser.RECIPROCAL, true),
                AttributeRule.newBooleanRule(ArbitraryBranchRatesParser.EXP, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "Construct BranchRateTransform object.";
    }

    @Override
    public Class getReturnType() {
        return ArbitraryBranchRates.BranchRateTransform.class;
    }

    @Override
    public String getParserName() {
        return BRANCH_RATE_TRANFORM;
    }
}
