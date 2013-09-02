/*
 * CountableMixtureBranchRatesParser.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.branchratemodel.CountableMixtureBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class CountableMixtureBranchRatesParser extends AbstractXMLObjectParser {

    public static final String COUNTABLE_CLOCK_BRANCH_RATES = "countableMixtureBranchRates";
    public static final String RATES = "rates";
    public static final String ALLOCATION = "rateCategories";

    public String getParserName() {
        return COUNTABLE_CLOCK_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);
        Parameter allocationParameter = (Parameter) xo.getElementFirstChild(ALLOCATION);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

//        final int length = treeModel.getC

        Logger.getLogger("dr.evomodel").info("Using a countable mixture molecular clock model.");

        return new CountableMixtureBranchRates(treeModel, ratesParameter, allocationParameter);
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
        return CountableMixtureBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(RATES, Parameter.class, "The molecular evolutionary rate parameter", false),
            new ElementRule(ALLOCATION, Parameter.class, "Allocation parameter", false),
    };
}
