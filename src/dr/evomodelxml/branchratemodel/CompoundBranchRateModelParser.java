/*
 * CompoundBranchRateModelParser.java
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

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.CompoundBranchRateModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class CompoundBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String COMPOUND_BRANCH_RATE_MODEL = "compoundBranchRateModel";

    public String getParserName() {
        return COMPOUND_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<BranchRateModel> branchRateModels = new ArrayList<BranchRateModel>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof BranchRateModel) {
                branchRateModels.add((BranchRateModel) xo.getChild(i));
            } else {

                Object rogueElement = xo.getChild(i);

                throw new XMLParseException("An element (" + rogueElement
                        + ") which is not a branchRateModel has been added to a " + COMPOUND_BRANCH_RATE_MODEL + " element");
            }

        }

        Logger.getLogger("dr.evomodel").info("Creating a compound branch rate model of " + branchRateModels.size() + " sub-models");

        return new CompoundBranchRateModel(branchRateModels);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides a strict clock model. " +
                "All branches have the same rate of molecular evolution.";
    }

    public Class getReturnType() {
        return CompoundBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(BranchRateModel.class, "The component branch rate models", 1, Integer.MAX_VALUE),
    };

}
