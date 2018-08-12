/*
 * BranchSpecificSubstitutionRateModelParser.java
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

package dr.evomodelxml.branchmodel;

import dr.evomodel.branchmodel.BranchSpecificRateBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchSpecificSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificRateBranchModelParser extends AbstractXMLObjectParser {

    public static final String BRANCH_SPECIFIC_SUBSTITUTION_RATE_MODEL="branchSpecificRateBranchModel";
    private static final String SINGLE_RATE="single_rate_subsitution_model";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Logger.getLogger("dr.evomodel").info("\nUsing branch-specific rate branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        ArbitraryBranchRates branchRates = (ArbitraryBranchRates) xo.getAttribute("arbitraryBranchRates", null);


        BranchSpecificSubstitutionModel branchSubstitutionModels = null;
        if (branchRates == null) {
            branchSubstitutionModels = new BranchSpecificSubstitutionModel.None(substitutionModel);
        }

        BranchSpecificRateBranchModel rateBranchModel = new BranchSpecificRateBranchModel(SINGLE_RATE, branchSubstitutionModels);
        return rateBranchModel;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SubstitutionModel.class, "The substitution model throughout the tree."),
                new ElementRule(TreeModel.class, "The tree."),
                new ElementRule("rateRule", ArbitraryBranchRates.class, "Branch-specific rates.", true)
        };
    }

    @Override
    public String getParserDescription() {
        return "This element represents a branch specific substitution rate model.";
    }

    @Override
    public Class getReturnType() {
        return BranchSpecificSubstitutionModel.class;
    }

    @Override
    public String getParserName() {
        return BRANCH_SPECIFIC_SUBSTITUTION_RATE_MODEL;
    }
}
