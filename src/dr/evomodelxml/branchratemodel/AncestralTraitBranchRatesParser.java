/*
 * DiscretizedBranchRatesParser.java
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

package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.AncestralTraitBranchRateModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.AncestralTraitTreeModel;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class AncestralTraitBranchRatesParser extends AbstractXMLObjectParser {

    private static final String ANCESTRAL_TRAIT_BRANCH_RATES = "ancestralTraitBranchRates";

    public String getParserName() {
        return ANCESTRAL_TRAIT_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
        AncestralTraitTreeModel treeModel = (AncestralTraitTreeModel) xo.getChild(AncestralTraitTreeModel.class);

        return new AncestralTraitBranchRateModel(branchRateModel, treeModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns an branch rate model for use with an AncestralTraitTree.";
    }

    public Class getReturnType() {
        return AncestralTraitBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BranchRateModel.class),
            new ElementRule(AncestralTraitTreeModel.class),
    };
}
