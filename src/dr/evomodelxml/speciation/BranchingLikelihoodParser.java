/*
 * BranchingLikelihoodParser.java
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

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.BranchingLikelihood;
import dr.evomodel.speciation.BranchingModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLSyntaxRule;

/**
 */
public class BranchingLikelihoodParser extends AbstractXMLObjectParser {

    public static final String BRANCHING_LIKELIHOOD = "branchingLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "branchingTree";

    public String getParserName() {
        return BRANCHING_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) {

        XMLObject cxo = xo.getChild(MODEL);
        BranchingModel branchingModel = (BranchingModel) cxo.getChild(BranchingModel.class);

        cxo = xo.getChild(TREE);
        TreeModel treeModel = (TreeModel) cxo.getChild(TreeModel.class);

        return new BranchingLikelihood(treeModel, branchingModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of the tree given the demographic function.";
    }

    public Class getReturnType() {
        return BranchingLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(BranchingModel.class)
            }),
            new ElementRule(TREE, new XMLSyntaxRule[]{
                    new ElementRule(TreeModel.class)
            }),
    };

}
