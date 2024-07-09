/*
 * AttributeBranchRateModelParser.java
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

import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.AttributeBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class AttributeBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String RATE_ATTRIBUTE_NAME = "rateAttribute";

    public String getParserName() {
        return AttributeBranchRateModel.ATTRIBUTE_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);


        final String rateAttributeName = (xo.hasAttribute(RATE_ATTRIBUTE_NAME) ?
                xo.getStringAttribute(RATE_ATTRIBUTE_NAME) : null);

        return new AttributeBranchRateModel(tree, rateAttributeName);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a branch rate model." +
                "The branch rates are specified by an attribute embedded in the nodes of the tree.";
    }

    public Class getReturnType() {
        return AttributeBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new StringAttributeRule(RATE_ATTRIBUTE_NAME,
                    "Optional name of a rate attribute to be read with the trees")
    };


}
