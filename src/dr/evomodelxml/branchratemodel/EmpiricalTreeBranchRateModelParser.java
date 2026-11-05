/*
 * BranchParameterParser.java
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

import dr.evomodel.branchratemodel.EmpiricalTreeBranchRateModel;
import dr.evomodel.tree.EmpiricalTreeDistributionModel;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Nidia Trovao
 */
public class EmpiricalTreeBranchRateModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "empiricalTreeBranchRates";
    private static final String ATTRIBUTE_NAME = "rateAttribute";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        EmpiricalTreeDistributionModel trees = (EmpiricalTreeDistributionModel)
                xo.getChild(EmpiricalTreeBranchRateModel.class);

        String rateAttributeName = xo.getStringAttribute(ATTRIBUTE_NAME);

        Object obj = trees.getNodeAttribute(trees.getExternalNode(0), rateAttributeName);
        if (!(obj instanceof Double)) {
            throw new XMLParseException("Unable to find double-attributed named '" + rateAttributeName + "'");
        }

        return new EmpiricalTreeBranchRateModel(xo.getId(), trees, rateAttributeName);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(EmpiricalTreeDistributionModel.class),
                AttributeRule.newStringRule(ATTRIBUTE_NAME),
        };
    }

    @Override
    public String getParserDescription() {
        return "Take branch rates from node-attribute in a collection of empirical trees.";
    }

    @Override
    public Class getReturnType() {
        return EmpiricalTreeBranchRateModel.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }
}
