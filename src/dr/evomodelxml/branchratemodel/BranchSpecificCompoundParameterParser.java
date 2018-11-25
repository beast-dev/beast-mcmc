/*
 * TreeSpecificCompoundParameterParser.java
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

import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.BranchParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificCompoundParameterParser extends AbstractXMLObjectParser {
    public static final String BRANCH_SPECIFIC_COMPOUND_PARAMETER = "branchSpecificCompoundParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CompoundParameter compoundParameter = new CompoundParameter((String) null);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        final int numNodes = treeModel.getNodeCount();
        for (int i = 0; i < numNodes; i++) {
            Parameter branchInnerParameter = new Parameter.Default((String) parameter.getId() + String.valueOf(i + 1), parameter.getParameterValue(0),
                    parameter.getBounds().getLowerLimit(0), parameter.getBounds().getUpperLimit(0));
            compoundParameter.addParameter(branchInnerParameter);
        }
        BranchParameter branchParameter = new BranchParameter(compoundParameter,
                treeModel,
                ArbitraryBranchRatesParser.parseTransform(xo),
                new TreeParameterModel(treeModel, compoundParameter, true, TreeTrait.Intent.BRANCH));

        CompoundParameter resultCompoundParameter = new CompoundParameter(null);
        for (int i = 0; i < numNodes; i++) {
            BranchParameter.IndividualBranchParameter individualBranchParameter = new BranchParameter.IndividualBranchParameter(branchParameter, i, compoundParameter.getParameter(i));
            resultCompoundParameter.addParameter(individualBranchParameter);
        }
        return resultCompoundParameter;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(TreeModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "A multidimensional parameter constructed for each branch of the tree.";
    }

    @Override
    public Class getReturnType() {
        return CompoundParameter.class;
    }

    @Override
    public String getParserName() {
        return BRANCH_SPECIFIC_COMPOUND_PARAMETER;
    }
}
