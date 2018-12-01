/*
 * BranchSpecificCompoundParameterParser.java
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
import dr.evomodel.tree.TreeModel;
import dr.inference.model.BranchParameter;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificCompoundParameterParser extends AbstractXMLObjectParser {
    public static final String BRANCH_SPECIFIC_COMPOUND_PARAMETER = "branchSpecificCompoundParameter";
    public static final String BRANCH_RATE_TRANFORM="branchRateTransform";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        CompoundParameter compoundParameter = (CompoundParameter) xo.getChild(CompoundParameter.class);
        final int numNodes = treeModel.getNodeCount();
        if (compoundParameter.getDimension() != numNodes) {
            throw new RuntimeException("Dimension mismatch!");
        }

        ArbitraryBranchRates.BranchRateTransform transform = (ArbitraryBranchRates.BranchRateTransform) xo.getChild(ArbitraryBranchRates.BranchRateTransform.class);
        if (transform == null) {
            transform = new ArbitraryBranchRates.BranchRateTransform.None();
        }

        BranchParameter branchParameter = new BranchParameter(
                xo.getName(),
                compoundParameter,
                treeModel,
                transform);

        List<BranchParameter.IndividualBranchParameter> individualBranchParameterList = new ArrayList<BranchParameter.IndividualBranchParameter>();
        for (int i = 0; i < numNodes; i++) {
            BranchParameter.IndividualBranchParameter individualBranchParameter =
                    new BranchParameter.IndividualBranchParameter(branchParameter, i, compoundParameter.getParameter(i));
            individualBranchParameterList.add(individualBranchParameter);
        }
        branchParameter.addTransformedParameterList(individualBranchParameterList);
        return branchParameter;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(CompoundParameter.class),
                new ElementRule(TreeModel.class),
                new ElementRule(ArbitraryBranchRates.BranchRateTransform.class, true),
        };
    }

    @Override
    public String getParserDescription() {
        return "A multidimensional parameter constructed for each branch of the tree.";
    }

    @Override
    public Class getReturnType() {
        return BranchParameter.class;
    }

    @Override
    public String getParserName() {
        return BRANCH_SPECIFIC_COMPOUND_PARAMETER;
    }
}
