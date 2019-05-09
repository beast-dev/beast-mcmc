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

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.BranchParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificCompoundParameterParser extends AbstractXMLObjectParser {

    public static final String BRANCH_SPECIFIC_COMPOUND_PARAMETER = "branchSpecificCompoundParameter";
    public static final String BRANCH_PARAMETER = "branchParameter";
    public static final String ROOT_PARAMETER = "rootParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        CompoundParameter compoundParameter = (CompoundParameter) xo.getChild(BRANCH_PARAMETER).getChild(CompoundParameter.class);
        Parameter rootParameter = (Parameter) xo.getChild(ROOT_PARAMETER).getChild(Parameter.class);
        final int numNodes = treeModel.getNodeCount();
        if (compoundParameter.getDimension() != numNodes - 1) {
            throw new RuntimeException("Dimension mismatch!");
        }
        if (rootParameter.getDimension() != 1) {
            throw new RuntimeException("Root parameter dimension should be one.");
        }
//
//        ArbitraryBranchRates.BranchRateTransform transform = (ArbitraryBranchRates.BranchRateTransform) xo.getChild(ArbitraryBranchRates.BranchRateTransform.class);
//        if (transform == null) {
//            transform = new ArbitraryBranchRates.BranchRateTransform.None();
//        }
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        BranchParameter branchParameter = new BranchParameter(
                xo.getName(),
                treeModel,
                branchRateModel,
                rootParameter
                );

        List<BranchParameter.IndividualBranchParameter> individualBranchParameterList = new ArrayList<BranchParameter.IndividualBranchParameter>();
        BranchParameter.IndividualBranchParameter rootBranchParameter = null;
        int v = 0;
        for (int i = 0; i < numNodes; i++) {
            NodeRef node = treeModel.getNode(i);
            BranchParameter.IndividualBranchParameter individualBranchParameter;
            if (treeModel.isRoot(node)) {
                rootBranchParameter =
//                        new BranchParameter.IndividualBranchParameter(branchParameter, numNodes - 1, rootParameter);
                        new BranchParameter.IndividualBranchParameter(branchParameter, numNodes - 1, rootParameter);
            } else {
                individualBranchParameter =
//                        new BranchParameter.IndividualBranchParameter(branchParameter, v, compoundParameter.getParameter(v));
                        new BranchParameter.IndividualBranchParameter(branchParameter, v, rootParameter);
                v++;
                individualBranchParameterList.add(individualBranchParameter);
            }
        }
        assert(rootBranchParameter != null);
        individualBranchParameterList.add(rootBranchParameter);
        branchParameter.addTransformedParameterList(individualBranchParameterList);
        return branchParameter;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(BRANCH_PARAMETER,
                        new XMLSyntaxRule[]{
                                new ElementRule(BranchRateModel.class, "The branch-specific substitution parameter.", 1, Integer.MAX_VALUE),
                        }),
                new ElementRule(TreeModel.class),
                new ElementRule(ROOT_PARAMETER,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class, "The root substitution parameter.", 1, 1),
                        }),
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
