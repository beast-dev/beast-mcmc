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

import dr.evolution.tree.NodeRef;
import dr.evomodel.branchmodel.ArbitrarySubstitutionParameterBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.BranchSpecificSubstitutionModelProvider;
import dr.evomodel.substmodel.ParameterReplaceableSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.BranchParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class ArbitrarySubstitutionParameterBranchModelParser extends AbstractXMLObjectParser {

    private static final String ARBITRARY_SUBSTITUTION_PARAMETER_BRANCH_MODEL ="arbitrarySubstitutionParameterBranchModel";
    private static final String SINGLE_RATE="replacedParameter";
    private static final String BRANCH_SPECIFIC_PARAMETER = "branchSpecificParameter";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Logger.getLogger("dr.evomodel").info("\nUsing branch-specific substitution parameter branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        if (!(substitutionModel instanceof ParameterReplaceableSubstitutionModel)) {
            throw new RuntimeException("The substitution model is not parameter replaceable!");
        }

        XMLObject cxo = xo.getChild(BRANCH_SPECIFIC_PARAMETER);
        XMLObject dxo = xo.getChild(SINGLE_RATE);


        BranchSpecificSubstitutionModelProvider substitutionModelProvider;
        ArbitrarySubstitutionParameterBranchModel branchParameterModel;

        assert (dxo.getChildCount() == cxo.getChildCount());

        List<BranchParameter> parameterList = new ArrayList<BranchParameter>();
        for (int i = 0; i < cxo.getChildCount(); i++) {

            Parameter rootParameter = (Parameter) dxo.getChild(i);
            BranchRateModel branchRateModel = (BranchRateModel) cxo.getChild(i);
            BranchParameter branchParameter = new BranchParameter("branchSpecific." + rootParameter.getId(),
                    tree,
                    branchRateModel,
                    rootParameter);
            parameterList.add(branchParameter);
        }

        List<SubstitutionModel> substitutionModelList = new ArrayList<SubstitutionModel>();

        ParameterReplaceableSubstitutionModel rootSubstitutionModel = (ParameterReplaceableSubstitutionModel) substitutionModel;
        List<Parameter> oldParameters = parseParameters(dxo);
//        List<Parameter> branchParameters = parseParameters(cxo);

        final int parameterCount = oldParameters.size();
        int v = 0;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum) {

            NodeRef node = tree.getNode(nodeNum);

            ParameterReplaceableSubstitutionModel branchSubstitutionModel = (ParameterReplaceableSubstitutionModel) substitutionModel;

            List<Parameter> newParameters = new ArrayList<Parameter>();

            if (tree.isRoot(node)) {
                for (int i = 0; i < parameterCount; i++) {
                    BranchParameter branchParameter = (BranchParameter) parameterList.get(i);
                    newParameters.add(branchParameter.getParameter(tree.getNodeCount() - 1));
                }
            } else {
                for (int i = 0; i < parameterCount; i++) {
                    BranchParameter branchParameter = (BranchParameter) parameterList.get(i);
                    newParameters.add(branchParameter.getParameter(v));
                }
                v++;
            }

            branchSubstitutionModel = branchSubstitutionModel.factory(
                    oldParameters,
                    newParameters);

            substitutionModelList.add(branchSubstitutionModel);
        }

        substitutionModelProvider = new BranchSpecificSubstitutionModelProvider.Default(substitutionModelList, tree);

        branchParameterModel = new ArbitrarySubstitutionParameterBranchModel(ARBITRARY_SUBSTITUTION_PARAMETER_BRANCH_MODEL,
                substitutionModelProvider, parameterList, tree);

        return branchParameterModel;
    }

    private List<Parameter> parseParameters(XMLObject cxo) {
        List<Parameter> parameters = new ArrayList<Parameter>();

        for (int i = 0; i < cxo.getChildCount(); i++) {
            parameters.add((Parameter) cxo.getChild(i));
        }

        return parameters;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SubstitutionModel.class, "The substitution model throughout the tree."),
                new ElementRule(BRANCH_SPECIFIC_PARAMETER,
                        new XMLSyntaxRule[]{
                                new ElementRule(BranchRateModel.class, "The branch-specific substitution parameter handled by BranchRateModels.", 1, Integer.MAX_VALUE),
                        }),
                new ElementRule(SINGLE_RATE,
                        new XMLSyntaxRule[]{
                                new ElementRule(Parameter.class, "The substitution parameter to be replaced.", 1, Integer.MAX_VALUE),
                        })
        };
    }

    @Override
    public String getParserDescription() {
        return "This element represents a branch specific substitution rate model.";
    }

    @Override
    public Class getReturnType() {
        return ArbitrarySubstitutionParameterBranchModel.class;
    }

    @Override
    public String getParserName() {
        return ARBITRARY_SUBSTITUTION_PARAMETER_BRANCH_MODEL;
    }
}
