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

import dr.evomodel.branchmodel.ArbitrarySubstitutionParameterBranchModel;
import dr.evomodel.substmodel.BranchSpecificSubstitutionModelProvider;
import dr.evomodel.substmodel.ParameterReplaceableSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
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

        List<CompoundParameter> parameterList = new ArrayList<CompoundParameter>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            parameterList.add((CompoundParameter) cxo.getChild(i));
        }

        //TODO: more generic SubstitutionModel construction
        List<SubstitutionModel> substitutionModelList = new ArrayList<SubstitutionModel>();

        ParameterReplaceableSubstitutionModel rootSubstitutionModel = (ParameterReplaceableSubstitutionModel) substitutionModel;
        for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum) {

            ParameterReplaceableSubstitutionModel branchSubstitutionModel = (ParameterReplaceableSubstitutionModel) substitutionModel;


            for (int i = 0; i < dxo.getChildCount(); i++) {
                Parameter oldParameter = (Parameter) dxo.getChild(i);
                Parameter branchParameter = (Parameter) cxo.getChild(i);


                if (!(branchParameter.getDimension() == tree.getNodeCount() && branchParameter instanceof CompoundParameter)) {
                    throw new RuntimeException("branchSubstitutionParameter miss-specified.");
                }

                CompoundParameter branchSpecificParameter = (CompoundParameter) branchParameter;

                if (tree.getRoot().getNumber() != tree.getNodeCount() - 1) {
                    throw new RuntimeException("Root node number is not the maximum.");
                }

                branchSubstitutionModel = branchSubstitutionModel.factory(
                        oldParameter,
                        branchSpecificParameter.getParameter(nodeNum));


            }
            substitutionModelList.add(branchSubstitutionModel);
        }
        for (int i = 0; i < dxo.getChildCount(); i++) {
            Parameter oldParameter = (Parameter) dxo.getChild(i);
            CompoundParameter branchSpecificParameter = (CompoundParameter) cxo.getChild(i);
            Parameter rootParameter = new Parameter.Default( "root.parameter",
                    branchSpecificParameter.getParameter(0).getParameterValue(0),
                    branchSpecificParameter.getParameter(0).getBounds().getLowerLimit(0),
                    branchSpecificParameter.getParameter(0).getBounds().getUpperLimit(0));
            rootSubstitutionModel = rootSubstitutionModel.factory(oldParameter, rootParameter);
        }
        substitutionModelList.add(rootSubstitutionModel);

        substitutionModelProvider = new BranchSpecificSubstitutionModelProvider.Default(substitutionModelList, tree);

        branchParameterModel = new ArbitrarySubstitutionParameterBranchModel(ARBITRARY_SUBSTITUTION_PARAMETER_BRANCH_MODEL,
                substitutionModelProvider, parameterList, tree);

        return branchParameterModel;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SubstitutionModel.class, "The substitution model throughout the tree."),
                new ElementRule(TreeModel.class, "The tree."),
//                new ElementRule(Parameter.class, "Substitution Parameters."),
                new ElementRule(BRANCH_SPECIFIC_PARAMETER,
                        new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, "The branch-specific substitution parameter.", 1, Integer.MAX_VALUE),
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
