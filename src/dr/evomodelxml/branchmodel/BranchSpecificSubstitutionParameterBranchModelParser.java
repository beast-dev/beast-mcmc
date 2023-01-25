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

import dr.evomodel.branchmodel.BranchSpecificSubstitutionParameterBranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.ParameterReplaceableSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificSubstitutionParameterBranchModelParser extends AbstractXMLObjectParser {
    private static final String OLD_PARAMETER = "substitutionParameters";
    private static final String NEW_PARAMETER = "branchRateModels";
    private static final String BRANCH_SPECIFIC_SUBSTITUTION_PARAMETER_BRANCH_MODEL = "branchSpecificSubstitutionParameterBranchModel";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Parameter> parameters = new ArrayList<>();
        XMLObject xoc = xo.getChild(OLD_PARAMETER);
        for (int i = 0; i < xoc.getChildCount(); ++i) {
            Parameter parameter = (Parameter) xoc.getChild(i);
            parameters.add(parameter);
        }

        List<BranchRateModel> branchRateModels = new ArrayList<>();
        xoc = xo.getChild(NEW_PARAMETER);
        for (int i = 0; i < xoc.getChildCount(); ++i) {
            BranchRateModel branchRateModel = (BranchRateModel) xoc.getChild(i);
            branchRateModels.add(branchRateModel);
        }

        ParameterReplaceableSubstitutionModel substitutionModel = (ParameterReplaceableSubstitutionModel) xo.getChild(ParameterReplaceableSubstitutionModel.class);
        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        return new BranchSpecificSubstitutionParameterBranchModel(BRANCH_SPECIFIC_SUBSTITUTION_PARAMETER_BRANCH_MODEL,
                parameters, branchRateModels, substitutionModel, tree);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ParameterReplaceableSubstitutionModel.class, "The substitution model throughout the tree."),
                new ElementRule(OLD_PARAMETER, Parameter.class, "Parameters of the substitution model to be replaced.",
                        0, Integer.MAX_VALUE),
                new ElementRule(NEW_PARAMETER, BranchRateModel.class, "BranchRateModels to store the branchSpecific parameters.",
                        0, Integer.MAX_VALUE),
                new ElementRule(TreeModel.class)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return BranchSpecificSubstitutionParameterBranchModel.class;
    }

    @Override
    public String getParserName() {
        return BRANCH_SPECIFIC_SUBSTITUTION_PARAMETER_BRANCH_MODEL;
    }
}
