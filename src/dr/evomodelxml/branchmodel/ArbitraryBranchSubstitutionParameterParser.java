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

import dr.evomodel.branchmodel.ArbitraryBranchSubstitutionParameter;
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
public class ArbitraryBranchSubstitutionParameterParser extends AbstractXMLObjectParser {

    public static final String ARBITRARY_BRANCH_SUBSTITUTION_PARAMETER_MODEL="branchSubstitutionParameterBranchModel";
    private static final String SINGLE_RATE="single_rate_subsitution_model";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Logger.getLogger("dr.evomodel").info("\nUsing branch-specific rate branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        if (!(substitutionModel instanceof ParameterReplaceableSubstitutionModel)) {
            throw new RuntimeException("The substitution model is not parameter replaceable!");
        }
        Parameter branchParameter = (Parameter) xo.getChild(Parameter.class);


        BranchSpecificSubstitutionModelProvider substitutionModelProvider = null;
        ArbitraryBranchSubstitutionParameter branchParameterModel = null;
        if (branchParameter == null || branchParameter.getDimension() == 1) {
            substitutionModelProvider = new BranchSpecificSubstitutionModelProvider.None(substitutionModel);
            branchParameterModel = new ArbitraryBranchSubstitutionParameter(SINGLE_RATE, substitutionModelProvider, branchParameter, tree);
        } else{
            final int numBranch = tree.getNodeCount() - 1;
            if (!(branchParameter.getDimension() == numBranch && branchParameter instanceof CompoundParameter)) {
                throw new RuntimeException("branchSubstitutionParameter miss-specified.");
            }
            //TODO: more generic SubstitutionModel construction
            List<SubstitutionModel> substitutionModelList = new ArrayList<SubstitutionModel>();
            int v = 0;
            for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum){
                if (!tree.isRoot(tree.getNode(nodeNum))) {
                    substitutionModelList.add(((ParameterReplaceableSubstitutionModel) substitutionModel).replaceParameter(
                            ((ParameterReplaceableSubstitutionModel) substitutionModel).getReplaceableParameter(),
                            ((CompoundParameter) branchParameter).getParameter(v)));
//                            new HKY(((CompoundParameter) branchRates.getRateParameter()).getParameter(v), substitutionModel.getFrequencyModel()));
                    v++;
                }
            }
            substitutionModelProvider = new BranchSpecificSubstitutionModelProvider.Default((CompoundParameter) branchParameter, substitutionModelList, tree);
            branchParameterModel = new ArbitraryBranchSubstitutionParameter(ARBITRARY_BRANCH_SUBSTITUTION_PARAMETER_MODEL,
                    substitutionModelProvider, branchParameter, tree);
        }

        return branchParameterModel;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SubstitutionModel.class, "The substitution model throughout the tree."),
                new ElementRule(TreeModel.class, "The tree."),
                new ElementRule(Parameter.class, "Substitution Parameters.")
        };
    }

    @Override
    public String getParserDescription() {
        return "This element represents a branch specific substitution rate model.";
    }

    @Override
    public Class getReturnType() {
        return BranchSpecificSubstitutionModelProvider.class;
    }

    @Override
    public String getParserName() {
        return ARBITRARY_BRANCH_SUBSTITUTION_PARAMETER_MODEL;
    }
}
