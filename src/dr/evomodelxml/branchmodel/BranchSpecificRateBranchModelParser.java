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

import dr.evomodel.branchmodel.BranchSpecificRateBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.substmodel.BranchSpecificSubstitutionModelProvider;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Xiang Ji
 */
public class BranchSpecificRateBranchModelParser extends AbstractXMLObjectParser {

    public static final String BRANCH_SPECIFIC_SUBSTITUTION_RATE_MODEL="branchSpecificRateBranchModel";
    private static final String SINGLE_RATE="single_rate_subsitution_model";
    private static final String BRANCH_SPECIFIC_RATE="branch_specific_rate_subsitution_model";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Logger.getLogger("dr.evomodel").info("\nUsing branch-specific rate branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        ArbitraryBranchRates branchRates = (ArbitraryBranchRates) xo.getChild(ArbitraryBranchRates.class);


        BranchSpecificSubstitutionModelProvider substitutionModelProvider = null;
        BranchSpecificRateBranchModel rateBranchModel = null;
        if (branchRates == null || branchRates.getRateParameter().getDimension() == 1) {
            substitutionModelProvider = new BranchSpecificSubstitutionModelProvider.None(substitutionModel);
            rateBranchModel = new BranchSpecificRateBranchModel(SINGLE_RATE, substitutionModelProvider);
        } else{
            final int numBranch = tree.getNodeCount() - 1;
            if (!(branchRates.getRateParameter().getDimension() == numBranch && branchRates.getRateParameter() instanceof CompoundParameter)) {
                throw new RuntimeException("rateParameter miss-specified.");
            }
            //TODO: more generic SubstitutionModel construction
            List<SubstitutionModel> substitutionModelList = new ArrayList<SubstitutionModel>();
            int v = 0;
            for (int nodeNum = 0; nodeNum < tree.getNodeCount(); ++nodeNum){
                if (!tree.isRoot(tree.getNode(nodeNum))) {
                    substitutionModelList.add(new HKY(((CompoundParameter) branchRates.getRateParameter()).getParameter(v), substitutionModel.getFrequencyModel()));
                    v++;
                }
            }
            substitutionModelProvider = new BranchSpecificSubstitutionModelProvider.Default(branchRates, substitutionModelList, tree);
            rateBranchModel = new BranchSpecificRateBranchModel(BRANCH_SPECIFIC_SUBSTITUTION_RATE_MODEL, substitutionModelProvider);
        }

        return rateBranchModel;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SubstitutionModel.class, "The substitution model throughout the tree."),
                new ElementRule(TreeModel.class, "The tree."),
                new ElementRule("rateRule", ArbitraryBranchRates.class, "Branch-specific rates.", true)
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
        return BRANCH_SPECIFIC_SUBSTITUTION_RATE_MODEL;
    }
}
