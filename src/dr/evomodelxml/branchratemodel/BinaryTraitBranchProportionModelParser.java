/*
 * BinaryTraitBranchProportionModelParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BinaryTraitBranchProportionModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.ContinuousTraitBranchRateModel;
import dr.evomodel.branchratemodel.DiscreteTraitBranchRateModel;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class BinaryTraitBranchProportionModelParser extends AbstractXMLObjectParser {

    public static final String RATE = "rate";
    public static final String RATES = "rates";
    public static final String RELATIVE_RATES = "relativeRates";

    public static final String INDICATORS = "indicators";
    public static final String TRAIT_INDEX = "traitIndex";
    public static final String TRAIT_NAME = "traitName";

    public String getParserName() {
        return BinaryTraitBranchProportionModel.BINARY_TRAIT_BRANCH_PROP_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);

        TreeTraitProvider traitProvider = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);
        DataType dataType = DataTypeUtils.getDataType(xo);

        Parameter rateParameter = new Parameter.Default("absoluteRates",new double[]{1.0, 0.0});
        Parameter indicatorsParameter = null;

        if (xo.getChild(INDICATORS) != null) {
            indicatorsParameter = (Parameter) xo.getElementFirstChild(INDICATORS);
        }

        int traitIndex = xo.getAttribute(TRAIT_INDEX, 1) - 1;
        String traitName = "states";

        Logger.getLogger("dr.evomodel").info("Using discrete trait branch rate model.\n" +
                "\tIf you use this model, please cite:\n" +
                "\t\tDrummond and Suchard (in preparation)");

        DiscreteTraitBranchRateModel traitRateModel;

        if (traitProvider == null) {
            // Use the version that reconstructs the trait using parsimony:
            traitRateModel = new DiscreteTraitBranchRateModel(treeModel, patternList, traitIndex, rateParameter);
        } else {
            if (traitName != null) {
                TreeTrait trait = traitProvider.getTreeTrait(traitName);
                if (trait == null) {
                    throw new XMLParseException("A trait called, " + traitName + ", was not available from the TreeTraitProvider supplied to " + getParserName() + ", with ID " + xo.getId());
                }
                traitRateModel = new DiscreteTraitBranchRateModel(traitProvider, dataType, treeModel, trait, traitIndex, rateParameter);
            } else {
                TreeTrait[] traits = new TreeTrait[dataType.getStateCount()];
                for (int i =  0; i < dataType.getStateCount(); i++) {
                    traits[i] = traitProvider.getTreeTrait(dataType.getCode(i));
                    if (traits[i] == null) {
                        throw new XMLParseException("A trait called, " + dataType.getCode(i) + ", was not available from the TreeTraitProvider supplied to " + getParserName() + ", with ID " + xo.getId());
                    }
                }
                traitRateModel = new DiscreteTraitBranchRateModel(traitProvider, traits, treeModel, rateParameter);
            }
        }
        return new BinaryTraitBranchProportionModel(traitRateModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This Branch Rate Model takes a binary discrete trait reconstruction (provided by a TreeTraitProvider) " +
                        "and gives the proportion of each branch of the tree spent in state 0 " +
                        "based on the child trait of that branch.";
    }

    public Class getReturnType() {
        return DiscreteTraitBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class, "The tree model"),
            new XORRule(
                    new AndRule(
                            new ElementRule(TreeTraitProvider.class, "The trait provider"),
                            new XORRule(
                                    new StringAttributeRule(DataType.DATA_TYPE, "The type of general data",
                                            DataType.getRegisteredDataTypeNames(), false),
                                    new ElementRule(DataType.class))
                    ),
                    new ElementRule(PatternList.class)),
            AttributeRule.newIntegerRule(TRAIT_INDEX, true),
            AttributeRule.newStringRule(TRAIT_NAME, true)
    };

}