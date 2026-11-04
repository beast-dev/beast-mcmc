/*
 * DiscreteTraitBranchRateModelParser.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.ContinuousTraitBranchRateModel;
import dr.evomodel.branchratemodel.DiscreteTraitBranchRateModel;
import dr.evomodel.continuous.SampledMultivariateTraitLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class DiscreteTraitBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String RATE = "rate";
    public static final String RATES = "rates";
    public static final String RELATIVE_RATES = "relativeRates";

    public static final String INDICATORS = "indicators";
    public static final String TRAIT_INDEX = "traitIndex";
    public static final String TRAIT_NAME = "traitName";

    public String getParserName() {
        return DiscreteTraitBranchRateModel.DISCRETE_TRAIT_BRANCH_RATE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);

        TreeTraitProvider traitProvider = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);
        DataType dataType = DataTypeUtils.getDataType(xo);

        Parameter rateParameter = null;
        Parameter relativeRatesParameter = null;
        Parameter indicatorsParameter = null;

        if (xo.getChild(RATE) != null) {
            rateParameter = (Parameter) xo.getElementFirstChild(RATE);
        }

        if (xo.getChild(RATES) != null) {
            rateParameter = (Parameter) xo.getElementFirstChild(RATES);
        }

        if (xo.getChild(RELATIVE_RATES) != null) {
            relativeRatesParameter = (Parameter) xo.getElementFirstChild(RELATIVE_RATES);
        }

        if (xo.getChild(INDICATORS) != null) {
            indicatorsParameter = (Parameter) xo.getElementFirstChild(INDICATORS);
        }

        int traitIndex = xo.getAttribute(TRAIT_INDEX, 1) - 1;
        String traitName = "states";

        Logger.getLogger("dr.evomodel").info("Using discrete trait branch rate model.\n" +
                "\tIf you use this model, please cite:\n" +
                "\t\tDrummond and Suchard (in preparation)");

        if (traitProvider == null) {
            // Use the version that reconstructs the trait using parsimony:
            return new DiscreteTraitBranchRateModel(treeModel, patternList, traitIndex, rateParameter);
        } else {
            if (traitName != null) {
                TreeTrait trait = traitProvider.getTreeTrait(traitName);
                if (trait == null) {
                    throw new XMLParseException("A trait called, " + traitName + ", was not available from the TreeTraitProvider supplied to " + getParserName() + ", with ID " + xo.getId());
                }
                if (relativeRatesParameter != null) {
                    return new DiscreteTraitBranchRateModel(traitProvider, dataType, treeModel, trait, traitIndex, rateParameter, relativeRatesParameter, indicatorsParameter);
                } else {
                    return new DiscreteTraitBranchRateModel(traitProvider, dataType, treeModel, trait, traitIndex, rateParameter);
                }
            } else {
                TreeTrait[] traits = new TreeTrait[dataType.getStateCount()];
                for (int i =  0; i < dataType.getStateCount(); i++) {
                    traits[i] = traitProvider.getTreeTrait(dataType.getCode(i));
                    if (traits[i] == null) {
                        throw new XMLParseException("A trait called, " + dataType.getCode(i) + ", was not available from the TreeTraitProvider supplied to " + getParserName() + ", with ID " + xo.getId());
                    }
                }
                return new DiscreteTraitBranchRateModel(traitProvider, traits, treeModel, rateParameter);
            }

        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This Branch Rate Model takes a discrete trait reconstruction (provided by a TreeTraitProvider) and " +
                        "gives the rate for each branch of the tree based on the child trait of " +
                        "that branch. The rates for each trait value are specified in a multidimensional parameter.";
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
            new XORRule(
                    new AndRule(
                            new ElementRule(RATE, Parameter.class, "The over-all rate"),
                            new ElementRule(RELATIVE_RATES, Parameter.class, "The state-specific relative rates")
                    ),
                    new ElementRule(RATES, Parameter.class, "The state-specific rates")),
            AttributeRule.newIntegerRule(TRAIT_INDEX, true),
            AttributeRule.newStringRule(TRAIT_NAME, true)
    };

}