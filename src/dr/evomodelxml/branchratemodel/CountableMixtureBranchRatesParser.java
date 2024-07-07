/*
 * CountableMixtureBranchRatesParser.java
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

import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.branchratemodel.CountableMixtureBranchRates;
import dr.evomodel.branchratemodel.CountableModelMixtureBranchRates;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import static dr.evomodelxml.branchratemodel.BranchCategoriesParser.parseCladeCategories;

/**
 */
public class CountableMixtureBranchRatesParser extends AbstractXMLObjectParser {

    public static final String COUNTABLE_CLOCK_BRANCH_RATES = "countableMixtureBranchRates";
    public static final String RATES = "rates";
    public static final String ALLOCATION = "rateCategories";
    public static final String CATEGORY = "category";
    public static final String RANDOMIZE = "randomize";
    public static final String RANDOM_EFFECTS = "randomEffects";
    public static final String FIXED_EFFECTS = "fixedEffects";
    public static final String IN_LOG_SPACE = "inLogSpace";
    public static final String TIME_COEFFICIENT = "branchTimeEffect";

    public String getParserName() {
        return COUNTABLE_CLOCK_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter = null;

        List<AbstractBranchRateModel> fixedEffects = null;
        if (xo.hasChildNamed(FIXED_EFFECTS)) {
            XMLObject cxo = xo.getChild(FIXED_EFFECTS);
            fixedEffects = new ArrayList<AbstractBranchRateModel>();
            for (int i = 0; i < cxo.getChildCount(); ++i) {
                fixedEffects.add((AbstractBranchRateModel)cxo.getChild(i));
            }
        } else {
            ratesParameter = (Parameter) xo.getElementFirstChild(RATES);
        }

        Parameter timeCoefficient;
        if (xo.hasChildNamed(TIME_COEFFICIENT)){
            timeCoefficient = (Parameter) xo.getElementFirstChild(TIME_COEFFICIENT);
        } else {
            timeCoefficient = null;
        }

        Parameter allocationParameter = (Parameter) xo.getElementFirstChild(ALLOCATION);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        List<AbstractBranchRateModel> randomEffects = null;
        if (xo.hasChildNamed(RANDOM_EFFECTS)) {
            XMLObject cxo = xo.getChild(RANDOM_EFFECTS);
            randomEffects = new ArrayList<AbstractBranchRateModel>();
            for (int i = 0; i < cxo.getChildCount(); ++i) {
                randomEffects.add((AbstractBranchRateModel)cxo.getChild(i));
            }
        }

        boolean inLogSpace = xo.getAttribute(IN_LOG_SPACE, false);

        Logger.getLogger("dr.evomodel").info("Using a countable mixture molecular clock model.");

        CountableBranchCategoryProvider.BranchCategoryModel cladeModel;

        if (!xo.getAttribute(RANDOMIZE, true)) {
            CountableBranchCategoryProvider.CladeBranchCategoryModel cm = new
                    CountableBranchCategoryProvider.CladeBranchCategoryModel(treeModel, allocationParameter);
            parseCladeCategories(xo, cm);
            cladeModel = cm;
        } else {
            CountableBranchCategoryProvider.IndependentBranchCategoryModel cm = new CountableBranchCategoryProvider.IndependentBranchCategoryModel(treeModel, allocationParameter);
            cm.randomize();
            cladeModel = cm;
        }

        if (fixedEffects != null) {
            return new CountableModelMixtureBranchRates(cladeModel, treeModel, fixedEffects, randomEffects, inLogSpace);
        } else {
            return new CountableMixtureBranchRates(cladeModel, treeModel, ratesParameter, timeCoefficient, randomEffects, inLogSpace);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a clock consisting of a mixture of fixed effects and random effects.";
    }

    public Class getReturnType() {
        return CountableMixtureBranchRates.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new XORRule(
                    new ElementRule(RATES, Parameter.class, "The molecular evolutionary rate parameter", false),
                    new ElementRule(FIXED_EFFECTS,
                            new XMLSyntaxRule[] {
                                    new ElementRule(AbstractBranchRateModel.class, 0, Integer.MAX_VALUE),
                            },
                            "Fixed effects", false)
            ),
            new ElementRule(ALLOCATION, Parameter.class, "Allocation parameter", false),
            new ElementRule(TIME_COEFFICIENT, Parameter.class, "The coefficient for the branch time-dependent rate function", true),
            new ElementRule(RANDOM_EFFECTS,
                    new XMLSyntaxRule[] {
                            new ElementRule(AbstractBranchRateModel.class, 0, Integer.MAX_VALUE),
                    },
                    "Possible random effects", true),
            AttributeRule.newBooleanRule(IN_LOG_SPACE, true),
            AttributeRule.newBooleanRule(RANDOMIZE, true),
            new ElementRule(LocalClockModelParser.CLADE,
                    new XMLSyntaxRule[]{
//                            AttributeRule.newBooleanRule(RELATIVE, true),
                            AttributeRule.newIntegerRule(CATEGORY, false),
                            AttributeRule.newBooleanRule(LocalClockModelParser.INCLUDE_STEM, true, "determines whether or not the stem branch above this clade is included in the siteModel (default false)."),
                            AttributeRule.newBooleanRule(LocalClockModelParser.EXCLUDE_CLADE, true, "determines whether to exclude actual branches of the clade from the siteModel (default false)."),
                            new ElementRule(Taxa.class, "A set of taxa which defines a clade to apply a different site model to"),
                    }, 0, Integer.MAX_VALUE),
    };
}
