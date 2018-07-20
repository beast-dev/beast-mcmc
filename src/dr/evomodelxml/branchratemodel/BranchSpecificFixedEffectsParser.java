/*
 * CountableMixtureBranchRatesParser.java
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

import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.branchratemodel.BranchCategoriesParser.ALLOCATION;
import static dr.evomodelxml.branchratemodel.BranchCategoriesParser.CATEGORY;
import static dr.evomodelxml.branchratemodel.BranchCategoriesParser.parseCladeCategories;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Xiang Ji
 */
public class BranchSpecificFixedEffectsParser extends AbstractXMLObjectParser {

    private static final String FIXED_EFFECTS = "fixedEffects";

    public String getParserName() {
        return FIXED_EFFECTS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter allocationParameter = (Parameter) xo.getElementFirstChild(ALLOCATION);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        Parameter coefficients = (Parameter) xo.getChild(Parameter.class);

        CountableBranchCategoryProvider.CladeBranchCategoryModel cladeModel =
                new CountableBranchCategoryProvider.CladeBranchCategoryModel(treeModel, allocationParameter);

            parseCladeCategories(xo, cladeModel);

            List<CountableBranchCategoryProvider> categories = new ArrayList<CountableBranchCategoryProvider>();
            categories.add(cladeModel);

            List<ContinuousBranchValueProvider> values = new ArrayList<ContinuousBranchValueProvider>();

            return new BranchSpecificFixedEffects.Default(
                    xo.getId(),
                    categories, values,
                    coefficients
        );
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides a clock consisting of a mixture of fixed effects and random effects.";
    }

    public Class getReturnType() {
        return BranchSpecificFixedEffects.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class),
            new ElementRule(Parameter.class),
            new ElementRule(LocalClockModelParser.CLADE,
                    new XMLSyntaxRule[]{
                            AttributeRule.newIntegerRule(CATEGORY, false),
                            AttributeRule.newBooleanRule(LocalClockModelParser.INCLUDE_STEM, true, "determines whether or not the stem branch above this clade is included in the siteModel (default false)."),
                            AttributeRule.newBooleanRule(LocalClockModelParser.EXCLUDE_CLADE, true, "determines whether to exclude actual branches of the clade from the siteModel (default false)."),
                            new ElementRule(Taxa.class, "A set of taxa which defines a clade to apply a different site model to"),
                    }, 0, Integer.MAX_VALUE),
    };
}
