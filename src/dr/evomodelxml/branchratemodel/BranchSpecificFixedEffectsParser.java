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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.branchratemodel.ContinuousBranchValueProvider;
import dr.evomodel.branchratemodel.CountableBranchCategoryProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import static dr.evomodelxml.branchratemodel.BranchCategoriesParser.*;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @author Xiang Ji
 */
public class BranchSpecificFixedEffectsParser extends AbstractXMLObjectParser {

    private static final String FIXED_EFFECTS = "fixedEffects";
    private static final String INCLUDE_INTERCEPT = "includeIntercept";

    public String getParserName() {
        return FIXED_EFFECTS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//        Parameter allocationParameter = (Parameter) xo.getElementFirstChild(ALLOCATION);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        Parameter coefficients = (Parameter) xo.getChild(Parameter.class);

        boolean includeIntercept = xo.getAttribute(INCLUDE_INTERCEPT, true);

        CountableBranchCategoryProvider.CladeBranchCategoryModel cladeModel =
                new CountableBranchCategoryProvider.CladeBranchCategoryModel(treeModel,
                        new Parameter.Default(treeModel.getNodeCount() -1));

        parseCladeCategories(xo, cladeModel);

        List<CountableBranchCategoryProvider> categories = new ArrayList<CountableBranchCategoryProvider>();
        categories.add(cladeModel);

        List<ContinuousBranchValueProvider> values = new ArrayList<ContinuousBranchValueProvider>();

        BranchSpecificFixedEffects.Default fixedEffects = new BranchSpecificFixedEffects.Default(
                xo.getId(),
                categories, values,
                coefficients,
                includeIntercept
        );

        double[][] designMatrix = fixedEffects.getDesignMatrix(treeModel);

        Logger.getLogger("dr.evomodel").info("Using a fixed effects model with initial design matrix:\n"
                + annotateDesignMatrix(designMatrix, treeModel));

        return fixedEffects;
    }

    private String annotateDesignMatrix(double[][] matrix, Tree tree) {
        StringBuilder sb = new StringBuilder();

        int offset = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (node != tree.getRoot()) {
                String row = new WrappedVector.Raw(matrix[offset]).toString();
                Taxon taxon = tree.getNodeTaxon(tree.getNode(i));
                String name = (taxon != null) ? taxon.getId() : "(not external)";
                sb.append("\t").append(row).append(" : ").append(name).append("\n");

                ++offset;
            }
        }

        return sb.toString();
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
            AttributeRule.newBooleanRule(INCLUDE_INTERCEPT, true),
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
