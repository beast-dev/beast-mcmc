/*
 * EstimableStemWeightBranchSpecificBranchModelParser.java
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

package dr.evomodelxml.branchmodel;

import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.BranchSpecificBranchModel;
import dr.evomodel.branchmodel.EstimableStemWeightBranchSpecificBranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andy Magee
 */
public class EstimableStemWeightBranchSpecificBranchModelParser extends AbstractXMLObjectParser {

    public static final String BRANCH_SPECIFIC_SUBSTITUTION_MODEL = "estimableStemWeightBranchSpecificSubstitutionModel";
    public static final String CLADE = "clade";
    public static final String EXTERNAL_BRANCHES = "externalBranches";
    public static final String BACKBONE = "backbone";
    public static final String STEM_WEIGHT = "stemWeight";
    private static final String ALLOW_SINGLETON = "allowSingleton";

    public String getParserName() {
        return BRANCH_SPECIFIC_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("\nUsing clade-specific branch model.");

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        EstimableStemWeightBranchSpecificBranchModel branchModel = new EstimableStemWeightBranchSpecificBranchModel(tree, substitutionModel);

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject xoc = (XMLObject) xo.getChild(i);
                if (xoc.getName().equals(CLADE)) {

                    double stemWeight = xoc.getAttribute(STEM_WEIGHT, 0.0);

                    substitutionModel = (SubstitutionModel) xoc.getChild(SubstitutionModel.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);

                    boolean allowSingleton = xoc.getAttribute(ALLOW_SINGLETON, false);
                    if (!allowSingleton && taxonList.getTaxonCount() == 1) {
                        throw new XMLParseException("A clade must be defined by at least two taxa");
                    }

                    try {
                        Parameter stemWeightParameter = null;
                        if (xoc.hasChildNamed(STEM_WEIGHT)) {
                            stemWeightParameter = (Parameter) xoc.getElementFirstChild(STEM_WEIGHT);
                            branchModel.addClade(taxonList, substitutionModel, stemWeightParameter);
                        } else {
                            throw new RuntimeException("stemWeight must be provided as a parameter.");
//                            branchModel.addClade(taxonList, substitutionModel, stemWeight);
                        }

                    } catch (TreeUtils.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                } else if (xoc.getName().equals(EXTERNAL_BRANCHES)) {

                    substitutionModel = (SubstitutionModel) xoc.getChild(SubstitutionModel.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);


                    try {
                        branchModel.addExternalBranches(taxonList, substitutionModel);

                    } catch (TreeUtils.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                } else if (xoc.getName().equals(BACKBONE)) {

                    substitutionModel = (SubstitutionModel) xoc.getChild(SubstitutionModel.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);

                    try {
                        branchModel.addBackbone(taxonList, substitutionModel);

                    } catch (TreeUtils.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                }

            }
        }

        branchModel.checkValidityOfBranchAssignments();

        return branchModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return
                "This element provides a branch model which allows different substitution models" +
                        "on different parts of the tree.";
    }

    public Class getReturnType() {
        return EstimableStemWeightBranchSpecificBranchModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TreeModel.class, "The tree"),
            new ElementRule(SubstitutionModel.class, "The substitution model for branches not explicitly included"),
            new ElementRule(EXTERNAL_BRANCHES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Taxa.class, "A substitution model will be applied to the external branches for these taxa"),
                            new ElementRule(SubstitutionModel.class, "The substitution model"),
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(CLADE,
                    new XMLSyntaxRule[]{
                            new ElementRule(STEM_WEIGHT, Parameter.class, "Parameter defining what proportion of the stem branch to include [0 <= w <= 1].", false),
                            AttributeRule.newBooleanRule(ALLOW_SINGLETON, true),
                            new ElementRule(Taxa.class, "A set of taxa which defines a clade to apply a different site model to"),
                            new ElementRule(SubstitutionModel.class, "The substitution model"),
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(BACKBONE,
                    new XMLSyntaxRule[]{
                            new ElementRule(Taxa.class, "A substitution model will be applied only to " +
                                    "the 'backbone' branches defined by these taxa."),
                            new ElementRule(SubstitutionModel.class, "The substitution model"),
                    }, 0, Integer.MAX_VALUE),
    };

}
