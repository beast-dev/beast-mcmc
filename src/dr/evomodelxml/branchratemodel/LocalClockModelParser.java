/*
 * LocalClockModelParser.java
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

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.LocalClockModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class LocalClockModelParser extends AbstractXMLObjectParser {

    public static final String LOCAL_CLOCK_MODEL = "localClockModel";
    public static final String RATE = BranchRateModel.RATE;
    public static final String RELATIVE = "relative";
    public static final String CLADE = "clade";
    public static final String INCLUDE_STEM = "includeStem";
    public static final String STEM_PROPORTION = "stemProportion";
    public static final String EXCLUDE_CLADE = "excludeClade";
    public static final String EXTERNAL_BRANCHES = "externalBranches";
    public static final String TRUNK = "trunk";
    public static final String INDEX = "index";


    public String getParserName() {
        return LOCAL_CLOCK_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

        Parameter globalRateParameter = (Parameter) xo.getElementFirstChild(RATE);
        LocalClockModel localClockModel = new LocalClockModel(tree, globalRateParameter);

        for (int i = 0; i < xo.getChildCount(); i++) {
            if (xo.getChild(i) instanceof XMLObject) {

                XMLObject xoc = (XMLObject) xo.getChild(i);
                if (xoc.getName().equals(CLADE)) {

                    boolean relative = xoc.getAttribute(RELATIVE, false);

                    Parameter rateParameter = (Parameter) xoc.getChild(Parameter.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);

                    if (taxonList.getTaxonCount() == 1) {
                        throw new XMLParseException("A local clock for a clade must be defined by at least two taxa");
                    }

                    boolean includeStem = false;
                    boolean excludeClade = false;
                    double stemProportion = 0.0;

                    if (xoc.hasAttribute(INCLUDE_STEM)) {
                        includeStem = xoc.getBooleanAttribute(INCLUDE_STEM);
                        // if includeStem=true then assume it is the whole stem
                        stemProportion = 1.0;
                    }

                    if (xoc.hasAttribute(STEM_PROPORTION)) {
                        stemProportion = xoc.getDoubleAttribute(STEM_PROPORTION);
                        if (stemProportion < 0.0 || stemProportion > 1.0) {
                            throw new XMLParseException("A stem proportion should be between 0, 1");
                        }
                    }

                    if (xoc.hasAttribute(EXCLUDE_CLADE)) {
                        excludeClade = xoc.getBooleanAttribute(EXCLUDE_CLADE);
                    }

                    try {
                        localClockModel.addCladeClock(taxonList, rateParameter, relative, stemProportion, excludeClade);

                    } catch (Tree.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                } else if (xoc.getName().equals(EXTERNAL_BRANCHES)) {

                    boolean relative = xoc.getAttribute(RELATIVE, false);

                    Parameter rateParameter = (Parameter) xoc.getChild(Parameter.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);


                    try {
                        localClockModel.addExternalBranchClock(taxonList, rateParameter, relative);

                    } catch (Tree.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                } else if (xoc.getName().equals(TRUNK)) {

                    boolean relative = xoc.getAttribute(RELATIVE, false);

                    Parameter indexParameter = null;
                    if (xoc.hasChildNamed(INDEX)) {
                        indexParameter = (Parameter) xoc.getElementFirstChild(INDEX);
                    }
                    Parameter rateParameter = (Parameter) xoc.getChild(Parameter.class);
                    TaxonList taxonList = (TaxonList) xoc.getChild(TaxonList.class);


                    try {
                        localClockModel.addTrunkClock(taxonList, rateParameter, indexParameter, relative);

                    } catch (Tree.MissingTaxonException mte) {
                        throw new XMLParseException("Taxon, " + mte + ", in " + getParserName() + " was not found in the tree.");
                    }
                }

            }
        }

        System.out.println("Using local clock branch rate model.");

        return localClockModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a branch rate model that adds a delta to each terminal branch length.";
    }

    public Class getReturnType() {
        return LocalClockModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(RATE, Parameter.class, "The molecular evolutionary rate parameter", false),
            new ElementRule(EXTERNAL_BRANCHES,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(RELATIVE, true),
                            new ElementRule(Taxa.class, "A local clock that will be applied only to " +
                                    "the external branches for these taxa"),
                            new ElementRule(Parameter.class, "The rate parameter"),
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(CLADE,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(RELATIVE, true),
                            AttributeRule.newBooleanRule(INCLUDE_STEM, true, "determines whether or not the stem branch above this clade is included in the siteModel (default false)."),
                            AttributeRule.newDoubleRule(STEM_PROPORTION, true, "proportion of stem to include in clade rate (default 0)."),
                            AttributeRule.newBooleanRule(EXCLUDE_CLADE, true, "determines whether to exclude actual branches of the clade from the siteModel (default false)."),
                            new ElementRule(Taxa.class, "A set of taxa which defines a clade to apply a different site model to"),
                            new ElementRule(Parameter.class, "The rate parameter")
                    }, 0, Integer.MAX_VALUE),
            new ElementRule(TRUNK,
                    new XMLSyntaxRule[]{
                            AttributeRule.newBooleanRule(RELATIVE, true),
                            new ElementRule(Taxa.class, "A local clock that will be applied only to " +
                                    "the 'trunk' branches defined by these taxa"),
                            new ElementRule(Parameter.class, "The rate parameter"),
                            new ElementRule(INDEX, Parameter.class, "The trunk taxon index", true),
                    }, 0, Integer.MAX_VALUE),
    };
}
