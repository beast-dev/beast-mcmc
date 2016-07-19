/*
 * AncestralStateTreeLikelihoodParser.java
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

package dr.oldevomodelxml.treelikelihood;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.xml.*;

/**
 */
@Deprecated // Switching to BEAGLE
public class AncestralStateTreeLikelihoodParser extends AbstractXMLObjectParser {

    public static final String RECONSTRUCTING_TREE_LIKELIHOOD = "ancestralTreeLikelihood";
    public static final String RECONSTRUCTION_TAG = AncestralStateTreeLikelihood.STATES_KEY;
    public static final String RECONSTRUCTION_TAG_NAME = "stateTagName";
    public static final String MAP_RECONSTRUCTION = "useMAP";
    public static final String MARGINAL_LIKELIHOOD = "useMarginalLikelihood";

    public String[] getParserNames() {
        return new String[]{
                getParserName(), "beast_" + getParserName()
        };
    }

    public String getParserName() {
        return RECONSTRUCTING_TREE_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        throw new XMLParseException("Ancestral state functionality is only support when using the BEAGLE library.\nAvailable from http://github.com/beagle-dev/beagle-lib/");

//        boolean useAmbiguities = xo.getAttribute(TreeLikelihoodParser.USE_AMBIGUITIES, false);
//        boolean storePartials = xo.getAttribute(TreeLikelihoodParser.STORE_PARTIALS, true);
//
//        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
//        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
//        SiteModel siteModel = (SiteModel) xo.getChild(SiteModel.class);
//
//        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
//
//        DataType dataType = ((SubstitutionModel) xo.getChild(SubstitutionModel.class)).getDataType();
//
//        boolean useMAP = xo.getAttribute(MAP_RECONSTRUCTION, false);
//        boolean useMarginalLogLikelihood = xo.getAttribute(MARGINAL_LIKELIHOOD, true);
//
//        // default tag is RECONSTRUCTION_TAG
//        String tag = xo.getAttribute(RECONSTRUCTION_TAG_NAME, RECONSTRUCTION_TAG);
//
//        boolean forceRescaling = xo.getAttribute(TreeLikelihoodParser.FORCE_RESCALING, false);
//
//        return new AncestralStateTreeLikelihood(patternList, treeModel, siteModel,
//                branchRateModel, useAmbiguities, storePartials, dataType, tag, forceRescaling,
//                useMAP, useMarginalLogLikelihood);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents the likelihood of a patternlist on a tree given the site model.";
    }

    public Class getReturnType() {
        return AncestralStateTreeLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            AttributeRule.newBooleanRule(TreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newBooleanRule(TreeLikelihoodParser.STORE_PARTIALS, true),
            AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
            AttributeRule.newBooleanRule(TreeLikelihoodParser.FORCE_RESCALING, true),
            AttributeRule.newBooleanRule(MAP_RECONSTRUCTION, true),
            AttributeRule.newBooleanRule(MARGINAL_LIKELIHOOD, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(SiteModel.class),
            new ElementRule(BranchRateModel.class, true)
    };
}
