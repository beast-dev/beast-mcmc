/*
 * StructuredCoalescentLikelihoodParser.java
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

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.BeagleTreeLikelihoodParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Guy Baele
 */
public class StructuredCoalescentLikelihoodParser extends AbstractXMLObjectParser {

    public static final String STRUCTURED_COALESCENT = "structuredCoalescent";
    public static final String TYPE = "type";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";
    public static final String SUBINTERVALS = "subIntervals";
    private static final String THREADS = "threads";
    public static final String USE_AMBIGUITIES = "useAmbiguities";

    public static final String MAP_RECONSTRUCTION = "useMAP";

    public static final Boolean USE_OLD_CODE = false;
    private static final boolean USE_DELEGATE = true;
    private static final boolean USE_BEAGLE = true;
    private static final boolean TRANSPOSE = true;

    public String getParserName() {
        return STRUCTURED_COALESCENT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TaxonList includeSubtree = null;

        if (xo.hasChildNamed(INCLUDE)) {
            includeSubtree = (TaxonList) xo.getElementFirstChild(INCLUDE);
        }

        List<TaxonList> excludeSubtrees = new ArrayList<TaxonList>();

        if (xo.hasChildNamed(EXCLUDE)) {
            XMLObject cxo = xo.getChild(EXCLUDE);
            for (int i = 0; i < cxo.getChildCount(); i++) {
                excludeSubtrees.add((TaxonList) cxo.getChild(i));
            }
        }

        int subIntervals = xo.getAttribute(SUBINTERVALS, 1);

        if (subIntervals != 1) {
            throw new XMLParseException("The number of sub-intervals currently has to be set to 1.");
        }

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, false);

        boolean useMAP = xo.getAttribute(MAP_RECONSTRUCTION, false);

        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        DataType dataType = patternList.getDataType();
        String tag = xo.getAttribute(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG);
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GeneralSubstitutionModel generalSubstitutionModel = (GeneralSubstitutionModel) xo.getChild(GeneralSubstitutionModel.class);
        Parameter popSizes = (Parameter) xo.getChild(Parameter.class);

        if (popSizes.getDimension() != 1) {
            if (generalSubstitutionModel.getDataType().getStateCount() != popSizes.getDimension()) {
                throw new XMLParseException("Mismatch between rate matrix and deme count.");
            }
        }

        int threads = xo.getAttribute(THREADS, 1);

        if (treeModel != null) {
            try {
                if (USE_OLD_CODE) {
                    return new OldStructuredCoalescentLikelihood(treeModel, branchRateModel, popSizes, patternList,
                            generalSubstitutionModel, subIntervals, includeSubtree, excludeSubtrees);
                } else {
                    if (USE_DELEGATE) {
                        final BastaLikelihoodDelegate delegate;
                        if (USE_BEAGLE) {
                            delegate = new BeagleBastaLikelihoodDelegate("name", treeModel,
                                    generalSubstitutionModel.getDataType().getStateCount(), TRANSPOSE);
                        } else {
                            delegate = (threads != 1) ?
                                    new ParallelBastaLikelihoodDelegate("name", treeModel,
                                            generalSubstitutionModel.getDataType().getStateCount(), threads, TRANSPOSE) :
                                    new GenericBastaLikelihoodDelegate("name", treeModel,
                                            generalSubstitutionModel.getDataType().getStateCount(), TRANSPOSE);
                        }
                        return new BastaLikelihood("name", treeModel, patternList, generalSubstitutionModel,
                                popSizes, branchRateModel, delegate, dataType, tag, useMAP, subIntervals, useAmbiguities);
                    } else {
                        return new FasterStructuredCoalescentLikelihood(treeModel, branchRateModel, popSizes, patternList,
                                dataType, tag, generalSubstitutionModel, subIntervals, includeSubtree, excludeSubtrees,
                                useMAP);
                    }
                }
            } catch (TreeUtils.MissingTaxonException mte) {
                throw new XMLParseException("treeModel missing a taxon from taxon list in " + getParserName() + " element");
            }
        } else {
            throw new XMLParseException("");
        }

    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A Bayesian structured coalescent approximation model.";
    }

    public Class getReturnType() {
        return OldStructuredCoalescentLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(BeagleTreeLikelihoodParser.USE_AMBIGUITIES, true),
            AttributeRule.newStringRule(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, true),
            AttributeRule.newIntegerRule(SUBINTERVALS, true),
            AttributeRule.newIntegerRule(THREADS, true),
            new ElementRule(PatternList.class),
            new ElementRule(TreeModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(Parameter.class, true)
    };

}
