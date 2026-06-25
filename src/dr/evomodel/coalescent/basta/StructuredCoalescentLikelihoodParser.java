/*
 * StructuredCoalescentLikelihoodParser.java
 *
 * Copyright (c) 2002-2025 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
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

    public static final String MAP_RECONSTRUCTION = "useMAP";
    public static final String MARGINAL_LIKELIHOOD = "useMarginalLikelihood";
    public static final String CONDITIONAL_PROBABILITIES_IN_LOG_SPACE = "conditionalProbabilitiesInLogSpace";
    public static final String RECONSTRUCTION_TAG = "states";
    public static final String RECONSTRUCTION_TAG_NAME = "stateTagName";
    public static final String GROWTH_RATES = "growthRates";
    public static final String POPSIZES = "popSizes";
    public static final String LOG_POP_SIZES = "logPopSizes";
    public static final String ANCHOR_TIME = "anchorTime";
    public static final String ROOT_ANCHORED = "rootAnchored";
    public static final String ANCESTRAL_GROWTH_RATE = "ancestralGrowthRate";
    public static final String ANCHOR_PROPORTION = "anchorProportion";
    public static final String GRID_POINTS = "gridPoints";
    public static final Boolean USE_OLD_CODE = false;
    private static final boolean USE_DELEGATE = true;
    private static final boolean USE_BEAGLE = true;
    public static final String USE_AMBIGUITIES = "useAmbiguities";
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

        boolean useMAP = xo.getAttribute(MAP_RECONSTRUCTION, false);
        boolean useMarginalLikelihood = xo.getAttribute(MARGINAL_LIKELIHOOD, true);
        boolean conditionalProbabilitiesInLogSpace = xo.getAttribute(CONDITIONAL_PROBABILITIES_IN_LOG_SPACE, false);
        String tag = xo.getAttribute(RECONSTRUCTION_TAG_NAME, RECONSTRUCTION_TAG);

        boolean useAmbiguities = xo.getAttribute(USE_AMBIGUITIES, true);
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        PatternList patternList = (PatternList) xo.getChild(PatternList.class);
        DataType dataType = patternList.getDataType();
        
        if (patternList.areUnique()) {
            throw new XMLParseException("Ancestral state reconstruction cannot be used with compressed (unique) patterns.");
        }

        MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
        GeneralSubstitutionModel generalSubstitutionModel = (GeneralSubstitutionModel) xo.getChild(GeneralSubstitutionModel.class);

        Parameter popSizes = xo.hasChildNamed(POPSIZES) ? (Parameter) xo.getElementFirstChild(POPSIZES) : null;
        Parameter logPopSizes = xo.hasChildNamed(LOG_POP_SIZES) ? (Parameter) xo.getElementFirstChild(LOG_POP_SIZES) : null;
        Parameter r = xo.hasChildNamed(GROWTH_RATES) ? (Parameter) xo.getElementFirstChild(GROWTH_RATES) : null;
        Parameter ancestralR = xo.hasChildNamed(ANCESTRAL_GROWTH_RATE) ? (Parameter) xo.getElementFirstChild(ANCESTRAL_GROWTH_RATE) : null;
        Parameter anchorProportion = xo.hasChildNamed(ANCHOR_PROPORTION) ? (Parameter) xo.getElementFirstChild(ANCHOR_PROPORTION) : null;
        Parameter gridPoints = xo.hasChildNamed(GRID_POINTS) ? (Parameter) xo.getElementFirstChild(GRID_POINTS) : null;

        if (popSizes == null && logPopSizes == null) {
            throw new XMLParseException("Either " + POPSIZES + " or " + LOG_POP_SIZES + " must be specified");
        }
        if (popSizes != null && logPopSizes != null) {
            throw new XMLParseException("Cannot specify both " + POPSIZES + " and " + LOG_POP_SIZES);
        }

        int stateCount = generalSubstitutionModel.getDataType().getStateCount();

        if (logPopSizes != null) {
            if (r == null) {
                throw new XMLParseException(LOG_POP_SIZES + " requires " + GROWTH_RATES + " to be specified");
            }
            int etaDim = logPopSizes.getDimension();
            if (!(etaDim == 1 || etaDim == stateCount)) {
                throw new XMLParseException(LOG_POP_SIZES + " dimension must be 1 or stateCount (" +
                        stateCount + "), got " + etaDim);
            }
            int rDim = r.getDimension();
            if (!(rDim == 1 || rDim == stateCount)) {
                throw new XMLParseException(GROWTH_RATES + " dimension must be 1 or stateCount (" +
                        stateCount + "), got " + rDim);
            }
        }

        if (popSizes != null) {
            if (r != null && ancestralR == null) {
                if (popSizes.getDimension() != r.getDimension()) {
                    throw new XMLParseException("Mismatch between popsizes and growth rates.");
                }
            }

            if (popSizes.getDimension() != 1) {
                if (gridPoints != null) {
                    int numGridSegments = gridPoints.getDimension() + 1;
                    int expectedDim = stateCount * numGridSegments;
                    if (popSizes.getDimension() != expectedDim) {
                        throw new XMLParseException("With grid points, popSizes dimension must be " +
                                expectedDim + " (stateCount=" + stateCount + " × (numGridPoints+1)=" +
                                numGridSegments + "), but got " + popSizes.getDimension());
                    }
                } else if (popSizes.getDimension() % stateCount != 0) {
                    throw new XMLParseException("popSizes dimension must be a multiple of stateCount (" +
                            stateCount + "), but got " + popSizes.getDimension());
                }
            }
        }

        int threads = xo.getAttribute(THREADS, 1);

        if (treeModel != null) {
            try {
                if (USE_OLD_CODE) {
                    if (popSizes == null) {
                        throw new XMLParseException(LOG_POP_SIZES + " is not supported with the old code path");
                    }
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

                        AbstractPopulationSizeModel populationSizeModel;

                        if (ancestralR != null) {
                            if (r == null) {
                                throw new XMLParseException(ANCESTRAL_GROWTH_RATE + " requires " + GROWTH_RATES +
                                        " (per-deme recent growth rates) to be specified");
                            }
                            if (ancestralR.getDimension() != 1) {
                                throw new XMLParseException(ANCESTRAL_GROWTH_RATE + " must have dimension 1");
                            }

                            boolean inputIsLog = logPopSizes != null;
                            Parameter anchorSize = inputIsLog ? logPopSizes : popSizes;
                            XMLObject anchorElement = xo.getChild(inputIsLog ? LOG_POP_SIZES : POPSIZES);
                            if (anchorSize.getDimension() != 1) {
                                throw new XMLParseException("The shared anchor population size (" +
                                        (inputIsLog ? LOG_POP_SIZES : POPSIZES) +
                                        ") must have dimension 1 for the shared-ancestral exponential model, got " +
                                        anchorSize.getDimension());
                            }

                            if (anchorProportion != null) {
                                if (!(anchorProportion.getDimension() == 1
                                        || anchorProportion.getDimension() == stateCount)) {
                                    throw new XMLParseException(ANCHOR_PROPORTION + " must have dimension 1 (shared anchor) or " +
                                            "stateCount (" + stateCount + ", per-deme anchor), got " + anchorProportion.getDimension());
                                }
                                populationSizeModel = new SharedAncestralExponentialGrowthPopulationSizeModel(
                                        "sharedAncestralExponential", anchorSize, ancestralR, r,
                                        anchorProportion, treeModel, stateCount, -1, inputIsLog);
                            } else if (anchorElement.hasAttribute(ANCHOR_TIME)) {
                                double anchorTime = anchorElement.getDoubleAttribute(ANCHOR_TIME);
                                populationSizeModel = new SharedAncestralExponentialGrowthPopulationSizeModel(
                                        "sharedAncestralExponential", anchorSize, ancestralR, r,
                                        anchorTime, stateCount, -1, inputIsLog);
                            } else {
                                throw new XMLParseException("The shared-ancestral exponential model requires either an " +
                                        ANCHOR_PROPORTION + " element (estimated anchor as a fraction of root height) or an " +
                                        ANCHOR_TIME + " attribute (fixed anchor time) on the " +
                                        (inputIsLog ? LOG_POP_SIZES : POPSIZES) + " element");
                            }
                        } else if (logPopSizes != null) {
                            XMLObject logPopSizesElement = xo.getChild(LOG_POP_SIZES);
                            if (logPopSizesElement.hasAttribute(ANCHOR_TIME)) {
                                double anchorTime = logPopSizesElement.getDoubleAttribute(ANCHOR_TIME);
                                populationSizeModel = new AnchoredExponentialGrowthPopulationSizeModel(
                                        "anchoredExponentialGrowth", logPopSizes, r,
                                        anchorTime, stateCount, -1);
                            } else {
                                populationSizeModel = new AnchoredExponentialGrowthPopulationSizeModel(
                                        "anchoredExponentialGrowth", logPopSizes, r,
                                        treeModel, stateCount, -1);
                            }
                        } else if (r != null) {
                            XMLObject popSizesElement = xo.getChild(POPSIZES);
                            if (popSizesElement.hasAttribute(ANCHOR_TIME)) {
                                double anchorTime = popSizesElement.getDoubleAttribute(ANCHOR_TIME);
                                populationSizeModel = new AnchoredExponentialGrowthPopulationSizeModel(
                                    "anchoredExponentialGrowth", popSizes, r,
                                    anchorTime, stateCount, -1, false);
                            } else if (popSizesElement.getAttribute(ROOT_ANCHORED, false)) {
                                populationSizeModel = new AnchoredExponentialGrowthPopulationSizeModel(
                                    "anchoredExponentialGrowth", popSizes, r,
                                    treeModel, stateCount, -1, false);
                            } else {
                                populationSizeModel = new ExponentialGrowthPopulationSizeModel(
                                    "exponentialGrowth", popSizes, r, stateCount, -1);
                            }
                        } else if (gridPoints != null) {
                            populationSizeModel = new PiecewiseConstantPopulationSizeModel(
                                "piecewiseConstant", popSizes, gridPoints, stateCount, -1);
                        } else if (popSizes.getDimension() == stateCount) {
                            populationSizeModel = new ConstantPopulationSizeModel(
                                "constant", popSizes, stateCount, -1);
                        } else if (popSizes.getDimension() % stateCount == 0 && popSizes.getDimension() > stateCount) {
                            populationSizeModel = new PiecewiseConstantPopulationSizeModel(
                                "piecewiseConstant", popSizes, stateCount, -1);
                        } else {
                            throw new XMLParseException(
                                "Unable to auto-detect population model. Parameter dimension " + popSizes.getDimension() + 
                                " should be " + stateCount + " (constant) or a multiple of " + stateCount + 
                                " (piecewise constant).");
                        }
                        
                        return new BastaLikelihood("name", treeModel, patternList, generalSubstitutionModel,
                                populationSizeModel, branchRateModel, delegate, subIntervals, useAmbiguities,
                                dataType, tag, useMAP, useMarginalLikelihood, conditionalProbabilitiesInLogSpace);
                    } else {
                        if (popSizes == null) {
                            throw new XMLParseException(LOG_POP_SIZES + " is not supported without the delegate code path");
                        }
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
        return BastaLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(SUBINTERVALS, true),
            AttributeRule.newIntegerRule(THREADS, true),
            AttributeRule.newBooleanRule(MAP_RECONSTRUCTION, true),
            AttributeRule.newBooleanRule(MARGINAL_LIKELIHOOD, true),
            AttributeRule.newBooleanRule(CONDITIONAL_PROBABILITIES_IN_LOG_SPACE, true),
            AttributeRule.newStringRule(RECONSTRUCTION_TAG_NAME, true),
            AttributeRule.newBooleanRule(USE_AMBIGUITIES, true),
            new ElementRule(PatternList.class),
            new ElementRule(MutableTreeModel.class),
            new ElementRule(BranchRateModel.class, true),
            new ElementRule(SubstitutionModel.class, true),
            new ElementRule(POPSIZES,
                    new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(ANCHOR_TIME, true),
                            AttributeRule.newBooleanRule(ROOT_ANCHORED, true),
                            new ElementRule(Parameter.class)}, true),
            new ElementRule(LOG_POP_SIZES,
                    new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(ANCHOR_TIME, true),
                            AttributeRule.newBooleanRule(ROOT_ANCHORED, true),
                            new ElementRule(Parameter.class)}, true),
            new ElementRule(GROWTH_RATES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "Growth rates (optional)", true),
            new ElementRule(ANCESTRAL_GROWTH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(ANCHOR_PROPORTION,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(GRID_POINTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };

}