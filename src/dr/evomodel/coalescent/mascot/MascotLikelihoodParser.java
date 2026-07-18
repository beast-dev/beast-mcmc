/*
 * MascotLikelihoodParser.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parses {@code <mascotLikelihood>}. Migration rates and population sizes are
 * declared as two separate elements ({@code <theta>} and {@code <popSizes>},
 * the latter matching BASTA's {@code <popSizes>} convention) rather than one
 * combined flat vector, and tip states are resolved from an
 * {@code <attributePatterns>}-style {@link PatternList} by taxon name (like
 * BASTA's own tip-state input) rather than a raw, tree-traversal-ordered
 * {@code <tipStates>} parameter.
 */
public final class MascotLikelihoodParser extends AbstractXMLObjectParser {

    public static final String STATE_COUNT = "stateCount";
    public static final String MAX_STEP = "maxStep";
    public static final String CHECK_PROBABILITIES = "checkProbabilities";
    public static final String THETA = "theta";
    public static final String POPULATION_SIZES = "popSizes";
    public static final String EPOCH_TIMES = "epochTimes";

    @Override
    public String getParserName() {
        return MascotLikelihood.MASCOT_LIKELIHOOD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        PatternList tipPatterns = (PatternList) xo.getChild(PatternList.class);
        Parameter migrationRates = (Parameter) xo.getElementFirstChild(THETA);
        Parameter popSizes = (Parameter) xo.getElementFirstChild(POPULATION_SIZES);
        Parameter epochTimes = xo.hasChildNamed(EPOCH_TIMES) ?
                (Parameter) xo.getElementFirstChild(EPOCH_TIMES) : null;
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        int stateCount = xo.getIntegerAttribute(STATE_COUNT);
        double maxStep = xo.getDoubleAttribute(MAX_STEP);
        boolean checkProbabilities = xo.getAttribute(CHECK_PROBABILITIES, false);

        return new MascotLikelihood(xo.getId(), treeModel, tipPatterns, migrationRates, popSizes, epochTimes,
                stateCount, maxStep, checkProbabilities, branchRateModel);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(STATE_COUNT),
            AttributeRule.newDoubleRule(MAX_STEP),
            AttributeRule.newBooleanRule(CHECK_PROBABILITIES, true),
            new ElementRule(TreeModel.class),
            new ElementRule(PatternList.class),
            new ElementRule(THETA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(BranchRateModel.class, true)
    };

    @Override
    public String getParserDescription() {
        return "A MASCOT-style marginal structured coalescent likelihood, with migration rates, " +
                "population sizes, and tip states declared as separate elements matching BASTA's " +
                "XML conventions.";
    }

    @Override
    public Class getReturnType() {
        return MascotLikelihood.class;
    }
}
