/*
 * MascotLikelihoodParser.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.PatternList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

/**
 * Parses {@code <mascotLikelihood>}. Migration rates are declared through a
 * single {@code <migrationModel>} element containing either native positive
 * rates in a {@link Parameter} or a BEAST substitution model whose realized
 * Q matrix supplies positive migration rates. Population sizes are declared separately in
 * {@code <popSizes>}, matching BASTA's convention, and tip states are resolved
 * from an
 * {@code <attributePatterns>}-style {@link PatternList} by taxon name (like
 * BASTA's own tip-state input) rather than a raw, tree-traversal-ordered
 * {@code <tipStates>} parameter.
 */
public final class MascotLikelihoodParser extends AbstractXMLObjectParser {

    public static final String STATE_COUNT = "stateCount";
    public static final String MAX_STEP = "maxStep";
    public static final String CHECK_PROBABILITIES = "checkProbabilities";
    public static final String MIGRATION_MODEL = "migrationModel";
    public static final String POPULATION_SIZES = "popSizes";
    public static final String EPOCH_TIMES = "epochTimes";
    // Alias for EPOCH_TIMES matching BASTA's PiecewiseConstantPopulationSizeModel
    // naming (same concept: strictly-increasing backward-time breakpoints).
    public static final String GRID_POINTS = "gridPoints";

    @Override
    public String getParserName() {
        return MascotLikelihood.MASCOT_LIKELIHOOD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        PatternList tipPatterns = (PatternList) xo.getChild(PatternList.class);
        Parameter migrationRates = null;
        SubstitutionModel migrationModel = null;
        Object migrationSource = xo.getElementFirstChild(MIGRATION_MODEL);
        if (migrationSource instanceof Parameter) {
            migrationRates = (Parameter) migrationSource;
        } else if (migrationSource instanceof SubstitutionModel) {
            migrationModel = (SubstitutionModel) migrationSource;
        } else {
            throw new XMLParseException("<" + MIGRATION_MODEL + "> must contain either a Parameter " +
                    "of native positive migration rates or a SubstitutionModel");
        }
        Parameter popSizes = (Parameter) xo.getElementFirstChild(POPULATION_SIZES);
        Parameter epochTimes = xo.hasChildNamed(EPOCH_TIMES) ?
                (Parameter) xo.getElementFirstChild(EPOCH_TIMES) : null;
        Parameter gridPoints = xo.hasChildNamed(GRID_POINTS) ?
                (Parameter) xo.getElementFirstChild(GRID_POINTS) : null;
        if (epochTimes != null && gridPoints != null) {
            throw new XMLParseException("Cannot specify both " + EPOCH_TIMES + " and " + GRID_POINTS +
                    " (" + GRID_POINTS + " is an alias for " + EPOCH_TIMES + ")");
        }
        if (gridPoints != null) {
            epochTimes = gridPoints;
        }
        BranchRateModel branchRateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);

        int stateCount = xo.getIntegerAttribute(STATE_COUNT);
        double maxStep = xo.getDoubleAttribute(MAX_STEP);
        boolean checkProbabilities = xo.getAttribute(CHECK_PROBABILITIES, false);

        if (migrationRates != null) {
            return new MascotLikelihood(xo.getId(), treeModel, tipPatterns, migrationRates, popSizes, epochTimes,
                    stateCount, maxStep, checkProbabilities, branchRateModel);
        }
        return new MascotLikelihood(xo.getId(), treeModel, tipPatterns, migrationModel, popSizes, epochTimes,
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
            new ElementRule(MIGRATION_MODEL, new XMLSyntaxRule[]{
                    new XORRule(
                            new ElementRule(Parameter.class),
                            new ElementRule(SubstitutionModel.class))}),
            new ElementRule(POPULATION_SIZES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(GRID_POINTS, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(BranchRateModel.class, true)
    };

    @Override
    public String getParserDescription() {
        return "A MASCOT-style marginal structured coalescent likelihood, with migration rates, " +
                "population sizes, and tip states declared as separate elements matching BASTA's XML conventions. " +
                "Migration rates may be supplied as a native positive-rate parameter or via a substitution model.";
    }

    @Override
    public Class getReturnType() {
        return MascotLikelihood.class;
    }
}
