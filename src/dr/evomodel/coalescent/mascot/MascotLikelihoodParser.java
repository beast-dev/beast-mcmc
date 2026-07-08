/*
 * MascotLikelihoodParser.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public final class MascotLikelihoodParser extends AbstractXMLObjectParser {

    public static final String STATE_COUNT = "stateCount";
    public static final String MAX_STEP = "maxStep";
    public static final String CHECK_PROBABILITIES = "checkProbabilities";
    public static final String TIP_STATES = "tipStates";
    public static final String THETA = "theta";
    public static final String EPOCH_TIMES = "epochTimes";

    @Override
    public String getParserName() {
        return MascotLikelihood.MASCOT_LIKELIHOOD;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        Parameter tipStates = (Parameter) xo.getElementFirstChild(TIP_STATES);
        Parameter theta = (Parameter) xo.getElementFirstChild(THETA);
        Parameter epochTimes = xo.hasChildNamed(EPOCH_TIMES) ?
                (Parameter) xo.getElementFirstChild(EPOCH_TIMES) : null;

        int stateCount = xo.getIntegerAttribute(STATE_COUNT);
        double maxStep = xo.getDoubleAttribute(MAX_STEP);
        boolean checkProbabilities = xo.getAttribute(CHECK_PROBABILITIES, false);

        return new MascotLikelihood(xo.getId(), treeModel, tipStates, theta, epochTimes,
                stateCount, maxStep, checkProbabilities);
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
            new ElementRule(TIP_STATES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(THETA, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true)
    };

    @Override
    public String getParserDescription() {
        return "A MASCOT-style marginal structured coalescent likelihood with a flat log-parameter vector.";
    }

    @Override
    public Class getReturnType() {
        return MascotLikelihood.class;
    }
}
