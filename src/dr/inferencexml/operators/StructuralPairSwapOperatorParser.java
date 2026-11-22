/*
 * StructuralPairSwapOperatorParser.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inferencexml.operators;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.StructuralPairSwapOperator;
import dr.xml.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author Marc A. Suchard
 */
public class StructuralPairSwapOperatorParser extends AbstractXMLObjectParser {

    private static final String RATES = "rates";
    private static final String LOCAL_SCALE = "localScale";
    private static final String ON_THRESHOLD = "onThreshold";
    private static final String ANCHOR_CAP = "anchorCap";
    private static final String REVERSE_PROBABILITY = "reverseProbability";
    private static final String SAME_TARGET_PROBABILITY = "sameTargetProbability";
    private static final String SAME_SOURCE_PROBABILITY = "sameSourceProbability";
    private static final String LINK_PROBABILITY = "linkProbability";
    private static final String TRACE = "trace";

    public String getParserName() {
        return StructuralPairSwapOperator.OPERATOR_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter rates = (Parameter) xo.getElementFirstChild(RATES);
        Parameter localScale = (Parameter) xo.getElementFirstChild(LOCAL_SCALE);

        int dim = rates.getDimension();
        int stateCount = (int) Math.round((1.0 + Math.sqrt(1.0 + 4.0 * dim)) / 2.0);
        if (stateCount < 3 || stateCount * (stateCount - 1) != dim) {
            throw new XMLParseException("Rates " + rates.getId() + " have dimension " + dim
                    + ", which is not K(K - 1) for an integer K >= 3");
        }
        if (localScale.getDimension() != dim) {
            throw new XMLParseException("Local scale " + localScale.getId() + " has dimension "
                    + localScale.getDimension() + ", not " + dim);
        }

        double onThreshold = xo.getAttribute(ON_THRESHOLD, 0.05);
        double anchorCap = xo.getAttribute(ANCHOR_CAP, 1.0);
        if (!(anchorCap > 0.0)) {
            throw new XMLParseException("'" + ANCHOR_CAP + "' must be positive");
        }

        double reverse = xo.getAttribute(REVERSE_PROBABILITY, 0.5);
        double sameTarget = xo.getAttribute(SAME_TARGET_PROBABILITY, 0.25);
        double sameSource = xo.getAttribute(SAME_SOURCE_PROBABILITY, 0.25);
        if (reverse < 0.0 || sameTarget < 0.0 || sameSource < 0.0
                || Math.abs(reverse + sameTarget + sameSource - 1.0) > 1E-10) {
            throw new XMLParseException("'" + REVERSE_PROBABILITY + "', '" + SAME_TARGET_PROBABILITY + "' and '"
                    + SAME_SOURCE_PROBABILITY + "' must be non-negative and sum to 1");
        }

        double link = xo.getAttribute(LINK_PROBABILITY, 0.8);
        if (link < 0.0 || link > 1.0) {
            throw new XMLParseException("'" + LINK_PROBABILITY + "' must be in [0, 1]");
        }

        PrintWriter trace = null;
        if (xo.hasAttribute(TRACE)) {
            File file = XMLParser.getLogFile(xo, TRACE);
            try {
                trace = new PrintWriter(new FileWriter(file), true);
            } catch (IOException e) {
                throw new XMLParseException("Cannot open '" + file.getAbsolutePath() + "' for " + getParserName());
            }
        }

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new StructuralPairSwapOperator(rates, localScale, stateCount, onThreshold, anchorCap,
                reverse, sameTarget, link, trace, weight);
    }

    public String getParserDescription() {
        return "Swaps an asymmetric CTMC rate (q = theta^2) and its local scale with its reverse, a rate into the " +
                "same target or a rate out of the same source, with the anchor drawn proportional to " +
                "min(q, anchorCap)^2 and partners favoured along triangles of on rates";
    }

    public Class getReturnType() {
        return StructuralPairSwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newDoubleRule(ON_THRESHOLD, true),
            AttributeRule.newDoubleRule(ANCHOR_CAP, true),
            AttributeRule.newDoubleRule(REVERSE_PROBABILITY, true),
            AttributeRule.newDoubleRule(SAME_TARGET_PROBABILITY, true),
            AttributeRule.newDoubleRule(SAME_SOURCE_PROBABILITY, true),
            AttributeRule.newDoubleRule(LINK_PROBABILITY, true),
            AttributeRule.newStringRule(TRACE, true),
            new ElementRule(RATES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(LOCAL_SCALE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
    };
}
