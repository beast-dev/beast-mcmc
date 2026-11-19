/*
 * ListedPairSwapOperatorParser.java
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
import dr.inference.operators.ListedPairSwapOperator;
import dr.inference.operators.MCMCOperator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class ListedPairSwapOperatorParser extends AbstractXMLObjectParser {

    private static final String PAIRS = "pairs";

    public String getParserName() {
        return ListedPairSwapOperator.OPERATOR_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<Parameter> parameters = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            parameters.add((Parameter) xo.getChild(i));
        }
        int dim = parameters.get(0).getDimension();
        for (Parameter p : parameters) {
            if (p.getDimension() != dim) {
                throw new XMLParseException("Parameter " + p.getId() + " has dimension " + p.getDimension()
                        + ", not " + dim);
            }
        }

        int[] pairs = xo.getIntegerArrayAttribute(PAIRS);
        if (pairs.length == 0 || pairs.length % 2 != 0) {
            throw new XMLParseException("'" + PAIRS + "' needs a non-empty, even number of indices");
        }
        for (int p = 0; p < pairs.length; p += 2) {
            if (pairs[p] == pairs[p + 1]) {
                throw new XMLParseException("Pair " + pairs[p] + " " + pairs[p + 1] + " swaps an index with itself");
            }
        }
        for (int i = 0; i < pairs.length; i++) {
            if (pairs[i] < 1 || pairs[i] > dim) {
                throw new XMLParseException("Index " + pairs[i] + " is outside 1.." + dim);
            }
            pairs[i] -= 1;
        }

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new ListedPairSwapOperator(parameters, pairs, weight);
    }

    public String getParserDescription() {
        return "Swaps the values at two indices drawn uniformly from a fixed list of one-based index pairs, " +
                "in every parameter given";
    }

    public Class getReturnType() {
        return ListedPairSwapOperator.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerArrayRule(PAIRS, false),
            new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
    };
}
