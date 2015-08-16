/*
 * LikelihoodBenchmarkerParser.java
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

package dr.inferencexml.model;

import dr.inference.model.Likelihood;
import dr.inference.model.LikelihoodBenchmarker;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 */
public class LikelihoodBenchmarkerParser extends AbstractXMLObjectParser {

    public static final String BENCHMARKER = "benchmarker";

    public String getParserName() {
        return BENCHMARKER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int iterationCount = 1000;

        if (xo.hasAttribute("iterationCount")) {
            iterationCount = xo.getIntegerAttribute("iterationCount");
        }

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            Object xco = xo.getChild(i);
            if (xco instanceof Likelihood) {
                likelihoods.add((Likelihood) xco);
            }
        }

        if (likelihoods.size() == 0) {
            throw new XMLParseException("No likelihoods for benchmarking");
        }

        return new LikelihoodBenchmarker(likelihoods, iterationCount);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element runs a benchmark on a series of likelihood calculators.";
    }

    public Class getReturnType() {
        return Likelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule("iterationCount"),
            new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE)
    };

}
