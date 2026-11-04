/*
 * BaselineIncrementFieldParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RandomField;
import dr.inference.model.Parameter;
import dr.math.distributions.BaselineIncrementField;
import dr.math.distributions.Distribution;
import dr.xml.*;

import static dr.inferencexml.distribution.RandomFieldParser.WEIGHTS_RULE;
import static dr.inferencexml.distribution.RandomFieldParser.parseWeightProvider;

public class BaselineIncrementFieldParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "baselineIncrementField";
    private static final String BASELINE = "baseline";
    private static final String INCREMENTS = "increments";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Distribution baseline = (Distribution) xo.getElementFirstChild(BASELINE);
        Distribution increments = (Distribution) xo.getElementFirstChild(INCREMENTS);

        RandomField.WeightProvider weights = parseWeightProvider(xo, 0);

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new BaselineIncrementField(id, baseline, increments, weights);
    }



    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BASELINE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(INCREMENTS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            WEIGHTS_RULE,
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return BaselineIncrementField.class; }
}
