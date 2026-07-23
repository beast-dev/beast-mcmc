/*
 * BayesianBridgeMarkovRandomFieldParser.java
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
import dr.inference.distribution.shrinkage.BayesianBridgeDistributionModel;
import dr.inference.distribution.shrinkage.JointBayesianBridgeDistributionModel;
import dr.inference.model.Parameter;
import dr.math.distributions.BayesianBridgeMarkovRandomField;
import dr.xml.*;

import static dr.inferencexml.distribution.RandomFieldParser.WEIGHTS_RULE;
import static dr.inferencexml.distribution.RandomFieldParser.parseWeightProvider;

public class BayesianBridgeMarkovRandomFieldParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "newBayesianBridgeMarkovRandomField";
    private static final String MEAN = "mean";
    private static final String MATCH_PSEUDO_DETERMINANT = "matchPseudoDeterminant";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {


        Parameter mean = xo.hasChildNamed(MEAN) ?
                (Parameter) xo.getElementFirstChild(MEAN) : null;

        JointBayesianBridgeDistributionModel bayesBridge = (JointBayesianBridgeDistributionModel)
                xo.getChild(JointBayesianBridgeDistributionModel.class);

        final int dim = bayesBridge.getLocalScale().getDimension();

        RandomField.WeightProvider weights = parseWeightProvider(xo, dim);

        boolean matchPseudoDeterminant = xo.getAttribute(MATCH_PSEUDO_DETERMINANT, false);

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        return new BayesianBridgeMarkovRandomField(id, bayesBridge, mean, weights, matchPseudoDeterminant);
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BayesianBridgeDistributionModel.class),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            AttributeRule.newBooleanRule(MATCH_PSEUDO_DETERMINANT, true),
            WEIGHTS_RULE,
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return RandomField.class; }
}
