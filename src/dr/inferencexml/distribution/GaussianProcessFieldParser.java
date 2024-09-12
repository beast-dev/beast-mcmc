/*
 * GaussianProcessFieldParser.java
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
import dr.inference.model.DesignMatrix;
import dr.inference.model.Parameter;
import dr.math.distributions.gp.AdditiveGaussianProcessDistribution;
import dr.math.distributions.gp.GaussianProcessKernel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.math.distributions.gp.AdditiveGaussianProcessDistribution.BasisDimension;
import static dr.inferencexml.distribution.RandomFieldParser.WEIGHTS_RULE;

public class GaussianProcessFieldParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "gaussianProcessField";
    private static final String DIMENSION = "dim";
    private static final String ORDER_VARIANCE = "orderVariance";
    private static final String MEAN = "mean";
    private static final String BASIS = "basis";
    private static final String NOISE = "gaussianNoise";

    public String getParserName() { return PARSER_NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        int dim = xo.getIntegerAttribute(DIMENSION);

        Parameter orderVariance = (Parameter) xo.getElementFirstChild(ORDER_VARIANCE);

        Parameter mean = xo.hasChildNamed(MEAN) ?
                (Parameter) xo.getElementFirstChild(MEAN) : null;

        Parameter noise = xo.hasChildNamed(NOISE) ? (Parameter) xo.getElementFirstChild(NOISE) : null;

        String id = xo.hasId() ? xo.getId() : PARSER_NAME;

        List<BasisDimension> bases = parseBases(xo);

        return new AdditiveGaussianProcessDistribution(id, dim, orderVariance, mean, noise, bases);
    }

    private List<BasisDimension> parseBases(XMLObject xo) {
        List<BasisDimension> bases = new ArrayList<>();
        for (XMLObject cxo : xo.getAllChildren(BASIS)) {
            GaussianProcessKernel kernel = (GaussianProcessKernel) cxo.getChild(GaussianProcessKernel.class);

            DesignMatrix design = (DesignMatrix) cxo.getChild(DesignMatrix.class);

            if (design != null) {
                bases.add(new BasisDimension(kernel, design));
            } else {
                RandomField.WeightProvider weights = (RandomField.WeightProvider)
                        cxo.getChild(RandomField.WeightProvider.class);

                bases.add(new BasisDimension(kernel, weights));
            }
        }

        return bases;
    }

    public XMLSyntaxRule[] getSyntaxRules() { return rules; }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(ORDER_VARIANCE, Parameter.class),
            new ElementRule(MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(BASIS, new XMLSyntaxRule[] {
                    new XORRule(
                            new ElementRule(DesignMatrix.class),
                            new ElementRule(RandomField.WeightProvider.class)),
                    new ElementRule(GaussianProcessKernel.class),
                    }, 1, Integer.MAX_VALUE),
            new ElementRule(NOISE, Parameter.class, "", true),
            WEIGHTS_RULE,
    };

    public String getParserDescription() { // TODO update
        return "Describes a normal distribution with a given mean and precision " +
                "that can be used in a distributionLikelihood element";
    }

    public Class getReturnType() { return RandomField.class; }
}
